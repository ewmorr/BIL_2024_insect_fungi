##############################################################################
# Fungal Community Seasonality
#
# The insect-fungal CoCA analysis (insect_fungal_coca_analysis.general_
# workflow.r) found that fungal abundance tracks the insect community
# mostly through a shared seasonal trend, not a direct relationship, and
# that testing was done on a pre-selected subset of taxa (top 200 by
# |correlation with date| or CoCA loading), which is appropriate for
# ranking specific taxa but not for estimating what fraction of the fungal
# community is seasonal.
#
# This script asks that question directly, two ways:
# Step 1: Per-taxon test of fungal abundance ~ date, run on EVERY
#         prevalence-filtered fungal taxon (no pre-selection), for an
#         unbiased estimate of the percentage of taxa that are
#         date-associated.
# Step 2: Whole-community view: rarefy the fungal ASV table to a common
#         depth (sequencing depth varies ~150x across samples), average the
#         resulting Bray-Curtis distance matrices, and run PERMANOVA + NMDS
#         -- mirroring the approach used for the insect community in
#         insect_exploratory/insect_ords.R (site + site:lure + site:date).
#
# Uses the same 69-sample matched insect/fungal dataset (Curculionidae +
# Latridiidae trap samples with paired ITS2 read data) as the CoCA script,
# for direct comparability, though only the fungal side is used here.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

required_pkgs <- c("vegan", "permute", "dplyr", "tidyr", "tibble", "ggplot2")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos = "https://cloud.r-project.org")

library(vegan)
library(permute)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(parallel)

source("library/library.R")

set.seed(1)

## ---- 1. Load and match insect + fungal data (same as the CoCA script) -----

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab %>% filter(Family %in% c("Curculionidae", "Latridiidae")) -> sp_tab.curcus_latris
sp_tab.curcus_latris.t <- t(sp_tab.curcus_latris %>% select(where(is.numeric)))
colnames(sp_tab.curcus_latris.t) <- sp_tab.curcus_latris$Finest.ID
sp_tab.curcus_latris.t[, colSums(sp_tab.curcus_latris.t) > 1] -> insect_full
insect_full <- insect_full[rowSums(insect_full) > 0, ]

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
cat(length(shared_ids), "of", nrow(insect_full), "insect samples have matching fungal ASV data.\n")
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

# Restrict to ASVs confirmed as Fungi -- ASV_tab.csv includes non-fungal
# (plant/animal/protist) and taxonomically unidentified ASVs; see
# ASVs_taxonomy.tsv (Kingdom column).
fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy[colnames(fungal), "Kingdom"] == "k__Fungi"
cat(sum(is_fungus), "of", ncol(fungal), "ASVs confirmed Kingdom == k__Fungi (dropping",
    sum(!is_fungus), "non-fungal/unidentified ASVs).\n")
fungal <- fungal[, is_fungus, drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(fungal) == meta$sample_id))
cat("Matched dataset:", nrow(fungal), "samples,", ncol(fungal), "fungal ASVs.\n")

## ---- 2. Per-taxon date test, unbiased (no pre-selection) -------------------

keep <- colSums(fungal > 0) >= 5              # present in >=5 of the samples, same filter as the CoCA script
fungal_f <- fungal[, keep, drop = FALSE]
cat(sum(keep), "of", ncol(fungal), "fungal taxa retained after prevalence filter.\n")

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}
fungal_clr <- clr_transform(fungal_f)

dat_base <- meta
stopifnot(all(dat_base$sample_id == rownames(fungal_clr)))

n_perm <- 999

test_date_taxon <- function(taxon_abund, dat_base, n_perm, block_var = "trap_id") {
  dat <- dat_base
  dat$y <- taxon_abund

  fit_stat <- function(d) {
    m <- lm(y ~ site + lure + date, data = d)
    coef(summary(m))["date", "t value"]
  }

  obs_t <- fit_stat(dat)

  ctrl <- how(within = Within(type = "free"), blocks = dat[[block_var]], nperm = n_perm)
  perm_ids <- shuffleSet(nrow(dat), control = ctrl)

  perm_t <- apply(perm_ids, 1, function(idx) {
    d2 <- dat
    d2$date <- dat$date[idx]
    fit_stat(d2)
  })

  p_perm <- (sum(abs(perm_t) >= abs(obs_t)) + 1) / (n_perm + 1)
  c(t_stat = obs_t, p_perm = p_perm)
}

# Testing every retained taxon (unlike the CoCA script's top-200
# pre-selected subset) is ~35x more work -- parallelized to keep it
# practical. mclapply's default mc.set.seed=TRUE reseeds each forked
# worker independently, so this isn't just replaying the same permutations
# across taxa.
n_cores <- max(1, min(8, detectCores() - 1))
cat("\nTesting all", ncol(fungal_clr), "prevalence-filtered fungal taxa (no pre-selection) for",
    "direct association with date, using", n_cores, "cores...\n")

t0 <- Sys.time()
fungal_date_all <- bind_rows(mclapply(colnames(fungal_clr), function(tax) {
  r <- test_date_taxon(fungal_clr[, tax], dat_base, n_perm)
  data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"])
}, mc.cores = n_cores))
fungal_date_all$q_value <- p.adjust(fungal_date_all$p_perm, method = "BH")
fungal_date_all <- fungal_date_all %>% arrange(q_value)
cat("Completed in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

pct_sig <- 100 * mean(fungal_date_all$q_value < 0.10)
cat(sum(fungal_date_all$q_value < 0.10), "of", nrow(fungal_date_all), "prevalence-filtered fungal taxa (",
    round(pct_sig, 1), "%) are significantly associated with date (q<0.10), unbiased",
    "(no pre-selection on date correlation)\n")
print(head(fungal_date_all, 15))
write.csv(fungal_date_all, "data/2024_fungi/fungal_taxa_date_association.all_taxa.csv", row.names = FALSE)

## ---- 3. Whole fungal community: rarefaction-based PERMANOVA (site/lure/date) ----
# Sequencing depth varies ~150x across these 69 samples, so the RAW
# (unfiltered) ASV counts are repeatedly rarefied to a common depth before
# computing distances, and the resulting distance matrices are averaged --
# the multiple_subsamples() + avg_matrix_list() approach from
# library/library.R, used elsewhere in this project's fungal community
# analyses.

rarefy_depth <- 5000   # all 69 samples exceed this (min depth = 5,120); none dropped
rarefy_iterations <- 100

cat("\nRarefying full fungal ASV table (", ncol(fungal), "ASVs,", nrow(fungal),
    "samples) to depth", rarefy_depth, "x", rarefy_iterations, "iterations...\n")
fungal_rare_list <- multiple_subsamples(x = fungal, depth = rarefy_depth, iterations = rarefy_iterations)
cat("Retained", nrow(fungal_rare_list[[1]]), "of", nrow(fungal), "samples after depth filtering.\n")

fungal_rare_dist_list <- lapply(fungal_rare_list, vegdist, method = "bray")
fungal_rare_dist_avg <- avg_matrix_list(fungal_rare_dist_list)

# meta rows must match the (possibly depth-filtered) rarefied sample set/order
meta_rare <- meta[match(rownames(fungal_rare_dist_avg), meta$sample_id), ]
stopifnot(all(rownames(fungal_rare_dist_avg) == meta_rare$sample_id))

cat("\n--- PERMANOVA: whole fungal community ~ site + site:lure + site:date (rarefied, averaged Bray-Curtis) ---\n")
fungal_permanova <- adonis2(as.dist(fungal_rare_dist_avg) ~ site + lure + date, data = meta_rare, by = "margin")
print(fungal_permanova)
write.csv(as.data.frame(fungal_permanova), "data/2024_fungi/fungal_community_permanova.csv")

cat("\nRunning NMDS on the averaged rarefied Bray-Curtis distance matrix...\n")
fungal_nmds <- metaMDS(as.dist(fungal_rare_dist_avg), try = 20, trymax = 100)
fungal_nmds.sites <- as.data.frame(scores(fungal_nmds, display = "sites"))
fungal_nmds.sites$sample_id <- rownames(fungal_nmds.sites)
fungal_nmds.sites <- fungal_nmds.sites %>% left_join(meta_rare, by = "sample_id")

fungal_nmds_plot <- ggplot(fungal_nmds.sites, aes(x = NMDS1, y = NMDS2, fill = as.numeric(date), shape = site)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_fill_gradient2(low = "#b2182b", high = "#2166ac", midpoint = median(as.numeric(fungal_nmds.sites$date))) +
  annotate(geom = "text", label = paste0("stress = ", round(fungal_nmds$stress, 2)),
           x = -Inf, y = -Inf, hjust = -0.1, vjust = -1) +
  labs(fill = "Collection date") +
  theme_bw()

pdf("figures/fungal_community_NMDS.date_and_site.pdf", width = 8, height = 6)
print(fungal_nmds_plot)
dev.off()

cat("\nDone. Outputs written to data/2024_fungi/fungal_taxa_date_association.all_taxa.csv,",
    "data/2024_fungi/fungal_community_permanova.csv,",
    "figures/fungal_community_NMDS.date_and_site.pdf\n")
