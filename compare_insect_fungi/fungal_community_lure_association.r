##############################################################################
# Fungal Community Lure Association (unbiased, all taxa)
#
# Companion to fungal_community_seasonality.r's date test: instead of asking
# whether individual fungal taxa track collection date, this asks whether
# they track trap lure -- run on EVERY prevalence-filtered fungal taxon (no
# pre-selection), for an unbiased ranking usable as a candidate-selection
# screen (see insect_fungal_coca_analysis.PCoA3_lure_axis.r), the same way
# fungal_taxa_date_association.all_taxa.csv is used to rank candidates for
# the PCoA1/date strategy.
#
# Lure is a trap-level constant (exactly one trap per site x lure
# combination, sampled repeatedly over dates), so it can't be tested with
# the within-trap date-permutation used for the date test -- and unlike
# date, lure is a multi-level factor, so the per-taxon statistic is an F-test
# (site+date vs. site+lure+date) rather than a single t-statistic. Both the
# F-test and the trap-within-site permutation scheme are the same ones used
# in insect_pcoa_lure_association.r to screen insect PCoA axes for lure
# association -- applied here to fungal taxa instead of insect ordination
# axes.
#
# Uses the same 69-sample matched insect/fungal dataset (Curculionidae +
# Latridiidae trap samples with paired ITS2 read data) as the CoCA script.
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

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(fungal) == meta$sample_id))
cat("Matched dataset:", nrow(fungal), "samples,", ncol(fungal), "fungal ASVs.\n")

# Design check: the lure permutation below relies on lure being constant
# within trap (one trap per site x lure combination) -- same check as
# insect_pcoa_lure_association.r.
trap_lure_map <- meta %>% distinct(trap_id, site, lure)
stopifnot(all(table(trap_lure_map$site, trap_lure_map$lure) == 1))
cat(nrow(trap_lure_map), "traps, one per site x lure combination, as expected.\n")

## ---- 2. Per-taxon lure test, unbiased (no pre-selection) -------------------

keep <- colSums(fungal > 0) >= 5              # present in >=5 of the samples, same filter as the CoCA/date scripts
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

test_lure_taxon <- function(taxon_abund, dat_base, trap_lure_map, n_perm) {
  dat <- dat_base
  dat$y <- taxon_abund

  fit_F <- function(d) {
    m0 <- lm(y ~ site + date, data = d)
    m1 <- lm(y ~ site + lure + date, data = d)
    anova(m0, m1)[2, "F"]
  }

  obs_F <- fit_F(dat)

  perm_F <- replicate(n_perm, {
    perm_map <- trap_lure_map %>%
      group_by(site) %>%
      mutate(lure = sample(lure)) %>%
      ungroup()
    lure_lookup <- setNames(as.character(perm_map$lure), perm_map$trap_id)
    d2 <- dat
    d2$lure <- lure_lookup[d2$trap_id]
    fit_F(d2)
  })

  p_perm <- (sum(perm_F >= obs_F) + 1) / (n_perm + 1)

  full_aov <- anova(lm(y ~ site + lure + date, data = dat))
  ss_resid <- full_aov["Residuals", "Sum Sq"]
  ss_lure <- full_aov["lure", "Sum Sq"]
  r2_lure <- ss_lure / (ss_lure + ss_resid)

  c(F_stat = obs_F, r2_lure = r2_lure, p_perm = p_perm)
}

# Testing every retained taxon (unlike the CoCA script's top-200 pre-selected
# subset) is substantial work -- parallelized to keep it practical, same
# approach as fungal_community_seasonality.r's date test. mclapply's default
# mc.set.seed=TRUE reseeds each forked worker independently.
n_cores <- max(1, min(8, detectCores() - 1))
cat("\nTesting all", ncol(fungal_clr), "prevalence-filtered fungal taxa (no pre-selection) for",
    "direct association with lure, using", n_cores, "cores...\n")

t0 <- Sys.time()
fungal_lure_all <- bind_rows(mclapply(colnames(fungal_clr), function(tax) {
  r <- test_lure_taxon(fungal_clr[, tax], dat_base, trap_lure_map, n_perm)
  data.frame(taxon = tax, F_stat = r["F_stat"], r2_lure = r["r2_lure"], p_perm = r["p_perm"])
}, mc.cores = n_cores))
fungal_lure_all$q_value <- p.adjust(fungal_lure_all$p_perm, method = "BH")
fungal_lure_all <- fungal_lure_all %>% arrange(q_value)
cat("Completed in", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")

pct_sig <- 100 * mean(fungal_lure_all$q_value < 0.10)
cat(sum(fungal_lure_all$q_value < 0.10), "of", nrow(fungal_lure_all), "prevalence-filtered fungal taxa (",
    round(pct_sig, 1), "%) are significantly associated with lure (q<0.10), unbiased",
    "(no pre-selection on lure association)\n")
print(head(fungal_lure_all, 15))
write.csv(fungal_lure_all, "data/2024_fungi/fungal_taxa_lure_association.all_taxa.csv", row.names = FALSE)

cat("\nDone. Output written to data/2024_fungi/fungal_taxa_lure_association.all_taxa.csv\n",
    "(rank by q_value / F_stat -- used as the candidate-selection screen for",
    "insect_fungal_coca_analysis.PCoA3_lure_axis.r, the lure-axis-only CoCA",
    "alternative structure, parallel to how fungal_taxa_date_association.all_taxa.csv",
    "feeds insect_fungal_coca_analysis.PCoA1_date_axis.r)\n")
