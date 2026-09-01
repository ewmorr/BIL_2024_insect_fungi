##############################################################################
# All-against-all: individual insect taxa x individual fungal taxa
# -- Prevalence-Filtered Insect Table Variant (Lineage B)
#
# Lineage-B counterpart of insect_fungal_pairwise_taxon_association.r. Method,
# model, permutation scheme, FDR strategy and outputs are unchanged from that
# script -- see it for the full rationale behind each choice. The ONLY change
# is the insect community table: instead of restricting to Curculionidae +
# Latridiidae and then singleton-filtering, the insect table is built with no
# family restriction and an individual-taxon >=5-sample PREVALENCE filter on
# the full trap-catch table (the same construction as
# insect_fungal_coca_analysis.prevalence_filtered_insect.r and the rest of
# Lineage B -- 153 taxa before fungal-sample matching).
#
# insect_fungal_coca_analysis.prevalence_filtered_insect.r tests fungal taxa
# against insect_PCoA1-3 -- composite axes summarizing OVERALL insect
# community variance. That's the right test for "does this fungal taxon track
# the insect community as a whole", but it can miss a fungal taxon that forms
# a tight symbiosis with ONE insect species whose signal doesn't dominate any
# top PCoA axis. This script instead tests every prevalence-filtered fungal
# taxon against every prevalence-filtered insect taxon individually.
#
# COMPUTATIONAL NOTE (unchanged from the parent script): a naive per-pair
# lm() + apply()-over-permutations loop over this grid is prohibitive. For a
# fixed insect-taxon predictor and a fixed permutation of it, EVERY fungal
# taxon's model shares the same design matrix X -- so the entire fungal
# response matrix is fit in one vectorized OLS solve
# (t(X) %*% X)^-1 %*% t(X) %*% Y rather than one lm() call per fungal taxon.
# Verified to reproduce lm()'s t-statistics exactly. The all-families insect
# table has ~3x more predictor taxa than the Curculionidae+Latridiidae one,
# so the full grid runs in roughly 3x the parent script's ~3 min single-core.
#
# SCOPE, FIRST PASS (unchanged from the parent script): one model per pair --
# fungal_taxon ~ site + lure + date + insect_taxon -- with both lure and date
# included as covariates throughout. Directly answers "does this fungal taxon
# track this specific insect taxon's abundance, beyond a shared seasonal
# trend and beyond the site/lure blocking structure". With/without
# partialling checks can be run on interesting candidates as a follow-up.
#
# Outputs go to their own _prevalence_filtered_insect-suffixed directory so
# the Lineage-A results are preserved, not overwritten.
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

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect"
out_fig_dir <- "figures/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and match insect + fungal data --------------------------------
# Loading/matching logic identical to
# insect_fungal_coca_analysis.prevalence_filtered_insect.r Step 1 (Lineage B),
# NOT the family-filtered parent -- see that script for column documentation.

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

# Individual-taxon prevalence filter (>=5 samples), all families -- in place
# of the Curculionidae+Latridiidae family restriction used in the parent.
sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]  # one all-zero row with no ID
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

keep_taxa <- colSums(sp_tab.t > 0) >= 5
insect_full <- sp_tab.t[, keep_taxa, drop = FALSE]
insect_full <- insect_full[rowSums(insect_full) > 0, ]
cat(sum(keep_taxa), "of", ncol(sp_tab.t),
    "insect taxa retained at >=5-sample prevalence (all families).\n")

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)

shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
cat(length(shared_ids), "of", nrow(insect_full),
    "insect samples have matching fungal ASV data.\n")

id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy[colnames(fungal), "Kingdom"] == "k__Fungi"
cat(sum(is_fungus), "of", ncol(fungal), "ASVs confirmed Kingdom == k__Fungi (dropping",
    sum(!is_fungus), "non-fungal/unidentified ASVs).\n")
fungal <- fungal[, is_fungus, drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))

cat("Final matched dataset:", nrow(insect), "samples,", ncol(insect),
    "insect taxa (>=5-sample prevalence filter, all families),", ncol(fungal), "fungal ASVs.\n")

## ---- 2. Filter and transform ---------------------------------------------

# Fungal: same prevalence filter as everywhere else in this project (>=5 of
# the matched samples). Fungal side is lineage-invariant -- this is the same
# ~6,342-taxon set the parent script tests.
keep_f <- colSums(fungal > 0) >= 5
fungal_f <- fungal[, keep_f, drop = FALSE]
cat(sum(keep_f), "of", ncol(fungal), "fungal taxa retained after prevalence filter.\n")

# Insect: re-apply the >=5-sample prevalence filter on the MATCHED sample set,
# exactly as the parent script does. The full-table filter above keeps a
# taxon present in >=5 of all ~200 trap-catch samples; here we additionally
# require >=5 of the fungal-matched samples, for the parent's stated reason --
# as an individual PREDICTOR (not a pooled community member), a taxon seen in
# only a handful of the matched samples can barely vary independently of
# site/lure/trap and would yield an unstable, uninterpretable coefficient
# rather than real signal. This is the Lineage-B analogue of the parent's
# post-match keep_i step (which took the family-filtered table from 52 to 39).
keep_i <- colSums(insect > 0) >= 5
cat(sum(keep_i), "of", ncol(insect),
    "insect taxa retained after re-applying the >=5-sample prevalence filter on the matched samples.\n")

# Hellinger-transform on the FULL matched insect table (all prevalence-
# filtered, all-families taxa), THEN subset to the matched-prevalence columns
# -- so each taxon's transformed value is still "relative abundance within the
# full (153-taxon) insect community", consistent with how insect_hel is built
# in insect_fungal_coca_analysis.prevalence_filtered_insect.r, not
# renormalized against a smaller basket of taxa.
insect_hel_full <- decostand(insect, method = "hellinger")
insect_hel <- insect_hel_full[, keep_i, drop = FALSE]

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}
fungal_clr <- clr_transform(fungal_f)

stopifnot(all(rownames(fungal_clr) == meta$sample_id))
stopifnot(all(rownames(insect_hel) == meta$sample_id))

n_insect_taxa <- ncol(insect_hel)
n_fungal_taxa <- ncol(fungal_clr)
cat("\nTesting all", n_fungal_taxa, "fungal taxa x", n_insect_taxa,
    "insect taxa =", n_fungal_taxa * n_insect_taxa, "pairs.\n")

## ---- 3. Vectorized per-pair permutation test ------------------------------
# Identical to the parent script. For a fixed insect-taxon predictor column
# (observed, or one permutation of it), every fungal taxon's model
# y ~ site + lure + date + insect_taxon shares the same design matrix X. Fit
# all n_fungal_taxa responses at once via OLS matrix algebra rather than one
# lm() call per fungal taxon.
#
# Permutation scheme is the general workflow's Step 5 scheme: shuffle the
# insect-taxon predictor only WITHIN trap (preserves the repeated-measures /
# trap structure under the null), same permutation index set reused across
# all insect taxa.

date_numeric <- as.numeric(meta$date)
date_centered <- date_numeric - mean(date_numeric)  # numerical conditioning only; doesn't affect t-stats

X_fixed <- model.matrix(~ site + lure + date_centered, data = meta)
p_fixed <- ncol(X_fixed)
focal_idx <- p_fixed + 1
n <- nrow(X_fixed)

Y <- as.matrix(fungal_clr)

vectorized_t <- function(X, Y, focal_idx) {
  XtX_inv <- solve(crossprod(X))
  beta <- XtX_inv %*% crossprod(X, Y)
  resid <- Y - X %*% beta
  df <- nrow(X) - ncol(X)
  sigma2 <- colSums(resid^2) / df
  se_focal <- sqrt(sigma2 * XtX_inv[focal_idx, focal_idx])
  beta[focal_idx, ] / se_focal
}

n_perm <- 999
ctrl <- how(within = Within(type = "free"), blocks = meta$trap_id, nperm = n_perm)
perm_ids <- shuffleSet(n, control = ctrl)

t0 <- Sys.time()
pairwise_list <- vector("list", n_insect_taxa)
names(pairwise_list) <- colnames(insect_hel)

for (itax in colnames(insect_hel)) {
  x_focal <- insect_hel[, itax]
  X_obs <- cbind(X_fixed, insect_taxon = x_focal)

  chk <- tryCatch(solve(crossprod(X_obs)), error = function(e) NULL)
  if (is.null(chk)) {
    cat("SKIPPING", itax, "-- design matrix is singular (perfectly collinear with site/lure/date).\n")
    next
  }

  obs_t <- vectorized_t(X_obs, Y, focal_idx)

  exceed <- integer(n_fungal_taxa)
  for (i in seq_len(n_perm)) {
    X_perm <- cbind(X_fixed, insect_taxon = x_focal[perm_ids[i, ]])
    t_perm <- vectorized_t(X_perm, Y, focal_idx)
    exceed <- exceed + (abs(t_perm) >= abs(obs_t))
  }
  p_perm <- (exceed + 1) / (n_perm + 1)

  pairwise_list[[itax]] <- data.frame(
    insect_taxon = itax,
    fungal_taxon = colnames(Y),
    t_stat = obs_t,
    p_perm = p_perm,
    row.names = NULL
  )
}
cat("Full grid permutation test completed in",
    round(as.numeric(Sys.time() - t0, units = "mins"), 2), "min\n")

full_grid <- bind_rows(pairwise_list)

# Primary FDR correction: WITHIN each insect taxon (n_fungal_taxa tests per
# taxon) -- consistent with how the general workflow corrects within each
# family of tests sharing the same predictor. A secondary GLOBAL correction
# (across every pair at once) is also reported as the stricter bar.
full_grid <- full_grid %>%
  group_by(insect_taxon) %>%
  mutate(q_value_within_insect = p.adjust(p_perm, method = "BH")) %>%
  ungroup() %>%
  mutate(q_value_global = p.adjust(p_perm, method = "BH"))

cat("\nTotal pairs tested:", nrow(full_grid), "\n")
cat("Significant at q_value_within_insect < 0.10:", sum(full_grid$q_value_within_insect < 0.10), "\n")
cat("Significant at q_value_global < 0.10:", sum(full_grid$q_value_global < 0.10), "\n")

write.csv(full_grid, file.path(out_data_dir, "fungal_insect_pairwise_full_grid.csv"), row.names = FALSE)

## ---- 4. Taxonomic annotation and significant-hit table -------------------

insect_taxonomy <- sp_tab %>%
  select(insect_taxon = Finest.ID, insect_order = Order, insect_family = Family,
         insect_subfamily = Subfamily, insect_genus = Genus, insect_species = Species) %>%
  distinct()

fungal_tax_df <- fungal_taxonomy %>%
  rownames_to_column("fungal_taxon") %>%
  select(fungal_taxon, fungal_kingdom = Kingdom, fungal_phylum = Phylum, fungal_class = Class,
         fungal_order = Order, fungal_family = Family, fungal_genus = Genus, fungal_species = Species)

sig_hits <- full_grid %>%
  filter(q_value_within_insect < 0.10) %>%
  left_join(insect_taxonomy, by = "insect_taxon") %>%
  left_join(fungal_tax_df, by = "fungal_taxon") %>%
  mutate(direction = ifelse(t_stat > 0, "positive", "negative")) %>%
  arrange(q_value_within_insect)

write.csv(sig_hits, file.path(out_data_dir, "fungal_insect_pairwise_significant_hits.csv"), row.names = FALSE)

cat("\n", nrow(sig_hits), "significant fungal-insect taxon pairs (q_value_within_insect < 0.10)\n")
cat(length(unique(sig_hits$fungal_taxon)), "unique fungal taxa involved (of", n_fungal_taxa, "tested)\n")
cat(length(unique(sig_hits$insect_taxon)), "unique insect taxa involved (of", n_insect_taxa, "tested)\n")

## ---- 5. Summary tables --------------------------------------------------

# Per insect taxon: how many fungal partners, split by direction
per_insect_summary <- sig_hits %>%
  count(insect_taxon, insect_order, insect_family, insect_subfamily, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0)
for (col in c("positive", "negative")) if (!col %in% names(per_insect_summary)) per_insect_summary[[col]] <- 0
per_insect_summary <- per_insect_summary %>%
  mutate(total = positive + negative) %>%
  arrange(desc(total))
cat("\n--- Fungal partners per insect taxon (top 15) ---\n")
print(head(as.data.frame(per_insect_summary), 15))
write.csv(per_insect_summary, file.path(out_data_dir, "fungal_partners_per_insect_taxon.csv"), row.names = FALSE)

# Per fungal taxon: how many insect taxa it's associated with (generalist vs.
# specialist signal -- a fungal taxon linked to many insect taxa is more
# likely tracking a shared community-level gradient than a tight symbiosis).
per_fungal_summary <- sig_hits %>%
  count(fungal_taxon, fungal_phylum, fungal_class, fungal_order, fungal_family, fungal_genus, fungal_species,
        name = "n_insect_partners") %>%
  arrange(desc(n_insect_partners))
cat("\n--- Insect partners per fungal taxon (top 15) ---\n")
print(head(as.data.frame(per_fungal_summary), 15))
write.csv(per_fungal_summary, file.path(out_data_dir, "insect_partners_per_fungal_taxon.csv"), row.names = FALSE)

cat("\nDistribution of n_insect_partners among significant fungal taxa:\n")
print(table(per_fungal_summary$n_insect_partners))

# Fungal-family-level tally of significant hits, by direction
family_tally <- sig_hits %>%
  count(fungal_family, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0)
for (col in c("positive", "negative")) if (!col %in% names(family_tally)) family_tally[[col]] <- 0
family_tally <- family_tally %>%
  mutate(total = positive + negative) %>%
  arrange(desc(total))
cat("\n--- Fungal family tally among significant hits ---\n")
print(head(as.data.frame(family_tally), 20))
write.csv(family_tally, file.path(out_data_dir, "fungal_family_tally.csv"), row.names = FALSE)

## ---- 6. Diagnostic plots ----------------------------------------------------

if (nrow(sig_hits) > 0) {
  # Fungal partner count per insect taxon
  partner_bar <- ggplot(per_insect_summary,
                         aes(x = reorder(insect_taxon, total), y = total, fill = insect_family)) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "Number of significantly associated fungal taxa (q<0.10)",
         title = "Fungal taxa associated with each insect taxon (Lineage B: prevalence-filtered, all families)",
         subtitle = "site + lure + date + insect_taxon model, within-insect-taxon FDR",
         fill = "Insect family") +
    theme_minimal()
  ggsave(file.path(out_fig_dir, "fungal_partners_per_insect_taxon.png"),
         partner_bar, width = 8, height = max(4, 0.22 * nrow(per_insect_summary)), dpi = 150, limitsize = FALSE)

  # Heatmap: top fungal families x insect taxa, count of significant hits
  top_families <- head(family_tally$fungal_family, 20)
  heat_df <- sig_hits %>%
    filter(fungal_family %in% top_families) %>%
    count(insect_taxon, fungal_family)

  heatmap_plot <- ggplot(heat_df, aes(x = insect_taxon, y = reorder(fungal_family, n, sum), fill = n)) +
    geom_tile() +
    scale_fill_gradient(low = "grey90", high = "#b2182b") +
    labs(x = NULL, y = NULL, fill = "n ASVs",
         title = "Significant fungal family x insect taxon associations (Lineage B)",
         subtitle = "Top 20 fungal families by total hit count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 8))
  ggsave(file.path(out_fig_dir, "family_by_insect_taxon_heatmap.png"),
         heatmap_plot, width = 12, height = 8, dpi = 150)

  # Distribution of specialist vs. generalist fungal signal
  specialist_hist <- ggplot(per_fungal_summary, aes(x = n_insect_partners)) +
    geom_bar(fill = "#2166ac") +
    labs(x = "Number of insect taxa a fungal ASV is significantly associated with",
         y = "Number of fungal ASVs",
         title = "Specialist (n=1) vs. generalist (n=many) fungal association signal (Lineage B)") +
    theme_minimal()
  ggsave(file.path(out_fig_dir, "fungal_specialist_generalist_distribution.png"),
         specialist_hist, width = 6, height = 5, dpi = 150)
} else {
  cat("\nNo significant pairs at q_value_within_insect < 0.10 -- skipping diagnostic plots.\n")
}

cat("\nDone. Outputs written to", out_data_dir, ":\n",
    "  fungal_insect_pairwise_full_grid.csv (all", nrow(full_grid), "pairs)\n",
    "  fungal_insect_pairwise_significant_hits.csv\n",
    "  fungal_partners_per_insect_taxon.csv\n",
    "  insect_partners_per_fungal_taxon.csv\n",
    "  fungal_family_tally.csv\n",
    "and", out_fig_dir, "\n")
