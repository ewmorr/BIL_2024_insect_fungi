##############################################################################
# All-against-all: individual insect taxa x individual fungal taxa
#
# insect_fungal_coca_analysis.general_workflow.r tests fungal taxa against
# insect_PCoA1-3 -- composite axes summarizing OVERALL insect community
# variance. That's the right test for "does this fungal taxon track the
# insect community as a whole", but it can miss a fungal taxon that forms a
# tight symbiosis with ONE insect species whose signal doesn't dominate any
# top PCoA axis. This script instead tests every prevalence-filtered fungal
# taxon against every prevalence-filtered insect taxon individually, i.e. a
# 6,342 (fungal ASVs) x 39 (insect taxa) = 247,338-pair grid.
#
# COMPUTATIONAL NOTE: a naive port of the general workflow's per-pair
# lm() + apply()-over-permutations loop to this grid was benchmarked at
# ~247M individual lm() refits -- multiple hours even parallelized (the
# closest precedent, fungal_community_seasonality.r's all-6,342-taxa x
# 1-predictor date screen, already needed 8 cores to be practical, and this
# grid is ~39x larger). Instead, for a fixed insect-taxon predictor and a
# fixed permutation of it, EVERY fungal taxon's model shares the same design
# matrix X -- so the entire 6,342-column response matrix can be fit in one
# vectorized OLS solve (t(X) %*% X)^-1 %*% t(X) %*% Y rather than one lm()
# call per fungal taxon. Benchmarked and verified to reproduce lm()'s
# t-statistics exactly (see prototype check during development); full grid
# (39 insect taxa x 999 permutations x 6,342 fungal responses each) runs in
# ~3 minutes single-core. This means the FULL grid can be tested with the
# same rigor (999-permutation, trap-blocked) as the general workflow's
# per-axis tests, with no need to pre-screen or subsample fungal candidates.
#
# SCOPE, FIRST PASS: per the project's usual practice of running the full
# with/without date and with/without lure partialling suite (see general
# workflow Step 5), that's 4x this already-large grid. As a first pass, we
# instead fit ONE model per pair -- fungal_taxon ~ site + lure + date +
# insect_taxon -- with both lure and date included as covariates throughout.
# This directly answers "does this fungal taxon track this specific insect
# taxon's abundance, beyond a shared seasonal trend and beyond the
# site/lure blocking structure". If particular pairs turn out to be
# interesting, the with/without partialling checks from the general
# workflow can be run on just those candidates as a follow-up.
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

## ---- 1. Load and match insect + fungal data --------------------------------
# Identical loading/matching logic to insect_fungal_coca_analysis.general_workflow.r
# Step 1 -- see that script for column documentation.

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
    "insect taxa (Curculionidae + Latridiidae),", ncol(fungal), "fungal ASVs.\n")

## ---- 2. Filter and transform -----------------------------------------------

# Fungal: same prevalence filter as the general workflow (>=5 of the matched
# samples) -- this is the full 6,342-taxon set, not a CoCA-loading subset.
keep_f <- colSums(fungal > 0) >= 5
fungal_f <- fungal[, keep_f, drop = FALSE]
cat(sum(keep_f), "of", ncol(fungal), "fungal taxa retained after prevalence filter.\n")

# Insect: apply the same >=5-sample prevalence filter, for the same reason
# (very rare taxa add multiple-testing burden with essentially no power --
# and, as individual PREDICTORS here rather than pooled community members,
# a taxon caught in only 1-4 samples can barely vary independently of
# site/lure/trap, which would show up as an unstable, uninterpretable
# coefficient rather than real signal).
keep_i <- colSums(insect > 0) >= 5
cat(sum(keep_i), "of", ncol(insect), "insect taxa retained after prevalence filter.\n")

# Hellinger-transform on the FULL matched insect table (all Curculionidae +
# Latridiidae taxa), THEN subset to the prevalence-filtered columns -- so
# each taxon's transformed value is still "relative abundance within the
# full insect community" (consistent with how insect_hel is used elsewhere
# in this project, e.g. as the basis for insect_PCoA1-3), not renormalized
# against a smaller basket of taxa.
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

## ---- 3. Vectorized per-pair permutation test --------------------------------
# For a fixed insect-taxon predictor column (observed, or one permutation of
# it), every fungal taxon's model y ~ site + lure + date + insect_taxon
# shares the same design matrix X. Fit all n_fungal_taxa responses at once
# via OLS matrix algebra rather than one lm() call per fungal taxon --
# verified during development to reproduce lm()'s t-statistics exactly.
#
# Permutation scheme is the same as the general workflow's Step 5: shuffle
# the insect-taxon predictor only WITHIN trap (preserves the repeated-
# measures/trap structure under the null), same permutation index set reused
# across all insect taxa (the block structure doesn't depend on which
# taxon's values are being permuted, so this is both valid and avoids
# redundant calls to shuffleSet).

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

# Primary FDR correction: WITHIN each insect taxon (6,342 tests per taxon) --
# consistent with how the general workflow corrects within each PCoA-axis
# spec, i.e. within each family of tests sharing the same predictor/
# hypothesis, rather than pooling every predictor's tests into one
# correction. A secondary GLOBAL correction (across all 247k tests at once)
# is also reported for reference/context, since it's the stricter, more
# conservative bar appropriate to "how many of these associations would
# survive testing everything against everything".
full_grid <- full_grid %>%
  group_by(insect_taxon) %>%
  mutate(q_value_within_insect = p.adjust(p_perm, method = "BH")) %>%
  ungroup() %>%
  mutate(q_value_global = p.adjust(p_perm, method = "BH"))

cat("\nTotal pairs tested:", nrow(full_grid), "\n")
cat("Significant at q_value_within_insect < 0.10:", sum(full_grid$q_value_within_insect < 0.10), "\n")
cat("Significant at q_value_global < 0.10:", sum(full_grid$q_value_global < 0.10), "\n")

write.csv(full_grid, "data/compare_insects_fungi_pairwise_taxa/fungal_insect_pairwise_full_grid.csv", row.names = FALSE)

## ---- 4. Taxonomic annotation and significant-hit table ---------------------

insect_taxonomy <- sp_tab.curcus_latris %>%
  select(insect_taxon = Finest.ID, insect_family = Family, insect_subfamily = Subfamily,
         insect_genus = Genus, insect_species = Species) %>%
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

write.csv(sig_hits, "data/compare_insects_fungi_pairwise_taxa/fungal_insect_pairwise_significant_hits.csv", row.names = FALSE)

cat("\n", nrow(sig_hits), "significant fungal-insect taxon pairs (q_value_within_insect < 0.10)\n")
cat(length(unique(sig_hits$fungal_taxon)), "unique fungal taxa involved (of", n_fungal_taxa, "tested)\n")
cat(length(unique(sig_hits$insect_taxon)), "unique insect taxa involved (of", n_insect_taxa, "tested)\n")

## ---- 5. Summary tables -------------------------------------------------------

# Per insect taxon: how many fungal partners, split by direction
per_insect_summary <- sig_hits %>%
  count(insect_taxon, insect_family, insect_subfamily, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(total = positive + negative) %>%
  arrange(desc(total))
cat("\n--- Fungal partners per insect taxon (top 15) ---\n")
print(head(as.data.frame(per_insect_summary), 15))
write.csv(per_insect_summary, "data/compare_insects_fungi_pairwise_taxa/fungal_partners_per_insect_taxon.csv", row.names = FALSE)

# Per fungal taxon: how many insect taxa it's associated with (generalist vs.
# specialist signal -- a fungal taxon linked to many insect taxa is more
# likely tracking a shared community-level gradient than a tight symbiosis;
# a fungal taxon linked to exactly one insect taxon is the more compelling
# candidate for a specific association).
per_fungal_summary <- sig_hits %>%
  count(fungal_taxon, fungal_phylum, fungal_class, fungal_order, fungal_family, fungal_genus, fungal_species,
        name = "n_insect_partners") %>%
  arrange(desc(n_insect_partners))
cat("\n--- Insect partners per fungal taxon (top 15) ---\n")
print(head(as.data.frame(per_fungal_summary), 15))
write.csv(per_fungal_summary, "data/compare_insects_fungi_pairwise_taxa/insect_partners_per_fungal_taxon.csv", row.names = FALSE)

cat("\nDistribution of n_insect_partners among significant fungal taxa:\n")
print(table(per_fungal_summary$n_insect_partners))

# Fungal-family-level tally of significant hits, by direction -- same style
# as the taxonomic breakdown in interpretation_update_08172026.md Section 2
family_tally <- sig_hits %>%
  count(fungal_family, direction) %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(total = positive + negative) %>%
  arrange(desc(total))
cat("\n--- Fungal family tally among significant hits ---\n")
print(head(as.data.frame(family_tally), 20))
write.csv(family_tally, "data/compare_insects_fungi_pairwise_taxa/fungal_family_tally.csv", row.names = FALSE)

## ---- 6. Diagnostic plots ----------------------------------------------------

# Fungal partner count per insect taxon
partner_bar <- ggplot(per_insect_summary,
                       aes(x = reorder(insect_taxon, total), y = total, fill = insect_subfamily)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Number of significantly associated fungal taxa (q<0.10)",
       title = "Fungal taxa associated with each insect taxon",
       subtitle = "site + lure + date + insect_taxon model, within-insect-taxon FDR",
       fill = "Insect subfamily") +
  theme_minimal()
ggsave("figures/compare_insects_fungi_pairwise_taxa/fungal_partners_per_insect_taxon.png",
       partner_bar, width = 8, height = 8, dpi = 150)

# Heatmap: top fungal families x insect taxa, count of significant hits
top_families <- head(family_tally$fungal_family, 20)
heat_df <- sig_hits %>%
  filter(fungal_family %in% top_families) %>%
  count(insect_taxon, fungal_family)

heatmap_plot <- ggplot(heat_df, aes(x = insect_taxon, y = reorder(fungal_family, n, sum), fill = n)) +
  geom_tile() +
  scale_fill_gradient(low = "grey90", high = "#b2182b") +
  labs(x = NULL, y = NULL, fill = "n ASVs",
       title = "Significant fungal family x insect taxon associations",
       subtitle = "Top 20 fungal families by total hit count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 8))
ggsave("figures/compare_insects_fungi_pairwise_taxa/family_by_insect_taxon_heatmap.png",
       heatmap_plot, width = 10, height = 8, dpi = 150)

# Distribution of specialist vs. generalist fungal signal
specialist_hist <- ggplot(per_fungal_summary, aes(x = n_insect_partners)) +
  geom_bar(fill = "#2166ac") +
  labs(x = "Number of insect taxa a fungal ASV is significantly associated with",
       y = "Number of fungal ASVs",
       title = "Specialist (n=1) vs. generalist (n=many) fungal association signal") +
  theme_minimal()
ggsave("figures/compare_insects_fungi_pairwise_taxa/fungal_specialist_generalist_distribution.png",
       specialist_hist, width = 6, height = 5, dpi = 150)

cat("\nDone. Outputs written to data/compare_insects_fungi_pairwise_taxa/:\n",
    "  fungal_insect_pairwise_full_grid.csv (all", nrow(full_grid), "pairs)\n",
    "  fungal_insect_pairwise_significant_hits.csv\n",
    "  fungal_partners_per_insect_taxon.csv\n",
    "  insect_partners_per_fungal_taxon.csv\n",
    "  fungal_family_tally.csv\n",
    "and figures/compare_insects_fungi_pairwise_taxa/:\n",
    "  fungal_partners_per_insect_taxon.png\n",
    "  family_by_insect_taxon_heatmap.png\n",
    "  fungal_specialist_generalist_distribution.png\n")
