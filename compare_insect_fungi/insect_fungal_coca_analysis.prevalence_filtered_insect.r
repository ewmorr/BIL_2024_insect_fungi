##############################################################################
# Insect-Fungal Community Association Analysis -- Prevalence-Filtered Insect
# Table Variant
#
# Identical workflow to insect_fungal_coca_analysis.general_workflow.r (CoCA
# -> per-taxon insect_PCoA1-3 association tests -> direct date-association
# screens), except the insect community table is built with an individual-
# taxon >=5-sample PREVALENCE filter (insect_exploratory/
# insect_ords.prevalence_filter.R) instead of restricting to Curculionidae +
# Latridiidae. That prevalence filter was already shown to give a low-stress,
# repeatable insect NMDS (stress 0.16) without discarding taxa outside those
# two families, so this script asks whether the CoCA/insect-fungal
# association results are similar under that broader insect table.
#
# Everything else (fungal loading/filtering, CoCA method, axis-selection
# diagnostics, per-taxon permutation scheme, lure-partialling checks, date
# screens) is unchanged from the general workflow script -- see that script
# for the full rationale behind each choice. Outputs go to their own
# directory so the original Curculionidae+Latridiidae results are preserved
# for reference, not overwritten.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

required_pkgs <- c("cocorresp", "vegan", "permute",
                    "dplyr", "tidyr", "tibble", "ggplot2")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, repos = "https://cloud.r-project.org")

library(cocorresp)
library(vegan)
library(permute)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_top3axes_prevalence_filtered_insect"
out_fig_dir <- "figures/compare_insects_fungi_top3axes_prevalence_filtered_insect"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and match insect + fungal data ---------------------------------
# insect_species_tab.csv : full trap-catch table, columns
#   Class..Finest.ID + one column per sample (col_names format, e.g.
#   "C1.E..5.1").
# insect_community_metadata.csv : SequenceID, sampleID, col_names, trapID,
#   Site, Lure, CollectionDate -- links insect sample columns to the fungal
#   SequenceID used in ASV_tab.csv.
# ASV_tab.csv : fungal ITS2 ASV read counts, rows = SequenceID.

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

# Individual-taxon prevalence filter (>=5 samples), all families -- see
# insect_exploratory/insect_ords.prevalence_filter.R, in place of the
# Curculionidae+Latridiidae family restriction used in the general workflow.
sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]  # one all-zero row with no ID
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

keep_taxa <- colSums(sp_tab.t > 0) >= 5
insect_full <- sp_tab.t[, keep_taxa, drop = FALSE]
insect_full <- insect_full[rowSums(insect_full) > 0, ]
cat(sum(keep_taxa), "of", ncol(sp_tab.t),
    "insect taxa retained at >=5-sample prevalence (all families).\n")

# map insect sample columns (col_names) -> fungal SequenceID, and restrict to
# samples present in both the insect and fungal tables
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

cat("Final matched dataset:", nrow(insect), "samples,", ncol(insect),
    "insect taxa (>=5-sample prevalence filter, all families),", ncol(fungal), "fungal ASVs.\n")

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))

## ---- 2. Filter and transform -----------------------------------------------

# Drop very rare fungal taxa before testing -- singletons/doubletons add
# multiple-testing burden without power to detect anything.
keep <- colSums(fungal > 0) >= 5              # present in >=5 of the samples
fungal_f <- fungal[, keep, drop = FALSE]
cat(sum(keep), "of", ncol(fungal), "fungal taxa retained after prevalence filter.\n")

# CoCA operates on Hellinger-transformed (or similarly scaled) community data
insect_hel <- decostand(insect, method = "hellinger")
fungal_hel <- decostand(fungal_f, method = "hellinger")

## ---- 3. Co-correspondence analysis -----------------------------------------
# Predictive CoCA: fungal composition is the response (y), insect composition
# is the predictor (x). (Symmetric CoCA -- method = "symmetric" -- is the
# alternative if you don't want to designate a predictor/response direction.)

coca_full <- coca(y = fungal_hel, x = insect_hel, method = "predictive")

# --- Decide how many PLS axes to retain ---
# cocorresp provides two complementary diagnostics for this: leave-one-out
# cross-validation (the package authors' preferred method, but slow) and a
# permutation test of each axis's significance. Inspect the printed output
# of both, then set n_axes_use by hand -- there isn't a single field to read
# this off of automatically.

max_axes_to_check <- min(3, ncol(insect_hel) - 1)

cat("\nLeave-one-out cross-validation (can take a while)...\n")
# With the broader (all-families) prevalence-filtered insect table, one LOO
# fold hits a LAPACK SVD convergence failure (La.svd error code 1) partway
# through -- a known occasional numerical flakiness of dense SVD on
# ill-conditioned leave-one-out submatrices, not something worth retrying.
# This diagnostic doesn't drive n_axes_use below (hardcoded to 3, same as the
# general workflow script), so skip gracefully on failure and continue.
coca_cv <- tryCatch(
  crossval(y = fungal_hel, x = insect_hel, n.axes = max_axes_to_check),
  error = function(e) {
    cat("crossval() failed (", conditionMessage(e),
        ") -- skipping LOO-CV diagnostic; n_axes_use is set by hand below regardless.\n")
    NULL
  }
)
if (!is.null(coca_cv)) print(summary(coca_cv))

cat("\nPermutation test of axis significance...\n")
coca_perm <- tryCatch(
  permutest(coca_full, permutations = 99, n.axes = max_axes_to_check),
  error = function(e) {
    cat("permutest() failed (", conditionMessage(e), ") -- skipping.\n")
    NULL
  }
)
if (!is.null(coca_perm)) print(coca_perm)

# Same choice as the general workflow script (n_axes_use = 3), kept identical
# for direct comparability between the two insect-table variants.
n_axes_use <- 3

coca_mod <- coca(y = fungal_hel, x = insect_hel, method = "predictive", n.axes = n_axes_use)

# Fungal (response, Y-block) species scores on the retained CoCA axes = how
# strongly each fungal taxon loads on the insect-explained structure.
sc <- scores(coca_mod, display = "species", choices = seq_len(n_axes_use))
fungal_scores <- as.data.frame(sc$species$Y)
colnames(fungal_scores) <- paste0("CoCA", seq_len(n_axes_use))
fungal_scores$taxon <- rownames(fungal_scores)

# Rank by vector length across the retained axes (overall association strength,
# axis-1-only if n_axes_use == 1)
score_cols <- paste0("CoCA", seq_len(n_axes_use))
fungal_scores$coca_strength <- sqrt(rowSums(fungal_scores[, score_cols, drop = FALSE]^2))
fungal_scores <- fungal_scores %>% arrange(desc(coca_strength))

cat("\nTop 15 fungal taxa by CoCA loading strength:\n")
print(head(fungal_scores[, c("taxon", "coca_strength")], 15))

# Quick biplot for visual inspection (response block = fungal community)
pdf(file.path(out_fig_dir, "coca_biplot.pdf"), width = 7, height = 7)
plot(coca_mod, which = "response", type = "text",
     main = "Predictive CoCA: fungal (response) ~ insect (predictor, prevalence-filtered)")
dev.off()

## ---- 4. Insect community axis for use as a predictor in Step 5 -------------
# Use an independent ordination of the insect data (NMDS or PCoA) rather than
# the CoCA axes themselves, so Step 5's test isn't circular with Step 3.

insect_dist <- vegdist(insect_hel, method = "euclidean")  # Hellinger + Euclidean ~ chord distance
insect_pcoa <- cmdscale(insect_dist, k = 3, eig = TRUE)
insect_axes <- as.data.frame(insect_pcoa$points)
colnames(insect_axes) <- c("insect_PCoA1", "insect_PCoA2", "insect_PCoA3")
insect_axes$sample_id <- rownames(insect)

## ---- 5. Per-taxon mixed-model / permutation test ---------------------------
# Test whether each fungal taxon's abundance tracks the insect PCoA axes,
# accounting for the repeated-measures structure (same trap sampled across
# dates) and for site/lure as blocking factors. Same rationale/notes as the
# general workflow script (site+lure as fixed blocking factors, trap-blocked
# permutation of the axis under test, date and lure covariate variants) --
# not repeated here in full; see that script's Step 5 header comment.

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}

fungal_clr <- clr_transform(fungal_f)

dat_base <- meta %>% left_join(insect_axes, by = "sample_id")
stopifnot(all(dat_base$sample_id == rownames(fungal_clr)))

n_candidates <- 200
candidate_taxa <- head(fungal_scores$taxon, n_candidates)
cat("\nRunning per-taxon permutation test on the top", length(candidate_taxa),
    "of", ncol(fungal_clr), "fungal taxa (by CoCA loading strength).\n")

n_perm <- 999

test_one_taxon <- function(taxon_abund, dat_base, n_perm, test_axis, include_date, include_lure = TRUE, block_var = "trap_id") {
  dat <- dat_base
  dat$y <- taxon_abund

  rhs_terms <- c("site", if (include_lure) "lure", if (include_date) "date", "insect_PCoA1", "insect_PCoA2", "insect_PCoA3")
  model_formula <- as.formula(paste("y ~", paste(rhs_terms, collapse = " + ")))

  fit_stat <- function(d) {
    m <- lm(model_formula, data = d)
    coef(summary(m))[test_axis, "t value"]
  }

  obs_t <- fit_stat(dat)

  ctrl <- how(within = Within(type = "free"), blocks = dat[[block_var]], nperm = n_perm)
  perm_ids <- shuffleSet(nrow(dat), control = ctrl)

  # permute only the axis under test; hold the other PCoA axis (and all other
  # covariates) fixed, so the null is specific to that axis's association
  perm_t <- apply(perm_ids, 1, function(idx) {
    d2 <- dat
    d2[[test_axis]] <- dat[[test_axis]][idx]
    fit_stat(d2)
  })

  p_perm <- (sum(abs(perm_t) >= abs(obs_t)) + 1) / (n_perm + 1)
  c(t_stat = obs_t, p_perm = p_perm)
}

run_pcoa_association_test <- function(test_axis, include_date, include_lure = TRUE) {
  res <- bind_rows(lapply(candidate_taxa, function(tax) {
    r <- test_one_taxon(fungal_clr[, tax], dat_base, n_perm, test_axis = test_axis,
                         include_date = include_date, include_lure = include_lure)
    data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"])
  }))
  res$q_value <- p.adjust(res$p_perm, method = "BH")
  res %>%
    left_join(fungal_scores %>% select(taxon, coca_strength), by = "taxon") %>%
    arrange(q_value)
}

pcoa_test_specs <- list(
  list(axis = "insect_PCoA1", include_date = FALSE, include_lure = TRUE,  file = "fungal_insect_association_results.PCoA1_only.csv"),
  list(axis = "insect_PCoA1", include_date = TRUE,  include_lure = TRUE,  file = "fungal_insect_association_results.PCoA1_plus_date.csv"),
  list(axis = "insect_PCoA2", include_date = FALSE, include_lure = TRUE,  file = "fungal_insect_association_results.PCoA2_only.csv"),
  list(axis = "insect_PCoA2", include_date = TRUE,  include_lure = TRUE,  file = "fungal_insect_association_results.PCoA2_plus_date.csv"),
  list(axis = "insect_PCoA3", include_date = FALSE, include_lure = TRUE,  file = "fungal_insect_association_results.PCoA3_only.csv"),
  list(axis = "insect_PCoA3", include_date = TRUE,  include_lure = TRUE,  file = "fungal_insect_association_results.PCoA3_plus_date.csv"),
  # Sanity check: lure is a strong driver of the INSECT community but there is
  # no biological reason to expect fungi respond to lure independently (all
  # the lures contain ethanol).
  list(axis = "insect_PCoA1", include_date = FALSE, include_lure = FALSE, file = "fungal_insect_association_results.PCoA1_only.no_lure.csv"),
  list(axis = "insect_PCoA1", include_date = TRUE,  include_lure = FALSE, file = "fungal_insect_association_results.PCoA1_plus_date.no_lure.csv"),
  list(axis = "insect_PCoA3", include_date = FALSE, include_lure = FALSE, file = "fungal_insect_association_results.PCoA3_only.no_lure.csv"),
  list(axis = "insect_PCoA3", include_date = TRUE,  include_lure = FALSE, file = "fungal_insect_association_results.PCoA3_plus_date.no_lure.csv")
)

pcoa_results <- lapply(pcoa_test_specs, function(spec) {
  cat("\n---", spec$axis, "association", if (spec$include_date) "WITH" else "WITHOUT", "date covariate,",
      if (spec$include_lure) "WITH" else "WITHOUT", "lure covariate ---\n")
  res <- run_pcoa_association_test(test_axis = spec$axis, include_date = spec$include_date, include_lure = spec$include_lure)
  cat(sum(res$q_value < 0.10), "of", nrow(res), "taxa significant at q<0.10\n")
  print(head(res, 15))
  write.csv(res, file.path(out_data_dir, spec$file), row.names = FALSE)
  res
})
names(pcoa_results) <- sapply(pcoa_test_specs, `[[`, "file")

results_with_date <- pcoa_results[["fungal_insect_association_results.PCoA1_plus_date.csv"]]
results_pcoa3_with_date <- pcoa_results[["fungal_insect_association_results.PCoA3_plus_date.csv"]]

lure_partialling_check <- function(axis_label) {
  key <- function(suffix) sprintf("fungal_insect_association_results.%s%s.csv", axis_label, suffix)
  summary_df <- data.frame(
    model = paste0(axis_label, c(", with lure, no date", ", without lure, no date",
                                  ", with lure + date", ", without lure + date")),
    n_sig = c(sum(pcoa_results[[key("_only")]]$q_value < 0.10),
              sum(pcoa_results[[key("_only.no_lure")]]$q_value < 0.10),
              sum(pcoa_results[[key("_plus_date")]]$q_value < 0.10),
              sum(pcoa_results[[key("_plus_date.no_lure")]]$q_value < 0.10))
  )
  print(summary_df, row.names = FALSE)

  t_cor_nodate <- cor(pcoa_results[[key("_only")]]$t_stat,
                       pcoa_results[[key("_only.no_lure")]]$t_stat[
                         match(pcoa_results[[key("_only")]]$taxon,
                               pcoa_results[[key("_only.no_lure")]]$taxon)])
  cat("Correlation of", axis_label, "t-statistics (no-date models), with vs without lure:", round(t_cor_nodate, 4), "\n")
  invisible(summary_df)
}

cat("\n--- Lure partialling check: does including lure change the insect_PCoA1 signal? ---\n")
lure_partialling_check("PCoA1")

cat("\n--- Lure partialling check: does including lure change the insect_PCoA3 signal? ---\n")
lure_partialling_check("PCoA3")

# Direct check: do the candidate fungal taxa respond to lure at all, independent
# of the insect axes? lure is constant within trap (trap = site x lure), so
# within-trap permutation can't test it meaningfully here -- a parametric,
# site-controlled F-test is used instead as a quick descriptive check (not
# permutation-corrected, unlike the tests above).
cat("\n--- Direct fungal ~ lure check (parametric, site-controlled; descriptive only) ---\n")
lure_direct_check <- bind_rows(lapply(candidate_taxa, function(tax) {
  d <- dat_base
  d$y <- fungal_clr[, tax]
  m0 <- lm(y ~ site, data = d)
  m1 <- lm(y ~ site + lure, data = d)
  a <- anova(m0, m1)
  data.frame(taxon = tax, p_lure = a[["Pr(>F)"]][2])
}))
lure_direct_check$q_lure <- p.adjust(lure_direct_check$p_lure, method = "BH")
cat(sum(lure_direct_check$q_lure < 0.10, na.rm = TRUE), "of", nrow(lure_direct_check),
    "candidate fungal taxa show a nominally significant lure effect (q<0.10, parametric, not permutation-corrected)\n")
write.csv(lure_direct_check, file.path(out_data_dir, "fungal_taxa_lure_direct_check.csv"), row.names = FALSE)

## ---- 6. Diagnostic plots (date-adjusted models, the stricter test) --------

make_volcano <- function(res, subtitle) {
  ggplot(res, aes(x = t_stat, y = -log10(q_value))) +
    geom_point(aes(color = q_value < 0.10)) +
    geom_hline(yintercept = -log10(0.10), linetype = "dashed") +
    labs(x = paste0("t-statistic (fungal abundance ~ ", subtitle, ")"),
         y = "-log10(FDR q-value)",
         title = paste("Fungal taxa associated with", subtitle),
         color = "q < 0.10") +
    theme_minimal()
}

volcano <- make_volcano(results_with_date, "insect PCoA1, date-adjusted")
ggsave(file.path(out_fig_dir, "fungal_insect_volcano.png"), volcano, width = 6, height = 5, dpi = 150)

volcano_pcoa3 <- make_volcano(results_pcoa3_with_date, "insect PCoA3, date-adjusted, with lure")
ggsave(file.path(out_fig_dir, "fungal_insect_volcano.PCoA3.with_lure.png"), volcano_pcoa3, width = 6, height = 5, dpi = 150)

volcano_pcoa3_nolure <- make_volcano(pcoa_results[["fungal_insect_association_results.PCoA3_plus_date.no_lure.csv"]],
                                      "insect PCoA3, date-adjusted, without lure")
ggsave(file.path(out_fig_dir, "fungal_insect_volcano.PCoA3.no_lure.png"), volcano_pcoa3_nolure, width = 6, height = 5, dpi = 150)

## ---- 7. Direct seasonal (date) association ----------------------------------
# Complementary to Step 5: which individual insect and fungal taxa vary with
# collection date at all.

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

# --- Insect taxa: all prevalence-filtered taxa (153, vs. 52 for the
# Curculionidae+Latridiidae table) -- still cheap enough to test all. ---
insect_hel_df <- as.data.frame(insect_hel)
stopifnot(all(rownames(insect_hel_df) == meta$sample_id))

cat("\nTesting", ncol(insect_hel_df), "insect taxa for direct association with collection date...\n")
insect_date_results <- bind_rows(lapply(colnames(insect_hel_df), function(tax) {
  r <- test_date_taxon(insect_hel_df[[tax]], meta, n_perm)
  data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"])
}))
insect_date_results$q_value <- p.adjust(insect_date_results$p_perm, method = "BH")
insect_date_results <- insect_date_results %>% arrange(q_value)
cat(sum(insect_date_results$q_value < 0.10), "of", nrow(insect_date_results),
    "insect taxa significantly associated with date (q<0.10)\n")
print(head(insect_date_results, 15))
write.csv(insect_date_results, file.path(out_data_dir, "insect_taxa_date_association.csv"), row.names = FALSE)

# --- Fungal taxa: screen by raw correlation with date first (cheap, all
# prevalence-filtered ASVs), then confirm the top candidates with the
# permutation test. ---
date_numeric <- as.numeric(dat_base$date)
fungal_date_cor <- apply(fungal_clr, 2, function(x) cor(x, date_numeric))
fungal_date_candidates <- names(sort(abs(fungal_date_cor), decreasing = TRUE))[seq_len(n_candidates)]

cat("\nTesting top", length(fungal_date_candidates), "of", ncol(fungal_clr),
    "fungal taxa (by |correlation with date|) for direct association with collection date...\n")
fungal_date_results <- bind_rows(lapply(fungal_date_candidates, function(tax) {
  r <- test_date_taxon(fungal_clr[, tax], dat_base, n_perm)
  data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"])
}))
fungal_date_results$q_value <- p.adjust(fungal_date_results$p_perm, method = "BH")
fungal_date_results <- fungal_date_results %>% arrange(q_value)
cat(sum(fungal_date_results$q_value < 0.10), "of", nrow(fungal_date_results),
    "tested fungal taxa significantly associated with date (q<0.10)\n")
print(head(fungal_date_results, 15))
write.csv(fungal_date_results, file.path(out_data_dir, "fungal_taxa_date_association.csv"), row.names = FALSE)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
