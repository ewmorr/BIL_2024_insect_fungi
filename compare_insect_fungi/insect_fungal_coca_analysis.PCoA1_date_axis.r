##############################################################################
# Insect-Fungal Association: PCoA1 (date-associated axis) candidate strategy
#
# Alternative to the 3-axis CoCA candidate screen in
# insect_fungal_coca_analysis.general_workflow.r (results moved to
# data/compare_insects_fungi_top3axes/). That script selected candidate
# fungal taxa by their loading strength on CoCA axes 1-3 jointly, then tested
# them against insect_PCoA1/2/3 individually. This script instead selects
# candidates using ONLY insect_PCoA1 -- the axis found to be overwhelmingly a
# seasonal/date gradient (r=0.78 with collection date, see interpretation.md)
# -- to see whether a single-axis-targeted candidate list gives a different
# (larger, smaller, or differently-composed) set of fungal hits than the
# 3-axis screen, and how much of that list survives date-adjustment.
#
# cocorresp::coca() requires non-negative "community-style" data for both
# blocks (row sums > 0), so the signed insect_PCoA1 axis scores can't be
# plugged in directly as a single-column CoCA predictor block. Candidate
# selection therefore reuses the already-computed unbiased per-taxon date
# test (fungal_community_seasonality.r ->
# data/2024_fungi/fungal_taxa_date_association.all_taxa.csv, every
# prevalence-filtered fungal taxon tested against date with no pre-
# selection) and ranks candidates by that result -- "top-loading on the
# date-associated axis" operationalized as "most strongly associated with
# date". Same logic mirrored for the lure-associated axis (PCoA3) in
# insect_fungal_coca_analysis.PCoA3_lure_axis.r, using a new parallel lure
# screen (fungal_community_lure_association.r ->
# fungal_taxa_lure_association.all_taxa.csv).
#
# From there, Step 5's per-taxon permutation-test machinery (fungal
# abundance ~ site + lure? + date? + insect_PCoA1/2/3, permuting only
# insect_PCoA1 within trap) is reused unchanged, restricted to the 4
# insect_PCoA1 variants (with/without date x with/without lure) -- targeted
# to this axis only, not cross-tested against PCoA2/PCoA3.
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

out_data_dir <- "data/compare_insects_fungi_PCoA1_date"
out_fig_dir <- "figures/compare_insects_fungi_PCoA1_date"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

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
cat(length(shared_ids), "of", nrow(insect_full),
    "insect samples have matching fungal ASV data.\n")

id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

cat("Final matched dataset:", nrow(insect), "samples,", ncol(insect),
    "insect taxa (Curculionidae + Latridiidae),", ncol(fungal), "fungal ASVs.\n")

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))

## ---- 2. Filter and transform -----------------------------------------------

keep <- colSums(fungal > 0) >= 5              # present in >=5 of the samples, same filter as the CoCA/date/lure screens
fungal_f <- fungal[, keep, drop = FALSE]
cat(sum(keep), "of", ncol(fungal), "fungal taxa retained after prevalence filter.\n")

insect_hel <- decostand(insect, method = "hellinger")

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}
fungal_clr <- clr_transform(fungal_f)

## ---- 3. Candidate selection: top 200 fungal taxa by date association ------
# Reuses the unbiased, all-taxa per-taxon date test (no pre-selection on
# date correlation) already run in fungal_community_seasonality.r, rather
# than recomputing it. That file is arranged by ascending q_value.

date_assoc <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv")
stopifnot(all(date_assoc$taxon %in% colnames(fungal_f)))

n_candidates <- 200
candidate_taxa <- head(date_assoc$taxon, n_candidates)
cat("\nCandidate set: top", length(candidate_taxa), "of", nrow(date_assoc),
    "fungal taxa by association with collection date (q_value rank),",
    "from data/2024_fungi/fungal_taxa_date_association.all_taxa.csv.\n")
cat(sum(date_assoc$q_value[seq_len(n_candidates)] < 0.10), "of the top", n_candidates,
    "are individually significant at q<0.10 in that unbiased screen.\n")

## ---- 4. Insect community PCoA (same Hellinger + Euclidean strategy) -------

insect_dist <- vegdist(insect_hel, method = "euclidean")  # Hellinger + Euclidean ~ chord distance
insect_pcoa <- cmdscale(insect_dist, k = 3, eig = TRUE)
insect_axes <- as.data.frame(insect_pcoa$points)
colnames(insect_axes) <- c("insect_PCoA1", "insect_PCoA2", "insect_PCoA3")
insect_axes$sample_id <- rownames(insect)

## ---- 5. Per-taxon permutation test, targeted to insect_PCoA1 --------------
# Same model/permutation machinery as insect_fungal_coca_analysis.general_
# workflow.r Step 5: site + lure? + date? + insect_PCoA1/2/3 (all three axes
# mutually controlled for), permuting only insect_PCoA1 within trap. Run in
# all 4 combinations of with/without date x with/without lure, but only for
# insect_PCoA1 -- this candidate set was selected specifically to answer
# "which fungal taxa track the date-associated axis," so it isn't
# cross-tested against PCoA2/PCoA3.

dat_base <- meta %>% left_join(insect_axes, by = "sample_id")
stopifnot(all(dat_base$sample_id == rownames(fungal_clr)))

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
    left_join(date_assoc %>% select(taxon, date_t_stat = t_stat, date_q_value = q_value), by = "taxon") %>%
    arrange(q_value)
}

pcoa_test_specs <- list(
  list(include_date = FALSE, include_lure = TRUE,  file = "fungal_insect_association_results.PCoA1_only.csv"),
  list(include_date = TRUE,  include_lure = TRUE,  file = "fungal_insect_association_results.PCoA1_plus_date.csv"),
  list(include_date = FALSE, include_lure = FALSE, file = "fungal_insect_association_results.PCoA1_only.no_lure.csv"),
  list(include_date = TRUE,  include_lure = FALSE, file = "fungal_insect_association_results.PCoA1_plus_date.no_lure.csv")
)

cat("\nRunning per-taxon permutation test on the", length(candidate_taxa),
    "date-ranked candidate taxa against insect_PCoA1 (4 variants)...\n")

pcoa_results <- lapply(pcoa_test_specs, function(spec) {
  cat("\n--- insect_PCoA1 association", if (spec$include_date) "WITH" else "WITHOUT", "date covariate,",
      if (spec$include_lure) "WITH" else "WITHOUT", "lure covariate ---\n")
  res <- run_pcoa_association_test(test_axis = "insect_PCoA1", include_date = spec$include_date, include_lure = spec$include_lure)
  cat(sum(res$q_value < 0.10), "of", nrow(res), "taxa significant at q<0.10\n")
  print(head(res, 15))
  write.csv(res, file.path(out_data_dir, spec$file), row.names = FALSE)
  res
})
names(pcoa_results) <- sapply(pcoa_test_specs, `[[`, "file")

results_with_date <- pcoa_results[["fungal_insect_association_results.PCoA1_plus_date.csv"]]

## ---- 6. Lure partialling check (same sanity check as the original script) --

cat("\n--- Lure partialling check: does including lure change the insect_PCoA1 signal? ---\n")
summary_df <- data.frame(
  model = c("PCoA1, with lure, no date", "PCoA1, without lure, no date",
            "PCoA1, with lure + date", "PCoA1, without lure + date"),
  n_sig = c(sum(pcoa_results[["fungal_insect_association_results.PCoA1_only.csv"]]$q_value < 0.10),
            sum(pcoa_results[["fungal_insect_association_results.PCoA1_only.no_lure.csv"]]$q_value < 0.10),
            sum(pcoa_results[["fungal_insect_association_results.PCoA1_plus_date.csv"]]$q_value < 0.10),
            sum(pcoa_results[["fungal_insect_association_results.PCoA1_plus_date.no_lure.csv"]]$q_value < 0.10))
)
print(summary_df, row.names = FALSE)

t_cor_nodate <- cor(pcoa_results[["fungal_insect_association_results.PCoA1_only.csv"]]$t_stat,
                     pcoa_results[["fungal_insect_association_results.PCoA1_only.no_lure.csv"]]$t_stat[
                       match(pcoa_results[["fungal_insect_association_results.PCoA1_only.csv"]]$taxon,
                             pcoa_results[["fungal_insect_association_results.PCoA1_only.no_lure.csv"]]$taxon)])
cat("Correlation of PCoA1 t-statistics (no-date models), with vs without lure:", round(t_cor_nodate, 4), "\n")
write.csv(summary_df, file.path(out_data_dir, "lure_partialling_check.summary.csv"), row.names = FALSE)

# Direct check: do the candidate fungal taxa respond to lure at all, independent
# of the insect axes? (parametric, site-controlled, descriptive only)
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

## ---- 7. Diagnostic plot (date-adjusted model, the stricter test) ----------

make_volcano <- function(res, subtitle) {
  ggplot(res, aes(x = t_stat, y = -log10(q_value))) +
    geom_point(aes(color = q_value < 0.10)) +
    geom_hline(yintercept = -log10(0.10), linetype = "dashed") +
    labs(x = paste0("t-statistic (fungal abundance ~ ", subtitle, ")"),
         y = "-log10(FDR q-value)",
         title = paste("Fungal taxa associated with", subtitle),
         subtitle = "Candidates selected by date association (insect_PCoA1 strategy)",
         color = "q < 0.10") +
    theme_minimal()
}

volcano <- make_volcano(results_with_date, "insect PCoA1, date-adjusted")
print(volcano)
ggsave(file.path(out_fig_dir, "fungal_insect_volcano.PCoA1.png"), volcano, width = 6, height = 5, dpi = 150)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
