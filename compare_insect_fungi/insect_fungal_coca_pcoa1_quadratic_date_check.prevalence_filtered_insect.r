##############################################################################
# insect_PCoA1 CoCA Hit Check -- Quadratic Date Adjustment
# (Prevalence-Filtered Insect Table Variant)
#
# insect_pcoa_date_lure_association.prevalence_filtered_insect.r found that
# ALL THREE insect PCoA axes have a significant quadratic (non-linear) date
# component under the prevalence-filtered insect table, not just insect_PCoA2
# -- insect_PCoA1 itself has F=26.6, p=0.001 for the quadratic term beyond
# its already-known strong linear date trend (r=0.84). The PCoA2 version of
# this check (insect_fungal_coca_pcoa2_quadratic_date_check.
# prevalence_filtered_insect.r) found that PCoA2's entire fungal signal
# (including Cytospora prunicola, ASV_616) collapsed to 0/200 once linear
# date was replaced with linear+quadratic date. This script runs the
# identical check for insect_PCoA1, to see whether ITS date-adjusted
# survivors -- in particular Cytospora prunicola (ASV_1905, the top PCoA1
# candidate and the taxon flagged throughout this project as the most robust
# cross-cutting signal) -- hold up the same way, or whether PCoA1's already-
# established linear-date adjustment was similarly incomplete.
#
# Same candidate set (top-200 by CoCA loading for PCoA1, read from
# insect_fungal_coca_analysis.prevalence_filtered_insect.r's saved output),
# same trap-blocked permutation of insect_PCoA1, same site+lure+insect_PCoA1/
# 2/3 covariates, only the date covariate changes (date_c + I(date_c^2)
# instead of plain date).
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(vegan)
library(permute)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_top3axes_prevalence_filtered_insect"

## ---- 1. Load and match insect + fungal data (identical to the prevalence- ----
## ---- filtered CoCA script) -------------------------------------------------

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

keep_taxa <- colSums(sp_tab.t > 0) >= 5
insect_full <- sp_tab.t[, keep_taxa, drop = FALSE]
insect_full <- insect_full[rowSums(insect_full) > 0, ]

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy[colnames(fungal), "Kingdom"] == "k__Fungi"
fungal <- fungal[, is_fungus, drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))
cat("Matched dataset:", nrow(insect), "samples,", ncol(insect), "insect taxa,", ncol(fungal), "fungal ASVs.\n")

keep <- colSums(fungal > 0) >= 5
fungal_f <- fungal[, keep, drop = FALSE]
cat(sum(keep), "of", ncol(fungal), "fungal taxa retained after prevalence filter.\n")

insect_hel <- decostand(insect, method = "hellinger")

## ---- 2. Insect PCoA1-3 (identical construction) ----------------------------

insect_dist <- vegdist(insect_hel, method = "euclidean")
insect_pcoa <- cmdscale(insect_dist, k = 3, eig = TRUE)
insect_axes <- as.data.frame(insect_pcoa$points)
colnames(insect_axes) <- c("insect_PCoA1", "insect_PCoA2", "insect_PCoA3")
insect_axes$sample_id <- rownames(insect)

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}
fungal_clr <- clr_transform(fungal_f)

dat_base <- meta %>% left_join(insect_axes, by = "sample_id")
dat_base$date_c <- as.numeric(dat_base$date) - mean(as.numeric(dat_base$date))
stopifnot(all(dat_base$sample_id == rownames(fungal_clr)))

## ---- 3. Candidate taxa -- reuse the already-identified top-200 PCoA1 set ---

pcoa1_only_prev <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA1_only.csv"))
candidate_taxa <- pcoa1_only_prev$taxon
cat("\nRe-testing", length(candidate_taxa), "candidate fungal taxa against insect_PCoA1,",
    "now adjusting for date quadratically.\n")

## ---- 4. Per-taxon permutation test, linear + quadratic date ----------------

n_perm <- 999

test_one_taxon_quaddate <- function(taxon_abund, dat_base, n_perm, test_axis = "insect_PCoA1", block_var = "trap_id") {
  dat <- dat_base
  dat$y <- taxon_abund

  model_formula <- y ~ site + lure + date_c + I(date_c^2) + insect_PCoA1 + insect_PCoA2 + insect_PCoA3

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

res_quaddate <- bind_rows(lapply(candidate_taxa, function(tax) {
  r <- test_one_taxon_quaddate(fungal_clr[, tax], dat_base, n_perm)
  data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"])
}))
res_quaddate$q_value <- p.adjust(res_quaddate$p_perm, method = "BH")
res_quaddate <- res_quaddate %>%
  left_join(pcoa1_only_prev %>% select(taxon, coca_strength), by = "taxon") %>%
  arrange(q_value)

cat(sum(res_quaddate$q_value < 0.10), "of", nrow(res_quaddate),
    "candidate fungal taxa significant vs. insect_PCoA1 at q<0.10 (linear + quadratic date adjusted)\n")
print(head(res_quaddate, 15))
write.csv(res_quaddate, file.path(out_data_dir, "fungal_insect_association_results.PCoA1_plus_quadratic_date.csv"),
          row.names = FALSE)

## ---- 5. Compare against the no-date and linear-date-adjusted versions, and ----
## ---- specifically check Cytospora (ASV_1905) ------------------------------

pcoa1_linear_date <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA1_plus_date.csv"))

comparison <- data.frame(
  model = c("PCoA1_only (no date)", "PCoA1_plus_date (linear)", "PCoA1_plus_quadratic_date (linear+quadratic)"),
  n_sig = c(sum(pcoa1_only_prev$q_value < 0.10),
            sum(pcoa1_linear_date$q_value < 0.10),
            sum(res_quaddate$q_value < 0.10))
)
cat("\n--- insect_PCoA1 CoCA hits across date-adjustment strategies ---\n")
print(comparison, row.names = FALSE)
write.csv(comparison, file.path(out_data_dir, "PCoA1_date_adjustment_comparison.csv"), row.names = FALSE)

cat("\n--- Cytospora prunicola (ASV_1905) across all three date-adjustment strategies ---\n")
cytospora_check <- bind_rows(
  pcoa1_only_prev %>% filter(taxon == "ASV_1905") %>% mutate(model = "PCoA1_only (no date)"),
  pcoa1_linear_date %>% filter(taxon == "ASV_1905") %>% mutate(model = "PCoA1_plus_date (linear)"),
  res_quaddate %>% filter(taxon == "ASV_1905") %>% mutate(model = "PCoA1_plus_quadratic_date")
)
print(cytospora_check %>% select(model, taxon, t_stat, p_perm, q_value), row.names = FALSE)

if ("ASV_1905" %in% res_quaddate$taxon[res_quaddate$q_value < 0.10]) {
  cat("\nASV_1905 (Cytospora prunicola) SURVIVES linear+quadratic date adjustment on insect_PCoA1.\n")
} else {
  cat("\nASV_1905 (Cytospora prunicola) does NOT survive linear+quadratic date adjustment on insect_PCoA1.\n")
}

fungal_taxonomy_clean <- fungal_taxonomy %>%
  rownames_to_column("taxon") %>%
  mutate(across(Kingdom:Species, ~ sub("^[a-z]__", "", .)))

quad_hits_annotated <- res_quaddate %>%
  filter(q_value < 0.10) %>%
  left_join(fungal_taxonomy_clean, by = "taxon")
cat("\n--- Taxonomy of all taxa surviving linear+quadratic date adjustment on insect_PCoA1 ---\n")
print(quad_hits_annotated %>% select(taxon, t_stat, q_value, Family, Genus, Species), row.names = FALSE)

cat("\nDone. Outputs written to:\n",
    " ", file.path(out_data_dir, "fungal_insect_association_results.PCoA1_plus_quadratic_date.csv"), "\n",
    " ", file.path(out_data_dir, "PCoA1_date_adjustment_comparison.csv"), "\n")
