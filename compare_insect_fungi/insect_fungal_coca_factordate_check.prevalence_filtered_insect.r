##############################################################################
# insect_PCoA1 and insect_PCoA2 CoCA Hit Check -- factor(date) Adjustment
# (Prevalence-Filtered Insect Table Variant)
#
# The linear (insect_fungal_coca_analysis.prevalence_filtered_insect.r) and
# quadratic (insect_fungal_coca_pcoa1/pcoa2_quadratic_date_check.
# prevalence_filtered_insect.r) date-adjustment checks each rule out one
# specific parametric functional form of date (a line, then a curve) as the
# source of a fungal taxon's apparent association with an insect PCoA axis.
# With only 6 distinct collection dates in this dataset, factor(date) is the
# maximally strict version of that check -- a 6-level factor absorbs ANY
# functional form of date (linear, quadratic, or otherwise), the same
# strictest-possible-deseasonalizing logic already used elsewhere in this
# project (insect_fungal_pairwise_taxon_association.residualized.r).
#
# Reruns both insect_PCoA1's and insect_PCoA2's top-200 CoCA candidate sets
# (read from insect_fungal_coca_analysis.prevalence_filtered_insect.r's saved
# output, not refit) with site + lure + factor(date) + insect_PCoA1/2/3 as
# covariates, same trap-blocked permutation of the axis under test as every
# other per-taxon test in this project.
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
cat(length(unique(meta$date)), "distinct collection dates:", paste(sort(unique(meta$date)), collapse = ", "), "\n")

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
dat_base$date_f <- factor(dat_base$date)
stopifnot(all(dat_base$sample_id == rownames(fungal_clr)))

## ---- 3. Per-taxon permutation test, factor(date) --------------------------
# Same trap-blocked permutation of the axis under test as every other
# per-taxon test in this project -- only the date covariate changes
# (factor(date), absorbing any functional form of date) instead of a linear
# or quadratic term.

n_perm <- 999

test_one_taxon_factordate <- function(taxon_abund, dat_base, n_perm, test_axis, block_var = "trap_id") {
  dat <- dat_base
  dat$y <- taxon_abund

  model_formula <- as.formula(paste("y ~ site + lure + date_f +",
                                     paste(paste0("insect_PCoA", 1:3), collapse = " + ")))

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

run_factordate_test <- function(axis_name, candidate_source_file) {
  candidates_prev <- read.csv(file.path(out_data_dir, candidate_source_file))
  candidate_taxa <- candidates_prev$taxon
  cat("\nRe-testing", length(candidate_taxa), "candidate fungal taxa against", axis_name,
      "with factor(date) adjustment...\n")

  res <- bind_rows(lapply(candidate_taxa, function(tax) {
    r <- test_one_taxon_factordate(fungal_clr[, tax], dat_base, n_perm, test_axis = axis_name)
    data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"])
  }))
  res$q_value <- p.adjust(res$p_perm, method = "BH")
  res <- res %>%
    left_join(candidates_prev %>% select(taxon, coca_strength), by = "taxon") %>%
    arrange(q_value)

  cat(sum(res$q_value < 0.10), "of", nrow(res), "candidate fungal taxa significant vs.", axis_name,
      "at q<0.10 (factor(date) adjusted)\n")
  print(head(res, 15))
  res
}

pcoa1_factordate <- run_factordate_test("insect_PCoA1", "fungal_insect_association_results.PCoA1_only.csv")
write.csv(pcoa1_factordate, file.path(out_data_dir, "fungal_insect_association_results.PCoA1_plus_factordate.csv"),
          row.names = FALSE)

pcoa2_factordate <- run_factordate_test("insect_PCoA2", "fungal_insect_association_results.PCoA2_only.csv")
write.csv(pcoa2_factordate, file.path(out_data_dir, "fungal_insect_association_results.PCoA2_plus_factordate.csv"),
          row.names = FALSE)

## ---- 4. Compare across all date-adjustment strategies for both axes -------

pcoa1_only <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA1_only.csv"))
pcoa1_linear <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA1_plus_date.csv"))
pcoa1_quad <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA1_plus_quadratic_date.csv"))

pcoa2_only <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA2_only.csv"))
pcoa2_linear <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA2_plus_date.csv"))
pcoa2_quad <- read.csv(file.path(out_data_dir, "fungal_insect_association_results.PCoA2_plus_quadratic_date.csv"))

comparison <- data.frame(
  axis = c(rep("PCoA1", 4), rep("PCoA2", 4)),
  model = rep(c("no date", "linear date", "linear+quadratic date", "factor(date)"), 2),
  n_sig = c(sum(pcoa1_only$q_value < 0.10), sum(pcoa1_linear$q_value < 0.10),
            sum(pcoa1_quad$q_value < 0.10), sum(pcoa1_factordate$q_value < 0.10),
            sum(pcoa2_only$q_value < 0.10), sum(pcoa2_linear$q_value < 0.10),
            sum(pcoa2_quad$q_value < 0.10), sum(pcoa2_factordate$q_value < 0.10))
)
cat("\n--- Hit count across all date-adjustment strategies, both axes ---\n")
print(comparison, row.names = FALSE)
write.csv(comparison, file.path(out_data_dir, "date_adjustment_comparison.all_strategies.csv"), row.names = FALSE)

## ---- 5. Specific taxa of interest: Cytospora (ASV_1905) and Tympanis (ASV_13242) ----

cat("\n--- ASV_1905 (Cytospora prunicola) across all PCoA1 date-adjustment strategies ---\n")
cytospora_check <- bind_rows(
  pcoa1_only %>% filter(taxon == "ASV_1905") %>% mutate(model = "no date"),
  pcoa1_linear %>% filter(taxon == "ASV_1905") %>% mutate(model = "linear date"),
  pcoa1_quad %>% filter(taxon == "ASV_1905") %>% mutate(model = "linear+quadratic date"),
  pcoa1_factordate %>% filter(taxon == "ASV_1905") %>% mutate(model = "factor(date)")
)
print(cytospora_check %>% select(model, t_stat, p_perm, q_value), row.names = FALSE)

cat("\n--- ASV_13242 (Tympanis sp.) across all PCoA1 date-adjustment strategies ---\n")
tympanis_check <- bind_rows(
  pcoa1_only %>% filter(taxon == "ASV_13242") %>% mutate(model = "no date"),
  pcoa1_linear %>% filter(taxon == "ASV_13242") %>% mutate(model = "linear date"),
  pcoa1_quad %>% filter(taxon == "ASV_13242") %>% mutate(model = "linear+quadratic date"),
  pcoa1_factordate %>% filter(taxon == "ASV_13242") %>% mutate(model = "factor(date)")
)
print(tympanis_check %>% select(model, t_stat, p_perm, q_value), row.names = FALSE)

## ---- 6. Taxonomy of the factor(date)-adjusted survivors --------------------

fungal_taxonomy_clean <- fungal_taxonomy %>%
  rownames_to_column("taxon") %>%
  mutate(across(Kingdom:Species, ~ sub("^[a-z]__", "", .)))

pcoa1_survivors <- pcoa1_factordate %>% filter(q_value < 0.10) %>% left_join(fungal_taxonomy_clean, by = "taxon")
pcoa2_survivors <- pcoa2_factordate %>% filter(q_value < 0.10) %>% left_join(fungal_taxonomy_clean, by = "taxon")

cat("\n--- Taxonomy: insect_PCoA1 survivors under factor(date) ---\n")
print(pcoa1_survivors %>% select(taxon, t_stat, q_value, Family, Genus, Species), row.names = FALSE)

cat("\n--- Taxonomy: insect_PCoA2 survivors under factor(date) ---\n")
if (nrow(pcoa2_survivors) > 0) {
  print(pcoa2_survivors %>% select(taxon, t_stat, q_value, Family, Genus, Species), row.names = FALSE)
} else {
  cat("(none)\n")
}

cat("\nDone. Outputs written to:\n",
    " ", file.path(out_data_dir, "fungal_insect_association_results.PCoA1_plus_factordate.csv"), "\n",
    " ", file.path(out_data_dir, "fungal_insect_association_results.PCoA2_plus_factordate.csv"), "\n",
    " ", file.path(out_data_dir, "date_adjustment_comparison.all_strategies.csv"), "\n")
