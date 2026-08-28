##############################################################################
# Insect Taxa ~ Date (Linear + Quadratic) -- Prevalence-Filtered Insect Table
#
# Standalone version of insect_fungal_coca_analysis.general_workflow.r's
# Step 7 insect-taxon date screen, broken out on its own (rather than
# rerunning the full prevalence-filtered CoCA workflow, which also fits the
# CoCA model, LOO-CV, and the ~10-variant per-taxon PCoA association grid --
# all unnecessary for just this test) so this specific comparison can be
# rerun quickly. Built on the same >=5-sample individual-taxon prevalence-
# filtered insect table (153 taxa, all families -- insect_exploratory/
# insect_ords.prevalence_filter.R) used throughout the
# *.prevalence_filtered_insect.r script family, instead of the
# Curculionidae+Latridiidae family restriction.
#
# Same linear+quadratic date test as the just-updated general_workflow.r and
# fungal_community_seasonality.r: tests both a linear date term (early-/
# late-season direction) and a quadratic term (hump/dip, non-monotonic
# seasonal pattern) in one model, permuting date within trap for both.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(vegan)
library(permute)
library(dplyr)
library(tidyr)
library(tibble)

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_top3axes_prevalence_filtered_insect"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and match insect + fungal data (identical to the prevalence- ----
## ---- filtered CoCA script -- fungal data used only to restrict the sample set) ----

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

keep_taxa <- colSums(sp_tab.t > 0) >= 5
insect_full <- sp_tab.t[, keep_taxa, drop = FALSE]
insect_full <- insect_full[rowSums(insect_full) > 0, ]
cat(sum(keep_taxa), "of", ncol(sp_tab.t), "insect taxa retained at >=5-sample prevalence (all families).\n")

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
cat(length(shared_ids), "of", nrow(insect_full), "insect samples have matching fungal ASV data.\n")
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(insect) == meta$sample_id))
cat("Matched dataset:", nrow(insect), "samples,", ncol(insect),
    "insect taxa (>=5-sample prevalence filter, all families).\n")

insect_hel <- decostand(insect, method = "hellinger")

## ---- 2. Linear + quadratic date test (identical logic to general_workflow.r) ----

n_perm <- 999

test_date_taxon <- function(taxon_abund, dat_base, n_perm, block_var = "trap_id") {
  dat <- dat_base
  dat$y <- taxon_abund

  fit_stats <- function(d) {
    d$date_c <- as.numeric(d$date) - mean(as.numeric(d$date))
    m <- lm(y ~ site + lure + date_c + I(date_c^2), data = d)
    cs <- coef(summary(m))
    c(t_linear = cs["date_c", "t value"], t_quad = cs["I(date_c^2)", "t value"])
  }

  obs <- fit_stats(dat)

  ctrl <- how(within = Within(type = "free"), blocks = dat[[block_var]], nperm = n_perm)
  perm_ids <- shuffleSet(nrow(dat), control = ctrl)

  perm_stats <- t(apply(perm_ids, 1, function(idx) {
    d2 <- dat
    d2$date <- dat$date[idx]
    fit_stats(d2)
  }))

  p_linear <- (sum(abs(perm_stats[, "t_linear"]) >= abs(obs["t_linear"])) + 1) / (n_perm + 1)
  p_quad <- (sum(abs(perm_stats[, "t_quad"]) >= abs(obs["t_quad"])) + 1) / (n_perm + 1)

  c(t_stat = unname(obs["t_linear"]), p_perm = p_linear,
    t_stat_quad = unname(obs["t_quad"]), p_perm_quad = p_quad)
}

classify_shape <- function(df) {
  df %>%
    mutate(shape = case_when(
      q_value_quad < 0.10 & t_stat_quad < 0 ~ "hump (peaks mid-season)",
      q_value_quad < 0.10 & t_stat_quad > 0 ~ "dip (troughs mid-season)",
      q_value < 0.10 & t_stat > 0 ~ "linear increase (late-season)",
      q_value < 0.10 & t_stat < 0 ~ "linear decrease (early-season)",
      TRUE ~ "no significant date pattern"
    )) %>%
    arrange(pmin(q_value, q_value_quad))
}

insect_hel_df <- as.data.frame(insect_hel)
stopifnot(all(rownames(insect_hel_df) == meta$sample_id))

cat("\nTesting", ncol(insect_hel_df), "insect taxa for direct association with collection date",
    "(linear + quadratic,", n_perm, "permutations each)...\n")
insect_date_results <- bind_rows(lapply(colnames(insect_hel_df), function(tax) {
  r <- test_date_taxon(insect_hel_df[[tax]], meta, n_perm)
  data.frame(taxon = tax, t_stat = r["t_stat"], p_perm = r["p_perm"],
             t_stat_quad = r["t_stat_quad"], p_perm_quad = r["p_perm_quad"])
}))
insect_date_results$q_value <- p.adjust(insect_date_results$p_perm, method = "BH")
insect_date_results$q_value_quad <- p.adjust(insect_date_results$p_perm_quad, method = "BH")
insect_date_results <- classify_shape(insect_date_results)

cat(sum(insect_date_results$q_value < 0.10), "of", nrow(insect_date_results),
    "insect taxa significantly associated with date (q<0.10, LINEAR term)\n")
cat(sum(insect_date_results$q_value_quad < 0.10), "of", nrow(insect_date_results),
    "show a significant QUADRATIC date term (q<0.10) -- real hump/dip shape.\n")
cat("\nShape breakdown:\n")
print(table(insect_date_results$shape))
print(head(insect_date_results, 15))

## ---- 3. Add taxonomy for context (Family/Subfamily -- all families now, not ----
## ---- just Curculionidae/Latridiidae) ---------------------------------------

taxon_family <- sp_tab %>% distinct(Finest.ID, Family, Subfamily, Order)
insect_date_results <- insect_date_results %>% left_join(taxon_family, by = c("taxon" = "Finest.ID"))

write.csv(insect_date_results, file.path(out_data_dir, "insect_taxa_date_association.csv"), row.names = FALSE)

cat("\nHump/dip taxa (quadratic q<0.10), with family context:\n")
print(insect_date_results %>% filter(q_value_quad < 0.10) %>%
        select(taxon, Family, Subfamily, t_stat_quad, q_value_quad, shape), row.names = FALSE)

cat("\nDone. Output written to", file.path(out_data_dir, "insect_taxa_date_association.csv"), "\n")
