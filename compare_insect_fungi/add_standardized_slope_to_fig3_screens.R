##############################################################################
# Retrofit: add the standardized-slope ("beta weight") effect-size column to
# the three fig3 input CSVs produced by scripts we deliberately do NOT rerun.
#
# WHY THIS EXISTS
# fig3 (both lineages) switched its volcano x-axis from the raw permutation
# t-statistic to a fully-standardized partial slope
#   beta_std = b * sd(x) / sd(y)
# (the QuantPsyc::lm.beta / effectsize::standardize_parameters(method =
# "basic") convention). That quantity is a property of the OBSERVED OLS fit
# only -- it involves no permutation -- so it can be added exactly without
# rerunning any permutation screen:
#   * insect_fungal_coca_analysis.general_workflow.r and
#     fungal_community_seasonality.r are seeded and (for the former) single-
#     threaded, but general_workflow.r is long and fungal_community_
#     seasonality.r uses mclapply with mc.set.seed = TRUE, so a full rerun
#     of the latter would also jitter the permutation q-values (Monte-Carlo
#     error at n_perm = 999) that fig4-fig8 / figS1 / responsive_taxa and
#     two interpretation docs all depend on. Not wanted for an x-axis change.
# So: the two source scripts were edited to emit beta_std / beta_std_quad
# natively for any FUTURE full rerun, and this script computes the identical
# column now from the same data + same model and merges it into the existing
# CSVs, touching nothing else.
#
# The two prevalence-filtered Lineage-B scripts behind fig3 panels a/d and
# c/f (insect_taxa_date_association.prevalence_filtered_insect.r,
# insect_fungal_coca_factordate_check.prevalence_filtered_insect.r) ARE
# seeded, single-threaded and fast, so those were edited and simply rerun --
# they are not handled here.
#
# CSVs updated by this script (new columns inserted, all existing columns and
# their values preserved byte-for-byte):
#   1. data/2024_insect_data/insect_taxa_date_association.csv
#        + beta_std, beta_std_quad   (insect Hellinger ~ site+lure+date_c+I(date_c^2))
#   2. data/2024_fungi/fungal_taxa_date_association.all_taxa.csv
#        + beta_std, beta_std_quad   (fungal CLR ~ site+lure+date_c+I(date_c^2))
#   3. data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv
#        + beta_std                  (fungal CLR ~ site+insect_PCoA1+insect_PCoA2+insect_PCoA3)
##############################################################################

library(vegan)
library(dplyr)

## ---- 1. Rebuild the Lineage-A matched dataset (identical to general_workflow.r ----
## ---- / fungal_community_seasonality.r sections 1-2) ----------------------------

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

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}
fungal_clr <- clr_transform(fungal_f)
insect_hel <- decostand(insect, method = "hellinger")

# Insect PCoA axes -- identical construction to general_workflow.r section 4
insect_dist <- vegdist(insect_hel, method = "euclidean")
insect_pcoa <- cmdscale(insect_dist, k = 3, eig = TRUE)
insect_axes <- as.data.frame(insect_pcoa$points)
colnames(insect_axes) <- c("insect_PCoA1", "insect_PCoA2", "insect_PCoA3")
insect_axes$sample_id <- rownames(insect)

## ---- 2. Effect-size helpers (OBSERVED fit only, no permutation) --------------
# beta_std = b * sd(x) / sd(y), fully-standardized partial slope.

std_slope_date <- function(y, meta) {
  d <- meta
  d$y <- y
  d$date_c <- as.numeric(d$date) - mean(as.numeric(d$date))
  cs <- coef(summary(lm(y ~ site + lure + date_c + I(date_c^2), data = d)))
  sd_y <- sd(d$y)
  c(beta_std      = unname(cs["date_c", "Estimate"])      * sd(d$date_c)   / sd_y,
    beta_std_quad = unname(cs["I(date_c^2)", "Estimate"]) * sd(d$date_c^2) / sd_y)
}

std_slope_axis <- function(y, dat_base, test_axis, rhs) {
  d <- dat_base
  d$y <- y
  cs <- coef(summary(lm(as.formula(paste("y ~", rhs)), data = d)))
  unname(cs[test_axis, "Estimate"]) * sd(d[[test_axis]]) / sd(d$y)
}

insert_after <- function(df, new_cols, after) {
  # place new_cols (a data.frame/named list, same nrow) immediately after
  # column `after`, keeping every other column in place
  stopifnot(after %in% names(df))
  keep_before <- names(df)[seq_len(match(after, names(df)))]
  keep_after  <- setdiff(names(df), keep_before)
  bind_cols(df[keep_before], as.data.frame(new_cols), df[keep_after])
}

## ---- 3. CSV 1: insect ~ date (Lineage A, general_workflow.r) -----------------

f1 <- "data/2024_insect_data/insect_taxa_date_association.csv"
csv1 <- read.csv(f1, check.names = FALSE)
stopifnot(all(csv1$taxon %in% colnames(insect_hel)))
b1 <- t(vapply(csv1$taxon, function(tx) std_slope_date(insect_hel[, tx], meta),
               numeric(2)))
csv1_out <- insert_after(csv1, list(beta_std = b1[, "beta_std"],
                                    beta_std_quad = b1[, "beta_std_quad"]),
                         after = "p_perm_quad")
write.csv(csv1_out, f1, row.names = FALSE)
cat("Updated", f1, "-- added beta_std, beta_std_quad for", nrow(csv1_out), "insect taxa.\n")

## ---- 4. CSV 2: fungal ~ date (lineage-invariant, fungal_community_seasonality.r) ----

f2 <- "data/2024_fungi/fungal_taxa_date_association.all_taxa.csv"
csv2 <- read.csv(f2, check.names = FALSE)
stopifnot(all(csv2$taxon %in% colnames(fungal_clr)))
b2 <- t(vapply(csv2$taxon, function(tx) std_slope_date(fungal_clr[, tx], meta),
               numeric(2)))
csv2_out <- insert_after(csv2, list(beta_std = b2[, "beta_std"],
                                    beta_std_quad = b2[, "beta_std_quad"]),
                         after = "p_perm_quad")
write.csv(csv2_out, f2, row.names = FALSE)
cat("Updated", f2, "-- added beta_std, beta_std_quad for", nrow(csv2_out), "fungal taxa.\n")

## ---- 5. CSV 3: fungal ~ insect_PCoA1, no lure/no date (general_workflow.r) ---

f3 <- "data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv"
csv3 <- read.csv(f3, check.names = FALSE)
dat_base <- meta %>% left_join(insect_axes, by = "sample_id")
stopifnot(all(dat_base$sample_id == rownames(fungal_clr)))
stopifnot(all(csv3$taxon %in% colnames(fungal_clr)))
rhs3 <- "site + insect_PCoA1 + insect_PCoA2 + insect_PCoA3"
b3 <- vapply(csv3$taxon,
             function(tx) std_slope_axis(fungal_clr[, tx], dat_base, "insect_PCoA1", rhs3),
             numeric(1))
csv3_out <- insert_after(csv3, list(beta_std = unname(b3)), after = "p_perm")
write.csv(csv3_out, f3, row.names = FALSE)
cat("Updated", f3, "-- added beta_std for", nrow(csv3_out), "fungal taxa.\n")

## ---- 6. Sanity check: sign(beta_std) must match sign(t_stat) ----------------
for (nm in list(c(f1, "csv1_out"), c(f2, "csv2_out"), c(f3, "csv3_out"))) {
  d <- get(nm[2])
  bad <- sum(sign(d$beta_std) != sign(d$t_stat) & d$t_stat != 0)
  cat(nm[1], ": sign(beta_std) != sign(t_stat) for", bad, "of", nrow(d), "rows",
      if (bad == 0) "(OK)\n" else "(INVESTIGATE)\n")
  if ("t_stat_quad" %in% names(d)) {
    badq <- sum(sign(d$beta_std_quad) != sign(d$t_stat_quad) & d$t_stat_quad != 0)
    cat(nm[1], ": sign(beta_std_quad) != sign(t_stat_quad) for", badq, "of", nrow(d), "rows",
        if (badq == 0) "(OK)\n" else "(INVESTIGATE)\n")
  }
}

cat("\nDone.\n")
