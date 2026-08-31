##############################################################################
# Insect-Fungal Community Procrustes Analysis -- Prevalence-Filtered Insect
# Table Variant (Lineage B)
#
# Identical whole-community Procrustes concordance test to
# insect_fungal_procrustes.r, except the insect community table is built with
# an individual-taxon >=5-sample PREVALENCE filter (all families -- see
# insect_exploratory/insect_ords.prevalence_filter.R) instead of restricting
# to Curculionidae + Latridiidae. Procrustes asks whether the two FULL
# community ordinations have congruent shape -- i.e. do samples close together
# in insect-space also tend to be close together in fungal-space.
#
# Two deliberate differences from the Lineage-A script beyond the insect
# table itself:
#   1. Adds the `Kingdom == "k__Fungi"` filter on the fungal ASV table. The
#      Lineage-A standalone procrustes script predates that fix (it is not one
#      of the "5 core scripts" the Kingdom filter was retrofitted into) and
#      still runs on the unfiltered ASV table; presentation_items/
#      fig2_nmds_procrustes.R already applies it and it is now the project
#      standard, so this port uses it too.
#   2. Handles the rarefaction sample drop gracefully. Once non-fungal reads
#      are stripped, one sample's total read count can fall below the common
#      rarefaction depth (the 69 -> 68 drop documented in
#      iterative_analysis_updates.md). The Procrustes fit is restricted to the
#      rarefaction survivors rather than asserting no sample is lost.
#
# Also computes and writes the matched-dataset insect community PERMANOVA
# (site + lure + date, marginal) so the Lineage-B fig2 port can annotate its
# insect panel without re-running an ordination -- the same role
# fungal_community_seasonality.r plays as the source of the (lineage-
# invariant) fungal PERMANOVA.
#
# Insect ordination: NMDS on raw Bray-Curtis, no autotransform -- the setting
# insect_ords.prevalence_filter.R confirmed gives a low, repeatable stress
# solution (0.16) for this prevalence-filtered table.
# Fungal ordination: NMDS on rarefaction-averaged Bray-Curtis, mirroring
# fungal_community_seasonality.r Step 3.
#
# Output goes to its own directory so the Lineage-A results are untouched.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(vegan)
library(permute)
library(dplyr)
library(ggplot2)
source("library/library.R")

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_procrustes_prevalence_filtered_insect"
out_fig_dir  <- "figures/compare_insects_fungi_procrustes_prevalence_filtered_insect"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load full insect table, apply individual-taxon prevalence filter ----
# All families, taxa present in >=5 samples -- the same threshold used for
# fungal ASVs throughout compare_insect_fungi/, in place of the
# Curculionidae+Latridiidae family restriction used in the Lineage-A script.

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

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
cat(length(shared_ids), "of", nrow(insect_full), "insect samples have matching fungal ASV data.\n")
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

# Restrict to ASVs confirmed as Fungi -- ASV_tab.csv includes non-fungal
# (plant/animal/protist) and taxonomically unidentified ASVs; see
# ASVs_taxonomy.tsv (Kingdom column). (Difference #1 from the Lineage-A script.)
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
cat("Matched dataset:", nrow(insect), "samples,", ncol(insect),
    "insect taxa (>=5-sample prevalence filter, all families),", ncol(fungal), "fungal ASVs.\n")

## ---- 2. Insect NMDS (raw Bray-Curtis, no autotransform) -------------------
# autotransform=FALSE is the setting insect_ords.prevalence_filter.R confirmed
# for this prevalence-filtered table (lower, more repeatable stress than
# vegan's default sqrt + Wisconsin autotransform).

insect_nmds <- metaMDS(insect, distance = "bray", binary = FALSE, try = 20, trymax = 100, autotransform = FALSE)
cat("Insect NMDS stress:", round(insect_nmds$stress, 3), "\n")

## ---- 3. Fungal NMDS (rarefaction-averaged Bray-Curtis) --------------------
# Mirrors fungal_community_seasonality.r Step 3.

rarefy_depth <- 5000
rarefy_iterations <- 100

fungal_rare_list <- multiple_subsamples(x = fungal, depth = rarefy_depth, iterations = rarefy_iterations)
n_dropped <- nrow(fungal) - nrow(fungal_rare_list[[1]])
if (n_dropped > 0) {
  dropped_ids <- setdiff(rownames(fungal), rownames(fungal_rare_list[[1]]))
  cat(n_dropped, "sample(s) dropped from the fungal rarefaction (total reads fell below depth",
      rarefy_depth, "after the taxonomy filter):", paste(dropped_ids, collapse = ", "), "\n")
}

fungal_rare_dist_list <- lapply(fungal_rare_list, vegdist, method = "bray")
fungal_rare_dist_avg <- avg_matrix_list(fungal_rare_dist_list)

fungal_nmds <- metaMDS(as.dist(fungal_rare_dist_avg), try = 20, trymax = 100)
cat("Fungal NMDS stress:", round(fungal_nmds$stress, 3), "\n")

## ---- 4. Insect community PERMANOVA (site / lure / date) ------------------
# Written out for the Lineage-B fig2 port to annotate its insect panel.
# Same single permutation scheme as insect_ords.prevalence_filter.R: permute
# freely within trap, blocked by site. As that script notes, this is a sound
# date/site test but only an approximate lure test (lure is trap-constant);
# see insect_ords.R's balanced re-test if a rigorous lure p-value is needed.

perm_insect <- how(
  blocks = meta$site,
  plots  = Plots(strata = meta$trap_id, type = "none"),
  within = Within(type = "free"),
  nperm  = 999
)
insect_permanova <- adonis2(insect ~ site + lure + date, data = meta,
                            permutations = perm_insect, by = "margin")
cat("\n--- PERMANOVA: prevalence-filtered insect community ~ site + lure + date (matched dataset) ---\n")
print(insect_permanova)
write.csv(as.data.frame(insect_permanova),
          file.path(out_data_dir, "insect_community_permanova.csv"))

## ---- 5. Procrustes analysis ---------------------------------------------
# Restrict both configurations to the fungal rarefaction survivors and align
# by sample_id (do not assume the two NMDS runs preserved input row order).
# (Difference #2 from the Lineage-A script: survivors, not stopifnot-equal.)

proc_sample_ids <- rownames(fungal_rare_dist_avg)
insect_scores <- scores(insect_nmds, display = "sites")
fungal_scores <- scores(fungal_nmds, display = "sites")
stopifnot(all(proc_sample_ids %in% rownames(insect_scores)))
stopifnot(all(proc_sample_ids %in% rownames(fungal_scores)))
insect_scores <- insect_scores[proc_sample_ids, ]
fungal_scores <- fungal_scores[proc_sample_ids, ]

# symmetric = TRUE: neither community is treated as the fixed target, so the
# fitted correlation isn't sensitive to an arbitrary choice of X vs Y
proc <- procrustes(X = insect_scores, Y = fungal_scores, symmetric = TRUE)
print(summary(proc))

proc_test <- protest(X = insect_scores, Y = fungal_scores, permutations = 999, symmetric = TRUE)
print(proc_test)

cat("\nProcrustes correlation (concordance between insect and fungal NMDS):",
    round(proc_test$t0, 3), ", permutation p =", proc_test$signif, "\n")

proc_summary <- data.frame(
  correlation = proc_test$t0,
  sum_of_squares_m12 = proc$ss,
  p_value = proc_test$signif,
  n_permutations = proc_test$permutations,
  n_samples = length(proc_sample_ids)
)
write.csv(proc_summary, file.path(out_data_dir, "insect_fungal_procrustes_summary.csv"), row.names = FALSE)

residuals_df <- data.frame(sample_id = names(residuals(proc)), residual = residuals(proc)) %>%
  left_join(meta, by = "sample_id") %>%
  arrange(desc(residual))
write.csv(residuals_df, file.path(out_data_dir, "insect_fungal_procrustes_residuals.csv"), row.names = FALSE)
cat("\nSamples with the largest insect-fungal procrustes residual (worst-matched):\n")
print(head(residuals_df, 10))

## ---- 6. Plots ---------------------------------------------------------------

pdf(file.path(out_fig_dir, "insect_fungal_procrustes.pdf"), width = 8, height = 8)
plot(proc, kind = 1, main = "Procrustes: insect vs. fungal community NMDS (prevalence-filtered insect table)",
     sub = paste0("correlation = ", round(proc_test$t0, 3), ", p = ", proc_test$signif))
plot(proc, kind = 2, main = "Procrustes residuals by sample")
dev.off()

# ggplot version: insect and fungal points for the same sample connected by an
# arrow after Procrustes rotation -- short arrows = concordant samples.
proc_plot_df <- data.frame(
  sample_id = rownames(proc$Yrot),
  fungal_x = proc$Yrot[, 1], fungal_y = proc$Yrot[, 2],
  insect_x = proc$X[, 1], insect_y = proc$X[, 2]
) %>% left_join(meta, by = "sample_id")

procrustes_overlay <- ggplot(proc_plot_df) +
  geom_segment(aes(x = insect_x, y = insect_y, xend = fungal_x, yend = fungal_y),
               arrow = arrow(length = unit(0.1, "cm")), color = "grey50", alpha = 0.6) +
  geom_point(aes(x = insect_x, y = insect_y, shape = site), size = 2.5, fill = "#2166ac", color = "black") +
  geom_point(aes(x = fungal_x, y = fungal_y, shape = site), size = 2.5, fill = "#b2182b", color = "black") +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  annotate(geom = "text", label = paste0("procrustes r = ", round(proc_test$t0, 2), ", p = ", proc_test$signif),
           x = -Inf, y = -Inf, hjust = -0.1, vjust = -1) +
  labs(x = "Dimension 1", y = "Dimension 2",
       title = "Insect (blue) vs. fungal (red) community NMDS after Procrustes rotation",
       subtitle = "insect table: individual-taxon >=5-sample prevalence filter, all families",
       shape = "Site") +
  theme_bw()

ggsave(file.path(out_fig_dir, "insect_fungal_procrustes.NMDS_overlay.png"), procrustes_overlay, width = 7, height = 6, dpi = 150)
print(procrustes_overlay)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
