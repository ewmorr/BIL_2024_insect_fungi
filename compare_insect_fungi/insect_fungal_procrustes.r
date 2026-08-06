##############################################################################
# Insect-Fungal Community Procrustes Analysis
#
# Direct test of whole-community concordance between the insect (Curculionidae
# + Latridiidae trap catch) and fungal (ITS2 ASV) community ordinations, on
# the same 69-sample matched dataset used throughout compare_insect_fungi/.
# Procrustes analysis asks a different question than the CoCA/per-taxon work
# in insect_fungal_coca_analysis.general_workflow.r and
# fungal_community_seasonality.r: not "do specific fungal taxa track specific
# insect axes/date", but "do the two FULL community ordinations, taken as a
# whole, have congruent shape" -- i.e., do samples that are close together in
# insect-space also tend to be close together in fungal-space.
#
# Insect ordination: NMDS on raw Bray-Curtis, no autotransform -- the
# approach insect_exploratory/insect_ords.R found gave a low, repeatable
# stress solution for this family-filtered table.
# Fungal ordination: NMDS on rarefaction-averaged Bray-Curtis, mirroring
# fungal_community_seasonality.r Step 3 -- sequencing depth varies ~150x
# across samples, so raw counts are repeatedly rarefied to a common depth
# and the resulting distance matrices averaged before ordination.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(vegan)
library(dplyr)
library(ggplot2)
source("library/library.R")

set.seed(1)

out_data_dir <- "data/2024_insect_data"
out_fig_dir <- "figures"

## ---- 1. Load and match insect + fungal data (same as compare_insect_fungi/*) ----

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab %>% filter(Family %in% c("Curculionidae", "Latridiidae")) -> sp_tab.curcus_latris
sp_tab.curcus_latris.t <- t(sp_tab.curcus_latris %>% select(where(is.numeric)))
colnames(sp_tab.curcus_latris.t) <- sp_tab.curcus_latris$Finest.ID

# drop singleton taxa (present as a single individual total) and any
# resulting all-zero sample rows, as in insect_exploratory/insect_ords.R
sp_tab.curcus_latris.t[, colSums(sp_tab.curcus_latris.t) > 1] -> insect_full
insect_full <- insect_full[rowSums(insect_full) > 0, ]

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
cat(length(shared_ids), "of", nrow(insect_full), "insect samples have matching fungal ASV data.\n")
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))
cat("Matched dataset:", nrow(insect), "samples,", ncol(insect), "insect taxa (Curculionidae + Latridiidae),",
    ncol(fungal), "fungal ASVs.\n")

## ---- 2. Insect NMDS (raw Bray-Curtis, no autotransform) -------------------
# Same settings insect_ords.R settled on for this family-filtered table --
# low stress with a solution that repeats, unlike the untransformed full
# insect table or any of the tried transforms.

insect_nmds <- metaMDS(insect, distance = "bray", binary = FALSE, try = 20, trymax = 100, autotransform = FALSE)
cat("Insect NMDS stress:", round(insect_nmds$stress, 3), "\n")

## ---- 3. Fungal NMDS (rarefaction-averaged Bray-Curtis) --------------------
# Mirrors fungal_community_seasonality.r Step 3.

rarefy_depth <- 5000       # all matched samples exceed this (min depth = 5,120)
rarefy_iterations <- 100

fungal_rare_list <- multiple_subsamples(x = fungal, depth = rarefy_depth, iterations = rarefy_iterations)
cat("Retained", nrow(fungal_rare_list[[1]]), "of", nrow(fungal), "samples after depth filtering.\n")
# procrustes requires the two configurations to have a 1:1 sample correspondence
stopifnot(nrow(fungal_rare_list[[1]]) == nrow(fungal))

fungal_rare_dist_list <- lapply(fungal_rare_list, vegdist, method = "bray")
fungal_rare_dist_avg <- avg_matrix_list(fungal_rare_dist_list)

fungal_nmds <- metaMDS(as.dist(fungal_rare_dist_avg), try = 20, trymax = 100)
cat("Fungal NMDS stress:", round(fungal_nmds$stress, 3), "\n")

## ---- 4. Procrustes analysis ------------------------------------------------
# Align sample order explicitly by sample_id rather than assuming the two
# NMDS runs preserved input row order identically.

insect_scores <- scores(insect_nmds, display = "sites")
fungal_scores <- scores(fungal_nmds, display = "sites")
stopifnot(all(meta$sample_id %in% rownames(insect_scores)))
stopifnot(all(meta$sample_id %in% rownames(fungal_scores)))
insect_scores <- insect_scores[meta$sample_id, ]
fungal_scores <- fungal_scores[meta$sample_id, ]

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
  n_permutations = proc_test$permutations
)
write.csv(proc_summary, file.path(out_data_dir, "insect_fungal_procrustes_summary.csv"), row.names = FALSE)

residuals_df <- data.frame(sample_id = names(residuals(proc)), residual = residuals(proc)) %>%
  left_join(meta, by = "sample_id") %>%
  arrange(desc(residual))
write.csv(residuals_df, file.path(out_data_dir, "insect_fungal_procrustes_residuals.csv"), row.names = FALSE)
cat("\nSamples with the largest insect-fungal procrustes residual (worst-matched):\n")
print(head(residuals_df, 10))

## ---- 5. Plots ---------------------------------------------------------------

pdf(file.path(out_fig_dir, "insect_fungal_procrustes.pdf"), width = 8, height = 8)
plot(proc, kind = 1, main = "Procrustes: insect vs. fungal community NMDS",
     sub = paste0("correlation = ", round(proc_test$t0, 3), ", p = ", proc_test$signif))
plot(proc, kind = 2, main = "Procrustes residuals by sample")
dev.off()

# ggplot version, annotated with site, for the polished figure: insect and
# fungal points for the same sample connected by an arrow after Procrustes
# rotation, so short arrows = concordant samples, long arrows = discordant
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
       shape = "Site") +
  theme_bw()

ggsave(file.path(out_fig_dir, "insect_fungal_procrustes.NMDS_overlay.png"), procrustes_overlay, width = 7, height = 6, dpi = 150)
print(procrustes_overlay)

cat("\nDone. Outputs written to", out_data_dir, "/insect_fungal_procrustes_summary.csv,",
    out_data_dir, "/insect_fungal_procrustes_residuals.csv, and\n",
    out_fig_dir, "/insect_fungal_procrustes.pdf (base procrustes + residual diagnostic plots),\n",
    out_fig_dir, "/insect_fungal_procrustes.NMDS_overlay.png (annotated ggplot overlay)\n")
