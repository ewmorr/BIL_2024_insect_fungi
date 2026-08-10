##############################################################################
# Presentation figure 2: community ordinations + whole-community concordance
#
# Panel a: insect NMDS (Curculionidae + Latridiidae trap catch, raw
#          Bray-Curtis, no autotransform -- insect_exploratory/insect_ords.R's
#          setting for this table), fill = collection date, shape = lure,
#          annotated with PERMANOVA R2/p for site, lure, date.
# Panel b: fungal NMDS (ITS2 ASVs, rarefied-to-5000 x 100 iterations,
#          Bray-Curtis averaged across iterations -- compare_insect_fungi/
#          fungal_community_seasonality.r's approach), same aesthetics,
#          annotated with the pre-computed PERMANOVA from
#          data/2024_fungi/fungal_community_permanova.csv (same matched
#          69-sample dataset and model).
# Panel c: Procrustes overlay of the two ordinations (compare_insect_fungi/
#          insect_fungal_procrustes.r), fill = community (insect/fungal),
#          shape = site, annotated with the Procrustes correlation.
#
# All three panels are built from ONE insect NMDS + ONE fungal NMDS run (the
# fungal NMDS/rarefaction is the slow part), so panel b's ordination and
# panel c's fungal configuration are guaranteed to be the identical run, not
# two independently-stochastic re-fits.
#
# The fungal rarefaction + NMDS + Procrustes fit is cached to
# data/presentation_items/fig2_ordination_cache.rds after the first run --
# delete that file (or set force_recompute <- TRUE below) to redo it from
# scratch, e.g. after the underlying data changes.
##############################################################################

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(scales)
source("library/library.R")

set.seed(1)

out_fig_dir <- "figures/presentation_items"
out_data_dir <- "data/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

cache_file <- file.path(out_data_dir, "fig2_ordination_cache.rds")
force_recompute <- T

## ---- 1. Load and match insect + fungal data (same as compare_insect_fungi/*) ----

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

stopifnot(all(rownames(insect) == meta$sample_id))
stopifnot(all(rownames(fungal) == meta$sample_id))
cat("Matched dataset:", nrow(insect), "samples,", ncol(insect), "insect taxa,", ncol(fungal), "fungal ASVs.\n")

## ---- 2. Ordinations + Procrustes (cached) -----------------------------------

if (!force_recompute && file.exists(cache_file)) {
  cat("Loading cached ordinations from", cache_file, "\n")
  cached <- readRDS(cache_file)
  insect_nmds <- cached$insect_nmds
  fungal_nmds <- cached$fungal_nmds
  proc <- cached$proc
  proc_test <- cached$proc_test
} else {
  cat("No cache found (or force_recompute=TRUE) -- computing ordinations...\n")

  insect_nmds <- metaMDS(insect, distance = "bray", binary = FALSE, try = 20, trymax = 100, autotransform = FALSE)
  cat("Insect NMDS stress:", round(insect_nmds$stress, 3), "\n")

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

  # Procrustes requires matching insect/fungal sample sets -- restrict both
  # to the fungal rarefaction survivors (panels a & b keep their own full
  # sample sets; only the Procrustes fit/panel c is affected by a drop).
  proc_sample_ids <- rownames(fungal_rare_dist_avg)
  insect_scores <- scores(insect_nmds, display = "sites")[proc_sample_ids, ]
  fungal_scores <- scores(fungal_nmds, display = "sites")[proc_sample_ids, ]
  proc <- procrustes(X = insect_scores, Y = fungal_scores, symmetric = TRUE)
  proc_test <- protest(X = insect_scores, Y = fungal_scores, permutations = 999, symmetric = TRUE)
  cat("Procrustes correlation:", round(proc_test$t0, 3), ", p =", proc_test$signif, "\n")

  saveRDS(list(insect_nmds = insect_nmds, fungal_nmds = fungal_nmds, proc = proc, proc_test = proc_test), cache_file)
  cat("Cached ordinations to", cache_file, "\n")
}

## ---- 3. PERMANOVA (site / lure / date) --------------------------------------
# Project convention throughout compare_insect_fungi/ and insect_exploratory/:
# Fit marginal effects for accurate r2 values
# note the p-values were hardcoded into the file as these were run in a separate
# script using permute::how() to derive accurate p-values

#insect_permanova <- adonis2(insect ~ site + lure + date, data = meta, by = "margin")
#print(insect_permanova)
#write.csv(as.data.frame(insect_permanova), file.path(out_data_dir, "insect_community_permanova.csv"))
insect_permanova <- read.csv(file.path(out_data_dir, "insect_community_permanova.csv"), row.names = 1, check.names = FALSE)

# fungal PERMANOVA: use the pre-computed result (compare_insect_fungi/
# fungal_community_seasonality.r), same matched dataset + model, rather than
# re-running the rarefied PERMANOVA here.
fungal_permanova <- read.csv("data/2024_fungi/fungal_community_permanova.csv", row.names = 1, check.names = FALSE)

extract_term <- function(tab, term) list(r2 = tab[term, "R2"], p = tab[term, "Pr(>F)"])

fmt_p <- function(p) ifelse(p < 0.001, "p<0.001", paste0("p=", sprintf("%.3f", p)))
fmt_stats <- function(tab, site_term = "site", lure_term = "lure", date_term = "date") {
  s <- extract_term(tab, site_term); l <- extract_term(tab, lure_term); d <- extract_term(tab, date_term)
  paste0(
    "Site  R²=", sprintf("%.2f", s$r2), ", ", fmt_p(s$p), "\n",
    "Lure  R²=", sprintf("%.2f", l$r2), ", ", fmt_p(l$p), "\n",
    "Date  R²=", sprintf("%.2f", d$r2), ", ", fmt_p(d$p)
  )
}

insect_stats_label <- fmt_stats(as.data.frame(insect_permanova))
fungal_stats_label <- fmt_stats(fungal_permanova)

## ---- 4. Panels a & b: NMDS plots ---------------------------------------------

nmds_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(0.4, "cm"),
    plot.tag = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 8, family = "mono", color = "grey20", lineheight = 1.15,
                                  margin = margin(t = 2, b = 4)),
    panel.grid = element_blank(),
    aspect.ratio = 1   # square plotted panel
  )

date_range <- range(as.numeric(meta$date))
date_mid <- median(as.numeric(meta$date))
date_fill_scale <- scale_fill_gradient2(
  low = "#b2182b", mid = "white", high = "#2166ac", midpoint = date_mid,
  breaks = pretty(date_range, n = 4),
  labels = function(b) format(as.Date(b, origin = "1970-01-01"), "%b %d"),
  name = "Collection\ndate"
)
lure_shape_scale <- scale_shape_manual(values = c(21, 22, 24), name = "Lure")

make_nmds_panel <- function(nmds_obj, meta_df, stats_label, title, tag, show_legend = TRUE) {
  plot_df <- as.data.frame(scores(nmds_obj, display = "sites")) %>%
    tibble::rownames_to_column("sample_id") %>%
    left_join(meta_df, by = "sample_id")

  label_text <- paste0("stress = ", round(nmds_obj$stress, 2), "\n", stats_label)

  p <- ggplot(plot_df, aes(x = NMDS1, y = NMDS2, fill = as.numeric(date), shape = lure)) +
    geom_point(size = 2.8, color = "black", stroke = 0.3) +
    date_fill_scale +
    lure_shape_scale +
    guides(shape = guide_legend(override.aes = list(fill = "grey40"), order = 1),
           fill = guide_colorbar(order = 2)) +
    labs(title = title, subtitle = label_text, tag = tag) +
    nmds_theme

  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

# Lure + Date is the same legend for both a and b (same meta, same scales) --
# shown once, to the right of panel b, rather than duplicated on both.
panel_a <- make_nmds_panel(insect_nmds, meta, insect_stats_label, "Insect community", "a", show_legend = FALSE)
panel_b <- make_nmds_panel(fungal_nmds, meta, fungal_stats_label, "Fungal community", "b")

## ---- 5. Panel c: Procrustes overlay ------------------------------------------

proc_plot_df <- data.frame(
  sample_id = rownames(proc$Yrot),
  fungal_x = proc$Yrot[, 1], fungal_y = proc$Yrot[, 2],
  insect_x = proc$X[, 1], insect_y = proc$X[, 2]
) %>% left_join(meta, by = "sample_id")

proc_points <- bind_rows(
  proc_plot_df %>% transmute(sample_id, x = insect_x, y = insect_y, site, organism = "Insect"),
  proc_plot_df %>% transmute(sample_id, x = fungal_x, y = fungal_y, site, organism = "Fungi")
) %>% mutate(organism = factor(organism, levels = c("Insect", "Fungi")))

proc_label <- paste0("Procrustes r = ", round(proc_test$t0, 2), "\np = ", proc_test$signif)

panel_c <- ggplot() +
  geom_segment(data = proc_plot_df, aes(x = insect_x, y = insect_y, xend = fungal_x, yend = fungal_y),
               arrow = arrow(length = unit(0.08, "cm")), color = "grey55", alpha = 0.6) +
  geom_point(data = proc_points, aes(x = x, y = y, fill = organism, shape = site),
             size = 2.8, color = "black", stroke = 0.3) +
  scale_fill_manual(values = c(Insect = "#2166ac", Fungi = "#b2182b"), name = "Community") +
  scale_shape_manual(values = c(21, 22, 23, 24), name = "Site") +
  guides(fill = guide_legend(override.aes = list(shape = 21), order = 1),
         shape = guide_legend(order = 2)) +
  labs(x = "Dimension 1", y = "Dimension 2", title = "Procrustes: insect vs. fungal",
       subtitle = proc_label, tag = "c") +
  nmds_theme

## ---- 6. Combine + save --------------------------------------------------------
# grid.arrange gives the three panels equal overall column width, but each
# one's internal split between subtitle/legend area and the actual plotted
# NMDS panel differs (panel c's subtitle is one line shorter, and its
# Community+Site legend is a different height than a/b's Lure+Date legend),
# so the bordered panel boxes don't line up top-to-bottom. gtable::cbind
# (size="max") fixes this the same way fig1 aligns panel widths: it matches
# row heights (title/subtitle/panel/axis) across the three grobs directly.

row_grob <- cbind(ggplotGrob(panel_a), ggplotGrob(panel_b), ggplotGrob(panel_c), size = "max")

ggsave(file.path(out_fig_dir, "fig2_nmds_procrustes.png"), row_grob, width = 15, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig2_nmds_procrustes.pdf"), row_grob, width = 15, height = 5.5)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig2_nmds_procrustes.png"), "and .pdf\n")
