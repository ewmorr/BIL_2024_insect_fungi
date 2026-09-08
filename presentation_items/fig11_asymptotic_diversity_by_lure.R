##############################################################################
# Presentation figure 11: insect + fungal ASYMPTOTIC (Chao-extrapolated)
# diversity vs. collection date, by lure. Lineage C, full insect table --
# same construction as fig9_alpha_diversity_by_lure.R, which this figure
# reproduces one-for-one except for the metric: observed diversity swapped
# for the asymptotic estimates from insect_fungal_asymptotic_richness_iNEXT.
# full_insect_table.r (data/compare_insects_fungi_asymptotic_richness_iNEXT_
# full_insect_table/{insect,fungal}_asymptotic_diversity.csv) -- see that
# script's header and project_organization.md's Lineage C section for why
# the fungal side is left unrarefied here (unlike fig9's rarefied-and-
# averaged fungal table) and why "simpson" here is Hill number q=2 diversity
# (higher = more even), the OPPOSITE direction from fig9's Simpson
# DOMINANCE (higher = more dominated).
#
# All design choices identical to fig9 (see that script's header for the
# full rationale): panel a = insect, panel b = fungal, each its own
# facet_grid block (not a single grid) because insect and fungal live on
# very different scales; rows = metric, columns = lure, cubic OLS trend line
# per panel; blocks combined with gtable::cbind(size="max") so panel heights
# line up despite fungal's wider legend.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)
source("library/library.R")

out_fig_dir <- "figures/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

asymp_data_dir <- "data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table"

## ---- Load insect + fungal asymptotic diversity (Lineage C, full insect table, iNEXT) ----

site_levels <- c("Durham", "Pease Airport", "Manchester Cedar Swamp", "Manchester Airport")
lure_levels <- c("Ethanol", "Alpha-pinene_EtOH", "Ips")
metrics <- c("shannon_asymptotic", "simpson_asymptotic", "richness_asymptotic")
# Short labels -- these sit in a switch="y" row strip (rotated text along the
# left edge, see make_block() below), which clips text longer than fig9's
# original labels ("Shannon diversity" etc.) allowed for. Full "Hill q=1/q=2"
# wording is kept in fig12's (horizontal, unconstrained) strip labels instead.
metric_labels <- c(shannon_asymptotic = "Shannon div. (q=1)",
                    simpson_asymptotic = "Simpson div. (q=2)",
                    richness_asymptotic = "Asymptotic richness")

load_asymp_long <- function(file) {
  read.csv(file.path(asymp_data_dir, file)) %>%
    mutate(date = as.Date(date),
           site = ifelse(site == "Pease", "Pease Airport", site),
           site = factor(site, levels = site_levels),
           lure = factor(lure, levels = lure_levels)) %>%
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric_labels[metric], levels = metric_labels))
}

insect_long <- load_asymp_long("insect_asymptotic_diversity.csv")
fungal_long <- load_asymp_long("fungal_asymptotic_diversity.csv")

cat("Insect asymptotic diversity (Lineage C, full table):", n_distinct(insect_long$sample_id), "samples\n")
cat("Fungal asymptotic diversity (unrarefied):", n_distinct(fungal_long$sample_id), "samples\n")

## ---- Shared look ------------------------------------------------------------

panel_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 9, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(size = 11),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.background.y = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    strip.placement = "outside",
    plot.tag = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8.5)
  )

site_palette <- setNames(twelvePaired[c(2, 4, 6, 8)], site_levels)

make_block <- function(df, tag, show_legend, show_row_strip) {
  p <- ggplot(df, aes(x = date, y = value)) +
    geom_point(aes(fill = site), shape = 21, size = 2.2, color = "black", stroke = 0.3) +
    geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = FALSE, color = "grey30", linewidth = 0.6) +
    facet_grid(rows = vars(metric), cols = vars(lure), scales = "free_y",
               switch = if (show_row_strip) "y" else NULL) +
    scale_fill_manual(values = site_palette, name = "Site") +
    scale_x_date(date_labels = "%b %d") +
    labs(x = "Collection date", y = NULL, tag = tag) +
    panel_theme
  if (!show_legend) p <- p + guides(fill = "none")
  if (!show_row_strip) p <- p + theme(strip.text.y = element_blank())
  p
}

## ---- Figure: insect block (left, tag a, row labels) + fungal block (right, tag b) ----

insect_block <- make_block(insect_long, "a", show_legend = FALSE, show_row_strip = TRUE)
fungal_block <- make_block(fungal_long, "b", show_legend = TRUE, show_row_strip = FALSE)

combined <- cbind(ggplotGrob(insect_block), ggplotGrob(fungal_block), size = "max")

ggsave(file.path(out_fig_dir, "fig11_asymptotic_diversity_by_lure.png"), combined,
       width = 18, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig11_asymptotic_diversity_by_lure.pdf"), combined,
       width = 18, height = 8)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig11_asymptotic_diversity_by_lure.png"), "and .pdf\n")
