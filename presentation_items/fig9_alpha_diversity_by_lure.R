##############################################################################
# Presentation figure 9: insect + fungal alpha diversity vs. collection date,
# by lure (Lineage C -- full insect table, singleton-filter only, NOT
# prevalence-filtered; see project_organization.md's "Lineage C" section and
# "When to prevalence-filter ... vs. not" for why alpha-diversity work uses
# this construction rather than Lineage A or B's insect table).
#
# Standalone script, independent of the Lineage A/B presentation_items split
# (fig1-fig8 live in presentation_items/ for Lineage A, presentation_items_
# prevalence_filtered_insect/ for Lineage B) -- alpha diversity only has a
# Lineage A/C split, no Lineage B variant, so this lives directly in
# presentation_items/ alongside fig1-fig8 rather than a dedicated dir.
#
# Reads the insect + fungal alpha-diversity tables already computed by
# insect_fungal_alpha_diversity.full_insect_table.r (data/compare_insects_
# fungi_alpha_diversity_full_insect_table/{insect,fungal}_alpha_diversity.
# csv) rather than recomputing them here -- same pattern fig4/fig5/fig7/fig8
# use for reading an upstream screen's output.
#
# Panel a: insect (full insect table) alpha diversity. Panel b: fungal
# (rarefied) alpha diversity. Each is its own facet_grid block, placed side
# by side -- NOT combined into one grid, because insect and fungal values
# live on very different scales -- e.g. richness ~4-70 taxa vs. ~270-950
# ASVs -- and facet_grid free_y only frees scales across ROWS, not columns,
# so a single grid crossing community x lure would force insect and fungal
# onto a shared per-metric y-axis. Each block: rows = alpha-diversity metric
# (Shannon, Simpson dominance, richness), columns = lure type, date on the
# x-axis, one point per sample colored by site, with a per-panel cubic
# (y ~ poly(date, 3)) OLS trend line -- matching the linear+quadratic+cubic
# date model insect_fungal_alpha_diversity.full_insect_table.r now tests via
# permutation (see alpha_diversity_covariate_tests.csv and alpha_diversity_
# date_shape.csv in the full_insect_table output dir; 4 of the 6 community x
# metric date tests have a significant cubic term, so this isn't just
# visually motivated). The trend line itself is still purely descriptive
# (the permutation test is what carries the actual significance call).
#
# The two blocks are combined with gtable::cbind(size="max") (same technique
# fig1_taxonomic_breakdown.R uses for its a/c vs. b/d columns) so panel/row
# heights -- which would otherwise drift apart given the two blocks'
# different legend widths -- line up between insect and fungal. The row strip
# (metric name) is shown once, on panel a (leftmost), switched via
# facet_grid(switch="y") to sit outside the y-axis text rather than between
# the axis and the panel; panel b's row strip is blanked.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)
source("library/library.R")

out_fig_dir <- "figures/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

alpha_data_dir <- "data/compare_insects_fungi_alpha_diversity_full_insect_table"

## ---- Load insect + fungal alpha diversity (Lineage C, full insect table) ----

site_levels <- c("Durham", "Pease Airport", "Manchester Cedar Swamp", "Manchester Airport")
lure_levels <- c("Ethanol", "Alpha-pinene_EtOH", "Ips")
metrics <- c("shannon", "simpson_dominance", "richness")
metric_labels <- c(shannon = "Shannon diversity", simpson_dominance = "Simpson dominance",
                    richness = "Richness (observed taxa)")

load_alpha_long <- function(file) {
  read.csv(file.path(alpha_data_dir, file)) %>%
    mutate(date = as.Date(date),
           site = ifelse(site == "Pease", "Pease Airport", site),
           site = factor(site, levels = site_levels),
           lure = factor(lure, levels = lure_levels)) %>%
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric_labels[metric], levels = metric_labels))
}

insect_long <- load_alpha_long("insect_alpha_diversity.csv")
fungal_long <- load_alpha_long("fungal_alpha_diversity.csv")

cat("Insect alpha diversity (Lineage C, full table):", n_distinct(insect_long$sample_id), "samples\n")
cat("Fungal alpha diversity (rarefied):", n_distinct(fungal_long$sample_id), "samples\n")

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

# cubic OLS trend line per panel, descriptive only. show_row_strip's block
# gets the metric row labels, switched to the left (outside the axis text)
# via facet_grid(switch = "y") -- the other block's row strip is blanked so
# the label isn't shown twice.
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

ggsave(file.path(out_fig_dir, "fig9_alpha_diversity_by_lure.png"), combined,
       width = 18, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig9_alpha_diversity_by_lure.pdf"), combined,
       width = 18, height = 8)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig9_alpha_diversity_by_lure.png"), "and .pdf\n")
