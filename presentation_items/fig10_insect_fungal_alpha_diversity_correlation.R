##############################################################################
# Presentation figure 10: per-sample insect vs. fungal alpha diversity
# correlation (Lineage C -- full insect table, singleton-filter only, NOT
# prevalence-filtered; see project_organization.md's "Lineage C" section).
# The alpha-diversity analog of the whole-community Procrustes/CoCA cross-
# community checks used elsewhere in this project.
#
# One panel per metric (Shannon, Simpson dominance, richness), NOT split by
# lure -- an earlier version of this figure faceted by lure (mirroring
# fig9_alpha_diversity_by_lure.R's metric x lure grid), but a follow-up test
# added to insect_fungal_alpha_diversity.full_insect_table.r (section 5b,
# fungal ~ site+date+insect*lure vs. fungal ~ site+date+insect+lure, 999-
# permutation F-test on the interaction term) found NO evidence lure
# modifies the insect~fungal relationship for any of the 3 metrics (p_perm
# 0.80-0.95, see alpha_diversity_insect_fungal_interaction_by_lure.csv) --
# so splitting the figure by lure wasn't visualizing a real effect, just
# cutting n roughly into thirds. Collapsed back to the pooled relationship.
#
# One point per sample plotted as insect value (x) vs. fungal value (y),
# colored by collection date (continuous gradient) and shaped by lure (the
# trap-type variable -- Ethanol, Alpha-pinene_EtOH, Ips -- constant within a
# trap; descriptive only, the trend line/Spearman annotation still pools
# across lure per the null interaction result above). facet_wrap(scales=
# "free") gives each metric panel its own independent x AND y range (Shannon
# ~0-5, Simpson dominance ~0-0.7, richness insect ~4-70 vs. fungal
# ~270-950) -- facet_wrap's free scales are already per-panel, unlike
# facet_grid's (which are tied to row/column position; see fig10's git
# history for that bug), so no ggh4x needed now that lure is no longer a
# second faceting dimension.
#
# Reads the insect + fungal alpha-diversity tables already computed by
# insect_fungal_alpha_diversity.full_insect_table.r (data/compare_insects_
# fungi_alpha_diversity_full_insect_table/{insect,fungal}_alpha_diversity.
# csv), inner-joined by sample_id, and the pooled per-metric Spearman
# correlation that same script's section 5 already computed and wrote to
# alpha_diversity_cross_community_correlation.csv (rather than recomputing
# it here) for the in-panel rho/p annotation.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
source("library/library.R")

out_fig_dir <- "figures/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

alpha_data_dir <- "data/compare_insects_fungi_alpha_diversity_full_insect_table"

## ---- Load + join insect + fungal alpha diversity (Lineage C, full insect table) ----

metrics <- c("shannon", "simpson_dominance", "richness")
metric_labels <- c(shannon = "Shannon diversity", simpson_dominance = "Simpson dominance",
                    richness = "Richness (observed taxa)")

lure_levels <- c("Ethanol", "Alpha-pinene_EtOH", "Ips")
lure_shapes <- c(Ethanol = 21, `Alpha-pinene_EtOH` = 24, Ips = 22)  # filled circle/triangle/square

insect_alpha <- read.csv(file.path(alpha_data_dir, "insect_alpha_diversity.csv")) %>%
  mutate(date = as.Date(date), lure = factor(lure, levels = lure_levels))
fungal_alpha <- read.csv(file.path(alpha_data_dir, "fungal_alpha_diversity.csv")) %>%
  mutate(date = as.Date(date))

combined <- insect_alpha %>%
  select(sample_id, date, lure, insect_shannon = shannon,
         insect_simpson_dominance = simpson_dominance, insect_richness = richness) %>%
  inner_join(
    fungal_alpha %>% select(sample_id, fungal_shannon = shannon,
                             fungal_simpson_dominance = simpson_dominance, fungal_richness = richness),
    by = "sample_id"
  )
cat(nrow(combined), "samples with both insect and fungal alpha diversity.\n")

combined_long <- bind_rows(lapply(metrics, function(m) {
  combined %>% transmute(sample_id, date, lure, metric = metric_labels[m],
                          insect_value = .data[[paste0("insect_", m)]],
                          fungal_value = .data[[paste0("fungal_", m)]])
})) %>% mutate(metric = factor(metric, levels = metric_labels))

## ---- Pooled per-metric Spearman correlation, for the in-panel label -------
# Read from the upstream script's output rather than recomputing -- same
# test, same number (n=68, pooled across lure).

cross_corr <- read.csv(file.path(alpha_data_dir, "alpha_diversity_cross_community_correlation.csv"))
corr_labels <- cross_corr %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         label = paste0("rho = ", round(rho, 2), "\np = ", signif(p_value, 2)))

cat("\n--- Pooled per-metric Spearman correlation (panel labels) ---\n")
print(corr_labels %>% select(metric, rho, p_value, n), row.names = FALSE)

## ---- Figure: one panel per metric, independent scales ---------------------

panel_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8.5)
  )

date_range_num <- range(as.numeric(combined_long$date))

fig10 <- ggplot(combined_long, aes(x = insect_value, y = fungal_value)) +
  geom_point(aes(fill = as.numeric(date), shape = lure), size = 2.4, color = "black", stroke = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30", linewidth = 0.6) +
  geom_text(data = corr_labels, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.3, size = 3, lineheight = 0.9) +
  facet_wrap(vars(metric), nrow = 1, scales = "free") +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac",
                        midpoint = mean(date_range_num), limits = date_range_num,
                        breaks = date_range_num, labels = format(range(combined_long$date), "%b %d"),
                        name = "Date") +
  scale_shape_manual(values = lure_shapes, name = "Lure") +
  guides(fill = guide_colorbar(order = 1), shape = guide_legend(order = 2, override.aes = list(fill = "grey50"))) +
  labs(x = "Insect alpha diversity", y = "Fungal alpha diversity") +
  panel_theme

ggsave(file.path(out_fig_dir, "fig10_insect_fungal_alpha_diversity_correlation.png"), fig10,
       width = 12, height = 4.5, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig10_insect_fungal_alpha_diversity_correlation.pdf"), fig10,
       width = 12, height = 4.5)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig10_insect_fungal_alpha_diversity_correlation.png"), "and .pdf\n")
