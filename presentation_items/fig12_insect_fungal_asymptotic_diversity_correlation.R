##############################################################################
# Presentation figure 12: per-sample insect vs. fungal ASYMPTOTIC (Chao-
# extrapolated) diversity correlation. Lineage C, full insect table -- same
# construction as fig10_insect_fungal_alpha_diversity_correlation.R, which
# this figure reproduces one-for-one except for the metric: observed
# diversity swapped for the asymptotic estimates from insect_fungal_
# asymptotic_richness_iNEXT.full_insect_table.r (data/compare_insects_fungi_
# asymptotic_richness_iNEXT_full_insect_table/). "Simpson" here is Hill
# number q=2 diversity (higher = more even), the OPPOSITE direction from
# fig10's Simpson DOMINANCE -- see that script's header and project_
# organization.md's Lineage C section.
#
# Not split by lure, same reasoning as fig10: the interaction test in
# insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r (fungal ~
# site+date+insect*lure vs. without the interaction, 999-perm F-test) found
# no evidence lure modifies the insect~fungal relationship for any of the 3
# asymptotic metrics (p_perm >= 0.075, see asymptotic_diversity_insect_
# fungal_interaction_by_lure.csv) -- so the pooled relationship is shown.
#
# One point per sample, insect value (x) vs. fungal value (y), colored by
# collection date and shaped by lure (the trap-type variable -- Ethanol,
# Alpha-pinene_EtOH, Ips -- constant within a trap, matching fig9/fig11's
# lure faceting even though fig12 itself pools lures for the correlation).
# facet_wrap(scales="free") gives each metric panel its own independent x/y
# range (asymptotic richness in particular spans a much wider range than
# observed richness, since the estimator's confidence interval upper bound
# can be large for poorly-covered samples).
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
source("library/library.R")

out_fig_dir <- "figures/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

asymp_data_dir <- "data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table"

## ---- Load + join insect + fungal asymptotic diversity (Lineage C, full insect table, iNEXT) ----

metrics <- c("shannon", "simpson", "richness")
metric_labels <- c(shannon = "Asymptotic Shannon diversity (Hill q=1)",
                    simpson = "Asymptotic Simpson diversity (Hill q=2)",
                    richness = "Asymptotic species richness")

lure_levels <- c("Ethanol", "Alpha-pinene_EtOH", "Ips")
lure_shapes <- c(Ethanol = 21, `Alpha-pinene_EtOH` = 24, Ips = 22)  # filled circle/triangle/square

insect_asymp <- read.csv(file.path(asymp_data_dir, "insect_asymptotic_diversity.csv")) %>%
  mutate(date = as.Date(date), lure = factor(lure, levels = lure_levels))
fungal_asymp <- read.csv(file.path(asymp_data_dir, "fungal_asymptotic_diversity.csv")) %>%
  mutate(date = as.Date(date))

combined <- insect_asymp %>%
  select(sample_id, date, lure, insect_shannon = shannon_asymptotic,
         insect_simpson = simpson_asymptotic, insect_richness = richness_asymptotic) %>%
  inner_join(
    fungal_asymp %>% select(sample_id, fungal_shannon = shannon_asymptotic,
                             fungal_simpson = simpson_asymptotic, fungal_richness = richness_asymptotic),
    by = "sample_id"
  )
cat(nrow(combined), "samples with both insect and fungal asymptotic diversity.\n")

combined_long <- bind_rows(lapply(metrics, function(m) {
  combined %>% transmute(sample_id, date, lure, metric = metric_labels[m],
                          insect_value = .data[[paste0("insect_", m)]],
                          fungal_value = .data[[paste0("fungal_", m)]])
})) %>% mutate(metric = factor(metric, levels = metric_labels))

## ---- Pooled per-metric Spearman correlation, for the in-panel label -------
# Read the ASYMPTOTIC rows from the upstream script's output rather than
# recomputing -- same test, same numbers.

cross_corr <- read.csv(file.path(asymp_data_dir, "asymptotic_diversity_cross_community_correlation.csv")) %>%
  filter(value_type == "asymptotic")
corr_labels <- cross_corr %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels),
         label = paste0("rho = ", round(rho, 2), "\np = ", signif(p_value, 2)))

cat("\n--- Pooled per-metric Spearman correlation, asymptotic values (panel labels) ---\n")
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

fig12 <- ggplot(combined_long, aes(x = insect_value, y = fungal_value)) +
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
  labs(x = "Insect asymptotic diversity", y = "Fungal asymptotic diversity") +
  panel_theme

ggsave(file.path(out_fig_dir, "fig12_insect_fungal_asymptotic_diversity_correlation.png"), fig12,
       width = 12, height = 4.5, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig12_insect_fungal_asymptotic_diversity_correlation.pdf"), fig12,
       width = 12, height = 4.5)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig12_insect_fungal_asymptotic_diversity_correlation.png"), "and .pdf\n")
