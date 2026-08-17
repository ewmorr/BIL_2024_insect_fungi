##############################################################################
# Insect family rank-abundance / rank-prevalence plots (diagnostic, not for
# publication).
#
# Shows every family caught (no prevalence/abundance filtering), ranked by
# (1) total individuals caught and (2) number of samples occurred in, to
# visualize how steeply abundance/prevalence drop off after the top
# families -- i.e. to justify restricting downstream community analyses to
# the top two families (Curculionidae + Latridiidae; see
# fig1_taxonomic_breakdown.R panel c and insect_taxonomic_groups_bars.R).
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

out_fig_dir <- "figures/insect_comps_2024"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

n_focal <- 2   # highlight the n_focal most-abundant families (the ones used downstream)

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
sample_cols <- sp_tab %>% select(where(is.numeric)) %>% colnames()
n_samples <- length(sample_cols)

insect_long <- sp_tab %>%
  mutate(Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "col_names", values_to = "count")

family_totals <- insect_long %>%
  group_by(Family) %>%
  summarize(total = sum(count), .groups = "drop") %>%
  filter(total > 0) %>%
  arrange(desc(total)) %>%
  mutate(rank = row_number(),
         focal = ifelse(rank <= n_focal, "Top families (used downstream)", "Other families"))

cat("All", nrow(family_totals), "families shown, ranked by total individuals caught.\n")
cat("Top", n_focal, "families:", paste(head(family_totals$Family, n_focal), collapse = ", "),
    "=", round(100 * sum(head(family_totals$total, n_focal)) / sum(family_totals$total), 1),
    "% of all individuals caught.\n")

# number of samples (out of n_samples) each family was detected in at all
# (count > 0), regardless of how many individuals -- same "focal" families
# (top n_focal by abundance, above) are highlighted here too for comparison
top_families <- head(family_totals$Family, n_focal)

family_prevalence <- insect_long %>%
  group_by(Family, col_names) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  filter(count > 0) %>%
  group_by(Family) %>%
  summarize(n_occ = n_distinct(col_names), .groups = "drop") %>%
  arrange(desc(n_occ)) %>%
  mutate(rank = row_number(),
         focal = ifelse(Family %in% top_families, "Top families (used downstream)", "Other families"))

cat("\nAll", nrow(family_prevalence), "families shown, ranked by number of samples occurred in (out of",
    n_samples, "samples total).\n")

base_plot <- function(data, y_label, title) {
  ggplot(data, aes(x = reorder(Family, rank), fill = focal)) +
    scale_fill_manual(values = c("Top families (used downstream)" = "#1f78b4",
                                  "Other families" = "grey70"),
                       name = NULL) +
    labs(x = NULL, y = y_label,
         title = title,
         subtitle = paste0("All ", nrow(data), " families caught; top ", n_focal,
                            " families used to filter downstream community analyses")) +
    theme_bw(base_size = 18) +
    theme(
      axis.text.x = element_text(size = 13, color = "black", angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(size = 24, face = "bold"),
      plot.subtitle = element_text(size = 17),
      legend.position = "top",
      legend.text = element_text(size = 16),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

p_log <- base_plot(family_totals, "Individuals caught (log scale, all samples/sites pooled)",
                    "Insect family rank abundance") +
  geom_col(aes(y = total), width = 0.7) +
  scale_y_log10(labels = scales::comma)

p_raw <- base_plot(family_totals, "Individuals caught (raw count, all samples/sites pooled)",
                    "Insect family rank abundance") +
  geom_col(aes(y = total), width = 0.7) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.03)))

p_prevalence <- base_plot(family_prevalence, paste0("Samples occurred in (of ", n_samples, " total)"),
                           "Insect family rank prevalence") +
  geom_col(aes(y = n_occ), width = 0.7) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.03)))

ggsave(file.path(out_fig_dir, "insect_family_rank_abundance_log.png"), p_log,
       width = 24, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "insect_family_rank_abundance_log.pdf"), p_log,
       width = 24, height = 9)

ggsave(file.path(out_fig_dir, "insect_family_rank_abundance_raw.png"), p_raw,
       width = 24, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "insect_family_rank_abundance_raw.pdf"), p_raw,
       width = 24, height = 9)

ggsave(file.path(out_fig_dir, "insect_family_rank_prevalence.png"), p_prevalence,
       width = 24, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "insect_family_rank_prevalence.pdf"), p_prevalence,
       width = 24, height = 9)

cat("\nDone. Wrote insect_family_rank_abundance_{log,raw}.{png,pdf} and insect_family_rank_prevalence.{png,pdf} to",
    out_fig_dir, "\n")
