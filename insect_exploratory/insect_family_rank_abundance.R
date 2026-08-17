##############################################################################
# Insect family rank-abundance plot (diagnostic, not for publication).
#
# Shows every family caught (no prevalence/abundance filtering), ranked by
# total individuals caught (summed across all samples/sites), to visualize
# how steeply abundance drops off after the top families -- i.e. to justify
# restricting downstream community analyses to the top two families
# (Curculionidae + Latridiidae; see fig1_taxonomic_breakdown.R panel c and
# insect_taxonomic_groups_bars.R).
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

out_fig_dir <- "figures/insect_comps_2024"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

n_focal <- 2   # highlight the n_focal most-abundant families (the ones used downstream)

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
sample_cols <- sp_tab %>% select(where(is.numeric)) %>% colnames()

family_totals <- sp_tab %>%
  mutate(Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "col_names", values_to = "count") %>%
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

base_plot <- function(y_label) {
  ggplot(family_totals, aes(x = reorder(Family, rank), y = total, fill = focal)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Top families (used downstream)" = "#1f78b4",
                                  "Other families" = "grey70"),
                       name = NULL) +
    labs(x = NULL, y = y_label,
         title = "Insect family rank abundance",
         subtitle = paste0("All ", nrow(family_totals), " families caught; top ", n_focal,
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

p_log <- base_plot("Individuals caught (log scale, all samples/sites pooled)") +
  scale_y_log10(labels = scales::comma)

p_raw <- base_plot("Individuals caught (raw count, all samples/sites pooled)") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.03)))

ggsave(file.path(out_fig_dir, "insect_family_rank_abundance_log.png"), p_log,
       width = 24, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "insect_family_rank_abundance_log.pdf"), p_log,
       width = 24, height = 9)

ggsave(file.path(out_fig_dir, "insect_family_rank_abundance_raw.png"), p_raw,
       width = 24, height = 9, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "insect_family_rank_abundance_raw.pdf"), p_raw,
       width = 24, height = 9)

cat("\nDone. Wrote insect_family_rank_abundance_{log,raw}.{png,pdf} to", out_fig_dir, "\n")
