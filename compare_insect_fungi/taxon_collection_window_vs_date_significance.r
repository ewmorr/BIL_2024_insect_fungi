##############################################################################
# Taxon Collection-Window vs. Date-Significance Diagnostic
#
# Motivating concern: the linear+quadratic date screen (insect_taxa_date_
# association.prevalence_filtered_insect.r for Lineage B insect taxa;
# fungal_community_seasonality.r for the shared/lineage-invariant fungal
# taxa) can only detect a real early-/late-season signal if a taxon has
# enough SPREAD across the 6 collection dates (5/1-7/10/24) to fit a slope
# against. A taxon confined to a narrow window near one edge of the season
# (e.g. caught only on the first collection date) could carry an obvious
# early-season signal that the model has little power to resolve, and would
# be misclassified as "no significant date pattern" purely from a lack of
# spread, not a lack of real seasonality.
#
# This is a diagnostic, not a new lineage or a change to any existing
# screen: for each taxon already tested in the two screens above, computes
# a seasonality-agnostic COLLECTION WINDOW (last collection date the taxon
# was PRESENT - first collection date present; raw presence/absence, no
# abundance weighting, no direction) and tests whether window length (or
# the closely related # of distinct dates present) predicts whether a
# taxon was called date-significant (q<0.10, linear OR quadratic term --
# the same "date-associated" union already used elsewhere in this project,
# e.g. fig6's date-associated definition).
#
# Each taxon set's window is computed by exactly reproducing that script's
# own load/match/filter logic, so the window is measured on the identical
# sample set and taxon list the significance screen was fit on.
##############################################################################

## ---- 0. Setup -------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

set.seed(1)

out_data_dir <- "data/compare_insects_fungi_collection_window"
out_fig_dir  <- "figures/compare_insects_fungi_collection_window"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

# Presence-based window/date-spread summary, shared by both taxon sets below.
compute_window <- function(presence_mat, meta) {
  stopifnot(all(rownames(presence_mat) == meta$sample_id))
  bind_rows(lapply(colnames(presence_mat), function(tax) {
    present <- presence_mat[, tax] > 0
    dts <- meta$date[present]
    data.frame(
      taxon = tax,
      n_samples_present = sum(present),
      n_dates_present = length(unique(dts)),
      first_date = min(dts),
      last_date = max(dts),
      window_days = as.numeric(max(dts) - min(dts))
    )
  }))
}

## ---- 1. Insect (Lineage B, 153-taxon prevalence-filtered table) -----------
## Reproduces insect_taxa_date_association.prevalence_filtered_insect.r's
## Step 1 load/match/filter exactly.

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab <- sp_tab[!is.na(sp_tab$Finest.ID), ]
stopifnot(!any(duplicated(sp_tab$Finest.ID)))

sp_tab.t <- t(sp_tab %>% select(where(is.numeric)))
colnames(sp_tab.t) <- sp_tab$Finest.ID

keep_taxa <- colSums(sp_tab.t > 0) >= 5
insect_full <- sp_tab.t[, keep_taxa, drop = FALSE]
insect_full <- insect_full[rowSums(insect_full) > 0, ]

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
id_map <- id_map %>% filter(SequenceID %in% shared_ids)

insect <- insect_full[id_map$col_names, , drop = FALSE]
rownames(insect) <- id_map$SequenceID

insect_meta <- id_map %>% transmute(sample_id = SequenceID, date = lubridate::mdy(CollectionDate))
stopifnot(all(rownames(insect) == insect_meta$sample_id))
cat("Insect (Lineage B): matched dataset", nrow(insect), "samples,", ncol(insect), "taxa.\n")

insect_window <- compute_window(insect, insect_meta)

insect_sig <- read.csv("data/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_taxa_date_association.csv")
insect_combined <- insect_window %>%
  left_join(insect_sig %>% select(taxon, t_stat, q_value, t_stat_quad, q_value_quad, shape, beta_std, beta_std_quad),
            by = "taxon") %>%
  mutate(date_significant = q_value < 0.10 | q_value_quad < 0.10)
stopifnot(nrow(insect_combined) == nrow(insect_window), !any(is.na(insect_combined$q_value)))

write.csv(insect_combined, file.path(out_data_dir, "insect_taxon_collection_window.csv"), row.names = FALSE)

## ---- 2. Fungal (shared/lineage-invariant prevalence-filtered screen) ------
## Reproduces fungal_community_seasonality.r's Step 1 load/match/filter
## exactly (the 69-sample Curculionidae+Latridiidae-matched set used to
## define the fungal sample match -- shared by both Lineage A and B fungal
## results).

sp_tab.curcus_latris <- sp_tab %>% filter(Family %in% c("Curculionidae", "Latridiidae"))
sp_tab.curcus_latris.t <- t(sp_tab.curcus_latris %>% select(where(is.numeric)))
colnames(sp_tab.curcus_latris.t) <- sp_tab.curcus_latris$Finest.ID
insect_full_a <- sp_tab.curcus_latris.t[, colSums(sp_tab.curcus_latris.t) > 1, drop = FALSE]
insect_full_a <- insect_full_a[rowSums(insect_full_a) > 0, ]

id_map_a <- insect_meta_full %>% filter(col_names %in% rownames(insect_full_a))
stopifnot(!any(duplicated(id_map_a$col_names)))
shared_ids_a <- intersect(id_map_a$SequenceID, rownames(fungal_full))
id_map_a <- id_map_a %>% filter(SequenceID %in% shared_ids_a)

fungal <- fungal_full[id_map_a$SequenceID, , drop = FALSE]

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy[colnames(fungal), "Kingdom"] == "k__Fungi"
fungal <- fungal[, is_fungus, drop = FALSE]

fungal_meta <- id_map_a %>% transmute(sample_id = SequenceID, date = lubridate::mdy(CollectionDate))
stopifnot(all(rownames(fungal) == fungal_meta$sample_id))

keep_f <- colSums(fungal > 0) >= 5
fungal_f <- fungal[, keep_f, drop = FALSE]
cat("Fungal (shared screen): matched dataset", nrow(fungal_f), "samples,", ncol(fungal_f), "taxa.\n")

fungal_window <- compute_window(fungal_f, fungal_meta)

fungal_sig <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv")
fungal_combined <- fungal_window %>%
  left_join(fungal_sig %>% select(taxon, t_stat, q_value, t_stat_quad, q_value_quad, shape, beta_std, beta_std_quad),
            by = "taxon") %>%
  mutate(date_significant = q_value < 0.10 | q_value_quad < 0.10)
stopifnot(nrow(fungal_combined) == nrow(fungal_window), !any(is.na(fungal_combined$q_value)))

write.csv(fungal_combined, file.path(out_data_dir, "fungal_taxon_collection_window.csv"), row.names = FALSE)

## ---- 3. Window length vs. date-significance --------------------------
## Two views per taxon set: (a) an unconditional Wilcoxon rank-sum test of
## window_days by date_significant, and (b) a logistic model of
## date_significant on window_days ADJUSTING for n_samples_present, since a
## taxon detected in more samples will mechanically tend to span more dates
## -- (b) isolates whether a narrow window predicts a missed call over and
## above simply being a rarely-detected taxon.

run_window_tests <- function(df, taxon_set_label) {
  wt <- wilcox.test(window_days ~ date_significant, data = df)
  glm_fit <- glm(date_significant ~ window_days + n_samples_present, data = df, family = binomial)
  glm_coef <- coef(summary(glm_fit))

  data.frame(
    taxon_set = taxon_set_label,
    n_taxa = nrow(df),
    n_date_significant = sum(df$date_significant),
    median_window_days.significant = median(df$window_days[df$date_significant]),
    median_window_days.not_significant = median(df$window_days[!df$date_significant]),
    wilcox_p = wt$p.value,
    glm_window_days_coef = glm_coef["window_days", "Estimate"],
    glm_window_days_p = glm_coef["window_days", "Pr(>|z|)"],
    glm_n_samples_present_coef = glm_coef["n_samples_present", "Estimate"],
    glm_n_samples_present_p = glm_coef["n_samples_present", "Pr(>|z|)"]
  )
}

test_results <- bind_rows(
  run_window_tests(insect_combined, "insect (Lineage B, 153-taxon prevalence-filtered)"),
  run_window_tests(fungal_combined, "fungal (shared prevalence-filtered screen)")
)
cat("\n--- Collection window vs. date-significance: test summary ---\n")
print(test_results)
write.csv(test_results, file.path(out_data_dir, "collection_window_vs_significance_tests.csv"), row.names = FALSE)

## ---- 4. Plots --------------------------------------------------------

make_window_plot <- function(df, title) {
  ggplot(df, aes(x = date_significant, y = window_days)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(aes(color = n_samples_present), width = 0.15, height = 1, alpha = 0.6) +
    scale_color_viridis_c(name = "# samples\npresent") +
    scale_x_discrete(labels = c("FALSE" = "not date-significant",
                                 "TRUE" = "date-significant\n(q<0.10, linear or quadratic)")) +
    labs(x = NULL, y = "Collection window (days; last - first date present)", title = title) +
    theme_bw()
}

p_insect <- make_window_plot(insect_combined, "Insect taxa (Lineage B, 153-taxon prevalence-filtered)")
p_fungal <- make_window_plot(fungal_combined, "Fungal taxa (shared prevalence-filtered screen)")
p_combined <- p_insect + p_fungal + plot_layout(ncol = 2)

ggsave(file.path(out_fig_dir, "collection_window_by_significance.png"), p_combined, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out_fig_dir, "collection_window_by_significance.pdf"), p_combined, width = 12, height = 5.5)

# Histogram of taxon counts (observations) per collection-window bin -- only
# 6 possible window_days values (0/14/28/42/56/70, 14-day collection
# spacing), so a stacked bar chart is the discrete-histogram equivalent.
# Stacked (not dodged) so the total bar height reads as "how many taxa land
# in this window bin", with the date-significant split visible inside each
# bar -- this is what the boxplot above can't show (equal medians say
# nothing about how many taxa sit at each window value).
make_window_hist <- function(df, title) {
  ggplot(df, aes(x = factor(window_days), fill = date_significant)) +
    geom_bar(color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c("FALSE" = "grey70", "TRUE" = "#0072B2"),
      labels = c("FALSE" = "not date-significant", "TRUE" = "date-significant"),
      name = NULL
    ) +
    labs(x = "Collection window (days; last - first date present)", y = "# taxa", title = title) +
    theme_bw() +
    theme(legend.position = "bottom")
}

p_insect_hist <- make_window_hist(insect_combined, "Insect taxa (Lineage B, 153-taxon prevalence-filtered)")
p_fungal_hist <- make_window_hist(fungal_combined, "Fungal taxa (shared prevalence-filtered screen)")
p_hist_combined <- p_insect_hist + p_fungal_hist + plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_fig_dir, "collection_window_histogram.png"), p_hist_combined, width = 12, height = 5.5, dpi = 300)
ggsave(file.path(out_fig_dir, "collection_window_histogram.pdf"), p_hist_combined, width = 12, height = 5.5)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
