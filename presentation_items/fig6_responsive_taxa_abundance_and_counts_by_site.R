##############################################################################
# Presentation figure 6: formalizes the finding from compare_insect_fungi/
# responsive_taxa_abundance_by_site.R (a small share of fungal ASVs carries
# a large share of sequence abundance) by showing relative abundance and
# ASV counts side by side, and splits that script's combined "responsive"
# category into its two component screens (date, insect_PCoA1) shown as
# separate panels instead of a single unioned category.
#
# Panel a: date-associated vs. not, relative sequence abundance, by site
# Panel b: insect_PCoA1-associated vs. not, relative sequence abundance, by site
# Panel c: date-associated vs. not, ASV count (richness), by site
# Panel d: insect_PCoA1-associated vs. not, ASV count (richness), by site
#
# "Date-associated"/"insect-associated" = q<0.10 in the same two screens as
# responsive_taxa_abundance_by_site.R (data/2024_fungi/fungal_taxa_date_
# association.all_taxa.csv; data/compare_insects_fungi_top3axes/fungal_
# insect_association_results.PCoA1_only.no_lure.csv -- tested on the top-200
# CoCA-loading candidates only, not all taxa). "Not associated" is
# everything else, including taxa tested-and-not-significant and taxa never
# tested (below the >=5-sample prevalence filter both screens use).
#
# Row 1 (a/b) is a proportion (relative abundance, sums to 100% per site);
# row 2 (c/d) is a raw count (number of distinct ASVs observed in at least
# one sample at that site) -- deliberately left unnormalized so the
# absolute rarity of the associated category is visible directly, not just
# its share.
#
# Axes/legends are de-duplicated across the 2x2 grid: row 1 (a/b) hides its
# x-axis text and legend (repeated in row 2 below); the right column (b/d)
# hides its y-axis title/text (repeated in the left column).
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)

out_fig_dir <- "figures/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

q_threshold <- 0.10

## ---- Load + match the same 69-sample insect/fungal dataset the association tests used --

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
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy[colnames(fungal), "Kingdom"] == "k__Fungi"
fungal <- fungal[, is_fungus, drop = FALSE]

cat("Matched dataset:", nrow(fungal), "samples,", ncol(fungal), "Kingdom==Fungi ASVs.\n")

site_levels <- c("Durham", "Pease Airport", "Manchester Cedar Swamp", "Manchester Airport")
meta <- id_map %>%
  transmute(SequenceID, Site = ifelse(Site == "Pease", "Pease Airport", Site)) %>%
  mutate(Site = factor(Site, levels = site_levels))

## ---- The two association screens (kept separate, not unioned) -------------

date_sig <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  filter(q_value < q_threshold) %>% pull(taxon)
insect_sig <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv") %>%
  filter(q_value < q_threshold) %>% pull(taxon)

cat(length(date_sig), "date-associated,", length(insect_sig), "insect-associated taxa,",
    "out of", ncol(fungal), "Kingdom==Fungi ASVs.\n")

## ---- Taxon-level long table (one row per sample x present taxon) ----------

taxon_long <- as.data.frame(fungal) %>%
  tibble::rownames_to_column("SequenceID") %>%
  pivot_longer(cols = -SequenceID, names_to = "taxon", values_to = "count") %>%
  filter(count > 0) %>%
  left_join(meta, by = "SequenceID")

## ---- Shared plotting logic --------------------------------------------------

pos_color <- "#0072B2"
neg_color <- "grey70"

# show_x_text: site labels (kept on row 2 only); show_legend: kept on row 2
# only; show_y: y-axis title + text (kept on the left column only).
panel_theme <- function(show_x_text = TRUE, show_legend = TRUE, show_y = TRUE) {
  theme_bw() +
    theme(
      axis.text = element_text(size = 10, color = "black"),
      axis.text.x = if (show_x_text) element_text(angle = 30, hjust = 1) else element_blank(),
      axis.ticks.x = if (show_x_text) element_line() else element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = if (show_y) element_text(size = 10, color = "black") else element_blank(),
      axis.title.y = if (show_y) element_text() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      plot.subtitle = element_text(size = 8.5, color = "grey30"),
      plot.tag = element_text(size = 14, face = "bold"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_blank(),
      panel.grid = element_blank()
    )
}

# Panel a/b: mean relative sequence abundance by site, stacked to 100%.
make_abundance_panel <- function(sig_taxa, pos_label, neg_label, title, tag, show_y = TRUE) {
  prop <- taxon_long %>%
    mutate(group = ifelse(taxon %in% sig_taxa, pos_label, neg_label)) %>%
    group_by(SequenceID, group) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    group_by(SequenceID) %>%
    mutate(sample_total = sum(count), prop = count / sample_total) %>%
    ungroup() %>%
    left_join(meta, by = "SequenceID")

  site_mean <- prop %>%
    group_by(Site, group) %>%
    summarise(mean_prop = mean(prop), .groups = "drop") %>%
    mutate(group = factor(group, levels = c(pos_label, neg_label)))

  ggplot(site_mean, aes(x = Site, y = mean_prop, fill = group)) +
    geom_col(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = setNames(c(pos_color, neg_color), c(pos_label, neg_label))) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    labs(y = "Mean relative sequence abundance", title = title, tag = tag) +
    panel_theme(show_x_text = FALSE, show_legend = FALSE, show_y = show_y)
}

# Panel c/d: number of distinct ASVs observed at each site, split by category.
make_count_panel <- function(sig_taxa, pos_label, neg_label, title, tag, show_y = TRUE) {
  site_counts <- taxon_long %>%
    distinct(Site, taxon) %>%
    mutate(group = ifelse(taxon %in% sig_taxa, pos_label, neg_label)) %>%
    count(Site, group) %>%
    mutate(group = factor(group, levels = c(pos_label, neg_label)))

  ggplot(site_counts, aes(x = Site, y = n, fill = group)) +
    geom_col(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = setNames(c(pos_color, neg_color), c(pos_label, neg_label))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y = "ASVs present (count)", title = title, tag = tag) +
    panel_theme(show_x_text = TRUE, show_legend = TRUE, show_y = show_y)
}

## ---- Build the four panels ---------------------------------------------------

panel_a <- make_abundance_panel(date_sig, "Date-associated", "Not date-associated",
                                 "Date-associated: relative abundance", "a", show_y = TRUE)
panel_b <- make_abundance_panel(insect_sig, "Insect-associated", "Not insect-associated",
                                 "Insect-associated: relative abundance", "b", show_y = FALSE)
panel_c <- make_count_panel(date_sig, "Date-associated", "Not date-associated",
                             "Date-associated: ASV count", "c", show_y = TRUE)
panel_d <- make_count_panel(insect_sig, "Insect-associated", "Not insect-associated",
                             "Insect-associated: ASV count", "d", show_y = FALSE)

## ---- Combine + save -----------------------------------------------------------
# All four panels share a flat (unfaceted) gtable structure, so gtable::
# rbind/cbind aligns axes cleanly -- same approach as fig1_taxonomic_
# breakdown.R's 2x2 layout.

ga <- ggplotGrob(panel_a); gb <- ggplotGrob(panel_b)
gc <- ggplotGrob(panel_c); gd <- ggplotGrob(panel_d)
left_col <- rbind(ga, gc, size = "max")
right_col <- rbind(gb, gd, size = "max")
combined <- cbind(left_col, right_col, size = "max")

ggsave(file.path(out_fig_dir, "fig6_responsive_taxa_abundance_and_counts_by_site.png"), combined,
       width = 9, height = 8, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig6_responsive_taxa_abundance_and_counts_by_site.pdf"), combined,
       width = 9, height = 8, bg = "white")

cat("\nDone. Wrote", file.path(out_fig_dir, "fig6_responsive_taxa_abundance_and_counts_by_site.png"), "and .pdf\n")
