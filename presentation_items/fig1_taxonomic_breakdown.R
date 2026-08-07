##############################################################################
# Presentation figure 1: taxonomic breakdown of the insect and fungal
# catch/communities, averaged across the four sites.
#
# Panel a: insect trap catch, Family level (Coleoptera makes up ~98% of
#          individuals caught -- see insect_exploratory/insect_ords.R -- so
#          Order level collapses to one bar; Family is the coarsest level
#          that still shows structure, and lines up with the Curculionidae/
#          Latridiidae focus of the rest of this project).
# Panel b: fungal ITS2 community, Class level, restricted to Kingdom==Fungi
#          (97.8% of reads; the remainder is co-amplified non-fungal
#          eukaryotic DNA -- plant, arthropod, etc. -- dropped here and
#          proportions renormalized over fungal reads only).
#
# For both panels: relative abundance is computed per sample, then averaged
# across all samples collected at each site (unweighted mean of proportions,
# not a pooled-count proportion), so no single high-read/high-catch sample
# dominates a site's bar.
#
# Standalone script -- edit n_top_insect_groups / n_top_fungal_groups below
# to change how many named categories are shown before folding the rest into
# "Other".
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("library/library.R")

out_fig_dir <- "figures/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

n_top_insect_groups <- 10   # + "Other" = 11 categories, ~90% of individuals
n_top_fungal_groups <- 10   # + "Other" = 11 categories, ~97% of fungal reads

## ---- shared look ------------------------------------------------------------

panel_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    plot.tag = element_text(size = 14, face = "bold"),
    panel.grid = element_blank()
  )

# ColorBrewer "Paired" 12-class palette (library/library.R::twelvePaired),
# already used elsewhere in this project. Assigned in fixed rank order
# (most-abundant group first) so a given taxon always gets the same color
# across panels/re-runs; "Other" is always grey.
stack_palette <- function(n_named) {
  pal <- c(twelvePaired[seq_len(n_named)], "grey65")
  names(pal) <- NULL
  pal
}

## ---- Panel a: insect Family-level breakdown by site ------------------------

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta <- read.csv("data/metadata/insect_community_metadata.csv")

sample_cols <- sp_tab %>% select(where(is.numeric)) %>% colnames()

insect_long <- sp_tab %>%
  mutate(Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "col_names", values_to = "count") %>%
  group_by(col_names, Family) %>%
  summarize(count = sum(count), .groups = "drop")

# rank families by total individuals caught (across all samples/sites)
family_rank <- insect_long %>%
  group_by(Family) %>%
  summarize(total = sum(count), .groups = "drop") %>%
  arrange(desc(total))

top_families <- head(family_rank$Family, n_top_insect_groups)

insect_prop <- insect_long %>%
  mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%
  group_by(col_names, Family) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  group_by(col_names) %>%
  mutate(sample_total = sum(count)) %>%
  ungroup() %>%
  filter(sample_total > 0) %>%    # drop the one zero-catch trap-date (P3.I..6.12);
  mutate(prop = count / sample_total) %>%  # a 0/0 sample carries no compositional info and would deflate its site's average
  left_join(insect_meta %>% select(col_names, Site), by = "col_names")

insect_site_mean <- insect_prop %>%
  group_by(Site, Family) %>%
  summarize(mean_prop = mean(prop), .groups = "drop") %>%
  mutate(Family = factor(Family, levels = c(top_families, "Other")))

cat("Insect panel: Family levels shown (rank order):\n")
print(levels(insect_site_mean$Family))

panel_a <- ggplot(insect_site_mean, aes(x = Site, y = mean_prop, fill = Family)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = stack_palette(n_top_insect_groups), name = "Family") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(y = "Mean relative abundance", tag = "a") +
  panel_theme +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

## ---- Panel b: fungal Class-level breakdown by site --------------------------

fungal_tab <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv") %>% rename(taxon = X)

# metadata table keyed by fungal SequenceID -> Site (same file used for
# insect col_names -> Site above; SequenceID is the ASV_tab row identifier)
fungal_meta <- insect_meta %>% select(SequenceID, Site) %>% distinct()

fungi_only_taxa <- taxonomy %>% filter(Kingdom == "k__Fungi") %>% pull(taxon)
fungal_tab_f <- fungal_tab[, colnames(fungal_tab) %in% fungi_only_taxa, drop = FALSE]
cat("\n", ncol(fungal_tab_f), "of", ncol(fungal_tab), "ASVs in this table are Kingdom==Fungi.\n")

class_map <- taxonomy %>%
  filter(Kingdom == "k__Fungi") %>%
  transmute(taxon, Class = ifelse(is.na(Class) | Class == "", "Unclassified", Class))

fungal_long <- as.data.frame(fungal_tab_f) %>%
  tibble::rownames_to_column("SequenceID") %>%
  pivot_longer(cols = -SequenceID, names_to = "taxon", values_to = "count") %>%
  filter(count > 0) %>%
  left_join(class_map, by = "taxon") %>%
  group_by(SequenceID, Class) %>%
  summarize(count = sum(count), .groups = "drop")

class_rank <- fungal_long %>%
  group_by(Class) %>%
  summarize(total = sum(count), .groups = "drop") %>%
  arrange(desc(total))

top_classes <- head(class_rank$Class, n_top_fungal_groups)

fungal_prop <- fungal_long %>%
  mutate(Class = ifelse(Class %in% top_classes, Class, "Other")) %>%
  group_by(SequenceID, Class) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  group_by(SequenceID) %>%
  mutate(sample_total = sum(count),
         prop = ifelse(sample_total > 0, count / sample_total, 0)) %>%
  ungroup() %>%
  left_join(fungal_meta, by = "SequenceID")

fungal_site_mean <- fungal_prop %>%
  group_by(Site, Class) %>%
  summarize(mean_prop = mean(prop), .groups = "drop") %>%
  mutate(Class = sub("^c__", "", Class),
         Class = factor(Class, levels = c(sub("^c__", "", top_classes), "Other")))

cat("\nFungal panel: Class levels shown (rank order):\n")
print(levels(fungal_site_mean$Class))

panel_b <- ggplot(fungal_site_mean, aes(x = Site, y = mean_prop, fill = Class)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = stack_palette(n_top_fungal_groups), name = "Class") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(y = "Mean relative abundance", tag = "b") +
  panel_theme +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

## ---- Combine + save ----------------------------------------------------------

combined <- gridExtra::arrangeGrob(panel_a, panel_b, ncol = 2)

ggsave(file.path(out_fig_dir, "fig1_taxonomic_breakdown.png"), combined, width = 11, height = 5, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig1_taxonomic_breakdown.pdf"), combined, width = 11, height = 5)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig1_taxonomic_breakdown.png"), "and .pdf\n")
