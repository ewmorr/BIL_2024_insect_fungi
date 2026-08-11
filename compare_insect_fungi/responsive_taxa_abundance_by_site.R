# Exploratory: how much of the fungal community's sequence abundance is
# made up of "responsive" taxa (significantly associated with either
# collection date or insect_PCoA1) vs. everything else, and does that split
# vary by site? Stacked bar, one bar per site, average relative abundance.
#
# "Responsive" = union of:
#  - date-associated (q<0.10 in the unbiased, all-taxa screen):
#    data/2024_fungi/fungal_taxa_date_association.all_taxa.csv
#  - insect-associated (q<0.10, no lure/date covariate -- same screen used
#    in fig3/fig4 panel c): data/compare_insects_fungi_top3axes/
#    fungal_insect_association_results.PCoA1_only.no_lure.csv (tested on the
#    top-200 CoCA-loading candidates only, not all taxa -- a taxon absent
#    from this file was simply never tested for this association, not shown
#    to be non-significant)
# "Unresponsive" = every other Kingdom==Fungi ASV, including both taxa
# tested-and-not-significant and taxa never tested (below the >=5-sample
# prevalence filter used by both screens, i.e. too rare to test) -- treated
# as one bucket since from this figure's perspective both are simply "not
# identified as responsive."
#
# Relative abundance is computed per sample as a fraction of that sample's
# total Kingdom==Fungi reads (same fungal-community denominator fig1 panel b
# uses), then averaged unweighted across samples within each site -- same
# convention as fig1_taxonomic_breakdown.R.

library(dplyr)
library(tidyr)
library(ggplot2)

fig_dir <- "figures"
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

## ---- Responsive taxon set (union of date- and insect_PCoA1-significant) ----

date_sig <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  filter(q_value < q_threshold) %>% pull(taxon)
insect_sig <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv") %>%
  filter(q_value < q_threshold) %>% pull(taxon)
responsive_taxa <- union(date_sig, insect_sig)

cat(length(date_sig), "date-associated,", length(insect_sig), "insect-associated,",
    length(responsive_taxa), "responsive taxa total (", length(intersect(date_sig, insect_sig)),
    "overlap ) out of", ncol(fungal), "Kingdom==Fungi ASVs in this matched dataset.\n")

## ---- Per-sample relative abundance, responsive vs. unresponsive ------------

long <- as.data.frame(fungal) %>%
  tibble::rownames_to_column("SequenceID") %>%
  pivot_longer(cols = -SequenceID, names_to = "taxon", values_to = "count") %>%
  filter(count > 0) %>%
  mutate(group = ifelse(taxon %in% responsive_taxa, "Responsive (date- or insect-associated)", "Unresponsive")) %>%
  group_by(SequenceID, group) %>%
  summarise(count = sum(count), .groups = "drop")

prop <- long %>%
  group_by(SequenceID) %>%
  mutate(sample_total = sum(count), prop = count / sample_total) %>%
  ungroup() %>%
  left_join(meta, by = "SequenceID")

site_mean <- prop %>%
  group_by(Site, group) %>%
  summarise(mean_prop = mean(prop), n_samples = n_distinct(SequenceID), .groups = "drop") %>%
  mutate(group = factor(group, levels = c("Responsive (date- or insect-associated)", "Unresponsive")))

cat("\nPer-site mean relative abundance:\n")
print(as.data.frame(site_mean %>% mutate(mean_prop = round(mean_prop, 4))))

## ---- Plot -------------------------------------------------------------------

group_colors <- c("Responsive (date- or insect-associated)" = "#0072B2", "Unresponsive" = "grey70")

p <- ggplot(site_mean, aes(x = Site, y = mean_prop, fill = group)) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(y = "Mean relative sequence abundance", x = NULL,
       title = "Relative abundance of date-/insect-associated fungal taxa by site",
       subtitle = paste0(length(responsive_taxa), " of ", ncol(fungal), " ASVs responsive (q<0.10): ",
                          length(date_sig), " date-assoc., ", length(insect_sig), " insect-assoc., ",
                          length(intersect(date_sig, insect_sig)), " both")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

print(p)
ggsave(file.path(fig_dir, "responsive_taxa_abundance_by_site.png"), p, width = 7.5, height = 6, dpi = 300, bg = "white")

cat("\nWrote", file.path(fig_dir, "responsive_taxa_abundance_by_site.png"), "\n")
