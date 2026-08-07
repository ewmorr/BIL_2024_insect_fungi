##############################################################################
# Presentation figure 1: taxonomic (and trait-based) breakdown of the insect
# and fungal catch/communities, averaged across the four sites.
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
# Panel c: insect trap catch, Genus level, restricted to Curculionidae +
#          Latridiidae -- the two focal families used throughout the rest of
#          this project's community analyses (69% of all individuals
#          caught). Proportions are renormalized within this two-family
#          subset (i.e. "genus composition of the focal insect community"),
#          not as a fraction of the total catch shown in panel a.
# Panel d: fungal ITS2 community (same Kingdom==Fungi set as panel b),
#          summarized by ecological trait instead of taxonomy -- ASV genus
#          is joined to the FungalTraits database (Polme et al. 2020) by
#          genus name, and taxa are grouped by primary_lifestyle x growth
#          form.
#
# For all four panels: relative abundance is computed per sample, then
# averaged across all samples collected at each site (unweighted mean of
# proportions, not a pooled-count proportion), so no single high-read/
# high-catch sample dominates a site's bar.
#
# Standalone script -- edit n_top_insect_groups / n_top_fungal_groups /
# n_top_insect_genera / n_top_trait_groups below to change how many named
# categories are shown per panel before folding the rest into "Other".
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(gtable)
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

# Display name + panel order for Site -- set once here so it's inherited by
# every panel (a-d), including fungal_meta below (derived from insect_meta).
site_levels <- c("Durham", "Pease Airport", "Manchester Cedar Swamp", "Manchester Airport")
insect_meta <- insect_meta %>%
  mutate(Site = ifelse(Site == "Pease", "Pease Airport", Site),
         Site = factor(Site, levels = site_levels))

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
  theme(axis.text.x = element_blank())   # site labels shown once, on panel c below

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
  theme(axis.text.x = element_blank(),   # site labels shown once, on panel d below
        axis.text.y = element_blank(), axis.title.y = element_blank())   # y-axis shown once, on panel a

## ---- Panel c: insect Genus-level breakdown (Curculionidae + Latridiidae only) ----

n_top_insect_genera <- 10   # + "Other" = 11 categories

insect_focal_long <- sp_tab %>%
  filter(Family %in% c("Curculionidae", "Latridiidae")) %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "col_names", values_to = "count") %>%
  group_by(col_names, Genus) %>%
  summarize(count = sum(count), .groups = "drop")

genus_rank <- insect_focal_long %>%
  group_by(Genus) %>%
  summarize(total = sum(count), .groups = "drop") %>%
  arrange(desc(total))

top_genera <- head(genus_rank$Genus, n_top_insect_genera)

# proportions renormalized within the Curculionidae+Latridiidae subset (the
# sample total below is computed AFTER filtering to just these two families)
insect_focal_prop <- insect_focal_long %>%
  mutate(Genus = ifelse(Genus %in% top_genera, Genus, "Other")) %>%
  group_by(col_names, Genus) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  group_by(col_names) %>%
  mutate(sample_total = sum(count)) %>%
  ungroup() %>%
  filter(sample_total > 0) %>%   # drop any sample with zero Curculionidae/Latridiidae catch
  mutate(prop = count / sample_total) %>%
  left_join(insect_meta %>% select(col_names, Site), by = "col_names")

insect_focal_site_mean <- insect_focal_prop %>%
  group_by(Site, Genus) %>%
  summarize(mean_prop = mean(prop), .groups = "drop") %>%
  mutate(Genus = factor(Genus, levels = c(top_genera, "Other")))

cat("\nInsect genus panel (Curculionidae+Latridiidae only): levels shown (rank order):\n")
print(levels(insect_focal_site_mean$Genus))

panel_c <- ggplot(insect_focal_site_mean, aes(x = Site, y = mean_prop, fill = Genus)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = stack_palette(n_top_insect_genera), name = "Genus") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(y = "Mean relative abundance", tag = "c") +
  panel_theme +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.text = element_text(face = "italic", size = 9))

## ---- Panel d: fungal community by primary_lifestyle x growth form (FungalTraits) ----
# Joins ASV genus-level taxonomy to the FungalTraits database (Polme et al.
# 2020, github.com/TartuNaturalHistoryMuseum/FungalTraits) by genus name.
# FungalTraits lives outside this repo -- update fungal_traits_path below if
# it moves.

fungal_traits_path <- "/Users/ericmorrison/repo/FungalTraits_Polme/Polme_FungalTraits_1.2_ver_16Dec_2020.csv"
n_top_trait_groups <- 10   # + "Other" = 11 categories

fungal_traits <- read.csv(fungal_traits_path, fileEncoding = "latin1") %>%
  distinct(GENUS, .keep_all = TRUE) %>%   # a handful of genera appear twice; keep first
  transmute(
    Genus = GENUS,
    primary_lifestyle = ifelse(is.na(primary_lifestyle) | primary_lifestyle == "", "unspecified", primary_lifestyle),
    growth_form = ifelse(is.na(Growth_form_template) | Growth_form_template == "", "unspecified", Growth_form_template)
  )

genus_map <- taxonomy %>%
  filter(Kingdom == "k__Fungi") %>%
  transmute(taxon, Genus = sub("^g__", "", Genus)) %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", NA, Genus)) %>%
  left_join(fungal_traits, by = "Genus") %>%
  mutate(
    # is.na(primary_lifestyle) after the join means the genus had no match in
    # FungalTraits (as opposed to matching a genus whose lifestyle field was
    # blank, which was already recoded to "unspecified" above)
    trait_group = case_when(
      is.na(Genus) ~ "Unclassified genus",
      is.na(primary_lifestyle) ~ "No FungalTraits match",
      TRUE ~ paste0(gsub("_", " ", primary_lifestyle), " : ", gsub("_", " ", growth_form))
    )
  )

genus_map_used <- genus_map %>% filter(taxon %in% colnames(fungal_tab_f))
n_matched <- sum(!genus_map_used$trait_group %in% c("Unclassified genus", "No FungalTraits match"))
cat("\n", n_matched, "of", nrow(genus_map_used), "fungal ASVs (in this table) matched to a FungalTraits genus entry.\n")

trait_long <- as.data.frame(fungal_tab_f) %>%
  tibble::rownames_to_column("SequenceID") %>%
  pivot_longer(cols = -SequenceID, names_to = "taxon", values_to = "count") %>%
  filter(count > 0) %>%
  left_join(genus_map %>% select(taxon, trait_group), by = "taxon") %>%
  group_by(SequenceID, trait_group) %>%
  summarize(count = sum(count), .groups = "drop")

trait_rank <- trait_long %>%
  group_by(trait_group) %>%
  summarize(total = sum(count), .groups = "drop") %>%
  arrange(desc(total))

top_traits <- head(trait_rank$trait_group, n_top_trait_groups)

trait_prop <- trait_long %>%
  mutate(trait_group = ifelse(trait_group %in% top_traits, trait_group, "Other")) %>%
  group_by(SequenceID, trait_group) %>%
  summarize(count = sum(count), .groups = "drop") %>%
  group_by(SequenceID) %>%
  mutate(sample_total = sum(count),
         prop = ifelse(sample_total > 0, count / sample_total, 0)) %>%
  ungroup() %>%
  left_join(fungal_meta, by = "SequenceID")

trait_site_mean <- trait_prop %>%
  group_by(Site, trait_group) %>%
  summarize(mean_prop = mean(prop), .groups = "drop") %>%
  mutate(trait_group = factor(trait_group, levels = c(top_traits, "Other")))

cat("\nFungal trait panel: primary_lifestyle : growth_form levels shown (rank order):\n")
print(levels(trait_site_mean$trait_group))

panel_d <- ggplot(trait_site_mean, aes(x = Site, y = mean_prop, fill = trait_group)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = stack_palette(n_top_trait_groups), name = "Lifestyle : growth form") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  labs(y = "Mean relative abundance", tag = "d") +
  panel_theme +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        axis.text.y = element_blank(), axis.title.y = element_blank(),   # y-axis shown once, on panel c
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.35, "cm"))

## ---- Combine + save ----------------------------------------------------------
# grid.arrange's null-unit columns size each panel's plot area independently,
# so a's/c's and b's/d's panel widths drift apart whenever their legends
# differ in text width (e.g. panel d's long "lifestyle : growth form"
# labels vs. panel b's short Class names) -- the axis-b (Site) labels then
# don't line up between rows. gtable::rbind/cbind fixes this directly: stack
# each column with rbind (size="max" forces both panels in a column to the
# wider of the two), then cbind the two columns (size="max" matches row
# heights left vs. right).

ga <- ggplotGrob(panel_a); gb <- ggplotGrob(panel_b)
gc <- ggplotGrob(panel_c); gd <- ggplotGrob(panel_d)

left_col <- rbind(ga, gc, size = "max")    # rbind/cbind dispatch to gtable's S3 methods once library(gtable) is loaded
right_col <- rbind(gb, gd, size = "max")
combined <- cbind(left_col, right_col, size = "max")

ggsave(file.path(out_fig_dir, "fig1_taxonomic_breakdown.png"), combined, width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig1_taxonomic_breakdown.pdf"), combined, width = 12, height = 10)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig1_taxonomic_breakdown.png"), "and .pdf\n")
