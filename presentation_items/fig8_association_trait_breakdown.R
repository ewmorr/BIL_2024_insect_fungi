##############################################################################
# Presentation figure 8: trophic mode breakdown of significant (q<0.10)
# fungal taxa from the association screens shown in fig3_volcano_plots.R /
# fig4_association_taxonomic_breakdown.R. Where fig4 asks "which taxonomic
# groups do the significant hits belong to", this figure asks the same
# question about ecological trait instead of taxonomy: ASV genus is joined
# to the FungalTraits database (Polme et al. 2020) by genus name, exactly as
# in fig1_taxonomic_breakdown.R panel d, and taxa are grouped by
# primary_lifestyle ("trophic mode"). growth_form is deliberately NOT part
# of the label/grouping here (unlike fig1 panel d and fig7, its
# growth-form-faceted companion) -- once bars are already split by
# taxonomic Class, growth form adds little beyond what Class x lifestyle
# already implies (e.g. almost everything in Saccharomycetes is some kind
# of yeast), so folding it in mostly just multiplied the number of bars
# without adding information.
# Same diverging-bar-chart technique as fig4 panels b/c, faceted by
# taxonomic Class exactly as fig4 does -- but unlike fig4 (where Order maps
# 1:1 onto Class), trait_group and Class are independent classifications, so
# this is a genuine cross-tab: a given trait_group (e.g. "plant pathogen")
# can and does have significant hits in several different classes, and shows
# up as a separate bar in each of those classes' facets. Bar counts are
# computed per (trait_group, Class, direction) rather than per trait_group
# alone. Two thresholds, both driven by the single per-panel min_facet_hits
# value (min_facet_hits_date / min_facet_hits_pcoa1 below): a Class facet
# with too few total hits is dropped from the plot entirely, and within a
# surviving Class facet, a lifestyle bar with too few hits *in that class*
# folds into that class's own "Other" bar -- this is a per-facet cutoff, not
# a global rank cutoff, so a lifestyle that's rare in one class but common
# in another keeps its own bar where it's common. Panel a (2116 hits) needs
# both cuts to stay legible; panel b (77 hits) has room to show every class
# and cell individually, so its threshold is effectively off.
#
# Taxa whose FungalTraits Secondary_lifestyle is "plant_pathogen" (e.g. a
# wood saprotroph that's also known to cause disease) get their own
# "<primary lifestyle> (+plant pathogen)" bar rather than being folded into
# the plain primary-lifestyle bar -- checked to add only a handful of bars
# (4 in panel a, 2 in panel b), all within Class facets already shown, so it
# doesn't blow up the category list.
#
# Panel a: fungal taxa ~ collection date -- data/2024_fungi/
#          fungal_taxa_date_association.all_taxa.csv (same unbiased,
#          every-prevalence-filtered-taxon input as fig4 panel b). 2116
#          significant taxa span 39 trait groups; lifestyle-within-class
#          cells with <6 hits fold into that class's "Other" bar (see
#          min_facet_hits_date above).
# Panel b: fungal taxa ~ insect_PCoA1, NO date or lure covariate (the more
#          permissive test -- see fig3 panel c) -- data/
#          compare_insects_fungi_top3axes/fungal_insect_association_results.
#          PCoA1_only.no_lure.csv, same top-200 CoCA-loading candidate taxa
#          as fig4 panel c. Only 14 trait groups appear among the 77 hits,
#          so all are shown individually -- no "Other" bucket needed here.
#
# ASVs with no genus-level ID ("Unclassified genus") or a genus not present
# in FungalTraits ("No FungalTraits match") are dropped entirely -- this
# figure only shows taxa with an actual trophic-mode call, not
# taxonomic-resolution artifacts. (They were a sizeable share of hits --
# ~539 of 2116 for the date screen -- so panel subtitles report the
# FungalTraits-matched taxon count, not the full significant-hit count from
# fig4.)
#
# Standalone script -- edit min_facet_hits_date / min_facet_hits_pcoa1 below
# to change the class-drop / cell-fold cutoffs. Also writes the per-group/
# per-facet/per-direction counts underlying each panel's bars to
# data/presentation_items/fig8_trait_breakdown.*.csv.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("library/library.R")

out_fig_dir <- "figures/presentation_items"
out_data_dir <- "data/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

q_threshold <- 0.10
min_facet_hits_date <- 6   # panel a (2116 taxa): Class facets with 5 or fewer total hits are dropped, AND lifestyle-within-class cells with 5 or fewer hits fold into that class's "Other" bar -- panel a is dense enough to need trimming
min_facet_hits_pcoa1 <- 1   # panel b (77 taxa): nothing dropped or folded -- far fewer hits, so there's room to show every class and cell individually

pos_color <- "#0072B2"   # Okabe-Ito blue -- same "significant" blue used in fig3/fig4
neg_color <- "#D55E00"   # Okabe-Ito vermillion

panel_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    plot.tag = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey30"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
    strip.placement = "outside",
    panel.spacing.y = unit(0.15, "lines"),
    plot.title.position = "plot"   # title/subtitle span the full plot width (incl. the Class strip column), not just the panel -- needed now that the strip/label columns are wide enough to otherwise crowd them off the right edge
  )

## ---- Shared plotting logic ----------------------------------------------------
# res: data.frame with columns `group` (trait label, lifestyle), `facet`
# (taxonomic Class), and `t_stat`. `facet` is used as-is per row -- NA/blank
# becomes "Unclassified". Two thresholds, both driven by the single
# min_facet_hits value (so "how big does something have to be to get its
# own bar/facet" means the same thing at both levels):
#   - a Class facet with fewer than min_facet_hits total significant hits
#     (summed across all its lifestyle bars) is dropped from the plot
#     entirely, not folded into anything;
#   - within a surviving Class facet, a (lifestyle, Class) cell with fewer
#     than min_facet_hits hits is folded into that facet's own "Other" bar
#     -- this is per-facet, not a global rank cutoff, so a lifestyle that's
#     small in one class but common in another keeps its own bar in the
#     class where it's common and folds only where it's rare.
#
# Unlike fig4's Order-within-Class facets (a 1:1 nesting), `group` and
# `facet` are independent here, so bars are counted per (group, facet,
# direction) triple -- the same trait group can appear as separate bars in
# several different Class facets. Each (group, facet) pair gets a unique
# axis position (axis_id = "<facet rank>|<group>") so ggplot doesn't force
# one shared vertical order for a label that legitimately repeats across
# facets; scale_x_discrete() then strips the facet-rank prefix back off for
# display, leaving just the lifestyle label.
#
# Returns list(plot, counts) -- counts is the per-group/per-facet/per-
# direction table underlying the bars, written to CSV below.

make_diverging_breakdown <- function(res, x_label, title, subtitle, tag,
                                      pos_label, neg_label, legend_title,
                                      min_facet_hits = 1) {
  res <- res %>%
    mutate(
      group = ifelse(is.na(group) | group == "", "Unclassified", group),
      facet = ifelse(is.na(facet) | facet == "", "Unclassified", facet),
      direction = ifelse(t_stat > 0, pos_label, neg_label)
    )

  all_facet_totals <- res %>% count(facet, name = "total") %>% arrange(desc(total))
  dropped_facets <- all_facet_totals %>% filter(total < min_facet_hits)
  if (nrow(dropped_facets) > 0) {
    cat("  [", title, "] dropping", nrow(dropped_facets), "classes below min_facet_hits =", min_facet_hits,
        "(", sum(dropped_facets$total), "taxa total):", paste(dropped_facets$facet, collapse = ", "), "\n")
  }
  res <- res %>% filter(facet %in% all_facet_totals$facet[all_facet_totals$total >= min_facet_hits])

  # Per-(lifestyle, Class) cell fold: a lifestyle bar with fewer than
  # min_facet_hits hits *within this class* becomes "Other" for just these
  # rows -- the same lifestyle can stay a named bar in a different class
  # where it clears the threshold.
  cell_totals <- res %>% count(group, facet, name = "cell_total")
  small_cells <- cell_totals %>% filter(cell_total < min_facet_hits)
  if (nrow(small_cells) > 0) {
    cat("  [", title, "] folding", nrow(small_cells), "lifestyle-within-class cells below min_facet_hits =", min_facet_hits,
        "(", sum(small_cells$cell_total), "taxa total) into each class's Other\n")
  }
  res <- res %>%
    left_join(cell_totals, by = c("group", "facet")) %>%
    mutate(group = ifelse(cell_total < min_facet_hits, "Other", group)) %>%
    select(-cell_total)

  row_totals <- res %>% count(group, facet, name = "total")

  facet_totals <- row_totals %>% group_by(facet) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
  # real Classes ranked by total hits; "Unclassified" Class pinned to the
  # bottom regardless of hit count, same convention fig4 uses.
  real_facets <- facet_totals$facet[facet_totals$facet != "Unclassified"]
  facet_levels <- c(real_facets, intersect("Unclassified", facet_totals$facet))

  row_totals <- row_totals %>%
    mutate(facet = factor(facet, levels = facet_levels)) %>%
    arrange(facet, desc(total)) %>%
    mutate(axis_id = paste0(sprintf("%03d", as.integer(facet)), "|", group))

  axis_order <- row_totals$axis_id
  axis_labels <- setNames(row_totals$group, row_totals$axis_id)

  counts <- res %>%
    count(group, facet, direction, name = "n") %>%
    left_join(row_totals %>% select(group, facet, axis_id), by = c("group", "facet")) %>%
    mutate(
      facet = factor(facet, levels = facet_levels),
      axis_id = factor(axis_id, levels = rev(axis_order)),
      signed_n = ifelse(direction == pos_label, n, -n)
    )

  p <- ggplot(counts, aes(x = axis_id, y = signed_n, fill = direction)) +
    geom_col(color = "white", linewidth = 0.2, width = 0.75) +
    geom_text(aes(label = n, hjust = ifelse(signed_n >= 0, -0.3, 1.3)), size = 2.8) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = axis_labels) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
    scale_fill_manual(values = setNames(c(pos_color, neg_color), c(pos_label, neg_label)),
                       name = legend_title) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = x_label, title = title, subtitle = subtitle, tag = tag) +
    panel_theme +
    facet_grid(rows = vars(facet), scales = "free_y", space = "free_y", switch = "y")

  out_counts <- counts %>%
    transmute(Class = as.character(facet), lifestyle = as.character(group), direction, n) %>%
    arrange(desc(n))

  list(plot = p, counts = out_counts)
}

## ---- Fungal trait + taxonomy lookup (shared by panels a/b) --------------------
## Trait side: same join as fig1_taxonomic_breakdown.R panel d (ASV genus ->
## FungalTraits genus -> primary_lifestyle x growth form). Taxonomy side:
## Class, same source/formatting as fig4_association_taxonomic_breakdown.R.

fungal_traits_path <- "/Users/ericmorrison/repo/FungalTraits_Polme/Polme_FungalTraits_1.2_ver_16Dec_2020.csv"

fungal_traits <- read.csv(fungal_traits_path, fileEncoding = "latin1") %>%
  distinct(GENUS, .keep_all = TRUE) %>%   # a handful of genera appear twice; keep first
  transmute(
    Genus = GENUS,
    primary_lifestyle = ifelse(is.na(primary_lifestyle) | primary_lifestyle == "", "unspecified", primary_lifestyle),
    secondary_plant_pathogen = !is.na(Secondary_lifestyle) & Secondary_lifestyle == "plant_pathogen"
  )

taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv") %>% rename(taxon = X)

genus_map <- taxonomy %>%
  filter(Kingdom == "k__Fungi") %>%
  transmute(taxon, Genus = sub("^g__", "", Genus), Class = sub("^c__", "", Class)) %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", NA, Genus)) %>%
  left_join(fungal_traits, by = "Genus") %>%
  mutate(
    # is.na(primary_lifestyle) after the join means the genus had no match in
    # FungalTraits (as opposed to matching a genus whose lifestyle field was
    # blank, which was already recoded to "unspecified" above)
    #
    # Taxa whose *secondary* lifestyle is plant_pathogen get their own
    # "<primary lifestyle> (+plant pathogen)" bucket rather than being
    # silently merged into the plain primary-lifestyle bucket -- guarded
    # against primary_lifestyle == "plant_pathogen" so a taxon already
    # primary-plant-pathogen never gets a redundant "(+plant pathogen)" tag.
    trait_group = case_when(
      is.na(Genus) ~ "Unclassified genus",
      is.na(primary_lifestyle) ~ "No FungalTraits match",
      secondary_plant_pathogen & primary_lifestyle != "plant_pathogen" ~
        paste0(gsub("_", " ", primary_lifestyle), " (+plant pathogen)"),
      TRUE ~ gsub("_", " ", primary_lifestyle)
    )
  ) %>%
  select(taxon, trait_group, Class)

## ---- Panel a: fungal taxa ~ collection date, by trait group x Class -----------

res_a <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(genus_map, by = "taxon") %>%
  filter(!trait_group %in% c("Unclassified genus", "No FungalTraits match")) %>%
  mutate(group = trait_group, facet = Class)

subtitle_a <- paste0(nrow(res_a), " taxa with a FungalTraits match across ", n_distinct(res_a$group),
                      " lifestyles; classes/cells with <", min_facet_hits_date, " hits omitted/folded, q<0.10")
out_a <- make_diverging_breakdown(
  res_a, x_label = "Number of significant taxa",
  title = "Fungal trait vs. collection date", subtitle = subtitle_a, tag = "a",
  pos_label = "Later season (t>0)", neg_label = "Earlier season (t<0)",
  legend_title = "Direction", min_facet_hits = min_facet_hits_date
)
panel_a <- out_a$plot

## ---- Panel b: fungal taxa ~ insect PCoA1, by trait group x Class (no lure, no date) ----

res_b <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(genus_map, by = "taxon") %>%
  filter(!trait_group %in% c("Unclassified genus", "No FungalTraits match")) %>%
  mutate(group = trait_group, facet = Class)

subtitle_b <- paste0(nrow(res_b), " taxa with a FungalTraits match across ", n_distinct(res_b$group),
                      " lifestyles, all classes and cells shown, q<0.10")
out_b <- make_diverging_breakdown(
  res_b, x_label = "Number of significant taxa",   # min_facet_hits_pcoa1 = 1: no class or cell folded, everything shown individually
  title = "Fungal trait vs. insect PCoA1", subtitle = subtitle_b, tag = "b",
  pos_label = "PCoA1+ / later-season insect community",
  neg_label = "PCoA1- / bark-beetle-dominated insect community",
  legend_title = "Direction", min_facet_hits = min_facet_hits_pcoa1
)
panel_b <- out_b$plot

## ---- Write per-group/per-Class counts underlying each panel -------------------

write.csv(out_a$counts, file.path(out_data_dir, "fig8_trait_breakdown.fungal_date.csv"), row.names = FALSE)
write.csv(out_b$counts, file.path(out_data_dir, "fig8_trait_breakdown.fungal_insectPCoA1.csv"), row.names = FALSE)

## ---- Combine + save --------------------------------------------------------------
# Both panels are faceted (different gtable row structures depending on how
# many trait-group x Class rows each has), so gridExtra::arrangeGrob places
# each panel as an opaque grob rather than trying to align them via
# gtable::cbind (same approach fig4 uses for its faceted panels b/c).

combined <- gridExtra::arrangeGrob(panel_a, panel_b, ncol = 2)

ggsave(file.path(out_fig_dir, "fig8_association_trait_breakdown.png"), combined,
       width = 15, height = 16, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig8_association_trait_breakdown.pdf"), combined,
       width = 15, height = 16, bg = "white")

cat("\nDone. Wrote", file.path(out_fig_dir, "fig8_association_trait_breakdown.png"), "and .pdf\n")
cat("Wrote per-group/per-Class counts to", out_data_dir, "(fig8_trait_breakdown.fungal_date, fig8_trait_breakdown.fungal_insectPCoA1)\n")
