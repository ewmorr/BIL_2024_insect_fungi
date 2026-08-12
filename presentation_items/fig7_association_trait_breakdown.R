##############################################################################
# Presentation figure 7: trophic mode / morphology breakdown of significant
# (q<0.10) fungal taxa from the association screens shown in
# fig3_volcano_plots.R / fig4_association_taxonomic_breakdown.R. Where fig4
# asks "which taxonomic groups do the significant hits belong to", this
# figure asks the same question about ecological trait instead of taxonomy:
# ASV genus is joined to the FungalTraits database (Polme et al. 2020) by
# genus name, exactly as in fig1_taxonomic_breakdown.R panel d, and taxa are
# grouped by primary_lifestyle x growth form ("trophic mode : morphology").
# Same diverging-bar-chart / facet-by-coarser-rank technique as fig4 panels
# b/c: bars are trait_group (primary_lifestyle : growth_form), faceted by
# growth_form -- trait_group maps 1:1 onto growth_form the same way an Order
# maps 1:1 onto a Class, so this is a direct swap of the facet variable, not
# a new technique. Bar labels show only the primary_lifestyle half of the
# combination (the growth_form half is already given by the facet strip, so
# repeating it on every bar would be redundant).
#
# fig8_association_trait_breakdown.R is this figure's companion: same
# question, same FungalTraits join, but faceted by taxonomic Class instead
# of growth_form and with growth_form dropped from the label entirely. Class
# facets read as "how does this trophic mode play out taxonomically", while
# this figure's growth_form facets read as "how does this trophic mode play
# out morphologically" -- genuinely different cuts of the same data, kept as
# separate figures rather than picking one.
#
# Panel a: fungal taxa ~ collection date -- data/2024_fungi/
#          fungal_taxa_date_association.all_taxa.csv (same unbiased,
#          every-prevalence-filtered-taxon input as fig4 panel b). 2116
#          significant taxa span 39 trait groups; only the top
#          n_top_trait_groups_date (by hit count) are shown individually,
#          long tail folded into "Other".
# Panel b: fungal taxa ~ insect_PCoA1, NO date or lure covariate (the more
#          permissive test -- see fig3 panel c) -- data/
#          compare_insects_fungi_top3axes/fungal_insect_association_results.
#          PCoA1_only.no_lure.csv, same top-200 CoCA-loading candidate taxa
#          as fig4 panel c. Only 14 trait groups appear among the 77 hits,
#          so all are shown -- no "Other" bucket needed here.
#
# "Unclassified genus" (ASV genus not assigned) and "No FungalTraits match"
# (genus assigned but absent from the FungalTraits database) are real,
# fairly large categories among the significant hits -- shown as their own
# bars/pseudo-facets rather than dropped, same treatment fig4 gives
# "Unclassified" taxonomy. (fig8 drops these instead -- with Class facets,
# they split thinly across many real classes rather than forming one
# legible bucket; here, with growth_form facets, they stay contained to a
# single pseudo-facet each, so there's no reason to drop them.)
#
# Taxa whose FungalTraits Secondary_lifestyle is "plant_pathogen" (e.g. a
# wood saprotroph that's also known to cause disease) get "(+plant
# pathogen)" appended to their primary_lifestyle label, same rule as fig8 --
# these still facet by growth_form as normal (the tag only affects the
# primary_lifestyle half of trait_group, not growth_form).
#
# Standalone script -- edit n_top_trait_groups_date below to change how many
# named trait groups are shown in panel a before folding the rest into
# "Other". Also writes the per-group/per-direction counts underlying each
# panel's bars to data/presentation_items/fig7_trait_breakdown.*.csv.
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
n_top_trait_groups_date <- 25   # + "Other" = 26 categories in panel a (of 41 trait groups present); panel b shows all 14 trait groups present, no folding needed

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
    panel.spacing.y = unit(0.15, "lines")
  )

## ---- Shared plotting logic (fig4_association_taxonomic_breakdown.R's ---------
## make_diverging_breakdown, generalized: facet variable here is a coarser
## trait rank (growth_form) instead of taxonomic Class. `facet` is allowed
## to be NA per-row (a group with no coarser-rank value of its own, e.g.
## "Unclassified genus"/"No FungalTraits match") -- those groups get a
## same-named pseudo-facet instead of a blank strip, same idea as fig4's
## "Unclassified" pseudo-Class but generalized to any/all such groups rather
## than one hardcoded name.
##
## res: data.frame with columns `group` (trait label, "lifestyle : growth
## form") and `t_stat`, plus `facet` when facet_by_rank = TRUE. Collapses
## `group` to the top n_top by significant-hit count (+ "Other") when there
## are more than n_top groups present; otherwise shows all groups as-is. Bar
## labels drop the " : growth form" half of `group` (already shown by the
## facet strip) -- pseudo-groups ("Other"/"Unclassified genus"/"No
## FungalTraits match"), which contain no " : ", are left as-is.
##
## Returns list(plot, counts) -- counts is the per-group/per-direction table
## underlying the bars, written to CSV below.

make_diverging_breakdown <- function(res, n_top, x_label, title, subtitle, tag,
                                      pos_label, neg_label, legend_title,
                                      facet_by_rank = FALSE) {
  res <- res %>% mutate(group = ifelse(is.na(group) | group == "", "Unclassified", group))

  group_rank <- res %>% count(group, name = "total") %>% arrange(desc(total))
  if (nrow(group_rank) > n_top) {
    top_groups <- head(group_rank$group, n_top)
    res <- res %>% mutate(group = ifelse(group %in% top_groups, group, "Other"))
    group_levels <- c(top_groups, "Other")
  } else {
    group_levels <- group_rank$group
  }

  if (facet_by_rank) {
    # group is already post-fold here (real trait groups, or "Other"); assign
    # each kept group's facet -- "Other" mixes multiple real facet values, so
    # it gets its own same-named pseudo-facet; a group whose facet value was
    # NA to begin with (no coarser-rank value of its own) uses its own name
    # as its pseudo-facet instead.
    group_facet <- res %>% distinct(group, facet) %>%
      mutate(facet = case_when(group == "Other" ~ "Other", is.na(facet) ~ group, TRUE ~ facet)) %>%
      distinct(group, facet)
    group_totals <- res %>% count(group, name = "total")   # post-fold totals, so "Other" gets its combined count

    facet_rank <- group_facet %>% left_join(group_totals, by = "group") %>%
      group_by(facet) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
    # "Other" and any no-coarser-rank pseudo-facets ("Unclassified genus",
    # "No FungalTraits match") are not real trait ranks -- pin them to the
    # bottom regardless of hit count rather than letting them rank among
    # real primary_lifestyle facets.
    real_facets <- facet_rank$facet[!facet_rank$facet %in% c("Other", "Unclassified genus", "No FungalTraits match")]
    facet_levels <- c(real_facets, intersect(c("Other", "No FungalTraits match", "Unclassified genus"), facet_rank$facet))

    group_order <- group_facet %>% left_join(group_totals, by = "group") %>%
      mutate(facet = factor(facet, levels = facet_levels)) %>% arrange(facet, desc(total))
    group_levels <- group_order$group
  }

  counts <- res %>%
    mutate(direction = ifelse(t_stat > 0, pos_label, neg_label)) %>%
    count(group, direction)
  if (facet_by_rank) counts <- counts %>% left_join(group_facet, by = "group") %>%
    mutate(facet = factor(facet, levels = facet_levels))
  counts <- counts %>%
    mutate(
      group = factor(group, levels = rev(group_levels)),
      signed_n = ifelse(direction == pos_label, n, -n)
    )

  # Display-only label shortening: "<lifestyle> : <growth form>" -> "<lifestyle>"
  # (the growth_form facet strip already shows that half). Pseudo-groups
  # ("Other"/"Unclassified genus"/"No FungalTraits match") contain no " : "
  # and are left untouched. The group/facet values themselves (used for
  # faceting and the counts table) are unaffected -- this just de-clutters
  # the printed axis label.
  group_labels <- setNames(sub(" : .*$", "", levels(counts$group)), levels(counts$group))

  p <- ggplot(counts, aes(x = group, y = signed_n, fill = direction)) +
    geom_col(color = "white", linewidth = 0.2, width = 0.75) +
    geom_text(aes(label = n, hjust = ifelse(signed_n >= 0, -0.3, 1.3)), size = 2.8) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = group_labels) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
    scale_fill_manual(values = setNames(c(pos_color, neg_color), c(pos_label, neg_label)),
                       name = legend_title) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = x_label, title = title, subtitle = subtitle, tag = tag) +
    panel_theme

  if (facet_by_rank) p <- p + facet_grid(rows = vars(facet), scales = "free_y", space = "free_y", switch = "y")

  out_counts <- if (facet_by_rank) {
    counts %>% transmute(facet = as.character(facet), group = as.character(group), direction, n) %>% arrange(desc(n))
  } else {
    counts %>% transmute(group = as.character(group), direction, n) %>% arrange(desc(n))
  }
  list(plot = p, counts = out_counts)
}

## ---- Fungal trait lookup (shared by panels a/b) --------------------------------
## Same join as fig1_taxonomic_breakdown.R panel d: ASV genus -> FungalTraits
## genus -> primary_lifestyle x growth form.

fungal_traits_path <- "/Users/ericmorrison/repo/FungalTraits_Polme/Polme_FungalTraits_1.2_ver_16Dec_2020.csv"

fungal_traits <- read.csv(fungal_traits_path, fileEncoding = "latin1") %>%
  distinct(GENUS, .keep_all = TRUE) %>%   # a handful of genera appear twice; keep first
  transmute(
    Genus = GENUS,
    primary_lifestyle = ifelse(is.na(primary_lifestyle) | primary_lifestyle == "", "unspecified", primary_lifestyle),
    growth_form = ifelse(is.na(Growth_form_template) | Growth_form_template == "", "unspecified", Growth_form_template),
    secondary_plant_pathogen = !is.na(Secondary_lifestyle) & Secondary_lifestyle == "plant_pathogen"
  )

taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv") %>% rename(taxon = X)

genus_map <- taxonomy %>%
  filter(Kingdom == "k__Fungi") %>%
  transmute(taxon, Genus = sub("^g__", "", Genus)) %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", NA, Genus)) %>%
  left_join(fungal_traits, by = "Genus") %>%
  mutate(
    # is.na(primary_lifestyle) after the join means the genus had no match in
    # FungalTraits (as opposed to matching a genus whose lifestyle field was
    # blank, which was already recoded to "unspecified" above)
    #
    # Taxa whose *secondary* lifestyle is plant_pathogen get "(+plant
    # pathogen)" appended to the primary_lifestyle half only -- guarded
    # against primary_lifestyle == "plant_pathogen" so a taxon already
    # primary-plant-pathogen never gets a redundant tag. growth_form is
    # unaffected, so these taxa still facet normally.
    lifestyle_label = case_when(
      secondary_plant_pathogen & primary_lifestyle != "plant_pathogen" ~
        paste0(gsub("_", " ", primary_lifestyle), " (+plant pathogen)"),
      TRUE ~ gsub("_", " ", primary_lifestyle)
    ),
    trait_group = case_when(
      is.na(Genus) ~ "Unclassified genus",
      is.na(primary_lifestyle) ~ "No FungalTraits match",
      TRUE ~ paste0(lifestyle_label, " : ", gsub("_", " ", growth_form))
    ),
    growth_form_facet = ifelse(is.na(primary_lifestyle), NA, gsub("_", " ", growth_form))
  ) %>%
  select(taxon, trait_group, growth_form_facet)

## ---- Panel a: fungal taxa ~ collection date, by trait group -------------------

res_a <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(genus_map, by = "taxon") %>%
  mutate(group = trait_group, facet = growth_form_facet)

subtitle_a <- paste0(nrow(res_a), " taxa across ", n_distinct(res_a$group), " trait groups (top ",
                      n_top_trait_groups_date, " shown), q<0.10")
out_a <- make_diverging_breakdown(
  res_a, n_top = n_top_trait_groups_date, x_label = "Number of significant taxa",
  title = "Fungal trait vs. collection date", subtitle = subtitle_a, tag = "a",
  pos_label = "Later season (t>0)", neg_label = "Earlier season (t<0)",
  legend_title = "Direction", facet_by_rank = TRUE
)
panel_a <- out_a$plot

## ---- Panel b: fungal taxa ~ insect PCoA1, by trait group (no lure, no date) ----

res_b <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(genus_map, by = "taxon") %>%
  mutate(group = trait_group, facet = growth_form_facet)

subtitle_b <- paste0(nrow(res_b), " taxa across ", n_distinct(res_b$group), " trait groups, q<0.10")
out_b <- make_diverging_breakdown(
  res_b, n_top = 999, x_label = "Number of significant taxa",   # 999: show all trait groups present, no folding
  title = "Fungal trait vs. insect PCoA1", subtitle = subtitle_b, tag = "b",
  pos_label = "PCoA1+ / later-season insect community",
  neg_label = "PCoA1- / bark-beetle-dominated insect community",
  legend_title = "Direction", facet_by_rank = TRUE
)
panel_b <- out_b$plot

## ---- Write per-group counts underlying each panel -------------------------------

write.csv(out_a$counts, file.path(out_data_dir, "fig7_trait_breakdown.fungal_date.csv"), row.names = FALSE)
write.csv(out_b$counts, file.path(out_data_dir, "fig7_trait_breakdown.fungal_insectPCoA1.csv"), row.names = FALSE)

## ---- Combine + save --------------------------------------------------------------
# Both panels are faceted (different gtable row structures depending on how
# many trait groups/facets each has), so gridExtra::arrangeGrob places each
# panel as an opaque grob rather than trying to align them via gtable::cbind
# (same approach fig4 uses for its faceted panels b/c).

combined <- gridExtra::arrangeGrob(panel_a, panel_b, ncol = 2)

ggsave(file.path(out_fig_dir, "fig7_association_trait_breakdown.png"), combined,
       width = 15, height = 10, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig7_association_trait_breakdown.pdf"), combined,
       width = 15, height = 10, bg = "white")

cat("\nDone. Wrote", file.path(out_fig_dir, "fig7_association_trait_breakdown.png"), "and .pdf\n")
cat("Wrote per-group counts to", out_data_dir, "(fig7_trait_breakdown.fungal_date, fig7_trait_breakdown.fungal_insectPCoA1)\n")
