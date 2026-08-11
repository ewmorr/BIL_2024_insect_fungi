##############################################################################
# Presentation figure 4: taxonomic breakdown of significant (q<0.10) taxa
# from the three per-taxon association screens shown as volcano plots in
# fig3_volcano_plots.R. Where fig3 shows every tested taxon as one point,
# this figure asks "of the significant hits, which taxonomic groups do they
# belong to, and which direction do they point?" -- diverging bar charts,
# taxonomic group on y, count of significant taxa on x (positive t-statistic
# to the right, negative to the left).
#
# Panel a: insect taxa ~ collection date, by Subfamily -- data/2024_insect_
#          data/insect_taxa_date_association.csv, taxonomy from insect_
#          species_tab.csv. Family collapses to 2 groups (98% Curculionidae),
#          so Subfamily is used instead -- it's also the level interpretation.
#          md uses to distinguish the Scolytinae (early-season bark/ambrosia
#          beetles) vs. Corticariinae/Cossoninae (later-season) pattern. Only
#          4 subfamilies appear among the 18 significant hits, so all are
#          shown (no "Other" bucket needed).
# Panel b: fungal taxa ~ collection date, by Order, faceted by Class --
#          data/2024_fungi/fungal_taxa_date_association.all_taxa.csv
#          (unbiased, every prevalence-filtered taxon, no pre-selection),
#          taxonomy from ASVs_taxonomy.tsv. 2116 significant taxa span 103
#          orders, so only the top n_top_fungal_date_orders (by hit count)
#          are shown individually and the long tail is folded into "Other".
#          Order (rather than Family, used in an earlier version of this
#          figure) was chosen to foreground the broad order-level trends
#          found in compare_insect_fungi/fungi_date_order_breakdown.R --
#          e.g. entire orders (Mycosphaerellales, Tremellales, Dothideales)
#          sit almost entirely on the later-season side, others
#          (Diaporthales, Helotiales) almost entirely earlier-season -- a
#          pattern the family-level view split across too many individual
#          bars to see clearly. Order-to-Class is a clean 1:1 mapping for
#          every real, individually-named order, so grouping the order bars
#          into Class facets (same facet_grid technique as fig5's
#          family-within-order panels) reads as one more level up the same
#          hierarchy. The "Other" fold and the "Unclassified" order both mix
#          taxa from several different classes, so each gets its own
#          pseudo-facet (labeled "Other"/"Unclassified") rather than a false
#          single-class assignment. Family-level detail within these orders
#          is the dedicated supplemental figure,
#          fig5_family_breakdown_by_order.R.
# Panel c: fungal taxa ~ insect_PCoA1, by Order, faceted by Class, NO date or
#          lure covariate (the more permissive test -- see fig3 panel c) --
#          data/compare_insects_fungi_top3axes/fungal_insect_association_
#          results.PCoA1_only.no_lure.csv, same top-200 CoCA-loading
#          candidate taxa as fig3 panel c. Only 31 orders appear among the
#          77 hits (vs. 259 families), so all are shown -- no "Other" bucket
#          needed here (though "Unclassified" order still gets its own
#          pseudo-facet, same reason as panel b). Direction here is the
#          taxon's association with insect_PCoA1 itself, not date directly
#          -- PCoA1 is ~78% collinear with date (interpretation.md), positive
#          = later-season/fungivorous-beetle-dominated end, negative =
#          earlier-season/bark-beetle-dominated end.
#
# Standalone script -- edit n_top_fungal_date_orders below to change how many
# named orders are shown in panel b before folding the rest into "Other".
# Also writes the per-group/per-direction counts underlying each panel's
# bars to data/presentation_items/fig4_taxonomic_breakdown.*.csv.
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
n_top_fungal_date_orders <- 15   # + "Other" = 16 categories in panel b; panel c shows all 31 orders present, no folding needed

pos_color <- "#0072B2"   # Okabe-Ito blue -- same "significant" blue used in fig3
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

## ---- Shared plotting logic ---------------------------------------------------
# res: data.frame with columns `group` (taxonomic label) and `t_stat`, plus
# `Class` when facet_by_class = TRUE. Collapses `group` to the top n_top by
# significant-hit count (+ "Other") when there are more than n_top groups
# present; otherwise shows all groups as-is.
#
# facet_by_class facets the bars by their subsuming Class (same facet_grid
# technique fig5 uses to facet families by order). Order-to-Class is a clean
# 1:1 mapping for every real, individually-named order in this dataset, so
# each kept group's Class is just its (unique) taxonomy value -- except the
# "Other" fold and "Unclassified" group, which each mix multiple classes and
# so get their own same-named pseudo-facet instead of a misleading single
# Class.
#
# Returns list(plot, counts) -- counts is the per-group/per-direction table
# underlying the bars, written to CSV below.

make_diverging_breakdown <- function(res, n_top, x_label, title, subtitle, tag,
                                      pos_label, neg_label, legend_title,
                                      facet_by_class = FALSE) {
  res <- res %>% mutate(group = ifelse(is.na(group) | group == "", "Unclassified", group))
  if (facet_by_class) {
    res <- res %>% mutate(Class = ifelse(is.na(Class) | Class == "", "Unclassified", Class))
  }

  group_rank <- res %>% count(group, name = "total") %>% arrange(desc(total))
  if (nrow(group_rank) > n_top) {
    top_groups <- head(group_rank$group, n_top)
    res <- res %>% mutate(group = ifelse(group %in% top_groups, group, "Other"))
    group_levels <- c(top_groups, "Other")
  } else {
    group_levels <- group_rank$group
  }

  if (facet_by_class) {
    # group is already post-fold here (real orders, or "Other"/"Unclassified");
    # collapse each to a single Class -- "Other"/"Unclassified" mix multiple
    # raw classes, so they're assigned their own same-named pseudo-Class here
    # rather than an arbitrary/misleading real one.
    group_class <- res %>% distinct(group, Class) %>%
      mutate(Class = case_when(group == "Other" ~ "Other", group == "Unclassified" ~ "Unclassified", TRUE ~ Class)) %>%
      distinct(group, Class)
    group_totals <- res %>% count(group, name = "total")   # post-fold totals, so "Other" gets its combined count

    facet_rank <- group_class %>% left_join(group_totals, by = "group") %>%
      group_by(Class) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
    # "Unclassified"/"Other" are pseudo-classes (each mixing several real
    # classes), not taxonomic ranks -- pin them to the bottom regardless of
    # their hit count rather than letting them rank among real classes.
    real_classes <- setdiff(facet_rank$Class, c("Unclassified", "Other"))
    facet_levels <- c(real_classes, intersect(c("Unclassified", "Other"), facet_rank$Class))

    group_order <- group_class %>% left_join(group_totals, by = "group") %>%
      mutate(Class = factor(Class, levels = facet_levels)) %>% arrange(Class, desc(total))
    group_levels <- group_order$group
  }

  counts <- res %>%
    mutate(direction = ifelse(t_stat > 0, pos_label, neg_label)) %>%
    count(group, direction)
  if (facet_by_class) counts <- counts %>% left_join(group_class, by = "group") %>%
    mutate(Class = factor(Class, levels = facet_levels))
  counts <- counts %>%
    mutate(
      group = factor(group, levels = rev(group_levels)),
      signed_n = ifelse(direction == pos_label, n, -n)
    )

  # Display-only label shortening: "<Class>_ord_Incertae_sedis" -> "Incertae_sedis".
  # The group/Class values themselves (used for faceting and the counts
  # table) are untouched -- an Incertae-sedis order still facets under its
  # own class, this just de-clutters the printed axis label.
  group_labels <- setNames(sub("^.*_ord_Incertae_sedis$", "Incertae_sedis", levels(counts$group)), levels(counts$group))

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

  if (facet_by_class) p <- p + facet_grid(rows = vars(Class), scales = "free_y", space = "free_y", switch = "y")

  out_counts <- if (facet_by_class) {
    counts %>% transmute(Class = as.character(Class), group = as.character(group), direction, n) %>% arrange(desc(n))
  } else {
    counts %>% transmute(group = as.character(group), direction, n) %>% arrange(desc(n))
  }
  list(plot = p, counts = out_counts)
}

## ---- Panel a: insect taxa ~ collection date, by Subfamily --------------------

insect_taxonomy <- read.csv("data/2024_insect_data/insect_species_tab.csv") %>%
  select(Family, Subfamily, Finest.ID) %>%
  distinct()

res_a <- read.csv("data/2024_insect_data/insect_taxa_date_association.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(insect_taxonomy, by = c("taxon" = "Finest.ID")) %>%
  mutate(group = ifelse(is.na(Subfamily) | Subfamily == "", Family, Subfamily))

subtitle_a <- paste0(nrow(res_a), " significant taxa (q<0.10) across ", n_distinct(res_a$group), " subfamilies")
out_a <- make_diverging_breakdown(
  res_a, n_top = 20, x_label = "Number of significant taxa",
  title = "Insect taxa vs. collection date", subtitle = subtitle_a, tag = "a",
  pos_label = "Later season (t>0)", neg_label = "Earlier season (t<0)",
  legend_title = "Direction"
)
panel_a <- out_a$plot

## ---- Fungal taxonomy lookup (shared by panels b/c) ----------------------------

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1) %>%
  tibble::rownames_to_column("taxon") %>%
  transmute(taxon, Order = sub("^o__", "", Order), Class = sub("^c__", "", Class))

## ---- Panel b: fungal taxa ~ collection date, by Order, faceted by Class -----

res_b <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(fungal_taxonomy, by = "taxon") %>%
  mutate(group = Order)

subtitle_b <- paste0(nrow(res_b), " taxa across ", n_distinct(res_b$group), " orders (top ",
                      n_top_fungal_date_orders, " shown), q<0.10")
out_b <- make_diverging_breakdown(
  res_b, n_top = n_top_fungal_date_orders, x_label = "Number of significant taxa",
  title = "Fungal taxa vs. collection date", subtitle = subtitle_b, tag = "b",
  pos_label = "Later season (t>0)", neg_label = "Earlier season (t<0)",
  legend_title = "Direction", facet_by_class = TRUE
)
panel_b <- out_b$plot

## ---- Panel c: fungal taxa ~ insect PCoA1, by Order, faceted by Class (no lure, no date) --

res_c <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(fungal_taxonomy, by = "taxon") %>%
  mutate(group = Order)

subtitle_c <- paste0(nrow(res_c), " taxa across ", n_distinct(res_c$group), " orders, q<0.10")
out_c <- make_diverging_breakdown(
  res_c, n_top = 999, x_label = "Number of significant taxa",   # 999: show all orders present (31), no folding
  title = "Fungal taxa vs. insect PCoA1", subtitle = subtitle_c, tag = "c",
  pos_label = "PCoA1+ / later-season insect community",
  neg_label = "PCoA1- / bark-beetle-dominated insect community",
  legend_title = "Direction", facet_by_class = TRUE
)
panel_c <- out_c$plot

## ---- Write per-group counts underlying each panel ------------------------------

write.csv(out_a$counts, file.path(out_data_dir, "fig4_taxonomic_breakdown.insect_date.by_subfamily.csv"), row.names = FALSE)
write.csv(out_b$counts, file.path(out_data_dir, "fig4_taxonomic_breakdown.fungal_date.by_order.csv"), row.names = FALSE)
write.csv(out_c$counts, file.path(out_data_dir, "fig4_taxonomic_breakdown.fungal_insectPCoA1.by_order.csv"), row.names = FALSE)

## ---- Combine + save -----------------------------------------------------------
# Panel a has no facets, while b/c are now faceted by Class (different
# gtable row structures), so gtable::cbind -- used when this figure's panels
# shared one flat row structure -- no longer applies (same issue fig5 hit
# combining its faceted panels). gridExtra::arrangeGrob places each panel as
# an opaque grob instead.

combined <- gridExtra::arrangeGrob(panel_a, panel_b, panel_c, ncol = 3)

# panel c (31 orders, unfolded) is taller than a/b, so the combined figure is
# taller than the family-level version of this figure was.
ggsave(file.path(out_fig_dir, "fig4_association_taxonomic_breakdown.png"), combined,
       width = 22, height = 10, dpi = 300, bg = "white")
ggsave(file.path(out_fig_dir, "fig4_association_taxonomic_breakdown.pdf"), combined,
       width = 22, height = 10, bg = "white")

cat("\nDone. Wrote", file.path(out_fig_dir, "fig4_association_taxonomic_breakdown.png"), "and .pdf\n")
cat("Wrote per-group counts to", out_data_dir, "(insect_date.by_subfamily, fungal_date.by_order, fungal_insectPCoA1.by_order)\n")
