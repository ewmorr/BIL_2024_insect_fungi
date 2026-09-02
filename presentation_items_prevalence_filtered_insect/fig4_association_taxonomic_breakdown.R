##############################################################################
# Presentation figure 4 -- Prevalence-Filtered Insect Table Variant (Lineage B)
#
# Lineage-B counterpart of presentation_items/fig4_association_taxonomic_
# breakdown.R. Where fig3 shows every tested taxon as one volcano point, this
# figure asks "of the significant hits, which taxonomic groups do they
# belong to, and which direction do they point?" -- diverging bar charts,
# taxonomic group on y, count of significant taxa on x.
#
# REVISION (2026-09-01, second pass): the first pass mirrored fig3's 2x3
# grid (rows = linear/quadratic term, in separate panels). Per user
# feedback, the linear and quadratic terms are now shown TOGETHER on one set
# of taxonomic-group axes per test -- each group gets up to 4 dodged
# (side-by-side) bars instead of two separate panels:
#   - Later season (linear, t>0)   / Earlier season (linear, t<0)
#   - Dip (quadratic, t>0)         / Hump (quadratic, t<0)
# using 4 Okabe-Ito colors (the original blue/vermillion pair for the linear
# term, plus sky-blue/orange for the quadratic term). This only applies
# where linear and quadratic are two terms of the SAME test (insect~date,
# fungal~date) -- panels c/d (fungal~insect_PCoA1, fungal~insect_PCoA2) are
# different predictor axes, not a linear/quadratic pair of one test, so they
# stay single-term, unchanged from the first pass.
#
# The result is a 2x2 grid instead of 2x3:
#   a  insect taxa ~ date (linear + quadratic combined), by Subfamily.
#   b  fungal taxa ~ date (linear + quadratic combined), by Order/Class.
#   c  fungal taxa ~ insect_PCoA1, factor(date)-adjusted, by Order/Class
#      (unchanged from the first pass).
#   d  fungal taxa ~ insect_PCoA2, factor(date)-adjusted, by Order/Class
#      (unchanged from the first pass; PCoA2 unstable across date-adjustment
#      strategies, caveat retained in the subtitle).
#
# Per user request, the top-n fold cutoff for a/b is now based on the
# COMBINED linear + quadratic hit count per group (a taxon significant in
# both terms counts twice toward its group's total), not ranked separately
# per term the way the two-panel version implicitly did.
#
# REVISION (2026-09-02, third pass): adopts the linear-term-takes-priority
# "peak timing" reading introduced in fig3_volcano_plots.R's 2026-09-02
# revision (see project_organization.md's "Peak-timing convention"). A taxon
# now gets exactly ONE of 4 mutually-exclusive categories instead of up to 2
# independent hit-instances:
#   - Peaks late season    -- significant linear term, t>0 (quadratic or not)
#   - Peaks early season   -- significant linear term, t<0 (quadratic or not)
#   - Peaks mid-season     -- NO significant linear term, significant
#                              quadratic hump (t_quad<0)
#   - Bimodal (early+late) -- NO significant linear term, significant
#                              quadratic dip (t_quad>0)
# Panels a/b's dodged bars, colors, and top-n fold ranking all switch to this
# scheme (a taxon now contributes to exactly one bar, so group totals are
# plain non-overlapping taxon counts, not hit-instance counts). Panels c/d
# (PCoA1/PCoA2, single-term) are unaffected -- different predictor axes, not
# part of this reading. fig5 (family-within-order breakdown of panel b) is
# revised the same way in the same pass; fig7/fig8 (trait breakdowns) are
# NOT yet revised -- still pending.
#
# Data sources unchanged from the first pass: insect_taxa_date_association.csv
# (153 taxa, Family/Subfamily/Order baked in), fungal_taxa_date_association.
# all_taxa.csv (lineage-invariant), fungal_insect_association_results.
# PCoA{1,2}_plus_factordate.csv.
#
# Standalone script -- edit n_top_fungal_date_orders / n_top_insect_date_
# subfamilies below to change the fold thresholds. Also writes the per-group/
# per-category counts underlying each panel's bars to data/presentation_
# items_prevalence_filtered_insect/fig4_taxonomic_breakdown.*.csv.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("library/library.R")

out_fig_dir <- "figures/presentation_items_prevalence_filtered_insect"
out_data_dir <- "data/presentation_items_prevalence_filtered_insect"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

lb_dir <- "data/compare_insects_fungi_top3axes_prevalence_filtered_insect"

q_threshold <- 0.10
n_top_fungal_date_orders <- 15        # matches Lineage A's fig4 panel b threshold
n_top_insect_date_subfamilies <- 15   # + "Other"; new for Lineage B -- Lineage A never needed folding here (only 4 subfamilies)

# Peak-timing colors (2026-09-02 revision): a taxon now gets exactly ONE of
# these 4 categories (linear term takes priority whenever it's significant --
# see header and project_organization.md's "Peak-timing convention"), not up
# to 2 independent hit-instances as in the first two passes. Colors carried
# over unchanged from the old linear/quadratic pairs for visual continuity.
lin_pos_color <- "#0072B2"    # Peaks late season (significant linear, t>0)
lin_neg_color <- "#D55E00"    # Peaks early season (significant linear, t<0)
quad_pos_color <- "#56B4E9"   # Bimodal, early+late (no sig linear; quadratic dip, t_quad>0)
quad_neg_color <- "#E69F00"   # Peaks mid-season (no sig linear; quadratic hump, t_quad<0)

# REVISION (2026-09-02, fourth pass): with 4 mutually-exclusive categories
# (not a signed t-statistic), a diverging left/right layout no longer maps
# to anything -- all 4 bars now run the same direction (0 -> positive).
# cat_levels controls the within-group dodge order; empirically, position_
# dodge2 + coord_flip renders the FIRST level at the BOTTOM of each group's
# band and the LAST level at the TOP (i.e. top-to-bottom is the REVERSE of
# this vector) -- so cat_levels is written back-to-front here to get the
# requested top-to-bottom reading order: early, bimodal, mid-season, late.
cat_levels <- c("Peaks late season", "Peaks mid-season",
                 "Bimodal (early + late)", "Peaks early season")
cat_colors <- c("Peaks late season" = lin_pos_color, "Peaks early season" = lin_neg_color,
                 "Bimodal (early + late)" = quad_pos_color, "Peaks mid-season" = quad_neg_color)[cat_levels]

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

## ---- Panel builder: peak-timing categories, dodged bars -----------------------
# res needs columns `group`, `t_stat`, `q_value`, `t_stat_quad`, `q_value_quad`,
# plus `Class` when facet_by_class = TRUE. REVISED 2026-09-02: a taxon now
# gets exactly ONE of the 4 cat_levels (linear term takes priority whenever
# it's significant -- see header), not up to 2 independent hit-instances as
# in the first two passes, so group ranking for the top-n fold and the bar
# counts are now plain non-overlapping taxon counts. Every kept group still
# gets all 4 cat levels (dodge slots) even when some are zero, via
# tidyr::complete(), so dodge geometry stays aligned across rows/facets.
# Returns list(plot, counts).

make_combined_breakdown <- function(res, n_top, x_label, title, subtitle, tag, facet_by_class = FALSE) {
  res <- res %>%
    mutate(group = ifelse(is.na(group) | group == "", "Unclassified", group),
           sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold,
           # Literal strings, NOT cat_levels[i] -- indexing into cat_levels
           # by position is fragile against reordering it for display (bit
           # us once already: the fourth-pass reorder below silently swapped
           # early/mid-season here until this fix).
           cat = case_when(
             sig_lin & t_stat > 0                  ~ "Peaks late season",
             sig_lin & t_stat < 0                  ~ "Peaks early season",
             !sig_lin & sig_quad & t_stat_quad > 0 ~ "Bimodal (early + late)",
             !sig_lin & sig_quad & t_stat_quad < 0 ~ "Peaks mid-season",
             TRUE ~ NA_character_
           )) %>%
    filter(!is.na(cat))
  if (facet_by_class) res <- res %>% mutate(Class = ifelse(is.na(Class) | Class == "", "Unclassified", Class))

  group_rank <- res %>% count(group, name = "total") %>% arrange(desc(total))

  if (nrow(group_rank) > n_top) {
    top_groups <- head(group_rank$group, n_top)
    res <- res %>% mutate(group = ifelse(group %in% top_groups, group, "Other"))
    group_levels <- c(top_groups, "Other")
  } else {
    group_levels <- group_rank$group
  }

  if (facet_by_class) {
    group_class <- res %>% distinct(group, Class) %>%
      mutate(Class = case_when(group == "Other" ~ "Other", group == "Unclassified" ~ "Unclassified", TRUE ~ Class)) %>%
      distinct(group, Class)
    group_totals <- res %>% count(group, name = "total")

    facet_rank <- group_class %>% left_join(group_totals, by = "group") %>%
      group_by(Class) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
    real_classes <- setdiff(facet_rank$Class, c("Unclassified", "Other"))
    facet_levels <- c(real_classes, intersect(c("Unclassified", "Other"), facet_rank$Class))

    group_order <- group_class %>% left_join(group_totals, by = "group") %>%
      mutate(Class = factor(Class, levels = facet_levels)) %>% arrange(Class, desc(total))
    group_levels <- group_order$group
  }

  counts <- res %>% count(group, cat, name = "n") %>%
    complete(group = group_levels, cat = cat_levels, fill = list(n = 0))
  # left_join must happen BEFORE `group` becomes a factor -- dplyr::left_join
  # silently coerces a factor join column back to character (dropping its
  # levels) when the other side is character, which upstream would have left
  # levels(counts$group) NULL and every scale_x_discrete label as "NA".
  if (facet_by_class) counts <- counts %>% left_join(group_class, by = "group") %>%
    mutate(Class = factor(Class, levels = facet_levels))
  counts <- counts %>%
    mutate(
      group = factor(group, levels = rev(group_levels)),
      cat = factor(cat, levels = cat_levels)
    )

  group_labels <- setNames(sub("^.*_ord_Incertae_sedis$", "Incertae_sedis", levels(counts$group)), levels(counts$group))

  # All 4 bars now run the same direction (0 -> positive) -- no more signed_n
  # diverging layout, since the 4 peak-timing categories aren't a signed
  # t-statistic (see cat_levels comment above).
  p <- ggplot(counts, aes(x = group, y = n, fill = cat)) +
    geom_col(color = "white", linewidth = 0.15, width = 0.8,
             position = position_dodge2(width = 0.8, padding = 0.1)) +
    geom_text(aes(label = ifelse(n == 0, "", n)), hjust = -0.3,
              position = position_dodge2(width = 0.8, padding = 0.1), size = 2.2) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = group_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    scale_fill_manual(values = cat_colors, name = "Peak timing", drop = FALSE) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = x_label, title = title, subtitle = subtitle, tag = tag) +
    panel_theme

  if (facet_by_class) p <- p + facet_grid(rows = vars(Class), scales = "free_y", space = "free_y", switch = "y")

  out_counts <- (if (facet_by_class) {
    counts %>% transmute(Class = as.character(Class), group = as.character(group), cat = as.character(cat), n)
  } else {
    counts %>% transmute(group = as.character(group), cat = as.character(cat), n)
  }) %>% filter(n > 0) %>% arrange(desc(n))

  list(plot = p, counts = out_counts)
}

## ---- Panel builder: SINGLE term (unchanged from the first pass) --------------
# Used for panels c/d (fungal ~ insect PCoA1/PCoA2) -- these are different
# predictor axes, not a linear/quadratic pair of one test, so there is no
# quadratic bar to dodge against.

make_single_term_breakdown <- function(res, n_top, x_label, title, subtitle, tag,
                                        pos_label, neg_label, legend_title,
                                        facet_by_class = FALSE) {
  res <- res %>% mutate(group = ifelse(is.na(group) | group == "", "Unclassified", group))
  if (facet_by_class) res <- res %>% mutate(Class = ifelse(is.na(Class) | Class == "", "Unclassified", Class))

  group_rank <- res %>% count(group, name = "total") %>% arrange(desc(total))
  if (nrow(group_rank) > n_top) {
    top_groups <- head(group_rank$group, n_top)
    res <- res %>% mutate(group = ifelse(group %in% top_groups, group, "Other"))
    group_levels <- c(top_groups, "Other")
  } else {
    group_levels <- group_rank$group
  }

  if (facet_by_class) {
    group_class <- res %>% distinct(group, Class) %>%
      mutate(Class = case_when(group == "Other" ~ "Other", group == "Unclassified" ~ "Unclassified", TRUE ~ Class)) %>%
      distinct(group, Class)
    group_totals <- res %>% count(group, name = "total")

    facet_rank <- group_class %>% left_join(group_totals, by = "group") %>%
      group_by(Class) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
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

  group_labels <- setNames(sub("^.*_ord_Incertae_sedis$", "Incertae_sedis", levels(counts$group)), levels(counts$group))

  p <- ggplot(counts, aes(x = group, y = signed_n, fill = direction)) +
    geom_col(color = "white", linewidth = 0.2, width = 0.75) +
    geom_text(aes(label = n, hjust = ifelse(signed_n >= 0, -0.3, 1.3)), size = 2.8) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = group_labels) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
    scale_fill_manual(values = setNames(c(lin_pos_color, lin_neg_color), c(pos_label, neg_label)),
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

## ---- Shared source data --------------------------------------------------------

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv") %>%
  rename(taxon = X) %>%
  transmute(taxon, Order = sub("^o__", "", Order), Class = sub("^c__", "", Class))

## ---- Panel a: insect taxa ~ date, LINEAR + QUADRATIC combined, by Subfamily ---
# insect_taxa_date_association.csv already carries Family/Subfamily/Order --
# no join needed. Subfamily uses two upstream fallback conventions: literal
# "NA" (unresolved, Family also "NA") and "NoSubfamily_<x>" (Family resolved,
# no subfamily rank) -- both collapse to Family here.

res_insect <- read.csv(file.path(lb_dir, "insect_taxa_date_association.csv")) %>%
  mutate(
    no_subfamily = is.na(Subfamily) | Subfamily %in% c("", "NA") | grepl("^NoSubfamily_", Subfamily),
    fallback_group = ifelse(is.na(Family) | Family %in% c("", "NA"), NA, Family),
    group = ifelse(no_subfamily, fallback_group, Subfamily)
  )

res_a <- res_insect %>%
  mutate(sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold,
         peak = case_when(
           sig_lin & t_stat > 0                  ~ "late",
           sig_lin & t_stat < 0                  ~ "early",
           !sig_lin & sig_quad & t_stat_quad > 0 ~ "bimodal",
           !sig_lin & sig_quad & t_stat_quad < 0 ~ "mid",
           TRUE ~ NA_character_
         )) %>%
  filter(!is.na(peak))
n_a_late <- sum(res_a$peak == "late"); n_a_early <- sum(res_a$peak == "early")
n_a_mid <- sum(res_a$peak == "mid"); n_a_bimodal <- sum(res_a$peak == "bimodal")
n_a_groups <- n_distinct(ifelse(is.na(res_a$group) | res_a$group == "", "Unclassified", res_a$group))
subtitle_a <- paste0(nrow(res_a), " taxa with a significant peak-timing call (", n_a_late, " late + ",
                      n_a_early, " early + ", n_a_mid, " mid-season + ", n_a_bimodal, " bimodal) across ",
                      n_a_groups, " subfamilies (top ", n_top_insect_date_subfamilies, " shown)")
out_a <- make_combined_breakdown(
  res_a, n_top = n_top_insect_date_subfamilies, x_label = "Number of significant taxa",
  title = "Insect taxa vs. date", subtitle = subtitle_a, tag = "a"
)
panel_a <- out_a$plot

## ---- Panel b: fungal taxa ~ date, LINEAR + QUADRATIC combined, by Order/Class ---

res_fungal_date <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  left_join(fungal_taxonomy, by = "taxon") %>%
  mutate(group = Order)

res_b <- res_fungal_date %>%
  mutate(sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold,
         peak = case_when(
           sig_lin & t_stat > 0                  ~ "late",
           sig_lin & t_stat < 0                  ~ "early",
           !sig_lin & sig_quad & t_stat_quad > 0 ~ "bimodal",
           !sig_lin & sig_quad & t_stat_quad < 0 ~ "mid",
           TRUE ~ NA_character_
         )) %>%
  filter(!is.na(peak))
n_b_late <- sum(res_b$peak == "late"); n_b_early <- sum(res_b$peak == "early")
n_b_mid <- sum(res_b$peak == "mid"); n_b_bimodal <- sum(res_b$peak == "bimodal")
n_b_groups <- n_distinct(ifelse(is.na(res_b$group) | res_b$group == "", "Unclassified", res_b$group))
subtitle_b <- paste0(nrow(res_b), " taxa with a significant peak-timing call (", n_b_late, " late + ",
                      n_b_early, " early + ", n_b_mid, " mid-season + ", n_b_bimodal, " bimodal) across ",
                      n_b_groups, " orders (top ", n_top_fungal_date_orders, " shown)")
out_b <- make_combined_breakdown(
  res_b, n_top = n_top_fungal_date_orders, x_label = "Number of significant taxa",
  title = "Fungal taxa vs. date", subtitle = subtitle_b, tag = "b", facet_by_class = TRUE
)
panel_b <- out_b$plot

## ---- Panel c: fungal taxa ~ insect PCoA1, factor(date)-adjusted, by Order/Class ----
# Unchanged from the first pass -- single term, no quadratic variant exists
# for this test.

res_c_all <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA1_plus_factordate.csv")) %>%
  left_join(fungal_taxonomy, by = "taxon")
res_c <- res_c_all %>%
  filter(q_value < q_threshold) %>%
  mutate(group = Order)

subtitle_c <- paste0(nrow(res_c), "/", nrow(res_c_all), " CoCA candidates across ", n_distinct(res_c$group),
                      " orders, q<0.10, factor(date)-adjusted")
out_c <- make_single_term_breakdown(
  res_c, n_top = 999, x_label = "Number of significant taxa",   # 999: only 10 orders present, no folding needed
  title = "Fungal taxa vs. insect PCoA1", subtitle = subtitle_c, tag = "c",
  pos_label = "PCoA1+ / later-season, fungivorous insect community",
  neg_label = "PCoA1- / earlier-season, bark-beetle-dominated insect community",
  legend_title = "Direction", facet_by_class = TRUE
)
panel_c <- out_c$plot

## ---- Panel d: fungal taxa ~ insect PCoA2, factor(date)-adjusted, by Order/Class ----
# Unchanged from the first pass.

res_d_all <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_factordate.csv")) %>%
  left_join(fungal_taxonomy, by = "taxon")
res_d <- res_d_all %>%
  filter(q_value < q_threshold) %>%
  mutate(group = Order)

n_d_unadj <- sum(read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_only.csv"))$q_value < q_threshold)
n_d_quad <- sum(read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_quadratic_date.csv"))$q_value < q_threshold)
subtitle_d <- paste0(nrow(res_d), "/", nrow(res_d_all), " CoCA candidates across ", n_distinct(res_d$group),
                      " orders, q<0.10,\nfactor(date)-adjusted -- but UNSTABLE across date controls:\n",
                      n_d_unadj, "/", nrow(res_d_all), " unadjusted, ", n_d_quad, "/", nrow(res_d_all),
                      " quadratic-date; small-sample-sensitive, see text")
out_d <- make_single_term_breakdown(
  res_d, n_top = 999, x_label = "Number of significant taxa",   # 999: only 4 orders present, no folding needed
  title = "Fungal taxa vs. insect PCoA2", subtitle = subtitle_d, tag = "d",
  pos_label = "PCoA2+", neg_label = "PCoA2-",
  legend_title = "Direction", facet_by_class = TRUE
)
panel_d <- out_d$plot

## ---- Write per-group counts underlying each panel ------------------------------

write.csv(out_a$counts, file.path(out_data_dir, "fig4_taxonomic_breakdown.insect_date_combined.by_subfamily.csv"), row.names = FALSE)
write.csv(out_b$counts, file.path(out_data_dir, "fig4_taxonomic_breakdown.fungal_date_combined.by_order.csv"), row.names = FALSE)
write.csv(out_c$counts, file.path(out_data_dir, "fig4_taxonomic_breakdown.fungal_insectPCoA1_factordate.by_order.csv"), row.names = FALSE)
write.csv(out_d$counts, file.path(out_data_dir, "fig4_taxonomic_breakdown.fungal_insectPCoA2_factordate.by_order.csv"), row.names = FALSE)

## ---- Combine + save -----------------------------------------------------------
# gridExtra::arrangeGrob, same reason as the first pass and Lineage A's fig5:
# panels are facet_grid'd with differing numbers of Class facets, so
# gtable::cbind's matching-row-count requirement doesn't apply.
# 3 columns, same layout idea as fig5 (busy panel(s) full-height, narrower
# single-term panels stacked in their own column): a and b (combined-term
# date tests, denser -- 4 dodged bars per group) each get a full-height
# column; c and d (single-term PCoA tests) stack in a third, narrower column.

combined <- gridExtra::arrangeGrob(panel_a, panel_b, gridExtra::arrangeGrob(panel_c, panel_d, ncol = 1),
                                    ncol = 3, widths = c(1, 1, 0.8))

ggsave(file.path(out_fig_dir, "fig4_association_taxonomic_breakdown.png"), combined,
       width = 26, height = 13, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_fig_dir, "fig4_association_taxonomic_breakdown.pdf"), combined,
       width = 26, height = 13, bg = "white", limitsize = FALSE)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig4_association_taxonomic_breakdown.png"), "and .pdf\n")
cat("Wrote per-group counts to", out_data_dir, "(insect_date_combined.by_subfamily, fungal_date_combined.by_order, fungal_insectPCoA{1,2}_factordate.by_order)\n")
