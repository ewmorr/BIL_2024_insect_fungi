##############################################################################
# Presentation figure 8 -- Prevalence-Filtered Insect Table Variant (Lineage B)
#
# Lineage-B counterpart of presentation_items/fig8_association_trait_
# breakdown.R. Same question as fig7 (trophic-mode breakdown of significant
# fungal hits via the FungalTraits/Polme et al. 2020 genus join) but faceted
# by taxonomic Class instead of growth_form, with growth_form dropped from
# the grouping entirely (same rationale as Lineage A: once split by Class,
# growth_form adds little). Because trait_group and Class are independent
# classifications (not a 1:1 nesting like fig4's Order-within-Class), a bar
# is drawn per (trait_group, Class, category) triple, with an
# axis_id = "<facet rank>|<group>" trick so the same lifestyle can occupy a
# separate row in each Class facet it appears in.
#
# Follows this lineage's fig4/fig5/fig7 design:
#   a  fungal trait vs. date -- LINEAR and QUADRATIC date terms shown
#      TOGETHER as 4 dodged bars per (trait_group, Class) row. The two
#      min_facet_hits thresholds (drop a Class facet entirely below the
#      threshold; fold a (lifestyle, Class) cell below the threshold into
#      that facet's own "Other") now use the COMBINED linear+quadratic hit
#      count, same convention as fig4/fig5/fig7's top-n fold.
#   b  fungal trait vs. insect_PCoA1, factor(date)-adjusted -- single term
#      (no quadratic counterpart to dodge against, same reasoning as fig4c/
#      fig7b).
#   c  fungal trait vs. insect_PCoA2, factor(date)-adjusted -- single term,
#      new addition (kept separate from PCoA1 rather than merged, matching
#      fig4's c/d split -- see fig7's header for the direction-semantics
#      rationale).
#
# Panel a (2490 combined-hit taxa, matches fig7/fig4) is the busiest and
# gets a full-height column; b/c (17 and 9 hits) stack in a narrower second
# column -- same layout as fig7/fig5.
#
# "Unclassified genus" / "No FungalTraits match" are DROPPED entirely here
# (not shown as their own bars), same as Lineage A's fig8 -- with Class
# facets they'd split thinly across many real classes rather than forming
# one legible bucket (fig7's growth_form facets keep them, for the opposite
# reason).
#
# Standalone script -- edit min_facet_hits_date below to change the panel-a
# drop/fold thresholds. Also writes the per-group/per-Class/per-category
# counts underlying each panel's bars to data/presentation_items_
# prevalence_filtered_insect/fig8_trait_breakdown.*.csv.
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
min_facet_hits_date <- 6   # panel a: Class facets / (lifestyle, Class) cells below this COMBINED (linear+quadratic) hit count are dropped/folded
min_facet_hits_pcoa <- 1   # panels b/c: far fewer hits, nothing dropped or folded

lin_pos_color <- "#0072B2"    # Later season (linear, t>0)
lin_neg_color <- "#D55E00"    # Earlier season (linear, t<0)
quad_pos_color <- "#56B4E9"   # Dip (quadratic, t>0)
quad_neg_color <- "#E69F00"   # Hump (quadratic, t<0)

cat_levels <- c("Later season (linear, t>0)", "Earlier season (linear, t<0)",
                 "Dip (quadratic, t>0)", "Hump (quadratic, t<0)")
pos_cats <- cat_levels[c(1, 3)]
cat_colors <- setNames(c(lin_pos_color, lin_neg_color, quad_pos_color, quad_neg_color), cat_levels)

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
    plot.title.position = "plot"
  )

## ---- Panel builder: COMBINED linear+quadratic, dodged bars, Class facets ------
# res needs `group` (trait_group), `facet` (Class), `t_stat`, `q_value`,
# `t_stat_quad`, `q_value_quad`. Unlike fig7's growth_form facets (1:1 with
# group after fold), `group`/`facet` here are independent, so bars key on
# (group, facet, cat) triples with an axis_id trick (same as Lineage A's
# fig8). Every kept (group,facet) row gets all 4 cat levels (dodge slots),
# zero-filled -- same alignment reasoning as fig4/fig7. NOTE: every join
# below happens on CHARACTER group/facet/axis_id columns, converted to
# factor only in the final mutate -- seefig4's note on the dplyr::left_join
# factor-coercion bug this ordering avoids.

make_combined_breakdown <- function(res, x_label, title, subtitle, tag, min_facet_hits = 1) {
  res <- res %>%
    mutate(
      group = ifelse(is.na(group) | group == "", "Unclassified", group),
      facet = ifelse(is.na(facet) | facet == "", "Unclassified", facet),
      sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold
    ) %>%
    filter(sig_lin | sig_quad)

  all_facet_totals <- res %>% group_by(facet) %>%
    summarise(total = sum(sig_lin) + sum(sig_quad), .groups = "drop") %>% arrange(desc(total))
  dropped_facets <- all_facet_totals %>% filter(total < min_facet_hits)
  if (nrow(dropped_facets) > 0) {
    cat("  [", title, "] dropping", nrow(dropped_facets), "classes below min_facet_hits =", min_facet_hits,
        "(", sum(dropped_facets$total), "hit-instances total):", paste(dropped_facets$facet, collapse = ", "), "\n")
  }
  res <- res %>% filter(facet %in% all_facet_totals$facet[all_facet_totals$total >= min_facet_hits])

  cell_totals <- res %>% group_by(group, facet) %>%
    summarise(cell_total = sum(sig_lin) + sum(sig_quad), .groups = "drop")
  small_cells <- cell_totals %>% filter(cell_total < min_facet_hits)
  if (nrow(small_cells) > 0) {
    cat("  [", title, "] folding", nrow(small_cells), "lifestyle-within-class cells below min_facet_hits =", min_facet_hits,
        "(", sum(small_cells$cell_total), "hit-instances total) into each class's Other\n")
  }
  res <- res %>%
    left_join(cell_totals, by = c("group", "facet")) %>%
    mutate(group = ifelse(cell_total < min_facet_hits, "Other", group)) %>%
    select(-cell_total)

  row_totals <- res %>% group_by(group, facet) %>%
    summarise(total = sum(sig_lin) + sum(sig_quad), .groups = "drop")

  facet_totals <- row_totals %>% group_by(facet) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
  real_facets <- facet_totals$facet[facet_totals$facet != "Unclassified"]
  facet_levels <- c(real_facets, intersect("Unclassified", facet_totals$facet))

  row_totals <- row_totals %>%
    mutate(facet_rank = match(facet, facet_levels)) %>%
    arrange(facet_rank, desc(total)) %>%
    mutate(axis_id = paste0(sprintf("%03d", facet_rank), "|", group))

  axis_order <- row_totals$axis_id
  axis_labels <- setNames(row_totals$group, row_totals$axis_id)

  bars <- bind_rows(
    res %>% filter(sig_lin) %>% transmute(group, facet, cat = ifelse(t_stat > 0, cat_levels[1], cat_levels[2])),
    res %>% filter(sig_quad) %>% transmute(group, facet, cat = ifelse(t_stat_quad > 0, cat_levels[3], cat_levels[4]))
  )

  full_grid <- row_totals %>% select(group, facet) %>% distinct() %>% crossing(cat = cat_levels)
  counts <- bars %>% count(group, facet, cat, name = "n") %>%
    right_join(full_grid, by = c("group", "facet", "cat")) %>%
    mutate(n = ifelse(is.na(n), 0, n))
  # Join axis_id on while group/facet are still character.
  counts <- counts %>% left_join(row_totals %>% select(group, facet, axis_id), by = c("group", "facet"))
  counts <- counts %>%
    mutate(
      facet = factor(facet, levels = facet_levels),
      axis_id = factor(axis_id, levels = rev(axis_order)),
      cat = factor(cat, levels = cat_levels),
      signed_n = ifelse(cat %in% pos_cats, n, -n)
    )

  p <- ggplot(counts, aes(x = axis_id, y = signed_n, fill = cat)) +
    geom_col(color = "white", linewidth = 0.15, width = 0.8,
             position = position_dodge2(width = 0.8, padding = 0.1)) +
    geom_text(aes(label = ifelse(n == 0, "", n), hjust = ifelse(signed_n >= 0, -0.3, 1.3)),
              position = position_dodge2(width = 0.8, padding = 0.1), size = 2.2) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = axis_labels) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
    scale_fill_manual(values = cat_colors, name = "Term / direction", drop = FALSE) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = x_label, title = title, subtitle = subtitle, tag = tag) +
    panel_theme +
    facet_grid(rows = vars(facet), scales = "free_y", space = "free_y", switch = "y")

  out_counts <- counts %>% filter(n > 0) %>%
    transmute(Class = as.character(facet), lifestyle = as.character(group), cat = as.character(cat), n) %>%
    arrange(desc(n))

  list(plot = p, counts = out_counts)
}

## ---- Panel builder: SINGLE term, Class facets (unchanged from Lineage A) ------

make_single_term_breakdown <- function(res, x_label, title, subtitle, tag,
                                        pos_label, neg_label, legend_title, min_facet_hits = 1) {
  res <- res %>%
    mutate(
      group = ifelse(is.na(group) | group == "", "Unclassified", group),
      facet = ifelse(is.na(facet) | facet == "", "Unclassified", facet),
      direction = ifelse(t_stat > 0, pos_label, neg_label)
    )

  all_facet_totals <- res %>% count(facet, name = "total") %>% arrange(desc(total))
  res <- res %>% filter(facet %in% all_facet_totals$facet[all_facet_totals$total >= min_facet_hits])

  cell_totals <- res %>% count(group, facet, name = "cell_total")
  res <- res %>%
    left_join(cell_totals, by = c("group", "facet")) %>%
    mutate(group = ifelse(cell_total < min_facet_hits, "Other", group)) %>%
    select(-cell_total)

  row_totals <- res %>% count(group, facet, name = "total")

  facet_totals <- row_totals %>% group_by(facet) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
  real_facets <- facet_totals$facet[facet_totals$facet != "Unclassified"]
  facet_levels <- c(real_facets, intersect("Unclassified", facet_totals$facet))

  row_totals <- row_totals %>%
    mutate(facet_rank = match(facet, facet_levels)) %>%
    arrange(facet_rank, desc(total)) %>%
    mutate(axis_id = paste0(sprintf("%03d", facet_rank), "|", group))

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
    scale_fill_manual(values = setNames(c(lin_pos_color, lin_neg_color), c(pos_label, neg_label)),
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

## ---- Fungal trait + taxonomy lookup (shared by all panels) --------------------

fungal_traits_path <- "/Users/ericmorrison/repo/FungalTraits_Polme/Polme_FungalTraits_1.2_ver_16Dec_2020.csv"

fungal_traits <- read.csv(fungal_traits_path, fileEncoding = "latin1") %>%
  distinct(GENUS, .keep_all = TRUE) %>%
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
    trait_group = case_when(
      is.na(Genus) ~ "Unclassified genus",
      is.na(primary_lifestyle) ~ "No FungalTraits match",
      secondary_plant_pathogen & primary_lifestyle != "plant_pathogen" ~
        paste0(gsub("_", " ", primary_lifestyle), " (+plant pathogen)"),
      TRUE ~ gsub("_", " ", primary_lifestyle)
    )
  ) %>%
  select(taxon, trait_group, Class)

## ---- Panel a: fungal trait vs. date, LINEAR + QUADRATIC combined, by Class -----

res_a <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  left_join(genus_map, by = "taxon") %>%
  filter(!trait_group %in% c("Unclassified genus", "No FungalTraits match")) %>%
  mutate(group = trait_group, facet = Class)

res_a_sig <- res_a %>% mutate(sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold) %>%
  filter(sig_lin | sig_quad)
n_a_lin <- sum(res_a_sig$sig_lin); n_a_quad <- sum(res_a_sig$sig_quad)
subtitle_a <- paste0(nrow(res_a_sig), " taxa with a FungalTraits match, significant in >=1 term (",
                      n_a_lin, " linear + ", n_a_quad, " quadratic = ", n_a_lin + n_a_quad, " hit-instances) across ",
                      n_distinct(res_a_sig$group), " lifestyles; classes/cells with <", min_facet_hits_date,
                      " hit-instances omitted/folded")
out_a <- make_combined_breakdown(
  res_a, x_label = "Number of significant taxa",
  title = "Fungal trait vs. date", subtitle = subtitle_a, tag = "a", min_facet_hits = min_facet_hits_date
)
panel_a <- out_a$plot

## ---- Panel b: fungal trait vs. insect PCoA1, factor(date)-adjusted, by Class ----

res_b <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA1_plus_factordate.csv")) %>%
  left_join(genus_map, by = "taxon") %>%
  filter(!trait_group %in% c("Unclassified genus", "No FungalTraits match")) %>%
  mutate(group = trait_group, facet = Class)
res_b_sig <- res_b %>% filter(q_value < q_threshold)

subtitle_b <- paste0(nrow(res_b_sig), " taxa with a FungalTraits match across ", n_distinct(res_b_sig$group),
                      " lifestyles, all classes/cells shown, q<0.10, factor(date)-adjusted")
out_b <- make_single_term_breakdown(
  res_b_sig, x_label = "Number of significant taxa",
  title = "Fungal trait vs. insect PCoA1", subtitle = subtitle_b, tag = "b",
  pos_label = "PCoA1+ / later-season, fungivorous insect community",
  neg_label = "PCoA1- / earlier-season, bark-beetle-dominated insect community",
  legend_title = "Direction", min_facet_hits = min_facet_hits_pcoa
)
panel_b <- out_b$plot

## ---- Panel c: fungal trait vs. insect PCoA2, factor(date)-adjusted, by Class ----

res_c <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_factordate.csv")) %>%
  left_join(genus_map, by = "taxon") %>%
  filter(!trait_group %in% c("Unclassified genus", "No FungalTraits match")) %>%
  mutate(group = trait_group, facet = Class)
res_c_sig <- res_c %>% filter(q_value < q_threshold)

n_c_unadj <- sum(read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_only.csv"))$q_value < q_threshold)
n_c_quad <- sum(read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_quadratic_date.csv"))$q_value < q_threshold)
subtitle_c <- paste0(nrow(res_c_sig), " taxa with a FungalTraits match across ", n_distinct(res_c_sig$group),
                      " lifestyles, q<0.10,\nfactor(date)-adjusted -- but UNSTABLE across date controls:\n",
                      n_c_unadj, "/200 unadjusted, ", n_c_quad, "/200 quadratic-date; small-sample-sensitive, see text")
out_c <- make_single_term_breakdown(
  res_c_sig, x_label = "Number of significant taxa",
  title = "Fungal trait vs. insect PCoA2", subtitle = subtitle_c, tag = "c",
  pos_label = "PCoA2+", neg_label = "PCoA2-",
  legend_title = "Direction", min_facet_hits = min_facet_hits_pcoa
)
panel_c <- out_c$plot

## ---- Write per-group/per-Class counts underlying each panel -------------------

write.csv(out_a$counts, file.path(out_data_dir, "fig8_trait_breakdown.fungal_date_combined.csv"), row.names = FALSE)
write.csv(out_b$counts, file.path(out_data_dir, "fig8_trait_breakdown.fungal_insectPCoA1_factordate.csv"), row.names = FALSE)
write.csv(out_c$counts, file.path(out_data_dir, "fig8_trait_breakdown.fungal_insectPCoA2_factordate.csv"), row.names = FALSE)

## ---- Combine + save --------------------------------------------------------------
# gridExtra::arrangeGrob, same reason as fig4/fig5/fig7. Panel a (combined
# dodge, busiest) gets a full-height column; b/c stack in a narrower column.

combined <- gridExtra::arrangeGrob(panel_a, gridExtra::arrangeGrob(panel_b, panel_c, ncol = 1),
                                    ncol = 2, widths = c(1.3, 1))

ggsave(file.path(out_fig_dir, "fig8_association_trait_breakdown.png"), combined,
       width = 18, height = 18, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_fig_dir, "fig8_association_trait_breakdown.pdf"), combined,
       width = 18, height = 18, bg = "white", limitsize = FALSE)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig8_association_trait_breakdown.png"), "and .pdf\n")
cat("Wrote per-group/per-Class counts to", out_data_dir, "(fig8_trait_breakdown.fungal_date_combined, fig8_trait_breakdown.fungal_insectPCoA{1,2}_factordate)\n")
