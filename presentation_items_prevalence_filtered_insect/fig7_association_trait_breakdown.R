##############################################################################
# Presentation figure 7 -- Prevalence-Filtered Insect Table Variant (Lineage B)
#
# Lineage-B counterpart of presentation_items/fig7_association_trait_
# breakdown.R. Same question (trophic mode x growth form breakdown of
# significant fungal hits, via the FungalTraits/Polme et al. 2020 genus
# join) and same growth_form-faceted diverging-bar technique, but following
# this lineage's fig4/fig5 design:
#
#   a  fungal trait vs. date -- LINEAR and QUADRATIC date terms shown
#      TOGETHER on one set of trait-group axes as 4 dodged bars (blue/
#      vermillion = linear later/earlier season, sky-blue/orange = quadratic
#      dip/hump), same mechanism as fig4 panel a/b. Top-n fold cutoff (top
#      25 trait groups + "Other") ranks by the COMBINED linear+quadratic hit
#      count, per the same user request that shaped fig4/fig5.
#   b  fungal trait vs. insect_PCoA1, factor(date)-adjusted -- single term
#      (unchanged design), same reasoning as fig4 panel c: PCoA1 has no
#      quadratic-term counterpart to dodge against (it's a different
#      predictor axis, not a linear/quadratic pair of one test).
#   c  fungal trait vs. insect_PCoA2, factor(date)-adjusted -- single term,
#      new addition (Lineage A's fig7 had only one PCoA panel; here PCoA1
#      and PCoA2 are kept as separate single-term panels, matching fig4's
#      c/d split rather than merged, since their direction semantics differ
#      -- PCoA1 is the guild-level date axis, PCoA2 a within-Scolytinae
#      species-composition axis, see the interpretation doc sec. 3).
#
# Panel a is the busiest (2490 combined-hit taxa across 42 trait groups) and
# gets a full-height column; b/c (17 and 9 hits respectively, few trait
# groups) stack in a narrower second column -- same layout idea as fig5.
#
# "Unclassified genus" / "No FungalTraits match" are kept as their own
# pseudo-facet bars (not dropped), same as Lineage A's fig7 (fig8's Class-
# faceted companion drops them instead -- see that script).
#
# Standalone script -- edit n_top_trait_groups_date below to change the
# panel-a fold threshold. Also writes the per-group/per-category counts
# underlying each panel's bars to data/presentation_items_prevalence_
# filtered_insect/fig7_trait_breakdown.*.csv.
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
n_top_trait_groups_date <- 25   # + "Other"; panel a only -- panels b/c show all trait groups present (few hits), no folding needed

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
    panel.spacing.y = unit(0.15, "lines")
  )

## ---- Panel builder: COMBINED linear+quadratic, dodged bars, growth-form facets ----
# res needs `group` (trait_group), `facet` (growth_form, NA allowed -- a
# group with no coarser-rank value of its own gets its OWN name as its
# pseudo-facet, same convention as Lineage A's fig7), `t_stat`, `q_value`,
# `t_stat_quad`, `q_value_quad`. Fold-cutoff ranking and dodge-slot
# mechanics are identical to fig4's make_combined_breakdown (join before
# factor-conversion -- see that script's note on the dplyr::left_join
# factor-coercion bug this avoids).

make_combined_breakdown <- function(res, n_top, x_label, title, subtitle, tag) {
  res <- res %>%
    mutate(group = ifelse(is.na(group) | group == "", "Unclassified", group),
           sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold) %>%
    filter(sig_lin | sig_quad)

  group_rank <- res %>% group_by(group) %>%
    summarise(total = sum(sig_lin) + sum(sig_quad), .groups = "drop") %>%
    arrange(desc(total))

  if (nrow(group_rank) > n_top) {
    top_groups <- head(group_rank$group, n_top)
    res <- res %>% mutate(group = ifelse(group %in% top_groups, group, "Other"))
    group_levels <- c(top_groups, "Other")
  } else {
    group_levels <- group_rank$group
  }

  group_facet <- res %>% distinct(group, facet) %>%
    mutate(facet = case_when(group == "Other" ~ "Other", is.na(facet) ~ group, TRUE ~ facet)) %>%
    distinct(group, facet)
  group_totals <- res %>% group_by(group) %>% summarise(total = sum(sig_lin) + sum(sig_quad), .groups = "drop")

  facet_rank <- group_facet %>% left_join(group_totals, by = "group") %>%
    group_by(facet) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
  real_facets <- facet_rank$facet[!facet_rank$facet %in% c("Other", "Unclassified genus", "No FungalTraits match")]
  facet_levels <- c(real_facets, intersect(c("Other", "No FungalTraits match", "Unclassified genus"), facet_rank$facet))

  group_order <- group_facet %>% left_join(group_totals, by = "group") %>%
    mutate(facet = factor(facet, levels = facet_levels)) %>% arrange(facet, desc(total))
  group_levels <- group_order$group

  bars <- bind_rows(
    res %>% filter(sig_lin) %>% transmute(group, cat = ifelse(t_stat > 0, cat_levels[1], cat_levels[2])),
    res %>% filter(sig_quad) %>% transmute(group, cat = ifelse(t_stat_quad > 0, cat_levels[3], cat_levels[4]))
  )

  counts <- bars %>% count(group, cat, name = "n") %>%
    complete(group = group_levels, cat = cat_levels, fill = list(n = 0))
  counts <- counts %>% left_join(group_facet, by = "group") %>%
    mutate(facet = factor(facet, levels = facet_levels))
  counts <- counts %>%
    mutate(
      group = factor(group, levels = rev(group_levels)),
      cat = factor(cat, levels = cat_levels),
      signed_n = ifelse(cat %in% pos_cats, n, -n)
    )

  group_labels <- setNames(sub(" : .*$", "", levels(counts$group)), levels(counts$group))

  p <- ggplot(counts, aes(x = group, y = signed_n, fill = cat)) +
    geom_col(color = "white", linewidth = 0.15, width = 0.8,
             position = position_dodge2(width = 0.8, padding = 0.1)) +
    geom_text(aes(label = ifelse(n == 0, "", n), hjust = ifelse(signed_n >= 0, -0.3, 1.3)),
              position = position_dodge2(width = 0.8, padding = 0.1), size = 2.2) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = group_labels) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
    scale_fill_manual(values = cat_colors, name = "Term / direction", drop = FALSE) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = x_label, title = title, subtitle = subtitle, tag = tag) +
    panel_theme +
    facet_grid(rows = vars(facet), scales = "free_y", space = "free_y", switch = "y")

  out_counts <- counts %>%
    transmute(facet = as.character(facet), group = as.character(group), cat = as.character(cat), n) %>%
    filter(n > 0) %>% arrange(desc(n))

  list(plot = p, counts = out_counts)
}

## ---- Panel builder: SINGLE term, growth-form facets (unchanged from Lineage A) ----

make_single_term_breakdown <- function(res, n_top, x_label, title, subtitle, tag,
                                        pos_label, neg_label, legend_title) {
  res <- res %>% mutate(group = ifelse(is.na(group) | group == "", "Unclassified", group))

  group_rank <- res %>% count(group, name = "total") %>% arrange(desc(total))
  if (nrow(group_rank) > n_top) {
    top_groups <- head(group_rank$group, n_top)
    res <- res %>% mutate(group = ifelse(group %in% top_groups, group, "Other"))
    group_levels <- c(top_groups, "Other")
  } else {
    group_levels <- group_rank$group
  }

  group_facet <- res %>% distinct(group, facet) %>%
    mutate(facet = case_when(group == "Other" ~ "Other", is.na(facet) ~ group, TRUE ~ facet)) %>%
    distinct(group, facet)
  group_totals <- res %>% count(group, name = "total")

  facet_rank <- group_facet %>% left_join(group_totals, by = "group") %>%
    group_by(facet) %>% summarise(total = sum(total), .groups = "drop") %>% arrange(desc(total))
  real_facets <- facet_rank$facet[!facet_rank$facet %in% c("Other", "Unclassified genus", "No FungalTraits match")]
  facet_levels <- c(real_facets, intersect(c("Other", "No FungalTraits match", "Unclassified genus"), facet_rank$facet))

  group_order <- group_facet %>% left_join(group_totals, by = "group") %>%
    mutate(facet = factor(facet, levels = facet_levels)) %>% arrange(facet, desc(total))
  group_levels <- group_order$group

  counts <- res %>%
    mutate(direction = ifelse(t_stat > 0, pos_label, neg_label)) %>%
    count(group, direction)
  counts <- counts %>% left_join(group_facet, by = "group") %>%
    mutate(facet = factor(facet, levels = facet_levels))
  counts <- counts %>%
    mutate(
      group = factor(group, levels = rev(group_levels)),
      signed_n = ifelse(direction == pos_label, n, -n)
    )

  group_labels <- setNames(sub(" : .*$", "", levels(counts$group)), levels(counts$group))

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
    panel_theme +
    facet_grid(rows = vars(facet), scales = "free_y", space = "free_y", switch = "y")

  out_counts <- counts %>%
    transmute(facet = as.character(facet), group = as.character(group), direction, n) %>% arrange(desc(n))
  list(plot = p, counts = out_counts)
}

## ---- Fungal trait lookup (shared by all panels) ---------------------------------

fungal_traits_path <- "/Users/ericmorrison/repo/FungalTraits_Polme/Polme_FungalTraits_1.2_ver_16Dec_2020.csv"

fungal_traits <- read.csv(fungal_traits_path, fileEncoding = "latin1") %>%
  distinct(GENUS, .keep_all = TRUE) %>%
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

## ---- Panel a: fungal trait vs. date, LINEAR + QUADRATIC combined ---------------

res_a <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  left_join(genus_map, by = "taxon") %>%
  mutate(group = trait_group, facet = growth_form_facet)

res_a_sig <- res_a %>% mutate(sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold) %>%
  filter(sig_lin | sig_quad)
n_a_lin <- sum(res_a_sig$sig_lin); n_a_quad <- sum(res_a_sig$sig_quad)
subtitle_a <- paste0(nrow(res_a_sig), " taxa significant in >=1 term (", n_a_lin, " linear + ", n_a_quad,
                      " quadratic = ", n_a_lin + n_a_quad, " hit-instances) across ",
                      n_distinct(ifelse(is.na(res_a_sig$group) | res_a_sig$group == "", "Unclassified", res_a_sig$group)),
                      " trait groups (top ", n_top_trait_groups_date, " shown)")
out_a <- make_combined_breakdown(
  res_a, n_top = n_top_trait_groups_date, x_label = "Number of significant taxa",
  title = "Fungal trait vs. date", subtitle = subtitle_a, tag = "a"
)
panel_a <- out_a$plot

## ---- Panel b: fungal trait vs. insect PCoA1, factor(date)-adjusted -------------

res_b <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA1_plus_factordate.csv")) %>%
  left_join(genus_map, by = "taxon") %>%
  mutate(group = trait_group, facet = growth_form_facet)
res_b_sig <- res_b %>% filter(q_value < q_threshold)

subtitle_b <- paste0(nrow(res_b_sig), "/", nrow(res_b), " CoCA candidates across ", n_distinct(res_b_sig$group),
                      " trait groups, q<0.10, factor(date)-adjusted")
out_b <- make_single_term_breakdown(
  res_b_sig, n_top = 999, x_label = "Number of significant taxa",
  title = "Fungal trait vs. insect PCoA1", subtitle = subtitle_b, tag = "b",
  pos_label = "PCoA1+ / later-season, fungivorous insect community",
  neg_label = "PCoA1- / earlier-season, bark-beetle-dominated insect community",
  legend_title = "Direction"
)
panel_b <- out_b$plot

## ---- Panel c: fungal trait vs. insect PCoA2, factor(date)-adjusted -------------

res_c <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_factordate.csv")) %>%
  left_join(genus_map, by = "taxon") %>%
  mutate(group = trait_group, facet = growth_form_facet)
res_c_sig <- res_c %>% filter(q_value < q_threshold)

n_c_unadj <- sum(read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_only.csv"))$q_value < q_threshold)
n_c_quad <- sum(read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_quadratic_date.csv"))$q_value < q_threshold)
subtitle_c <- paste0(nrow(res_c_sig), "/", nrow(res_c), " CoCA candidates across ", n_distinct(res_c_sig$group),
                      " trait groups, q<0.10,\nfactor(date)-adjusted -- but UNSTABLE across date controls:\n",
                      n_c_unadj, "/", nrow(res_c), " unadjusted, ", n_c_quad, "/", nrow(res_c),
                      " quadratic-date; small-sample-sensitive, see text")
out_c <- make_single_term_breakdown(
  res_c_sig, n_top = 999, x_label = "Number of significant taxa",
  title = "Fungal trait vs. insect PCoA2", subtitle = subtitle_c, tag = "c",
  pos_label = "PCoA2+", neg_label = "PCoA2-",
  legend_title = "Direction"
)
panel_c <- out_c$plot

## ---- Write per-group counts underlying each panel -------------------------------

write.csv(out_a$counts, file.path(out_data_dir, "fig7_trait_breakdown.fungal_date_combined.csv"), row.names = FALSE)
write.csv(out_b$counts, file.path(out_data_dir, "fig7_trait_breakdown.fungal_insectPCoA1_factordate.csv"), row.names = FALSE)
write.csv(out_c$counts, file.path(out_data_dir, "fig7_trait_breakdown.fungal_insectPCoA2_factordate.csv"), row.names = FALSE)

## ---- Combine + save --------------------------------------------------------------
# gridExtra::arrangeGrob, same reason as fig4/fig5: facet_grid'd panels have
# differing row structures, so gtable::cbind's matching-row-count
# requirement doesn't apply. Panel a (combined dodge, busiest) gets a full-
# height column; b/c (single-term, few hits) stack in a narrower column --
# same layout idea as fig5.

combined <- gridExtra::arrangeGrob(panel_a, gridExtra::arrangeGrob(panel_b, panel_c, ncol = 1),
                                    ncol = 2, widths = c(1.3, 1))

ggsave(file.path(out_fig_dir, "fig7_association_trait_breakdown.png"), combined,
       width = 18, height = 14, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_fig_dir, "fig7_association_trait_breakdown.pdf"), combined,
       width = 18, height = 14, bg = "white", limitsize = FALSE)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig7_association_trait_breakdown.png"), "and .pdf\n")
cat("Wrote per-group counts to", out_data_dir, "(fig7_trait_breakdown.fungal_date_combined, fig7_trait_breakdown.fungal_insectPCoA{1,2}_factordate)\n")
