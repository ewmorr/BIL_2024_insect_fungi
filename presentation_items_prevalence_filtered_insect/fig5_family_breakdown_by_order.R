##############################################################################
# Presentation figure 5 (supplemental) -- Prevalence-Filtered Insect Table
# Variant (Lineage B). Family-level detail behind this lineage's fig4 order-
# level fungal panels -- "which families make up each major order?" Bars are
# grouped/labeled by their subsuming Order (facet strips), families within
# each order ranked by significant-hit count.
#
# REVISION (2026-09-01, second pass): matches fig4's second-pass redesign.
# fig4's date panel (a) now shows the linear and quadratic date terms
# TOGETHER as 4 dodged (side-by-side) bars per taxonomic group instead of two
# separate panels -- this figure's family-within-order breakdown of that
# panel follows the same design. fig4's PCoA1/PCoA2 panels (c/d) stay
# single-term (different predictor axes, not a linear/quadratic pair of one
# test), so their family-within-order breakdowns (b/c here) are unchanged
# from the first pass.
#
# Three panels (was four in the first pass, which had separate linear/
# quadratic date panels):
#   a  fungal taxa ~ date (linear + quadratic combined), family within order
#      -- mirrors fig4 panel a's dodge design; restricted to the top n_top_
#      fungal_date_orders orders by COMBINED (linear+quadratic) hit count,
#      matching fig4 panel a's fold rule, then to the top n_top_families_
#      per_order families within each of those orders (also by combined
#      count).
#   b  fungal taxa ~ insect_PCoA1, factor(date)-adjusted, family within order
#      -- unchanged from the first pass (single term).
#   c  fungal taxa ~ insect_PCoA2, factor(date)-adjusted, family within order
#      -- unchanged from the first pass (single term). Small n (9 taxa
#      across 4 orders) but genuinely informative: Helotiales alone splits
#      into 3 different families here (incl. Tympanidaceae) that fig4's
#      order-level bar collapses into one, and Diaporthales resolves to
#      Valsaceae specifically -- the "Tympanis-as-a-group" / "Valsaceae/
#      Cytospora-as-a-group" framing the interpretation doc recommends for
#      this axis (sec. 6).
#
# Order==NA / Order=="Unclassified" facets are dropped from every panel, same
# reasoning as the first pass: taxa unresolved at Order are also unresolved
# at Family, so that facet would only ever repeat fig4's order-level
# "Unclassified" bar with no new detail.
#
# Standalone script -- edit n_top_fungal_date_orders / n_top_families_per_
# order below to change the fold thresholds (n_top_fungal_date_orders must
# match fig4_association_taxonomic_breakdown.R panel a/b).
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
n_top_fungal_date_orders <- 15   # must match fig4_association_taxonomic_breakdown.R panel a/b
n_top_families_per_order <- 6    # + "Other" within any order spanning more families than this

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
    axis.text = element_text(size = 7.5, color = "black"),
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

## ---- Panel builder: COMBINED linear+quadratic, dodged bars, family within order --
# res needs `Order`, `Family`, `t_stat`, `q_value`, `t_stat_quad`,
# `q_value_quad`. order_levels fixes which orders appear/facet order (caller
# computes this so it matches fig4 panel a/b). Families within each order
# fold to the top n_top_families_per_order (+"Other") by COMBINED hit count.
# Every kept family gets all 4 cat levels (dodge slots), zero-filled, via
# complete() -- same alignment reasoning as fig4's make_combined_breakdown.
# NOTE: any left_join on a `fam_id`/`group`-type column must happen BEFORE
# that column is converted to a factor -- dplyr::left_join silently coerces
# a factor join column back to character (dropping its levels) against a
# character partner, which upstream produced blank/"NA" axis labels in
# fig4's first draft of this redesign; fixed there and avoided here by
# joining Order/Family back on while fam_id is still character.

make_combined_family_by_order_breakdown <- function(res, order_levels, x_label, title, subtitle, tag) {
  res <- res %>%
    filter(Order %in% order_levels) %>%
    mutate(Order = factor(Order, levels = order_levels),
           Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family),
           sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold) %>%
    filter(sig_lin | sig_quad)

  fam_rank <- res %>% group_by(Order, Family) %>%
    summarise(total = sum(sig_lin) + sum(sig_quad), .groups = "drop") %>%
    group_by(Order) %>% arrange(desc(total), .by_group = TRUE) %>%
    mutate(rk = row_number()) %>% ungroup()

  res <- res %>%
    left_join(fam_rank %>% select(Order, Family, rk), by = c("Order", "Family")) %>%
    mutate(Family = ifelse(rk <= n_top_families_per_order, Family, "Other"))

  fam_order <- res %>% group_by(Order, Family) %>%
    summarise(total = sum(sig_lin) + sum(sig_quad), .groups = "drop") %>%
    mutate(is_other = Family == "Other") %>%
    arrange(Order, is_other, desc(total)) %>%
    mutate(fam_id = paste(Order, Family, sep = "___"))

  res <- res %>% mutate(fam_id = paste(Order, Family, sep = "___"))

  bars <- bind_rows(
    res %>% filter(sig_lin) %>% transmute(fam_id, cat = ifelse(t_stat > 0, cat_levels[1], cat_levels[2])),
    res %>% filter(sig_quad) %>% transmute(fam_id, cat = ifelse(t_stat_quad > 0, cat_levels[3], cat_levels[4]))
  )

  counts <- bars %>% count(fam_id, cat, name = "n") %>%
    complete(fam_id = fam_order$fam_id, cat = cat_levels, fill = list(n = 0))
  # Join Order/Family back on while fam_id is still character (see note above).
  counts <- counts %>% left_join(fam_order %>% select(fam_id, Order, Family), by = "fam_id")
  counts <- counts %>%
    mutate(
      fam_id = factor(fam_id, levels = rev(fam_order$fam_id)),
      cat = factor(cat, levels = cat_levels),
      signed_n = ifelse(cat %in% pos_cats, n, -n)
    )

  fam_labels <- setNames(fam_order$Family, fam_order$fam_id)

  p <- ggplot(counts, aes(x = fam_id, y = signed_n, fill = cat)) +
    geom_col(color = "white", linewidth = 0.12, width = 0.8,
             position = position_dodge2(width = 0.8, padding = 0.1)) +
    geom_text(aes(label = ifelse(n == 0, "", n), hjust = ifelse(signed_n >= 0, -0.3, 1.3)),
              position = position_dodge2(width = 0.8, padding = 0.1), size = 1.9) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    facet_grid(rows = vars(Order), scales = "free_y", space = "free_y", switch = "y") +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = fam_labels) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
    scale_fill_manual(values = cat_colors, name = "Term / direction", drop = FALSE) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = x_label, title = title, subtitle = subtitle, tag = tag) +
    panel_theme

  list(plot = p, counts = counts %>% filter(n > 0) %>%
         transmute(Order = as.character(Order), Family, direction = as.character(cat), n) %>%
         arrange(Order, desc(n)))
}

## ---- Panel builder: SINGLE term, family within order (unchanged from first pass) --

make_family_by_order_breakdown <- function(res, order_levels, x_label, title, subtitle, tag,
                                            pos_label, neg_label, legend_title) {
  res <- res %>%
    filter(Order %in% order_levels) %>%
    mutate(Order = factor(Order, levels = order_levels),
           Family = ifelse(is.na(Family) | Family == "", "Unclassified", Family))

  fam_rank <- res %>% count(Order, Family, name = "total") %>%
    group_by(Order) %>% arrange(desc(total), .by_group = TRUE) %>%
    mutate(rk = row_number()) %>% ungroup()

  res <- res %>%
    left_join(fam_rank %>% select(Order, Family, rk), by = c("Order", "Family")) %>%
    mutate(Family = ifelse(rk <= n_top_families_per_order, Family, "Other"))

  fam_order <- res %>% count(Order, Family, name = "total") %>%
    mutate(is_other = Family == "Other") %>%
    arrange(Order, is_other, desc(total)) %>%
    mutate(fam_id = paste(Order, Family, sep = "___"))

  counts <- res %>%
    mutate(direction = ifelse(t_stat > 0, pos_label, neg_label),
           fam_id = paste(Order, Family, sep = "___")) %>%
    count(Order, Family, fam_id, direction) %>%
    mutate(
      fam_id = factor(fam_id, levels = rev(fam_order$fam_id)),
      signed_n = ifelse(direction == pos_label, n, -n)
    )

  fam_labels <- setNames(fam_order$Family, fam_order$fam_id)

  p <- ggplot(counts, aes(x = fam_id, y = signed_n, fill = direction)) +
    geom_col(color = "white", linewidth = 0.15, width = 0.8) +
    geom_text(aes(label = n, hjust = ifelse(signed_n >= 0, -0.3, 1.3)), size = 2.4) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    facet_grid(rows = vars(Order), scales = "free_y", space = "free_y", switch = "y") +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = fam_labels) +
    scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
    scale_fill_manual(values = setNames(c(lin_pos_color, lin_neg_color), c(pos_label, neg_label)),
                       name = legend_title) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = x_label, title = title, subtitle = subtitle, tag = tag) +
    panel_theme

  list(plot = p, counts = counts %>% transmute(Order = as.character(Order), Family, direction, n) %>%
         arrange(Order, desc(n)))
}

## ---- Fungal taxonomy lookup (shared by all panels) -----------------------------

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv") %>%
  rename(taxon = X) %>%
  transmute(taxon, Order = sub("^o__", "", Order), Family = sub("^f__", "", Family))

res_fungal_date <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  left_join(fungal_taxonomy, by = "taxon")

## ---- Panel a: fungal taxa ~ date, LINEAR + QUADRATIC combined, family within order ----
# Order selection matches fig4 panel a: top n_top_fungal_date_orders orders
# by COMBINED (linear+quadratic) hit count.

res_a_sig <- res_fungal_date %>%
  mutate(Order = ifelse(is.na(Order) | Order == "", "Unclassified", Order),
         sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold) %>%
  filter(sig_lin | sig_quad)

order_rank_a <- res_a_sig %>% group_by(Order) %>%
  summarise(total = sum(sig_lin) + sum(sig_quad), .groups = "drop") %>% arrange(desc(total))
order_levels_a <- head(order_rank_a$Order, n_top_fungal_date_orders)
order_levels_a <- setdiff(order_levels_a, "Unclassified")

res_a_kept <- res_a_sig %>% filter(Order %in% order_levels_a)
n_a_lin <- sum(res_a_kept$sig_lin); n_a_quad <- sum(res_a_kept$sig_quad)
subtitle_a <- paste0(nrow(res_a_kept), " taxa significant in >=1 term (", n_a_lin, " linear + ", n_a_quad,
                      " quadratic = ", n_a_lin + n_a_quad, " hit-instances) in ", length(order_levels_a),
                      " of the orders from fig4a/b, by family")
out_a <- make_combined_family_by_order_breakdown(
  res_a_sig, order_levels_a, x_label = "Number of significant taxa",
  title = "Fungal taxa vs. date -- family within order", subtitle = subtitle_a, tag = "a"
)
panel_a <- out_a$plot

## ---- Panel b: fungal taxa ~ insect PCoA1, factor(date)-adjusted, family within order (unchanged) ----

res_b_all <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA1_plus_factordate.csv")) %>%
  left_join(fungal_taxonomy, by = "taxon") %>%
  filter(q_value < q_threshold) %>%
  mutate(Order = ifelse(is.na(Order) | Order == "", "Unclassified", Order))

order_rank_b <- res_b_all %>% count(Order, name = "total") %>% arrange(desc(total))
order_levels_b <- setdiff(order_rank_b$Order, "Unclassified")

subtitle_b <- paste0(sum(res_b_all$Order %in% order_levels_b), " taxa across ", length(order_levels_b),
                      " of the orders from fig4c, by family (q<0.10, factor(date)-adjusted)")
out_b <- make_family_by_order_breakdown(
  res_b_all, order_levels_b, x_label = "Number of significant taxa",
  title = "Fungal taxa vs. insect PCoA1 -- family within order", subtitle = subtitle_b, tag = "b",
  pos_label = "PCoA1+ / later-season, fungivorous insect community",
  neg_label = "PCoA1- / earlier-season, bark-beetle-dominated insect community",
  legend_title = "Direction"
)
panel_b <- out_b$plot

## ---- Panel c: fungal taxa ~ insect PCoA2, factor(date)-adjusted, family within order (unchanged) ----

res_c_all <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_factordate.csv")) %>%
  left_join(fungal_taxonomy, by = "taxon") %>%
  filter(q_value < q_threshold) %>%
  mutate(Order = ifelse(is.na(Order) | Order == "", "Unclassified", Order))

order_rank_c <- res_c_all %>% count(Order, name = "total") %>% arrange(desc(total))
order_levels_c <- setdiff(order_rank_c$Order, "Unclassified")

subtitle_c <- paste0(sum(res_c_all$Order %in% order_levels_c), " taxa across ", length(order_levels_c),
                      " of the orders from fig4d, by family (q<0.10, factor(date)-adjusted) --\nsmall n, unstable across date controls, see fig4d")
out_c <- make_family_by_order_breakdown(
  res_c_all, order_levels_c, x_label = "Number of significant taxa",
  title = "Fungal taxa vs. insect PCoA2 -- family within order", subtitle = subtitle_c, tag = "c",
  pos_label = "PCoA2+", neg_label = "PCoA2-",
  legend_title = "Direction"
)
panel_c <- out_c$plot

## ---- Write per-group counts underlying each panel ------------------------------

write.csv(out_a$counts, file.path(out_data_dir, "fig5_family_breakdown_by_order.fungal_date_combined.csv"), row.names = FALSE)
write.csv(out_b$counts, file.path(out_data_dir, "fig5_family_breakdown_by_order.fungal_insectPCoA1_factordate.csv"), row.names = FALSE)
write.csv(out_c$counts, file.path(out_data_dir, "fig5_family_breakdown_by_order.fungal_insectPCoA2_factordate.csv"), row.names = FALSE)

## ---- Combine + save -----------------------------------------------------------
# gridExtra::arrangeGrob, same reason as fig4 and the first pass: panels'
# facet_grid row counts (number of orders shown) differ, so gtable::cbind's
# matching-row-count requirement doesn't apply. Panel a (combined dodge,
# many families/order) is the busiest and widest; b/c are narrower single-
# term panels stacked below it.

combined <- gridExtra::arrangeGrob(panel_a, gridExtra::arrangeGrob(panel_b, panel_c, ncol = 1),
                                    ncol = 2, widths = c(1.3, 1))

ggsave(file.path(out_fig_dir, "fig5_family_breakdown_by_order.png"), combined,
       width = 20, height = 20, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_fig_dir, "fig5_family_breakdown_by_order.pdf"), combined,
       width = 20, height = 20, bg = "white", limitsize = FALSE)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig5_family_breakdown_by_order.png"), "and .pdf\n")
cat("Wrote per-group counts to", out_data_dir, "(fungal_date_combined, fungal_insectPCoA{1,2}_factordate)\n")
