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
# REVISION (2026-09-02, third pass): matches fig4's third-pass redesign --
# panel a now uses the linear-term-takes-priority "peak timing" reading (see
# fig4_association_taxonomic_breakdown.R's header and project_organization.
# md's "Peak-timing convention"): each taxon gets exactly ONE of 4 mutually-
# exclusive categories (Peaks late/early season, Peaks mid-season, Bimodal),
# not up to 2 independent hit-instances. Order selection (top n_top_fungal_
# date_orders) and family-within-order ranking both switch to plain
# non-overlapping taxon counts to match. Panels b/c (PCoA1/PCoA2) unchanged.
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

# Peak-timing colors (2026-09-02 revision) -- see fig4_association_
# taxonomic_breakdown.R's header. Colors unchanged from the old linear/
# quadratic pairs for visual continuity.
lin_pos_color <- "#0072B2"    # Peaks late season (significant linear, t>0)
lin_neg_color <- "#D55E00"    # Peaks early season (significant linear, t<0)
quad_pos_color <- "#56B4E9"   # Bimodal, early+late (no sig linear; quadratic dip, t_quad>0)
quad_neg_color <- "#E69F00"   # Peaks mid-season (no sig linear; quadratic hump, t_quad<0)

# REVISION (2026-09-02, fourth pass, matches fig4): with 4 mutually-
# exclusive categories (not a signed t-statistic), a diverging left/right
# layout no longer maps to anything -- all 4 bars now run the same
# direction (0 -> positive). cat_levels controls the within-group dodge
# order; empirically, position_dodge2 + coord_flip renders the FIRST level
# at the BOTTOM of each group's band and the LAST level at the TOP (i.e.
# top-to-bottom is the REVERSE of this vector) -- so cat_levels is written
# back-to-front here to get the requested top-to-bottom reading order:
# early, bimodal, mid-season, late.
cat_levels <- c("Peaks late season", "Peaks mid-season",
                 "Bimodal (early + late)", "Peaks early season")
cat_colors <- c("Peaks late season" = lin_pos_color, "Peaks early season" = lin_neg_color,
                 "Bimodal (early + late)" = quad_pos_color, "Peaks mid-season" = quad_neg_color)[cat_levels]

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

## ---- Panel builder: peak-timing categories, dodged bars, family within order --
# res needs `Order`, `Family`, `t_stat`, `q_value`, `t_stat_quad`,
# `q_value_quad`. order_levels fixes which orders appear/facet order (caller
# computes this so it matches fig4 panel a/b). REVISED 2026-09-02: a taxon
# gets exactly ONE of the 4 cat_levels (linear term takes priority whenever
# it's significant -- see fig4's header), not up to 2 independent hit-
# instances, so family ranking within each order and the bar counts are
# plain non-overlapping taxon counts. Families within each order fold to the
# top n_top_families_per_order (+"Other"). Every kept family gets all 4 cat
# levels (dodge slots), zero-filled, via complete() -- same alignment
# reasoning as fig4's make_combined_breakdown.
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
           sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold,
           # Literal strings, NOT cat_levels[i] -- see fig4's note on why
           # positional indexing into cat_levels is fragile.
           cat = case_when(
             sig_lin & t_stat > 0                  ~ "Peaks late season",
             sig_lin & t_stat < 0                  ~ "Peaks early season",
             !sig_lin & sig_quad & t_stat_quad > 0 ~ "Bimodal (early + late)",
             !sig_lin & sig_quad & t_stat_quad < 0 ~ "Peaks mid-season",
             TRUE ~ NA_character_
           )) %>%
    filter(!is.na(cat))

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

  res <- res %>% mutate(fam_id = paste(Order, Family, sep = "___"))

  counts <- res %>% count(fam_id, cat, name = "n") %>%
    complete(fam_id = fam_order$fam_id, cat = cat_levels, fill = list(n = 0))
  # Join Order/Family back on while fam_id is still character (see note above).
  counts <- counts %>% left_join(fam_order %>% select(fam_id, Order, Family), by = "fam_id")
  counts <- counts %>%
    mutate(
      fam_id = factor(fam_id, levels = rev(fam_order$fam_id)),
      cat = factor(cat, levels = cat_levels)
    )

  fam_labels <- setNames(fam_order$Family, fam_order$fam_id)

  # All 4 bars now run the same direction (0 -> positive) -- no more
  # signed_n diverging layout, matching fig4's fourth-pass revision.
  p <- ggplot(counts, aes(x = fam_id, y = n, fill = cat)) +
    geom_col(color = "white", linewidth = 0.12, width = 0.8,
             position = position_dodge2(width = 0.8, padding = 0.1)) +
    geom_text(aes(label = ifelse(n == 0, "", n)), hjust = -0.3,
              position = position_dodge2(width = 0.8, padding = 0.1), size = 1.9) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
    facet_grid(rows = vars(Order), scales = "free_y", space = "free_y", switch = "y") +
    coord_flip(clip = "off") +
    scale_x_discrete(labels = fam_labels) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    scale_fill_manual(values = cat_colors, name = "Peak timing", drop = FALSE) +
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
         sig_lin = q_value < q_threshold, sig_quad = q_value_quad < q_threshold,
         peak = case_when(
           sig_lin & t_stat > 0                  ~ "late",
           sig_lin & t_stat < 0                  ~ "early",
           !sig_lin & sig_quad & t_stat_quad > 0 ~ "bimodal",
           !sig_lin & sig_quad & t_stat_quad < 0 ~ "mid",
           TRUE ~ NA_character_
         )) %>%
  filter(!is.na(peak))

order_rank_a <- res_a_sig %>% count(Order, name = "total") %>% arrange(desc(total))
order_levels_a <- head(order_rank_a$Order, n_top_fungal_date_orders)
order_levels_a <- setdiff(order_levels_a, "Unclassified")

res_a_kept <- res_a_sig %>% filter(Order %in% order_levels_a)
n_a_late <- sum(res_a_kept$peak == "late"); n_a_early <- sum(res_a_kept$peak == "early")
n_a_mid <- sum(res_a_kept$peak == "mid"); n_a_bimodal <- sum(res_a_kept$peak == "bimodal")
subtitle_a <- paste0(nrow(res_a_kept), " taxa with a significant peak-timing call (", n_a_late, " late + ",
                      n_a_early, " early + ", n_a_mid, " mid-season + ", n_a_bimodal, " bimodal) in ",
                      length(order_levels_a), " of the orders from fig4a/b, by family")
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
