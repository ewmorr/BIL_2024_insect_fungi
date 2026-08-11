##############################################################################
# Presentation figure 5 (supplemental): family-level detail behind fig4's
# order-level panels b/c -- "which families make up each major order?"
# Bars are grouped and labeled by their subsuming Order (facet strips), and
# within each order, families are ranked by significant-hit count. As a
# supplement (rather than the main summary figure), this shows more detail/
# more taxa than fig4 -- every named order fig4 shows individually is broken
# down here, with per-order family lists folded to the top few + "Other"
# only where a single order's family list is itself long-tailed. The
# Order=="Unclassified" facet is dropped in both panels (see below) -- this
# figure's own panels are tagged a/b, not fig4's b/c, since it's a standalone
# figure.
#
# Panel a: fungal taxa ~ collection date, by Family within Order -- same
#          data as fig4 panel b (data/2024_fungi/fungal_taxa_date_
#          association.all_taxa.csv, unbiased screen), restricted to the
#          same 15 major orders fig4 panel b names individually (its
#          "Other" bucket is not broken down further here), minus
#          Order=="Unclassified" -- taxa with no resolved Order also have no
#          resolved Family, so that facet only ever contained a single
#          "Unclassified" family bar identical to fig4's order-level bar,
#          adding no detail. Each remaining order's family list is folded to
#          the top n_top_families_per_order + an "Other" row where an order
#          spans more families than that (e.g. Pleosporales alone spans 24
#          families) -- 87 family-level bars total across the 14 orders.
# Panel b: fungal taxa ~ insect_PCoA1, by Family within Order -- same data
#          as fig4 panel c (data/compare_insects_fungi_top3axes/fungal_
#          insect_association_results.PCoA1_only.no_lure.csv). All 31
#          orders from fig4 panel c are shown except Order=="Unclassified"
#          (same reason as panel a); no remaining order spans more than 3
#          families, so no per-order folding is needed -- 40 family-level
#          bars total.
#
# Both panels share fig4's diverging-bar/direction-color design; the count
# (x) axis is shared across an entire panel's facets (not free per facet),
# so bar lengths stay comparable across orders.
#
# Standalone script -- edit n_top_families_per_order below to change how
# many named families are shown per order before folding the rest into that
# order's "Other".
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
n_top_fungal_date_orders <- 15   # must match fig4_association_taxonomic_breakdown.R panel b
n_top_families_per_order <- 6    # + "Other" within any order spanning more families than this

pos_color <- "#0072B2"
neg_color <- "#D55E00"

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

## ---- Shared plotting logic ---------------------------------------------------
# res: data.frame with columns `Order`, `Family`, `t_stat`. `order_levels`
# fixes which orders appear and their top-to-bottom facet order (rank by
# hit count, descending -- computed by the caller so it can match fig4).
# Families within each order are folded to the top n_top_families_per_order
# (+ "Other") when an order spans more than that many families.
# Returns list(plot, counts).

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

  # row order within each order facet: by total hit count, descending (ties
  # broken by keeping "Other" last within its order)
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
    scale_fill_manual(values = setNames(c(pos_color, neg_color), c(pos_label, neg_label)),
                       name = legend_title) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = NULL, y = x_label, title = title, subtitle = subtitle, tag = tag) +
    panel_theme

  list(plot = p, counts = counts %>% transmute(Order = as.character(Order), Family, direction, n) %>%
         arrange(Order, desc(n)))
}

## ---- Fungal taxonomy lookup (shared by panels b/c) ----------------------------

fungal_taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1) %>%
  tibble::rownames_to_column("taxon") %>%
  transmute(taxon, Order = sub("^o__", "", Order), Family = sub("^f__", "", Family))

## ---- Panel b: fungal taxa ~ collection date, family within order -------------

res_b_all <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(fungal_taxonomy, by = "taxon") %>%
  mutate(Order = ifelse(is.na(Order) | Order == "", "Unclassified", Order))

order_rank_b <- res_b_all %>% count(Order, name = "total") %>% arrange(desc(total))
order_levels_b <- head(order_rank_b$Order, n_top_fungal_date_orders)
order_levels_b <- setdiff(order_levels_b, "Unclassified")   # dropped -- families within Order=="Unclassified" are themselves entirely Family=="Unclassified", so this facet never showed anything the order-level fig4 bar didn't already

subtitle_b <- paste0(sum(res_b_all$Order %in% order_levels_b), " taxa in ", length(order_levels_b),
                      " of the orders from fig4b, by family (q<0.10)")
out_b <- make_family_by_order_breakdown(
  res_b_all, order_levels_b, x_label = "Number of significant taxa",
  title = "Fungal taxa vs. date -- family within order", subtitle = subtitle_b, tag = "a",
  pos_label = "Later season (t>0)", neg_label = "Earlier season (t<0)",
  legend_title = "Direction"
)
panel_b <- out_b$plot

## ---- Panel c: fungal taxa ~ insect PCoA1, family within order (no lure/date) --

res_c_all <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv") %>%
  filter(q_value < q_threshold) %>%
  left_join(fungal_taxonomy, by = "taxon") %>%
  mutate(Order = ifelse(is.na(Order) | Order == "", "Unclassified", Order))

order_rank_c <- res_c_all %>% count(Order, name = "total") %>% arrange(desc(total))
order_levels_c <- setdiff(order_rank_c$Order, "Unclassified")   # dropped -- see panel a's note above; same 31 orders fig4 panel c shows individually, minus this one

subtitle_c <- paste0(sum(res_c_all$Order %in% order_levels_c), " taxa across ", length(order_levels_c),
                      " of the orders from fig4c, by family (q<0.10)")
out_c <- make_family_by_order_breakdown(
  res_c_all, order_levels_c, x_label = "Number of significant taxa",
  title = "Fungal taxa vs. insect PCoA1", subtitle = subtitle_c, tag = "b",
  pos_label = "PCoA1+ / later-season insect community",
  neg_label = "PCoA1- / bark-beetle-dominated insect community",
  legend_title = "Direction"
)
panel_c <- out_c$plot

## ---- Write per-group counts underlying each panel ------------------------------

write.csv(out_b$counts, file.path(out_data_dir, "fig5_family_breakdown_by_order.fungal_date.csv"), row.names = FALSE)
write.csv(out_c$counts, file.path(out_data_dir, "fig5_family_breakdown_by_order.fungal_insectPCoA1.csv"), row.names = FALSE)

## ---- Combine + save -----------------------------------------------------------
# Panel b (15 order facets, 88 family bars) and panel c (31 order facets, 41
# family bars) have different facet_grid row structures, so gtable::cbind
# (which requires matching gtable row counts, used elsewhere in this
# project when panels share a row structure) doesn't apply here --
# gridExtra::arrangeGrob places each panel as an opaque grob instead, so
# panel c's facets simply end up with more vertical breathing room per
# family than panel b's, filling the same total figure height. Tall figure
# -- this is a supplement, meant for close reading / a separate page, not a
# slide.

combined <- gridExtra::arrangeGrob(panel_b, panel_c, ncol = 2, widths = c(1, 1))

ggsave(file.path(out_fig_dir, "fig5_family_breakdown_by_order.png"), combined,
       width = 16, height = 20, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_fig_dir, "fig5_family_breakdown_by_order.pdf"), combined,
       width = 16, height = 20, bg = "white", limitsize = FALSE)

cat("\nDone. Wrote", file.path(out_fig_dir, "fig5_family_breakdown_by_order.png"), "and .pdf\n")
cat("Wrote per-group counts to", out_data_dir, "(fungal_date, fungal_insectPCoA1)\n")
