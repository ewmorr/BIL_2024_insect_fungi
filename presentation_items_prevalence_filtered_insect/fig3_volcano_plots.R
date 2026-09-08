##############################################################################
# Presentation figure 3 -- Prevalence-Filtered Insect Table Variant (Lineage B)
#
# Lineage-B counterpart of presentation_items/fig3_volcano_plots.R. The
# Lineage-A figure was a 1x3 row of volcano plots on the LINEAR date term
# (panels a/b) plus one fungal~insect_PCoA1 panel (c). The Lineage-B direct
# date screens now also carry a QUADRATIC date term and a hump/dip `shape`
# classification (see iterative_analysis_updates.md, 2026-08-31 entries), so
# this figure is a 2x3 grid with rows = functional form of the seasonal test:
#
#   Row 1 (LINEAR term):
#     a  insect abundance ~ collection date, linear term, direct test, all
#        153 prevalence-filtered insect taxa, no pre-selection
#        (compare_insect_fungi/insect_taxa_date_association.prevalence_
#        filtered_insect.r).
#     b  fungal abundance ~ collection date, linear term, every prevalence-
#        filtered fungal taxon, no pre-selection (compare_insect_fungi/
#        fungal_community_seasonality.r -- FUNGAL-ONLY, lineage-invariant, so
#        this is the SAME file the Lineage-A fig3 reads).
#     c  fungal abundance ~ insect_PCoA1, top-200 CoCA-loading candidates,
#        **factor(date)-adjusted** (insect_fungal_coca_factordate_check.
#        prevalence_filtered_insect.r, PCoA1_plus_factordate.csv). factor(date)
#        is a 6-level term absorbing ANY functional form of the season, so
#        this is the strictest available deseasonalizing of the axis -- the
#        DESEASONALIZED insect-community association, not the permissive
#        (date-collinear) version. PCoA1's hit count is stable once any
#        non-linear date flexibility is allowed (9 linear -> 19 quadratic ->
#        17-18 factor(date)); see interpretation_update_prevalence_filtered_
#        insect_08282026.md sec. 5.
#
#   Row 2 (QUADRATIC term / non-monotonic seasonality):
#     d  insect ~ collection date, QUADRATIC term. x-axis sign maps to shape:
#        negative beta_std_quad = concave-down = hump (mid-season peak);
#        positive = dip.
#     e  fungal ~ collection date, QUADRATIC term, same convention.
#     f  fungal abundance ~ insect_PCoA2, top-200 CoCA candidates,
#        **factor(date)-adjusted** (PCoA2_plus_factordate.csv). insect_PCoA2
#        is the axis whose fungal hits first surfaced the non-linear date
#        trend, and it is itself hump-shaped in date (hence its place in
#        row 2). Its factor(date)-adjusted hit count (9/200) is the value
#        carried here, but PCoA2 is NOT stable across date controls -- 6/200
#        unadjusted, 0/200 under a quadratic date term, 9/200 under
#        factor(date) -- so the subtitle reports all three and flags it as
#        small-sample-sensitive (interpretation_update_..._08282026.md sec. 5;
#        iterative_analysis_updates.md CORRECTION entries). The genus/family-
#        level reading (Valsaceae/Cytospora-as-a-group, Tympanis-as-a-group
#        reappearing via different ASVs) is the defensible one, not any single
#        ASV -- and specifically not Cytospora prunicola / ASV_1905 by name.
#
# Points are colored by FDR significance on the term being plotted (q<0.10 on
# q_value; for row 2 panels d/e that is q_value_quad). The few taxa with the
# largest |beta_std| per side (standardized slope of that term -- see the
# 2026-09-08 note below) are labeled -- by Genus/species (or Family) from
# ASVs_taxonomy.tsv for the fungal panels (b/c/e/f), or the already-resolved
# insect Genus/species name for a/d.
#
# REVISION (2026-09-02): after examining example taxa in each linear x
# quadratic combination (figS1_date_shape_bin_examples.R/figS1_date_shape_
# bin_counts.csv), the non-linearity itself wasn't judged particularly
# surprising -- what matters for interpretation is WHEN a taxon peaks: early
# season, late season, mid-season (hump with no significant linear trend),
# or bimodal/early+late (dip with no significant linear trend). This adopts
# a "linear term takes priority" reading (the opposite priority from the
# `shape` column fungal_community_seasonality.r writes, which gives the
# quadratic term priority when both are significant) -- a significant linear
# trend is read as early/late-season regardless of whether the quadratic
# term is ALSO significant, and the quadratic term only gets its own
# "mid-season"/"bimodal" interpretation for taxa with NO significant linear
# trend. Two changes follow for this figure specifically (fig4/fig5/fig7/fig8
# are NOT yet updated to match -- to be revisited):
#   - Panels a/b (linear-term volcanoes): point color no longer marks linear
#     significance (that's what the dashed horizontal line already shows) --
#     it now marks whether the SAME taxon also carries a significant
#     quadratic term, so a reader can see at a glance how much of the
#     linear-significant set (above the dashed line) is "pure" linear vs.
#     linear+quadratic, and how many quadratic-significant taxa sit below
#     the line (linear n.s. -- these are the taxa row 2 focuses on).
#   - Panels d/e (quadratic-term volcanoes): now EXCLUDE any taxon with a
#     significant linear term before plotting -- under the priority rule
#     above, a taxon with a significant linear trend is classified by that
#     trend regardless of its quadratic term, so it no longer belongs in the
#     "hump = mid-season / dip = bimodal" reading these panels are making.
#     Panels d/e are therefore a smaller, cleaner subset than before.
##############################################################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(gtable)
library(grid)
source("library/library.R")

out_fig_dir <- "figures/presentation_items_prevalence_filtered_insect"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

lb_dir <- "data/compare_insects_fungi_top3axes_prevalence_filtered_insect"

q_threshold <- 0.10
n_label_per_side <- 4   # top labeled hits per direction

sig_colors <- c("TRUE" = "#0072B2", "FALSE" = "grey70")   # Okabe-Ito blue vs. neutral grey

# 4-way combo palette for panels a/b (linear x quadratic significance,
# Okabe-Ito): grey = neither term significant, blue = linear-only (matches
# sig_colors' "significant" blue for visual continuity with the other
# panels), orange = quadratic-only, reddish purple = both.
combo_colors <- c("Not significant" = "grey70", "Linear only" = "#0072B2",
                   "Quadratic only" = "#E69F00", "Both significant" = "#CC79A7")
combo_levels <- names(combo_colors)

## ---- Taxonomy lookup for labeling -------------------------------------------

taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv") %>%
  rename(taxon = X) %>%
  mutate(
    Genus = sub("^g__", "", Genus), Species = sub("^s__", "", Species),
    Family = sub("^f__", "", Family)
  )

build_taxon_label <- function(genus, species, family) {
  case_when(
    !is.na(genus) & genus != "" & !is.na(species) & species != "" ~ paste(genus, species),
    !is.na(genus) & genus != "" ~ paste(genus, "sp."),
    !is.na(family) & family != "" ~ paste0(family, " (unclassified)"),
    TRUE ~ "Unclassified"
  )
}
taxonomy$label <- build_taxon_label(taxonomy$Genus, taxonomy$Species, taxonomy$Family)

## ---- Shared plotting logic ------------------------------------------------
# x_col / q_col let the same builder draw the linear-term volcano (beta_std,
# q_value) and the quadratic-term volcano (beta_std_quad, q_value_quad).
#
# The volcano x-axis is a fully-standardized partial slope ("beta weight",
# beta_std = b * sd(x) / sd(y); the QuantPsyc::lm.beta /
# effectsize::standardize_parameters(method = "basic") convention), NOT the
# raw permutation t-statistic -- it is a deterministic property of the
# observed OLS fit, added to the screen CSVs without rerunning the
# permutation tests (the p/q-values are unchanged). Chosen over t or a
# t-derived Cohen's d because the predictor (collection date / insect PCoA
# axis) is continuous, so a standardized regression slope is the meaningful
# effect size, not a two-group standardized mean difference. sign(beta_std)
# == sign(t_stat) throughout, so every direction/shape reading below is
# unchanged.

volcano_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    legend.position = "none",
    plot.tag = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.grid = element_blank()
  )

make_volcano <- function(res, x_label, title, subtitle, tag,
                         x_col = "beta_std", q_col = "q_value",
                         color_col = q_col, color_legend_title = NULL,
                         color_labels = c("FALSE" = "Not significant",
                                          "TRUE" = paste0("Significant (q<", q_threshold, ")")),
                         combo_color = FALSE) {
  res <- res %>%
    mutate(.t = .data[[x_col]], .q = .data[[q_col]], sig = .q < q_threshold,
           .color_sig = .data[[color_col]] < q_threshold)

  # combo_color (panels a/b): 4-way category crossing the plotted term's own
  # significance (`sig`, linear here) with the OTHER term's significance
  # (`.color_sig`, quadratic here) -- see combo_colors above. Non-combo
  # panels (c/d/e/f) keep the original binary TRUE/FALSE coloring on
  # whichever single term color_col points at (defaults to q_col, the term
  # being plotted).
  if (combo_color) {
    res <- res %>% mutate(.plot_color = factor(case_when(
      sig & .color_sig   ~ "Both significant",
      sig & !.color_sig  ~ "Linear only",
      !sig & .color_sig  ~ "Quadratic only",
      TRUE ~ "Not significant"
    ), levels = combo_levels))
  } else {
    res$.plot_color <- res$.color_sig
  }

  # geom_point draws in row order, so the (huge) "Not significant"/FALSE
  # majority can visually bury the sparser significant categories wherever
  # they happen to interleave in the source CSV's row order -- panel b's
  # orange ("Quadratic only") points were invisible under the grey cloud
  # before this fix. Sort so grey is always drawn first (bottom layer) and
  # every significant category draws on top of it, in combo_levels' order
  # (Not significant < Linear only < Quadratic only < Both significant) --
  # the exact level order doesn't matter much among the significant
  # categories, only that none of them sort before "Not significant."
  res <- res %>% arrange(if (combo_color) as.integer(.plot_color) else .plot_color)

  if (!"label" %in% names(res)) {
    # Fungal panels: taxon is an ASV ID -> resolve to a readable name.
    res <- res %>% left_join(taxonomy %>% select(taxon, label), by = "taxon")
  }
  # Insect panels set label = taxon directly (already a resolved name).

  # Rank labeled hits by |beta_std| of the plotted term rather than q-value --
  # with n_perm=999 many taxa tie at the minimum achievable q, which would
  # pick an arbitrary/overlapping cluster of "top" hits. Labeling is still
  # keyed off `sig` (significance of the term being volcano-plotted,
  # x_col/q_col) even when the point COLOR (color_col) marks something else
  # (panels a/b) -- this figure still highlights the top |beta_std| hits on
  # the plotted axis.
  to_label <- bind_rows(
    res %>% filter(sig, .t > 0) %>% arrange(desc(.t)) %>% head(n_label_per_side),
    res %>% filter(sig, .t < 0) %>% arrange(.t) %>% head(n_label_per_side)
  )

  # geom_text_repel only repels labels away from each other and from the
  # x/y position of the rows it's given -- it has no awareness of the
  # geom_point layer's marks. With many taxa tied at the same minimum
  # q-value (n_perm=999), the labeled points often sit in a dense horizontal
  # band of unlabeled ties, so a label can land right on top of nearby
  # points with nothing to repel it away. Pass every significant point as a
  # repel "obstacle" (blank label for the ones we're not annotating) so
  # their empty boxes still occupy space and push the real labels clear.
  label_obstacles <- res %>%
    filter(sig) %>%
    mutate(label = if_else(taxon %in% to_label$taxon, label, ""))

  p <- ggplot(res, aes(x = .t, y = -log10(.q))) +
    geom_hline(yintercept = -log10(q_threshold), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey85") +
    geom_point(aes(color = .plot_color), size = 1.6, alpha = 0.75) +
    ggrepel::geom_text_repel(
      data = label_obstacles, aes(label = label), size = 3.3, fontface = "italic",
      max.overlaps = Inf, segment.size = 0.25, segment.color = "grey50",
      min.segment.length = 0, box.padding = 0.5, force = 12, force_pull = 0.15,
      max.time = 3, max.iter = 20000, direction = "both",
      nudge_y = 0.2
    ) +
    { if (combo_color) scale_color_manual(values = combo_colors, name = color_legend_title, drop = FALSE)
      else scale_color_manual(values = sig_colors, name = color_legend_title, labels = color_labels) } +
    scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.3))) +
    labs(x = x_label, y = expression(-log[10](italic(q)*"-value")),
         title = title, subtitle = subtitle, tag = tag) +
    volcano_theme

  if (is.null(color_legend_title)) {
    p <- p + guides(color = "none")
  } else {
    p <- p + theme(legend.position = "bottom", legend.title = element_text(size = 9, face = "bold"),
                    legend.text = element_text(size = 8)) +
      guides(color = guide_legend(nrow = 1))
  }
  p
}

## ---- Row 1, panel a: insect ~ date, LINEAR term (153 taxa, unbiased) -----
# Color is now the 4-way linear x quadratic combo (combo_color = TRUE): grey
# = neither term significant, blue = linear only, orange = quadratic only,
# purple = both. The dashed horizontal line still shows linear significance
# via y-position -- points below it colored "Quadratic only" are exactly the
# taxa row 2 (panels d/e) focuses on; under the priority rule (see header) a
# taxon's peak-timing call is "early/late season" whenever linear is
# significant, quadratic or not, so "Linear only" and "Both significant"
# get the SAME peak-timing reading despite the different color.

res_ad <- read.csv(file.path(lb_dir, "insect_taxa_date_association.csv")) %>%
  mutate(label = taxon)   # resolved Genus/species names already

a_combo_n <- c(
  neither = sum(res_ad$q_value >= q_threshold & res_ad$q_value_quad >= q_threshold),
  lin_only = sum(res_ad$q_value < q_threshold & res_ad$q_value_quad >= q_threshold),
  quad_only = sum(res_ad$q_value >= q_threshold & res_ad$q_value_quad < q_threshold),
  both = sum(res_ad$q_value < q_threshold & res_ad$q_value_quad < q_threshold)
)
subtitle_a <- paste0("n=", nrow(res_ad), " taxa (unbiased): ", a_combo_n["neither"], " neither, ",
                     a_combo_n["lin_only"], " linear only, ", a_combo_n["quad_only"], " quadratic only,\n",
                     a_combo_n["both"], " both (q<0.10)")
panel_a <- make_volcano(res_ad,
                        expression(paste("standardized slope (", beta, ", insect ~ collection date, linear term)")),
                        "Insect taxa vs. date (linear)", subtitle_a, "a",
                        color_col = "q_value_quad", color_legend_title = "Term significance",
                        combo_color = TRUE)

## ---- Row 1, panel b: fungal ~ date, LINEAR term (all taxa; lineage-invariant) ----
# Same 4-way combo coloring as panel a, but no legend of its own -- panel a's
# single-row legend (same categories/colors, shared color scale) covers both.

res_be <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv")
b_combo_n <- c(
  neither = sum(res_be$q_value >= q_threshold & res_be$q_value_quad >= q_threshold),
  lin_only = sum(res_be$q_value < q_threshold & res_be$q_value_quad >= q_threshold),
  quad_only = sum(res_be$q_value >= q_threshold & res_be$q_value_quad < q_threshold),
  both = sum(res_be$q_value < q_threshold & res_be$q_value_quad < q_threshold)
)
subtitle_b <- paste0("n=", nrow(res_be), " prevalence-filtered taxa (unbiased): ", b_combo_n["neither"],
                     " neither, ", b_combo_n["lin_only"], " linear only, ", b_combo_n["quad_only"],
                     " quadratic only,\n", b_combo_n["both"], " both (q<0.10)")
panel_b <- make_volcano(res_be,
                        expression(paste("standardized slope (", beta, ", fungal ~ collection date, linear term)")),
                        "Fungal taxa vs. date (linear)", subtitle_b, "b",
                        color_col = "q_value_quad", combo_color = TRUE)

## ---- Row 1, panel c: fungal ~ insect_PCoA1, factor(date)-adjusted (CoCA candidates) ----

res_c <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA1_plus_factordate.csv"))
subtitle_c <- paste0("n=", nrow(res_c), " CoCA candidates, ",
                     sum(res_c$q_value < q_threshold),
                     " sig (q<0.10), factor(date)-adjusted")
panel_c <- make_volcano(res_c,
                        expression(paste("standardized slope (", beta, ", fungal ~ insect PCoA1 | factor(date))")),
                        "Fungal taxa vs. insect community (PCoA1)", subtitle_c, "c")

## ---- Row 2, panel d: insect ~ date, QUADRATIC term, non-linear-trending taxa only ----
# Under the linear-priority peak-timing rule (see header), a taxon with a
# significant linear term is called early/late-season regardless of its
# quadratic term, so it's excluded here -- hump = mid-season peak / dip =
# bimodal (early+late) peak is only the right reading for taxa with no
# significant linear trend in the first place.

res_d <- res_ad %>% filter(q_value >= q_threshold)
subtitle_d <- paste0("n=", nrow(res_d), " taxa with NO significant linear term (of ", nrow(res_ad),
                     " total, ", nrow(res_ad) - nrow(res_d), " excluded),\n",
                     sum(res_d$q_value_quad < q_threshold),
                     " significant (q<0.10) -- mid-season peak (hump) or bimodal (dip)")
quad_x_label <- expression(paste("quadratic-term standardized slope ", beta,
                                 "   (", beta, " < 0: hump = mid-season peak,  ",
                                 beta, " > 0: dip = bimodal)"))
panel_d <- make_volcano(res_d, quad_x_label,
                        "Insect taxa vs. date (quadratic, non-linear-trending taxa)", subtitle_d, "d",
                        x_col = "beta_std_quad", q_col = "q_value_quad")

## ---- Row 2, panel e: fungal ~ date, QUADRATIC term, non-linear-trending taxa only ----

res_e <- res_be %>% filter(q_value >= q_threshold)
subtitle_e <- paste0("n=", nrow(res_e), " taxa with NO significant linear term (of ", nrow(res_be),
                     " total, ", nrow(res_be) - nrow(res_e), " excluded),\n",
                     sum(res_e$q_value_quad < q_threshold),
                     " significant (q<0.10) -- mid-season peak (hump) or bimodal (dip)")
panel_e <- make_volcano(res_e, quad_x_label,
                        "Fungal taxa vs. date (quadratic, non-linear-trending taxa)", subtitle_e, "e",
                        x_col = "beta_std_quad", q_col = "q_value_quad")

## ---- Row 2, panel f: fungal ~ insect_PCoA2, factor(date)-adjusted (CoCA candidates) ----
# Headline count is the factor(date)-adjusted one, matching panel c. But PCoA2
# is unstable across date controls, so the subtitle also reports the
# unadjusted and quadratic-date counts and flags the instability (see header).

res_f <- read.csv(file.path(lb_dir, "fungal_insect_association_results.PCoA2_plus_factordate.csv"))
n_f_unadj <- sum(read.csv(file.path(lb_dir,
  "fungal_insect_association_results.PCoA2_only.csv"))$q_value < q_threshold)
n_f_quad <- sum(read.csv(file.path(lb_dir,
  "fungal_insect_association_results.PCoA2_plus_quadratic_date.csv"))$q_value < q_threshold)
subtitle_f <- paste0(
  sum(res_f$q_value < q_threshold), "/", nrow(res_f),
  " sig (q<0.10), factor(date)-adjusted -- but UNSTABLE across date controls:\n",
  n_f_unadj, "/", nrow(res_f), " unadjusted, ", n_f_quad, "/", nrow(res_f),
  " quadratic-date; small-sample-sensitive, see text")
panel_f <- make_volcano(res_f,
                        expression(paste("standardized slope (", beta, ", fungal ~ insect PCoA2 | factor(date))")),
                        "Fungal taxa vs. insect community (PCoA2)", subtitle_f, "f")

## ---- Combine + save -----------------------------------------------------
# Align panel axes across the 2x3 grid the same way fig1/fig2 do, rather than
# gridExtra::arrangeGrob's independent null-unit sizing: build each COLUMN by
# rbind-ing its two panels (size="max" -> both get the wider column's
# y-axis-title/y-axis-text/panel widths, so the panel boxes and x-axis tick
# labels line up top-to-bottom), then cbind the three columns (size="max" ->
# corresponding rows -- title, subtitle, panel, x-axis -- are matched in
# height left-to-right, so panel f's 2-line subtitle doesn't shove its panel
# box out of line with c's). rbind/cbind dispatch to gtable's S3 methods once
# library(gtable) is loaded.

col1 <- rbind(ggplotGrob(panel_a), ggplotGrob(panel_d), size = "max")
col2 <- rbind(ggplotGrob(panel_b), ggplotGrob(panel_e), size = "max")
col3 <- rbind(ggplotGrob(panel_c), ggplotGrob(panel_f), size = "max")
combined <- cbind(col1, col2, col3, size = "max")

# ggrepel computes label positions per-device using that device's font
# metrics, so labels can overlap in one output format but not another from
# the same plot object (cairo_pdf would fix this but is unavailable here --
# no X11/cairo libs). Render once to PNG, grab the fully-drawn (post-repel)
# grob, and replay that frozen layout to the PDF device so both formats show
# identical label placement.
png_path <- file.path(out_fig_dir, "fig3_volcano_plots.png")
pdf_path <- file.path(out_fig_dir, "fig3_volcano_plots.pdf")

grDevices::png(png_path, width = 16.5, height = 10.4, units = "in", res = 300, bg = "white")
grid::grid.draw(combined)
frozen <- grid::grid.grab()
grDevices::dev.off()

grDevices::pdf(pdf_path, width = 16.5, height = 10.4, bg = "white")
grid::grid.draw(frozen)
grDevices::dev.off()

cat("\nDone. Wrote", png_path, "and .pdf\n")
