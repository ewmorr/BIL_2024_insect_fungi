##############################################################################
# Presentation figure 3: volcano plots for the three main per-taxon
# insect/fungal screens.
#
# Panel a: insect abundance ~ collection date, direct test, run on all 52
#          insect taxa (Curculionidae + Latridiidae) with no pre-selection --
#          compare_insect_fungi/insect_fungal_coca_analysis.general_
#          workflow.r, file data/2024_insect_data/insect_taxa_date_
#          association.csv. Confirms the spring bark/ambrosia-beetle vs.
#          later-season fungivorous-beetle turnover described in
#          interpretation.md.
# Panel b: fungal abundance ~ collection date, direct test, run on every
#          prevalence-filtered fungal taxon with NO pre-selection (unlike
#          panel c's CoCA-loading screen) -- compare_insect_fungi/
#          fungal_community_seasonality.r, file data/2024_fungi/
#          fungal_taxa_date_association.all_taxa.csv. The unbiased estimate
#          of how much of the fungal community is seasonal (~33%).
# Panel c: fungal abundance ~ insect_PCoA1, no lure/date covariates -- the
#          top-200 CoCA-loading candidate taxa (compare_insect_fungi/
#          insect_fungal_coca_analysis.general_workflow.r), file
#          data/compare_insects_fungi_top3axes/fungal_insect_association_
#          results.PCoA1_only.no_lure.csv. Answers "does fungal abundance
#          track insect community composition" before any adjustment for the
#          shared seasonal trend (insect_PCoA1 is ~78% collinear with date --
#          see interpretation.md -- so this is the more permissive test;
#          the date-adjusted version collapses to 0/200 significant).
#
# Points are colored by FDR significance (q<0.10); the several
# most-significant taxa on each side (positive/negative t-statistic) are
# labeled -- by Genus/species (or Family, if unclassified at genus level)
# from ASVs_taxonomy.tsv for the fungal panels (b/c), or directly by the
# already-resolved insect Genus/species name for panel a.
##############################################################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
source("library/library.R")

out_fig_dir <- "figures/presentation_items"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

q_threshold <- 0.10
n_label_per_side <- 4   # top labeled hits per direction (positive/negative t-stat)

sig_colors <- c("TRUE" = "#0072B2", "FALSE" = "grey70")   # Okabe-Ito blue vs. neutral grey

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

## ---- Shared plotting logic ---------------------------------------------------

volcano_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11),
    legend.position = "none",
    plot.tag = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    panel.grid = element_blank()
  )

make_volcano <- function(res, x_label, title, subtitle, tag) {
  res <- res %>% mutate(sig = q_value < q_threshold)
  if (!"label" %in% names(res)) {
    # Fungal panels: taxon is an ASV ID, resolve to a readable Genus/species
    # (or Family) via the taxonomy lookup built above.
    res <- res %>% left_join(taxonomy %>% select(taxon, label), by = "taxon")
  }
  # Insect panels: taxon is already a resolved Genus/species name (see
  # panel a below, which sets label = taxon directly, no join needed).

  # Rank labeled hits by |t-statistic| (effect size) rather than q-value --
  # with n_perm=999 many taxa tie at the minimum achievable q-value, which
  # would otherwise pick an arbitrary/overlapping cluster of "top" hits.
  to_label <- bind_rows(
    res %>% filter(sig, t_stat > 0) %>% arrange(desc(t_stat)) %>% head(n_label_per_side),
    res %>% filter(sig, t_stat < 0) %>% arrange(t_stat) %>% head(n_label_per_side)
  )

  ggplot(res, aes(x = t_stat, y = -log10(q_value))) +
    geom_hline(yintercept = -log10(q_threshold), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey85") +
    geom_point(aes(color = sig), size = 1.6, alpha = 0.75) +
    ggrepel::geom_text_repel(
      data = to_label, aes(label = label), size = 3.3, fontface = "italic",
      max.overlaps = Inf, segment.size = 0.25, segment.color = "grey50",
      min.segment.length = 0, box.padding = 0.5, force = 12, force_pull = 0.15,
      max.time = 3, max.iter = 20000, direction = "both",
      nudge_y = 0.2
    ) +
    scale_color_manual(values = sig_colors) +
    scale_x_continuous(expand = expansion(mult = c(0.14, 0.14))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.3))) +
    labs(x = x_label, y = expression(-log[10](italic(q)*"-value")),
         title = title, subtitle = subtitle, tag = tag) +
    volcano_theme
}

## ---- Panel a: insect ~ collection date (all 52 taxa, unbiased) -------------------

res_a <- read.csv("data/2024_insect_data/insect_taxa_date_association.csv") %>%
  mutate(label = taxon)   # already a resolved Genus/species name, no taxonomy join needed
subtitle_a <- paste0("n=", nrow(res_a), " taxa (unbiased), ",
                      sum(res_a$q_value < q_threshold), " significant (q<0.10)")
panel_a <- make_volcano(res_a, "t-statistic (insect abundance ~ collection date)",
                         "Insect taxa vs. collection date", subtitle_a, "a")

## ---- Panel b: fungal ~ collection date (all prevalence-filtered taxa) ------------

res_b <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv")
subtitle_b <- paste0("n=", nrow(res_b), " prevalence-filtered taxa (unbiased), ",
                      sum(res_b$q_value < q_threshold), " significant (q<0.10)")
panel_b <- make_volcano(res_b, "t-statistic (fungal abundance ~ collection date)",
                         "Fungal taxa vs. collection date", subtitle_b, "b")

## ---- Panel c: fungal ~ insect PCoA1 (no lure, no date; CoCA-selected candidates) ----

res_c <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv")
subtitle_c <- paste0("n=", nrow(res_c), " CoCA-selected candidates, ",
                      sum(res_c$q_value < q_threshold), " significant (q<0.10)")
panel_c <- make_volcano(res_c, "t-statistic (fungal abundance ~ insect PCoA1)",
                         "Fungal taxa vs. insect community", subtitle_c, "c")

## ---- Combine + save -----------------------------------------------------------

combined <- gridExtra::arrangeGrob(panel_a, panel_b, panel_c, ncol = 3)

# ggrepel computes label positions per-device using that device's font
# metrics, which can make labels overlap in one output format but not
# another even from the same plot object (cairo_pdf would normally fix
# this, but it's unavailable on this machine -- no X11/cairo libs). Instead,
# render once to the PNG device, grab the fully-drawn (post-repel) grob,
# and replay that exact frozen layout to the PDF device so both formats
# show byte-identical label placement.
png_path <- file.path(out_fig_dir, "fig3_volcano_plots.png")
pdf_path <- file.path(out_fig_dir, "fig3_volcano_plots.pdf")

grDevices::png(png_path, width = 16.5, height = 5.2, units = "in", res = 300, bg = "white")
grid::grid.draw(combined)
frozen <- grid::grid.grab()
grDevices::dev.off()

grDevices::pdf(pdf_path, width = 16.5, height = 5.2, bg = "white")
grid::grid.draw(frozen)
grDevices::dev.off()

cat("\nDone. Wrote", file.path(out_fig_dir, "fig3_volcano_plots.png"), "and .pdf\n")
