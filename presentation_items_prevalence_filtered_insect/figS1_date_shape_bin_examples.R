##############################################################################
# Supplemental figure S1 -- fungal date-shape bin examples
# (Prevalence-Filtered Insect Table Variant / Lineage B presentation set --
# see README.md in this directory. Underlying data is fungal-only and
# lineage-invariant, but this figure exists to support the fig3/fig4/fig5
# linear+quadratic redesign discussion, so it's filed here as a supplemental
# item alongside fig1-fig8 rather than in compare_insect_fungi/.)
#
# Diagnostic/illustrative companion to fungal_community_seasonality.r's
# linear+quadratic date screen (fungal_taxa_date_association.all_taxa.csv,
# lineage-invariant -- read by both Lineage A and Lineage B's fig4/fig5).
# Prompted by a question while reviewing lineage_B_prevalence_filtered_
# insect_top_level_interpretation.md sec.7: the fig4-panel-b hit-instance
# counts (e.g. Dothideomycetes' 401 linear-later + 138 quadratic-hump) are
# NOT mutually exclusive -- a taxon can be significant in both the linear
# and quadratic terms at once (see fig4's make_combined_breakdown() and its
# "a taxon can contribute up to 2 hit-instances" comment). This script picks
# a couple of example taxa from each linear x quadratic combination and
# plots their actual CLR abundance vs. date, to make the distinction --
# and what a "later-season but also hump-shaped" taxon literally looks like
# -- visually concrete. It also writes out the bin-count crosstab itself as
# a small supplemental table (figS1_date_shape_bin_counts.csv), since that
# was the number that prompted the fig3 redesign (see fig3_volcano_plots.R's
# 2026-09-02 revision note and iterative_analysis_updates.md).
#
# Bins (q<0.10 threshold, matching project convention): the full 3x3 cross
# of linear (later / earlier / not sig) x quadratic (hump / dip / not sig),
# minus the "neither sig" cell (that's just "no significant date pattern",
# not a shape to illustrate) -- 8 bins total:
#   - late_only        linear later-season (t>0) sig, quadratic NOT sig
#   - early_only        linear earlier-season (t<0) sig, quadratic NOT sig
#   - late_and_hump      linear later-season (t>0) sig AND quadratic hump (t_quad<0) sig
#   - late_and_dip       linear later-season (t>0) sig AND quadratic dip (t_quad>0) sig
#   - early_and_hump     linear earlier-season (t<0) sig AND quadratic hump (t_quad<0) sig
#   - early_and_dip      linear earlier-season (t<0) sig AND quadratic dip (t_quad>0) sig
#   - hump_only          quadratic hump (t_quad<0) sig, linear NOT sig
#   - dip_only           quadratic dip (t_quad>0) sig, linear NOT sig
#
# NOTE: this is the raw linear x quadratic significance crosstab, kept for
# reference -- it is NOT the same as the "peak timing" classification
# (early/late/mid-season/bimodal) adopted for fig3 onward on 2026-09-02,
# which gives the LINEAR term priority whenever it's significant (see
# fig3_volcano_plots.R header). The two disagree exactly on the
# late_and_dip/early_and_hump/late_and_hump/early_and_dip bins above, where
# this table still lets the quadratic sign be read as an independent
# "shape" -- that's deliberate here (it's what motivated the redesign), but
# don't use these 8 bins as if they were the adopted peak-timing labels.
#
# Rebuilds the same matched 69-sample fungal dataset + prevalence filter +
# CLR transform that fungal_community_seasonality.r used to produce the
# association CSV (identical code), so the plotted points are the exact
# values the permutation test was run on -- checked against the CSV's
# taxon list via stopifnot() below. The overlaid trend line is a plain
# y ~ poly(date, 2) OLS fit (site/lure NOT included) purely for visual
# guidance, matching the linear+quadratic date_c + I(date_c^2) structure
# of the tested model -- descriptive only, like fig9's cubic trend lines;
# the permutation test (which DOES include site+lure) is what carries the
# actual significance call.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
source("library/library.R")

out_fig_dir <- "figures/presentation_items_prevalence_filtered_insect"
out_data_dir <- "data/presentation_items_prevalence_filtered_insect"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1)

## ---- Rebuild the matched fungal dataset (identical steps to fungal_community_seasonality.r) ----

sp_tab <- read.csv("data/2024_insect_data/insect_species_tab.csv")
insect_meta_full <- read.csv("data/metadata/insect_community_metadata.csv")

sp_tab %>% filter(Family %in% c("Curculionidae", "Latridiidae")) -> sp_tab.curcus_latris
sp_tab.curcus_latris.t <- t(sp_tab.curcus_latris %>% select(where(is.numeric)))
colnames(sp_tab.curcus_latris.t) <- sp_tab.curcus_latris$Finest.ID
sp_tab.curcus_latris.t[, colSums(sp_tab.curcus_latris.t) > 1] -> insect_full
insect_full <- insect_full[rowSums(insect_full) > 0, ]

id_map <- insect_meta_full %>% filter(col_names %in% rownames(insect_full))
stopifnot(!any(duplicated(id_map$col_names)))

fungal_full <- read.csv("data/2024_fungi/ASV_tab.csv", row.names = 1)
shared_ids <- intersect(id_map$SequenceID, rownames(fungal_full))
id_map <- id_map %>% filter(SequenceID %in% shared_ids)
fungal <- fungal_full[id_map$SequenceID, , drop = FALSE]

fungal_taxonomy_raw <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv", row.names = 1, check.names = FALSE)
is_fungus <- !is.na(fungal_taxonomy_raw[colnames(fungal), "Kingdom"]) &
  fungal_taxonomy_raw[colnames(fungal), "Kingdom"] == "k__Fungi"
fungal <- fungal[, is_fungus, drop = FALSE]

meta <- id_map %>%
  transmute(sample_id = SequenceID, site = Site, lure = Lure, trap_id = trapID,
            date = lubridate::mdy(CollectionDate))
stopifnot(all(rownames(fungal) == meta$sample_id))

keep <- colSums(fungal > 0) >= 5              # same prevalence filter as fungal_community_seasonality.r
fungal_f <- fungal[, keep, drop = FALSE]

clr_transform <- function(mat, pseudocount = 1) {
  mat <- as.matrix(mat) + pseudocount
  props <- sweep(mat, 1, rowSums(mat), "/")
  log_props <- log(props)
  sweep(log_props, 1, rowMeans(log_props), "-")
}
fungal_clr <- clr_transform(fungal_f)
stopifnot(all(rownames(fungal_clr) == meta$sample_id))

cat("Rebuilt matched dataset:", nrow(fungal_clr), "samples,", ncol(fungal_clr),
    "prevalence-filtered fungal taxa (should match fungal_community_seasonality.r's screen).\n")

## ---- Load the date-association results + assign each taxon to a bin --------

q_threshold <- 0.10
res <- read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv")
stopifnot(setequal(res$taxon, colnames(fungal_clr)))  # confirms we rebuilt the exact same tested taxon set

res <- res %>%
  mutate(
    sig_lin = q_value < q_threshold,
    sig_quad = q_value_quad < q_threshold,
    bin = case_when(
      sig_lin & t_stat > 0 & !sig_quad                  ~ "late_only",
      sig_lin & t_stat < 0 & !sig_quad                  ~ "early_only",
      sig_lin & t_stat > 0 & sig_quad & t_stat_quad < 0 ~ "late_and_hump",
      sig_lin & t_stat > 0 & sig_quad & t_stat_quad > 0 ~ "late_and_dip",
      sig_lin & t_stat < 0 & sig_quad & t_stat_quad < 0 ~ "early_and_hump",
      sig_lin & t_stat < 0 & sig_quad & t_stat_quad > 0 ~ "early_and_dip",
      !sig_lin & sig_quad & t_stat_quad < 0             ~ "hump_only",
      !sig_lin & sig_quad & t_stat_quad > 0             ~ "dip_only",
      TRUE ~ NA_character_
    )
  )

bin_levels <- c("late_only", "early_only", "late_and_hump", "late_and_dip",
                 "early_and_hump", "early_and_dip", "hump_only", "dip_only")
bin_labels <- c(
  late_only      = "Later season only\n(linear t>0, quad n.s.)",
  early_only     = "Earlier season only\n(linear t<0, quad n.s.)",
  late_and_hump  = "Later season + hump\n(linear t>0, quad hump)",
  late_and_dip   = "Later season + dip\n(linear t>0, quad dip)",
  early_and_hump = "Earlier season + hump\n(linear t<0, quad hump)",
  early_and_dip  = "Earlier season + dip\n(linear t<0, quad dip)",
  hump_only      = "Hump only, no linear trend\n(linear n.s., quad hump)",
  dip_only       = "Dip only, no linear trend\n(linear n.s., quad dip)"
)

bin_counts <- as.data.frame(table(bin = factor(res$bin, levels = bin_levels), useNA = "ifany"))
names(bin_counts) <- c("bin", "n_taxa")
bin_counts$bin <- as.character(bin_counts$bin)
bin_counts$bin[is.na(bin_counts$bin) | bin_counts$bin == ""] <- "no_significant_pattern"
bin_counts$label <- c(gsub("\n", " ", bin_labels[bin_levels]), "No significant linear or quadratic term")[
  match(bin_counts$bin, c(bin_levels, "no_significant_pattern"))]
bin_counts <- bin_counts[c("bin", "label", "n_taxa")]

cat("\nBin sizes (of", nrow(res), "prevalence-filtered fungal taxa):\n")
print(bin_counts)
write.csv(bin_counts, file.path(out_data_dir, "figS1_date_shape_bin_counts.csv"), row.names = FALSE)

## ---- Pick n_examples per bin, strongest by whichever term(s) define the bin ----
# "Strongest" = smallest permutation p-value on the defining term(s) (min
# possible p_perm given 999 permutations is 1/1000, so several taxa may tie
# at the floor within a bin -- ties broken arbitrarily by row order, which
# is fine for illustrative examples).

n_examples <- 2

res <- res %>%
  mutate(strength = case_when(
    bin %in% c("late_and_hump", "late_and_dip", "early_and_hump", "early_and_dip") ~
      pmax(-log10(p_perm), -log10(p_perm_quad)),
    bin %in% c("hump_only", "dip_only")   ~ -log10(p_perm_quad),
    bin %in% c("late_only", "early_only") ~ -log10(p_perm),
    TRUE ~ NA_real_
  ))

examples <- res %>%
  filter(!is.na(bin)) %>%
  arrange(desc(strength)) %>%
  group_by(bin) %>%
  slice_head(n = n_examples) %>%
  mutate(example_rank = row_number()) %>%
  ungroup()

n_found <- table(factor(examples$bin, levels = bin_levels))
if (any(n_found < n_examples)) {
  cat("\nNote: fewer than", n_examples, "example(s) available for:",
      paste(names(n_found)[n_found < n_examples], collapse = ", "),
      "-- see bin sizes above (that combination may be genuinely rare/absent).\n")
}
if (nrow(examples) == 0) stop("No taxa matched any of the 8 bins -- nothing to plot.")

## ---- Taxon labels (Genus species where resolved, else Family/Order + ASV id) ----

taxon_label <- function(tax) {
  tx <- fungal_taxonomy_raw[tax, ]
  genus <- sub("^g__", "", tx[["Genus"]]);   species <- sub("^s__", "", tx[["Species"]])
  family <- sub("^f__", "", tx[["Family"]]); order   <- sub("^o__", "", tx[["Order"]])
  if (!is.na(genus) && genus != "" && !is.na(species) && species != "" && species != "sp") {
    paste0(genus, " ", species, " (", tax, ")")
  } else if (!is.na(genus) && genus != "") {
    paste0(genus, " sp. (", tax, ")")
  } else if (!is.na(family) && family != "") {
    paste0(family, " (", tax, ")")
  } else if (!is.na(order) && order != "") {
    paste0(order, " (", tax, ")")
  } else {
    tax
  }
}
examples$label <- vapply(examples$taxon, taxon_label, character(1))
examples$stat_text <- sprintf("linear t=%.2f, q=%.3f | quadratic t=%.2f, q=%.3f",
                               examples$t_stat, examples$q_value, examples$t_stat_quad, examples$q_value_quad)

examples <- examples %>%
  mutate(bin = factor(bin, levels = bin_levels, labels = bin_labels),
         example_rank = factor(paste("Example", example_rank), levels = paste("Example", seq_len(n_examples))))

## ---- Assemble the long plotting data: CLR abundance vs. date per example taxon ----

site_levels <- c("Durham", "Pease Airport", "Manchester Cedar Swamp", "Manchester Airport")

plot_dat <- lapply(seq_len(nrow(examples)), function(i) {
  ex <- examples[i, ]
  data.frame(
    bin = ex$bin, example_rank = ex$example_rank, taxon = ex$taxon,
    sample_id = meta$sample_id, site = meta$site, date = meta$date,
    clr_abund = fungal_clr[meta$sample_id, ex$taxon]
  )
}) %>% bind_rows() %>%
  mutate(site = ifelse(site == "Pease", "Pease Airport", site),
         site = factor(site, levels = site_levels))

site_palette <- setNames(twelvePaired[c(2, 4, 6, 8)], site_levels)

## ---- Plot: bin x example_rank grid, one panel per example taxon -------------
# facet_grid keeps a clean bin (row) x example-rank (col) grid even when a
# bin has fewer than n_examples taxa (that cell is simply blank). Per-panel
# taxon identity + test stats can't live in the shared row/col strips (the
# taxon differs between the two example columns of the same row), so
# they're placed as in-panel text at the top-left corner instead -- the
# same annotate(x=-Inf,...) technique fig9 uses for its stress-value label,
# which still positions correctly per panel under facet scales="free_y".

p <- ggplot(plot_dat, aes(x = date, y = clr_abund)) +
  geom_point(aes(fill = site), shape = 21, size = 2.2, color = "black", stroke = 0.3) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color = "grey30", linewidth = 0.6) +
  geom_text(data = examples, aes(x = -Inf, y = Inf, label = paste0(label, "\n", stat_text)),
            inherit.aes = FALSE, hjust = -0.02, vjust = 1.2, size = 2.5, lineheight = 0.9) +
  facet_grid(rows = vars(bin), cols = vars(example_rank), scales = "free_y", switch = "y") +
  scale_fill_manual(values = site_palette, name = "Site") +
  scale_x_date(date_labels = "%b %d") +
  labs(x = "Collection date", y = NULL,
       title = "Example fungal taxa by linear x quadratic date-term bin",
       subtitle = paste0("q<", q_threshold, " (BH-FDR) on fungal_taxa_date_association.all_taxa.csv; ",
                          "quadratic OLS trend line is descriptive only -- significance is from the\n",
                          "permutation test (site+lure+date_c+I(date_c^2))")) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
    strip.placement = "outside",
    axis.text = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    plot.subtitle = element_text(size = 8.5, color = "grey30")
  )

ggsave(file.path(out_fig_dir, "figS1_date_shape_bin_examples.png"), p,
       width = 11, height = 22, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_fig_dir, "figS1_date_shape_bin_examples.pdf"), p,
       width = 11, height = 22, bg = "white", limitsize = FALSE)

cat("\nDone. Wrote", file.path(out_fig_dir, "figS1_date_shape_bin_examples.png"), "and .pdf,",
    "and", file.path(out_data_dir, "figS1_date_shape_bin_counts.csv"), "\n")
