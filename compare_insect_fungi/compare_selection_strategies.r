##############################################################################
# Compare fungal-taxon candidate-selection strategies
#
# Three ways of picking a "top-loading" set of candidate fungal taxa to test
# against the insect community have now been run:
#   1. top3axes_CoCA -- top 200 by co-correspondence-analysis (CoCA) loading
#      strength across insect CoCA axes 1-3 jointly (insect_fungal_coca_
#      analysis.general_workflow.r, results in
#      data/compare_insects_fungi_top3axes/)
#   2. PCoA1_date     -- top 200 by an unbiased per-taxon test of association
#      with collection date, the single axis (insect_PCoA1) found to be
#      overwhelmingly a seasonal gradient (insect_fungal_coca_analysis.
#      PCoA1_date_axis.r, results in data/compare_insects_fungi_PCoA1_date/)
#   3. PCoA3_lure     -- top 200 by an unbiased per-taxon test of association
#      with lure, the single axis (insect_PCoA3) found to be by far the
#      strongest lure-associated insect ordination axis (insect_fungal_coca_
#      analysis.PCoA3_lure_axis.r, results in
#      data/compare_insects_fungi_PCoA3_lure/)
#
# This script asks how much these three candidate-selection strategies agree
# -- both at the level of WHICH taxa get selected as candidates, and at the
# level of WHICH candidates come out significant in each strategy's targeted
# downstream test.
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

out_data_dir <- "data/compare_selection_strategies"
out_fig_dir <- "figures/compare_selection_strategies"
dir.create(out_data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load the three candidate sets (top 200 taxa each) -----------------

top3axes_candidates <- read.csv("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.csv")$taxon
pcoa1_date_candidates <- head(read.csv("data/2024_fungi/fungal_taxa_date_association.all_taxa.csv")$taxon, 200)
pcoa3_lure_candidates <- head(read.csv("data/2024_fungi/fungal_taxa_lure_association.all_taxa.csv")$taxon, 200)

stopifnot(length(top3axes_candidates) == 200, length(pcoa1_date_candidates) == 200, length(pcoa3_lure_candidates) == 200)

candidate_sets <- list(top3axes_CoCA = top3axes_candidates,
                        PCoA1_date = pcoa1_date_candidates,
                        PCoA3_lure = pcoa3_lure_candidates)

## ---- 2. Taxonomy for annotation -------------------------------------------

taxonomy <- read.delim("data/2024_fungi/ASVs_taxonomy.tsv") %>% rename(taxon = X)

## ---- 3. Candidate-set membership + overlap --------------------------------

all_taxa <- Reduce(union, candidate_sets)
membership <- data.frame(taxon = all_taxa) %>%
  mutate(in_top3axes_CoCA = taxon %in% candidate_sets$top3axes_CoCA,
         in_PCoA1_date = taxon %in% candidate_sets$PCoA1_date,
         in_PCoA3_lure = taxon %in% candidate_sets$PCoA3_lure) %>%
  left_join(taxonomy, by = "taxon") %>%
  arrange(desc(in_top3axes_CoCA + in_PCoA1_date + in_PCoA3_lure))
write.csv(membership, file.path(out_data_dir, "candidate_taxa_membership.csv"), row.names = FALSE)

pair_overlap <- function(a, b) {
  ia <- length(intersect(a, b))
  list(n_a = length(a), n_b = length(b), n_intersect = ia,
       jaccard = ia / length(union(a, b)))
}

overlap_summary <- bind_rows(
  data.frame(comparison = "top3axes_CoCA vs PCoA1_date", pair_overlap(candidate_sets$top3axes_CoCA, candidate_sets$PCoA1_date)),
  data.frame(comparison = "top3axes_CoCA vs PCoA3_lure", pair_overlap(candidate_sets$top3axes_CoCA, candidate_sets$PCoA3_lure)),
  data.frame(comparison = "PCoA1_date vs PCoA3_lure", pair_overlap(candidate_sets$PCoA1_date, candidate_sets$PCoA3_lure)),
  data.frame(comparison = "all three sets", n_a = NA, n_b = NA,
             n_intersect = length(Reduce(intersect, candidate_sets)),
             jaccard = length(Reduce(intersect, candidate_sets)) / length(all_taxa))
)
cat("\n--- Candidate-set overlap (top 200 taxa per strategy) ---\n")
print(overlap_summary, row.names = FALSE)
write.csv(overlap_summary, file.path(out_data_dir, "candidate_set_overlap_summary.csv"), row.names = FALSE)

## ---- 4. UpSet-style region plot (ggplot2 only, no extra deps) -------------

region_counts <- membership %>%
  count(in_top3axes_CoCA, in_PCoA1_date, in_PCoA3_lure, name = "n") %>%
  filter(in_top3axes_CoCA | in_PCoA1_date | in_PCoA3_lure) %>%
  rowwise() %>%
  mutate(region = paste(c("top3axes_CoCA", "PCoA1_date", "PCoA3_lure")[c(in_top3axes_CoCA, in_PCoA1_date, in_PCoA3_lure)], collapse = " & ")) %>%
  ungroup() %>%
  arrange(desc(n))

region_plot <- ggplot(region_counts, aes(x = reorder(region, n), y = n)) +
  geom_col(fill = "#2166ac") +
  geom_text(aes(label = n), hjust = -0.2) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "Number of fungal taxa",
       title = "Overlap among top-200 candidate sets",
       subtitle = "3-axis CoCA loading vs. date-only screen (PCoA1) vs. lure-only screen (PCoA3)") +
  theme_minimal() +
  theme(plot.margin = margin(5.5, 30, 5.5, 5.5))
print(region_plot)
ggsave(file.path(out_fig_dir, "candidate_overlap_regions.png"), region_plot, width = 7.5, height = 5, dpi = 150)

## ---- 5. Overlap among SIGNIFICANT hits within each strategy's own test ----
# Not just "were the same taxa considered as candidates" but "did the same
# taxa come out significant under each strategy's targeted downstream test."
# Compared at two adjustment levels (no-date and date-adjusted) since the
# date-adjusted level collapsed to ~0 hits for PCoA1 in the original 3-axis
# screen (see interpretation.md) -- the no-date comparison is where any
# overlap is actually likely to be visible.

sig_taxa <- function(path, q_col = "q_value", threshold = 0.10) {
  if (!file.exists(path)) return(character(0))
  d <- read.csv(path)
  d$taxon[d[[q_col]] < threshold]
}

sig_sets <- list(
  top3axes_PCoA1_only = sig_taxa("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.csv"),
  PCoA1_date_PCoA1_only = sig_taxa("data/compare_insects_fungi_PCoA1_date/fungal_insect_association_results.PCoA1_only.csv"),
  top3axes_PCoA1_plus_date = sig_taxa("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_plus_date.csv"),
  PCoA1_date_PCoA1_plus_date = sig_taxa("data/compare_insects_fungi_PCoA1_date/fungal_insect_association_results.PCoA1_plus_date.csv"),
  top3axes_PCoA3_only_no_lure = sig_taxa("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA3_only.no_lure.csv"),
  PCoA3_lure_PCoA3_only_no_lure = sig_taxa("data/compare_insects_fungi_PCoA3_lure/fungal_insect_association_results.PCoA3_only.no_lure.csv"),
  top3axes_PCoA3_plus_date_no_lure = sig_taxa("data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA3_plus_date.no_lure.csv"),
  PCoA3_lure_PCoA3_plus_date_no_lure = sig_taxa("data/compare_insects_fungi_PCoA3_lure/fungal_insect_association_results.PCoA3_plus_date.no_lure.csv")
)
cat("\nSignificant-hit counts per test (q<0.10):\n")
print(sapply(sig_sets, length))

sig_comparisons <- list(
  list(label = "PCoA1, no date: top3axes vs PCoA1_date-screen", a = sig_sets$top3axes_PCoA1_only, b = sig_sets$PCoA1_date_PCoA1_only),
  list(label = "PCoA1, plus date: top3axes vs PCoA1_date-screen", a = sig_sets$top3axes_PCoA1_plus_date, b = sig_sets$PCoA1_date_PCoA1_plus_date),
  list(label = "PCoA3, no lure/no date: top3axes vs PCoA3_lure-screen", a = sig_sets$top3axes_PCoA3_only_no_lure, b = sig_sets$PCoA3_lure_PCoA3_only_no_lure),
  list(label = "PCoA3, no lure/plus date: top3axes vs PCoA3_lure-screen", a = sig_sets$top3axes_PCoA3_plus_date_no_lure, b = sig_sets$PCoA3_lure_PCoA3_plus_date_no_lure)
)

sig_overlap_summary <- bind_rows(lapply(sig_comparisons, function(cmp) {
  ov <- pair_overlap(cmp$a, cmp$b)
  data.frame(comparison = cmp$label, n_a = ov$n_a, n_b = ov$n_b, n_intersect = ov$n_intersect,
             jaccard = ov$jaccard, shared_taxa = paste(intersect(cmp$a, cmp$b), collapse = ";"))
}))
cat("\n--- Overlap among significant hits (q<0.10) between strategies, same axis ---\n")
print(sig_overlap_summary %>% select(-shared_taxa), row.names = FALSE)
write.csv(sig_overlap_summary, file.path(out_data_dir, "significant_hits_overlap.csv"), row.names = FALSE)

cat("\nDone. Outputs written to", out_data_dir, "and", out_fig_dir, "\n")
