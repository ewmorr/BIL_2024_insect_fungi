# Project organization

This document is the **structural map** of the repo: what lives where, and
why. `iterative_analysis_updates.md` is the companion **chronological log**
of findings, decisions, and the reasoning behind them, written chat-style and
appended to as work happens. Use `iterative_analysis_updates.md` to
understand *why* a choice was made or *what a result was*; use this document
to find *where a script or output lives* and *which lineage it belongs to*.
Keep both updated going forward -- this one in place (it describes current
state, not history), the other by appending new dated sections.

Last updated: 2026-09-01.

## The core fact to understand before touching anything

**There are two parallel insect-table constructions running through this
entire project**, and almost every script, data file, and figure belongs to
exactly one of them. Confusing the two is the single easiest way to produce
numbers that don't match anything already reported.

### Lineage A -- "family-filtered" (original, kept as the reference/backup)

The insect trap-catch table (`data/2024_insect_data/insect_species_tab.csv`)
restricted to **Curculionidae + Latridiidae** (the two dominant families,
~69% of individuals), then a singleton-taxon filter (`colSums > 1`) and an
empty-sample drop. Ends up **~52-58 taxa** depending on which script's own
downstream filter is applied. This is the table every original analysis in
this project was built on: `insect_ords.R`, `insect_fungal_coca_analysis.
general_workflow.r`, all the `insect_pcoa_*`, `insect_fungal_procrustes.r`,
`insect_fungal_pairwise_taxon_association*.r`, and every `presentation_items/
fig*.R` script. **This lineage is being kept as-is, untouched, as the
backup/reference set of results.**

### Lineage B -- "prevalence-filtered" (new, expected to become the headline)

The same insect table, but with **no family restriction** -- all ~98
families -- filtered instead by an **individual-taxon >=5-sample prevalence
filter** (`colSums(x > 0) >= 5`), the same threshold and logic already used
for fungal ASVs everywhere in this project. Ends up **153 taxa**. This
construction was validated in `insect_exploratory/insect_ords.
prevalence_filter.R` (low-stress, repeatable NMDS -- see
`iterative_analysis_updates.md`) and is being ported through the rest of the
pipeline. As of this writing it is expected to become the headline lineage,
with Lineage A retained for comparison/robustness-checking, not replaced.

### Lineage C -- "full insect table" (deliberate, alpha-diversity-only -- NOT the same as Lineage B, and not meant to be)

`compare_insect_fungi/insect_fungal_alpha_diversity.full_insect_table.r`
builds the insect table as all families with **only** the singleton-taxon
filter -- **no** >=5-sample prevalence filter. That gives **277 taxa**, a
different number from Lineage B's 153. This is intentional and settled, not
an inconsistency to fix: **the >=5-sample prevalence filter is a
noise-reduction step for ordination/covariance-based methods** (it keeps a
handful of very rare taxa from dominating a distance matrix or CoCA loading),
**and is not appropriate for alpha-diversity metrics**, where the raw
taxon-presence signal (richness, Shannon, Simpson dominance) is exactly what
is being measured -- prevalence-filtering before computing richness would
mean discarding real taxa and artificially deflating it. So for
alpha-diversity work specifically, "all families, singleton-filter only" (no
prevalence filter) is the CORRECT all-families construction, not a
lesser/orphaned stand-in for Lineage B. See "When to prevalence-filter vs.
not" below for the general rule this follows. Output lives in
`data/compare_insects_fungi_alpha_diversity_full_insect_table/` /
`figures/compare_insects_fungi_alpha_diversity_full_insect_table/` --
directory name kept as-is (predates the Lineage A/B/C naming introduced in
this document) but no longer flagged as needing reconciliation with Lineage
B; there will not be a prevalence-filtered alpha-diversity variant.

As of 2026-09-01, this script also tests the date effect with linear +
quadratic + cubic terms in one model (extending this project's linear+
quadratic convention one order further, added specifically because fig9's
trend lines showed more structure than a hump/dip could capture -- see
`alpha_diversity_date_shape.csv`), and tests whether **lure modifies the
insect~fungal correlation** (`fungal ~ site+date+insect*lure`, permutation
F-test on the interaction term, `alpha_diversity_insect_fungal_interaction_
by_lure.csv`) -- no metric showed a significant interaction (p_perm
0.80-0.95), which is why fig10 below shows the pooled (lure-blind)
correlation rather than splitting by lure.

Two presentation figures cover this lineage, both living directly in
`presentation_items/` alongside the Lineage A fig1-fig8 (not a dedicated
Lineage C folder) since alpha diversity only has a Lineage A/C split, no
Lineage B counterpart to disambiguate from:
- `fig9_alpha_diversity_by_lure.R` -- insect (panel a) and fungal (panel b)
  alpha diversity vs. collection date, faceted metric x lure, cubic OLS
  trend lines matching the permutation-tested cubic date model above.
- `fig10_insect_fungal_alpha_diversity_correlation.R` -- per-sample insect
  vs. fungal alpha diversity, one panel per metric, colored by date, pooled
  Spearman correlation annotated (an earlier lure-faceted version was
  dropped once the interaction test above found no lure effect).

## Directory and naming conventions

Two different conventions are in play, depending on whether original and
variant scripts share a folder or not:

1. **Shared script folders** (`compare_insect_fungi/`, `insect_exploratory/`)
   -- original and Lineage-B scripts live side by side, so the Lineage-B
   script's *filename* carries a `.prevalence_filtered_insect.r` (or
   `.prevalence_filter.R` for the one ordination script) suffix, and its
   *data/figure output directories* carry a matching `_prevalence_filtered_
   insect` suffix appended to the name of the directory the original
   analysis would have used. Example: `insect_fungal_coca_analysis.
   general_workflow.r` (Lineage A) writes to `data/compare_insects_fungi_
   top3axes/`; its Lineage-B counterpart `insect_fungal_coca_analysis.
   prevalence_filtered_insect.r` writes to `data/compare_insects_fungi_
   top3axes_prevalence_filtered_insect/`.

2. **Dedicated per-lineage top-level folders** (the new `presentation_items_
   prevalence_filtered_insect/`, mirroring `presentation_items/`) -- since
   the folder itself already disambiguates the lineage, scripts inside keep
   the SAME filename as their Lineage-A counterpart (`fig1_taxonomic_
   breakdown.R`, not `fig1_taxonomic_breakdown.prevalence_filtered_insect.R`).
   Outputs go to `data/presentation_items_prevalence_filtered_insect/` and
   `figures/presentation_items_prevalence_filtered_insect/`.

**Never overwrite a Lineage-A output.** Every Lineage-B script must write to
a `_prevalence_filtered_insect`-suffixed (or dedicated Lineage-B) directory,
never back into the original `data/2024_insect_data/`, `data/compare_
insects_fungi_top3axes/`, `figures/insect_exploratory/`, etc. (One
inconsistency of this kind was found and fixed on 2026-08-28: `insect_ords.
prevalence_filter.R` was writing into `data/2024_insect_data/` and `figures/
insect_exploratory/` with an `insect_prevalence_filter.`-prefixed filename
instead of its own directory. Moved to `data/insect_exploratory_prevalence_
filtered_insect/` and `figures/insect_exploratory_prevalence_filtered_
insect/`; the script's `out_data_dir`/`out_fig_dir` variables were updated
to match. No other script referenced the old paths.)

## Script + output index

Every table below lists: script -> what it does -> where its outputs live.
"Data" and "Figures" columns give the directory (files inside follow the
script's own naming, described in `iterative_analysis_updates.md` when not
obvious).

### Lineage A (family-filtered, Curculionidae+Latridiidae) -- reference/backup

| Script | What it does | Data output | Figure output |
|---|---|---|---|
| `insect_exploratory/insect_ords.R` | Original exploratory NMDS work; messy/iterative, left as-is | `data/2024_insect_data/` (various) | `figures/FEDRR_all_2024/`, `figures/insect_comps_2024/` |
| `compare_insect_fungi/insect_fungal_coca_analysis.general_workflow.r` | CoCA + per-taxon insect_PCoA1-3 association tests + direct taxon~date screens (now linear+quadratic, see below) | `data/compare_insects_fungi_top3axes/`, `data/2024_insect_data/insect_taxa_date_association.csv` | `figures/compare_insects_fungi_top3axes/` |
| `compare_insect_fungi/insect_fungal_coca_analysis.PCoA1_date_axis.r` | Single-axis (date-screen-selected candidates) CoCA variant | `data/compare_insects_fungi_PCoA1_date/` | `figures/compare_insects_fungi_PCoA1_date/` |
| `compare_insect_fungi/insect_fungal_coca_analysis.PCoA3_lure_axis.r` | Single-axis (lure-screen-selected candidates) CoCA variant | `data/compare_insects_fungi_PCoA3_lure/` | `figures/compare_insects_fungi_PCoA3_lure/` |
| `compare_insect_fungi/compare_selection_strategies.r` | Compares the three candidate-selection strategies above | `data/compare_selection_strategies/` | `figures/compare_selection_strategies/` |
| `compare_insect_fungi/insect_pcoa_lure_association.r` | Screens all insect PCoA axes for lure association | `data/2024_insect_data/insect_pcoa_axes.lure_association.csv` | `figures/insect_pcoa_scree.lure_highlighted.png`, `figures/insect_pcoa_top_lure_axis.by_lure.png` |
| `compare_insect_fungi/insect_pcoa_date_lure_association.r` | insect_PCoA1-3 vs. date/site/lure (linear only) | `data/2024_insect_data/insect_pcoa_axes.date_site_lure_association.csv` | `figures/insect_pcoa_top3axes.by_date.png`, `by_lure.png` |
| `compare_insect_fungi/fungal_community_seasonality.r` | Unbiased per-taxon fungal~date screen (ALL 6,342 prevalence-filtered fungal taxa; now linear+quadratic) + rarefied whole-community PERMANOVA/NMDS | `data/2024_fungi/fungal_taxa_date_association.all_taxa.csv`, `fungal_community_permanova.csv` | `figures/fungal_community_NMDS.date_and_site.pdf` |
| `compare_insect_fungi/fungal_community_lure_association.r` | Unbiased per-taxon fungal~lure screen (all taxa) | `data/2024_fungi/fungal_taxa_lure_association.all_taxa.csv` | -- |
| `compare_insect_fungi/insect_fungal_procrustes.r` | Whole-community Procrustes (insect NMDS vs. fungal rarefied NMDS) | `data/2024_insect_data/insect_fungal_procrustes_*.csv` | `figures/insect_fungal_procrustes.*` |
| `compare_insect_fungi/insect_fungal_pairwise_taxon_association.r` + `.residualized.r` | All-against-all individual fungal x insect taxon screen | `data/compare_insects_fungi_pairwise_taxa/` | `figures/compare_insects_fungi_pairwise_taxa/` |
| `compare_insect_fungi/responsive_taxa_abundance_by_site.R` | Abundance share of date/insect-responsive taxa by site | -- | `figures/responsive_taxa_abundance_by_site.png` |
| `compare_insect_fungi/fungi_date_order_breakdown.R`, `fungi_insect_order_breakdown.R`, `fungi_top_taxa_date_insectPCoA1.R` | Exploratory taxonomic drill-downs, superseded by fig4/fig5 | -- | `figures/fungi_*_breakdown.png` |
| `compare_insect_fungi/insect_fungal_alpha_diversity.r` | Shannon/Simpson-dominance/richness, insect (raw Curc.+Latr. table) vs. fungal (rarefy-and-average); cross-community correlation + site/lure/date covariate tests | `data/compare_insects_fungi_alpha_diversity/` | `figures/compare_insects_fungi_alpha_diversity/` |
| `presentation_items/fig1_taxonomic_breakdown.R` ... `fig8_association_trait_breakdown.R` | Polished manuscript/presentation figures (8 total) | `data/presentation_items/` | `figures/presentation_items/` |

### Lineage B (prevalence-filtered, all families, >=5-sample taxon filter)

| Script | What it does | Data output | Figure output |
|---|---|---|---|
| `insect_exploratory/insect_ords.prevalence_filter.R` | Validates the prevalence-filter approach: NMDS (stress 0.16, repeatable) + PERMANOVA | `data/insect_exploratory_prevalence_filtered_insect/` | `figures/insect_exploratory_prevalence_filtered_insect/` |
| `compare_insect_fungi/insect_fungal_coca_analysis.prevalence_filtered_insect.r` | Full CoCA workflow rebuilt on the 153-taxon table (mirrors `general_workflow.r`) | `data/compare_insects_fungi_top3axes_prevalence_filtered_insect/` | `figures/compare_insects_fungi_top3axes_prevalence_filtered_insect/` |
| `compare_insect_fungi/insect_pcoa_date_lure_association.prevalence_filtered_insect.r` | insect_PCoA1-3 vs. date/site/lure (linear) + quadratic-date test + Pearson taxon-driver correlations for all 3 axes | same dir as above | same dir as above |
| `compare_insect_fungi/insect_fungal_coca_pcoa1_quadratic_date_check.prevalence_filtered_insect.r` | Re-tests PCoA1's top-200 CoCA candidates with linear+quadratic date | same dir as above (`fungal_insect_association_results.PCoA1_plus_quadratic_date.csv`) | -- |
| `compare_insect_fungi/insect_fungal_coca_pcoa2_quadratic_date_check.prevalence_filtered_insect.r` | Same for PCoA2 | same dir (`..PCoA2_plus_quadratic_date.csv`) | -- |
| `compare_insect_fungi/insect_fungal_coca_factordate_check.prevalence_filtered_insect.r` | Strictest check: re-tests PCoA1 AND PCoA2 with `factor(date)` (6-level, absorbs any functional form) | same dir (`..PCoA{1,2}_plus_factordate.csv`, `date_adjustment_comparison.all_strategies.csv`) | -- |
| `compare_insect_fungi/insect_taxa_date_association.prevalence_filtered_insect.r` | Standalone direct taxon~date screen (linear+quadratic) for all 153 insect taxa -- broken out so it doesn't require rerunning the full CoCA script | `data/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_taxa_date_association.csv` | -- |
| `compare_insect_fungi/insect_fungal_procrustes.prevalence_filtered_insect.r` | Whole-community Procrustes (insect NMDS vs. fungal rarefied NMDS) rebuilt on the 153-taxon table; also writes the matched-69-sample insect PERMANOVA for the fig2 port. Adds the `Kingdom=="k__Fungi"` filter + graceful rarefaction-drop handling that the Lineage-A standalone `insect_fungal_procrustes.r` lacks (both already in `fig2_nmds_procrustes.R`) | `data/compare_insects_fungi_procrustes_prevalence_filtered_insect/` | `figures/compare_insects_fungi_procrustes_prevalence_filtered_insect/` |
| `compare_insect_fungi/insect_fungal_pairwise_taxon_association.prevalence_filtered_insect.r` | All-against-all individual fungal x insect taxon screen rebuilt on the 153-taxon table (mirrors `insect_fungal_pairwise_taxon_association.r`; same vectorized-OLS trap-blocked 999-perm test, `fungal~site+lure+date+insect_taxon`, within-insect-taxon + global BH). Insect predictors: the 153-taxon table re-filtered to >=5 prevalence in the 69 matched samples (still 153); fungal side lineage-invariant (6,342). 153 x 6,342 = 970,326 pairs, ~11 min single-core | `data/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect/` | `figures/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect/` |
| `compare_insect_fungi/insect_fungal_pairwise_taxon_association.residualized.prevalence_filtered_insect.r` | Lineage-B counterpart of `insect_fungal_pairwise_taxon_association.residualized.r`: re-runs the grid with each insect predictor residualized against `site + factor(date)` (strict, any-functional-form seasonal control), fungal-side model `site + lure + factor(date) + insect_taxon_resid`. Reads the main grid above for the side-by-side comparison; the flagged-taxon list is data-driven (top 5 insect taxa by main-grid hit count) rather than hard-coded | same dir as above (`*.date_site_residualized.csv`, `original_vs_residualized_comparison.csv`) | same dir as above (`original_vs_residualized_comparison.png`) |
| `presentation_items_prevalence_filtered_insect/fig2_nmds_procrustes.R` | Lineage-B counterpart of `presentation_items/fig2_nmds_procrustes.R` (same 3-panel NMDS + Procrustes layout). Reads its insect PERMANOVA from the procrustes script above; reuses the shared (lineage-invariant) `data/2024_fungi/fungal_community_permanova.csv` | `data/presentation_items_prevalence_filtered_insect/` (ordination cache) | `figures/presentation_items_prevalence_filtered_insect/` |
| `presentation_items_prevalence_filtered_insect/fig3_volcano_plots.R` | Lineage-B counterpart of `presentation_items/fig3_volcano_plots.R`, expanded from 1x3 to **2x3, rows = functional form**: row 1 = LINEAR-term date volcanoes (a insect~date, b fungal~date) + c fungal~PCoA1; row 2 = QUADRATIC-term date volcanoes (d insect~date, e fungal~date) + f fungal~PCoA2. Panels c/f use the **factor(date)-adjusted** results (`PCoA{1,2}_plus_factordate.csv`) -- the deseasonalized association, not the permissive unadjusted one; f's subtitle flags PCoA2 as unstable across date controls (6/0/9 of 200 unadjusted/quadratic/factor(date)). Panels b/e read the shared `fungal_taxa_date_association.all_taxa.csv` | -- (reads existing CSVs) | `figures/presentation_items_prevalence_filtered_insect/fig3_volcano_plots.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig4_association_taxonomic_breakdown.R` | Lineage-B counterpart of `presentation_items/fig4_association_taxonomic_breakdown.R`. **3 columns**: a insect~date and b fungal~date show the LINEAR and QUADRATIC date terms TOGETHER on one set of taxonomic-group axes (by subfamily / order+class respectively) as 4 dodged bars per group (blue/vermillion = linear later/earlier season, sky-blue/orange = quadratic dip/hump), each a full-height column; the top-n fold cutoff for a/b uses the COMBINED linear+quadratic hit count. c (fungal~PCoA1) and d (fungal~PCoA2) stack in a third, narrower column, staying single-term (factor(date)-adjusted, unchanged design from fig3) since PCoA1/PCoA2 are different predictor axes, not a linear/quadratic pair of one test | `data/presentation_items_prevalence_filtered_insect/fig4_taxonomic_breakdown.*.csv` | `figures/presentation_items_prevalence_filtered_insect/fig4_association_taxonomic_breakdown.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig5_family_breakdown_by_order.R` | Lineage-B counterpart of `presentation_items/fig5_family_breakdown_by_order.R`. **3 panels**: a = family-within-order breakdown of fig4 panel a/b's fungal~date COMBINED (dodged linear+quadratic) result; b/c = family-within-order breakdown of fig4 panel c/d's fungal~PCoA1/PCoA2 (single-term, unchanged). Revised 2026-09-01 alongside fig4, dropping from 4 panels to 3 | `data/presentation_items_prevalence_filtered_insect/fig5_family_breakdown_by_order.*.csv` | `figures/presentation_items_prevalence_filtered_insect/fig5_family_breakdown_by_order.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig6_responsive_taxa_abundance_and_counts_by_site.R` | Lineage-B counterpart of `presentation_items/fig6_responsive_taxa_abundance_and_counts_by_site.R` (same 2x2 layout: a/b relative sequence abundance by site, c/d ASV count by site, date- vs. insect-associated). Rebuilt on the Lineage-B matched dataset (all-families, >=5-sample prevalence-filtered insect table, same 69-sample match as fig2). "Date-associated" and "insect-associated" are each a COMBINED (union) definition instead of Lineage A's single-test one: date = q<0.10 on linear OR quadratic date term (2490/15422 taxa); insect = q<0.10 on PCoA1 OR PCoA2, factor(date)-adjusted (26/15422 taxa, no overlap between the two hit sets) | -- (reads existing CSVs) | `figures/presentation_items_prevalence_filtered_insect/fig6_responsive_taxa_abundance_and_counts_by_site.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig7_association_trait_breakdown.R` | Lineage-B counterpart of `presentation_items/fig7_association_trait_breakdown.R` (FungalTraits/Polme et al. 2020 genus-join trophic-mode x growth-form breakdown). **2 columns**: a fungal trait vs. date, LINEAR+QUADRATIC combined dodged bars (top 25 trait groups + "Other" by combined hit count), faceted by growth_form, full-height; b fungal trait vs. insect PCoA1 and c vs. PCoA2 (both single-term, factor(date)-adjusted) stack in a narrower column. "Unclassified genus"/"No FungalTraits match" kept as their own pseudo-facet bars (unchanged from Lineage A) | `data/presentation_items_prevalence_filtered_insect/fig7_trait_breakdown.*.csv` | `figures/presentation_items_prevalence_filtered_insect/fig7_association_trait_breakdown.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig8_association_trait_breakdown.R` | Lineage-B counterpart of `presentation_items/fig8_association_trait_breakdown.R` (same trait join, faceted by taxonomic Class instead of growth_form, `axis_id`-keyed since trait_group x Class is a genuine cross-tab, not a 1:1 nesting). **2 columns**, same layout as fig7: a fungal trait vs. date, LINEAR+QUADRATIC combined dodged bars, with the Class-facet-drop / lifestyle-within-class-cell-fold thresholds (`min_facet_hits_date`) now applied to the COMBINED hit count; b/c fungal trait vs. PCoA1/PCoA2 (single-term, factor(date)-adjusted) stack narrower. "Unclassified genus"/"No FungalTraits match" dropped entirely (unchanged from Lineage A) | `data/presentation_items_prevalence_filtered_insect/fig8_trait_breakdown.*.csv` | `figures/presentation_items_prevalence_filtered_insect/fig8_association_trait_breakdown.{png,pdf}` |

### Lineage C (full insect table, alpha-diversity-only by design -- see above)

| Script | What it does | Data output | Figure output |
|---|---|---|---|
| `compare_insect_fungi/insect_fungal_alpha_diversity.full_insect_table.r` | Same alpha-diversity comparison as Lineage A's version, but insect table = all families, singleton-filter only (277 taxa, NOT the 153-taxon Lineage B table). Date effects tested linear+quadratic+cubic (`alpha_diversity_date_shape.csv`); also tests an insect x lure interaction on the insect~fungal correlation (`alpha_diversity_insect_fungal_interaction_by_lure.csv`) -- not significant for any metric | `data/compare_insects_fungi_alpha_diversity_full_insect_table/` | `figures/compare_insects_fungi_alpha_diversity_full_insect_table/` |
| `presentation_items/fig9_alpha_diversity_by_lure.R` | Presentation figure: insect (panel a) + fungal (panel b) alpha diversity vs. collection date, faceted metric x lure, cubic OLS trend lines. Lives directly in `presentation_items/` (not a dedicated Lineage C folder) since alpha diversity has no Lineage B counterpart to disambiguate from | -- (reads existing CSVs) | `figures/presentation_items/fig9_alpha_diversity_by_lure.{png,pdf}` |
| `presentation_items/fig10_insect_fungal_alpha_diversity_correlation.R` | Presentation figure: per-sample insect vs. fungal alpha diversity, one panel per metric, colored by date, pooled Spearman correlation annotated. An earlier lure-faceted version was dropped after the interaction test above found no lure effect | -- (reads existing CSVs) | `figures/presentation_items/fig10_insect_fungal_alpha_diversity_correlation.{png,pdf}` |

## Known gaps / open items (update as these get resolved)

As of 2026-09-01, ported to Lineage B: insect NMDS/PERMANOVA, the full CoCA
workflow, PCoA1-3 vs. date/site/lure (incl. quadratic and factor(date)
robustness checks), taxon-level PCoA drivers, the direct insect~date screen,
the whole-community Procrustes, the all-against-all pairwise taxon screen
(main + site/factor(date)-residualized), and all 8 presentation figures
(fig1 lineage-invariant; fig2 NMDS+Procrustes; fig3 volcano plots, linear-
row/quadratic-row 2x3; fig4/fig5 taxonomic breakdowns and fig7/fig8 trait
breakdowns, all four using the same combined-term design -- linear+quadratic
date terms as dodged bars in one panel, PCoA1/PCoA2 kept as separate single-
term panels stacked in a narrower column; fig6 responsive-taxa abundance/
count by site, using a combined linear-or-quadratic date definition and a
combined PCoA1-or-PCoA2 insect definition). **Lineage B port is now
feature-complete with Lineage A** other than the items below, which are
either deliberately out of scope or minor variants not yet rebuilt:

- **Alpha diversity** -- NOT a gap: settled as Lineage A (family-filtered,
  backup) + Lineage C (full all-families table, headline), by deliberate
  choice (see Lineage C section above and "When to prevalence-filter vs.
  not" below). There will not be a Lineage-B (prevalence-filtered)
  alpha-diversity variant -- don't add one without revisiting that decision
  first.
- **All-against-all pairwise taxon screen** -- DONE (2026-09-01):
  `insect_fungal_pairwise_taxon_association.prevalence_filtered_insect.r` +
  `.residualized.prevalence_filtered_insect.r`. Result tracks Lineage A
  closely -- only ~5 of 153 insect taxa carry any q<0.10 fungal partners
  (0 pairs survive a global BH correction), dominated by early-season
  Scolytinae (*Pityogenes hopkinsi* alone ~1,590), and after residualizing
  the insect predictor against `site + factor(date)` only *P. hopkinsi*
  retains partners (~38% of its hits). See
  `lineage_B_prevalence_filtered_insect_top_level_interpretation.md` sec. 6.
- **`presentation_items_prevalence_filtered_insect/`** -- NOT a gap anymore:
  all 8 figures are done (fig1 lineage-invariant, no port needed; fig2/fig3
  2026-08-31; fig4-fig8 2026-09-01). fig1's panel c was changed on
  2026-08-31 to show the top 10 insect genera across all families rather
  than within Curculionidae+Latridiidae -- see `iterative_analysis_updates.md`.
- **PCoA1_date_axis / PCoA3_lure_axis / compare_selection_strategies**
  variants -- not yet rebuilt for Lineage B.
- **Direct fungal~date and fungal~lure unbiased screens**
  (`fungal_community_seasonality.r`, `fungal_community_lure_association.r`)
  don't depend on the insect table at all (fungal-only), so they do NOT need
  a Lineage B counterpart -- their existing Lineage-A output already applies
  to both lineages. Only re-run these if the FUNGAL-side filtering/table
  changes, not the insect side.

## When to prevalence-filter an all-families insect table, vs. not

This is a general rule, not a one-off decision, worth applying to any future
all-families insect analysis:

- **Ordination / covariance / community-composition methods** (NMDS, PCoA,
  CoCA, PERMANOVA, per-taxon association tests) -- use the >=5-sample
  prevalence filter (Lineage B). Rare, sparsely-observed taxa contribute
  disproportionate noise to a distance matrix or loading and add
  multiple-testing burden without power; filtering them out is standard
  practice and is exactly what's already done on the fungal side throughout
  this project.
- **Alpha-diversity metrics** (richness, Shannon, Simpson dominance) -- do
  NOT prevalence-filter (Lineage C's "full table" construction, singleton
  filter only). These metrics are direct summaries of how many/how evenly
  distributed taxa are observed per sample; discarding taxa that happen to
  occur in fewer than 5 samples across the whole dataset would mean
  systematically undercounting richness and distorting diversity for
  reasons unrelated to the sample being described. The singleton filter
  (drop taxa with a single individual across the ENTIRE dataset) is kept
  because that's noise-floor cleanup, not the same kind of aggressive
  per-taxon filtering the >=5-sample rule applies.

If a future analysis doesn't obviously fall into either bucket, ask rather
than defaulting to whichever table is already loaded in the script you're
extending.

## Cross-lineage methodological conventions (apply to BOTH A and B)

These don't change between lineages -- they're properties of the fungal
side, the statistical approach, or shared infrastructure:

- **Fungal ASV table**: always Kingdom==Fungi filtered
  (`ASVs_taxonomy.tsv`), then a >=5-sample prevalence filter for per-taxon
  tests (`fungal_f`). Community-level analyses (rarefied PERMANOVA/NMDS,
  alpha diversity) use the Kingdom-filtered table WITHOUT the prevalence
  filter, rarefied to depth 5000 x 100 iterations
  (`multiple_subsamples()`/`avg_matrix_list()` in `library/library.R`).
- **CLR transform** (`clr_transform()`, pseudocount=1) for per-taxon fungal
  abundance in every permutation-based lm() test.
- **Insect PCoA construction**: Hellinger transform (`decostand(..,
  "hellinger")`) + Euclidean distance (`vegdist(.., "euclidean")`) +
  classical PCoA (`cmdscale()`), kept to the first 3 axes
  (insect_PCoA1/2/3). This combination is mathematically equivalent to
  Hellinger distance and avoids the negative-eigenvalue problem Bray-Curtis
  causes in classical (metric) PCoA -- see `iterative_analysis_updates.md`'s
  "Explain Euclidean distance usage" entry for the full rationale. NMDS
  panels (fig2-style) use raw Bray-Curtis directly instead, since NMDS
  doesn't need a Euclidean-embeddable distance.
- **Permutation schemes for trap-level fixed effects** (999 permutations
  throughout):
  - *Date*: varies within trap (repeated sampling over the season) -->
    permute date labels WITHIN trap, holding site/lure fixed.
  - *Lure*: constant within trap (one trap per site x lure combination) -->
    permute lure labels among the traps WITHIN each site.
  - *Site*: also trap-constant, no grouping level above it --> permute the
    site label freely among ALL traps (same logic as lure, one level up).
  - *Quadratic/factor(date)*: same within-trap date permutation as linear
    date, just with a richer model on each side of the comparison.
- **Significance convention**: q<0.10 (BH-FDR) for per-taxon/per-axis
  permutation tests, uncorrected p-values for the small, fixed set of
  community-level PERMANOVA terms (site/lure/date) -- these are NOT run
  through the same multiple-testing correction as a taxon screen.
- **Diversity metrics** (`compare_insect_fungi/insect_fungal_alpha_
  diversity*.r`): Shannon (`vegan::diversity(index="shannon")`), Simpson
  **dominance** (D = sum(p_i^2), computed directly -- NOT vegan's
  `index="simpson"`, which returns 1-D, a diversity not a dominance
  measure), richness (`vegan::specnumber()`, observed taxa, no
  extrapolation).
- **Linear+quadratic date test pattern** (as of 2026-08-28, added to
  `fungal_community_seasonality.r` and `insect_fungal_coca_analysis.
  general_workflow.r`, and used natively in the Lineage-B scripts): fits
  `y ~ site + lure + date_c + I(date_c^2)` (date centered to avoid a
  huge-magnitude collinear squared term), extracts both terms' t-statistics
  from ONE model per permutation. Output columns: `t_stat`/`p_perm`/
  `q_value` (linear -- UNCHANGED semantics from the original linear-only
  version, since a dozen+ downstream scripts filter on `q_value`) plus
  `t_stat_quad`/`p_perm_quad`/`q_value_quad`/`shape` (new). `shape` is
  "hump (peaks mid-season)" / "dip (troughs mid-season)" if the quadratic
  term is significant (quadratic takes priority over linear when both are),
  else "linear increase (late-season)" / "linear decrease (early-season)" if
  only the linear term is, else "no significant date pattern". **Gotcha
  already hit once**: `c(t_stat = obs["t_linear"], ...)` in R mangles the
  name to `"t_stat.t_linear"` when the right-hand side is already a named
  scalar -- always wrap in `unname()` before combining into a new named
  vector, or downstream `r["t_stat"]`-style extraction silently returns NA.

## Top-level docs in this repo

- `project_organization.md` (this file) -- structural map, current state.
- `iterative_analysis_updates.md` -- chronological log of analyses run and
  findings, in the order they happened. Append new dated `###` sections
  here; don't edit old ones except to add an explicit correction entry (see
  its own "CORRECTION:" sections for the house style on how to handle a
  finding that later turns out to be wrong -- don't silently rewrite, add a
  new entry that supersedes it).
- `lineage_A_family_filtered_insect_top_level_interepretation.md` -- a point-in-time re-derivation of
  top-level taxonomic/ecological summaries after the Kingdom==Fungi filter
  fix, explicitly noting which earlier numbers it supersedes. Not a living
  document in the same way as the other two -- treat as a dated snapshot
  because this line of analysis is considered complete for now.
  Covers Lineage A (family-filtered) only.
- `lineage_B_prevalence_filtered_insect_top_level_interpretation.md` -- the
  Lineage B (prevalence-filtered insect table) counterpart to the above:
  top-line numbers for everything ported to Lineage B so far, with explicit
  comparison to Lineage A throughout. When more Lineage B analyses are ported 
  (see "Known gaps"), edit this document to add new interpretation
  and summaries, matching the convention of 
  `lineage_A_family_filtered_insect_top_level_interepretation.md` already established. 
  If new analyses modify previous findings update this document accordingly
  rather than retaining outdated numbers and interpretation. We will keep this document up to date
  with the current lineage B analyses until this analytical line is considered complete.
- `lineage_C_alpha_div_full_insect_table_top_level_interpretation.md` -- the
  Lineage C (full insect table, alpha-diversity-only) counterpart to the
  above two: top-line numbers for the alpha-diversity comparison (value
  ranges, insect~fungal correlation, the lure-interaction test, within-
  community site/lure/date effects, and the linear+quadratic+cubic date
  shape classification), matching the same convention. Since Lineage C has
  no ordination/CoCA/Procrustes/pairwise-screen analog (alpha-diversity-only
  by design, see the Lineage C section above), this document only covers
  `insect_fungal_alpha_diversity.full_insect_table.r` and its two
  presentation figures (fig9, fig10) -- update it if that script or those
  figures change, rather than retaining outdated numbers.
