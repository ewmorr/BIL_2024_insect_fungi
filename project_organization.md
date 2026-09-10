# Project organization

This document is the **structural map** of the repo: what lives where, and
why. `iterative_analysis_updates.md` is the companion **chronological log**
of findings, decisions, and the reasoning behind them, written chat-style and
appended to as work happens. Use `iterative_analysis_updates.md` to
understand *why* a choice was made or *what a result was*; use this document
to find *where a script or output lives* and *which lineage it belongs to*.
Keep both updated going forward -- this one in place (it describes current
state, not history), the other by appending new dated sections.

Last updated: 2026-09-10.

## The core fact to understand before touching anything

**There are two parallel insect-table constructions running through this
entire project**, and almost every script, data file, and figure belongs to
exactly one of them. Confusing the two is the single easiest way to produce
numbers that don't match anything already reported. (Two later lineages,
C and D, do not add a third and fourth construction -- C uses a raw
all-families table, and D uses a taxonomic *slice* of that same raw table.
See their subsections below.)

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
builds the insect table as all families, **raw counts** -- **no** >=5-sample
prevalence filter, and (as of 2026-09-09) **no global singleton filter
either**. That gives **426 taxa**, a different number from Lineage B's 153.
This is intentional and settled, not an inconsistency to fix: **the
>=5-sample prevalence filter is a noise-reduction step for ordination/
covariance-based methods** (it keeps a handful of very rare taxa from
dominating a distance matrix or CoCA loading), **and is not appropriate for
alpha-diversity metrics**, where the raw taxon-presence signal (richness,
Shannon, Simpson dominance) is exactly what is being measured --
prevalence-filtering before computing richness would mean discarding real
taxa and artificially deflating it. The milder `colSums > 1` global
singleton filter this script used through 2026-09-08 was dropped on
2026-09-09 for consistency with the asymptotic (`iNEXT`) script, which
*cannot* use it (Chao estimators are built on the observed singleton count);
with both Lineage C scripts now on the raw table, "observed richness" means
the same thing in both. Effect was modest -- observed insect richness rose a
median of ~2 taxa/sample (Shannon/Simpson dominance barely moved). See "When
to prevalence-filter vs. not" below and the 2026-09-09 entry in
`iterative_analysis_updates.md`. (Lineage A's family-filtered
`insect_fungal_alpha_diversity.r` is frozen as the backup and keeps the
singleton filter -- unchanged.) Output lives in
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

A third Lineage C script, `insect_abundance_fungal_alpha_diversity.full_
insect_table.r` (added 2026-09-10), tests a distinct hypothesis on the same
table: if insects act as *non-targeted* dispersal vectors (picking up fungal
propagules at random as they move through the environment) then total insect
**abundance** -- individuals per sample, `rowSums()` of the raw 426-taxon
table -- rather than insect diversity should predict fungal alpha diversity,
most clearly fungal richness. It reads the fungal alpha-diversity CSVs the
main script writes (observed + asymptotic) rather than recomputing them, and
tests each fungal metric both as a plain pooled Spearman and season-
residualized (`log1p(abundance)` and the fungal metric each residualized on
`site + factor(date)`, Spearman on residuals, free-permutation p). Result is
a flat null on every metric under both analyses and both value types (all
|rho| <= 0.11, p_perm >= 0.38); total insect abundance is itself near-flat
with date (rho 0.01), so this is not a seasonality artifact. Outputs land in
the existing Lineage C alpha-diversity dirs with an `insect_abundance_`
filename prefix. Written up in
`lineage_C_alpha_div_full_insect_table_top_level_interpretation.md` §13 and
the 2026-09-10 entry in `iterative_analysis_updates.md`.

Four presentation figures cover this lineage, all living directly in
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
- `fig11_asymptotic_diversity_by_lure.R` / `fig12_insect_fungal_asymptotic_
  diversity_correlation.R` (added 2026-09-08) -- fig9/fig10 reproduced
  one-for-one on the ASYMPTOTIC (Chao-extrapolated, via `iNEXT`) diversity
  estimates from `insect_fungal_asymptotic_richness_iNEXT.full_insect_
  table.r` instead of observed diversity; see that script's entry below and
  its header for why the fungal side is left unrarefied there (unlike fig9/
  fig10's rarefied-and-averaged fungal table) and why "Simpson" in fig11/
  fig12 is Hill number q=2 diversity (higher = more even), the OPPOSITE
  direction from fig9/fig10's Simpson dominance.

### Lineage D -- "bark-beetle x Ophiostomatales co-occurrence" (targeted, presence/absence-only)

A targeted, hypothesis-driven line, NOT a whole-community method and NOT a new
insect-table construction: it takes the full all-families trap-catch table
(the Lineage C raw table), restricts it **taxonomically** to Subfamily ==
Scolytinae (the bark/ambrosia beetles the EDRR traps are baited for), reduces
to **presence/absence**, and screens each beetle species against each
**Ophiostomatales** ("blue-stain") fungal taxon for co-occurrence using
**Jaccard similarity** (`J = a/(a+b+c)`, shared absence excluded) with a
restricted-permutation null. Beetle x fungal-taxon pairs only -- no
beetle-vs-beetle, no fungus-vs-fungus. The insect prevalence floor here
(>=3 samples) only decides which beetles carry enough presences to test; it
is not a Lineage-B-style statistical filter of the table.

Fungal side: raw (un-rarefied) ASV table filtered to Order ==
`o__Ophiostomatales`, presence/absence, rolled up three ways (named species /
genus / individual ASV), each grid run separately. Matched dataset = the same
69 insect+fungal samples the `compare_insect_fungi/` scripts use.

Two permutation nulls matter (see the Lineage D interpretation doc for why):
- `within_trap` (blocks = site, plots = trap, within = free) -- controls
  site/trap/lure but NOT season; both Scolytinae and Ophiostomatales are
  spring-weighted, so this null inflates co-occurrence from shared phenology.
- `within_date` / `within_site_x_date` (blocks = collection date, or site x
  date) -- holds each date's beetle prevalence fixed, so shared seasonality
  cannot generate signal. This is the presence/absence analogue of the
  `.residualized.r` scripts' "residualize on site + factor(date)".

Scripts + outputs live in `insect_exploratory/` and
`data/`+`figures/insect_exploratory/{target_genera,scolytinae}/` (see the
Lineage D index table below). Findings are written up in
`lineage_D_insect_fungus_cooccurence_patterns.md`.

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
| `compare_insect_fungi/add_standardized_slope_to_fig3_screens.R` | One-off retrofit (2026-09-08): recomputes the OBSERVED-fit standardized partial slope `beta_std`(`_quad`) = `b*sd(x)/sd(y)` for the taxa already in three fig3 input CSVs and merges the column(s) in place -- no permutation, all other values byte-identical. Exists because `fungal_community_seasonality.r` uses `mclapply`/`mc.set.seed` so a full rerun would jitter q-values that ripple into fig4-fig8/figS1; the source scripts were still edited to emit the column natively for future full reruns | rewrites in place: `data/2024_insect_data/insect_taxa_date_association.csv`, `data/2024_fungi/fungal_taxa_date_association.all_taxa.csv`, `data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv` | -- |
| `presentation_items/fig1_taxonomic_breakdown.R` ... `fig8_association_trait_breakdown.R` | Polished manuscript/presentation figures (8 total). fig3 volcano x-axis is the standardized slope `beta_std` (not the t-statistic) as of 2026-09-08 -- see Lineage-B fig3 row and the 2026-09-08 log entry | `data/presentation_items/` | `figures/presentation_items/` |

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
| `presentation_items_prevalence_filtered_insect/fig3_volcano_plots.R` | Lineage-B counterpart of `presentation_items/fig3_volcano_plots.R`, expanded from 1x3 to **2x3, rows = functional form**: row 1 = LINEAR-term date volcanoes (a insect~date, b fungal~date) + c fungal~PCoA1; row 2 = QUADRATIC-term date volcanoes (d insect~date, e fungal~date) + f fungal~PCoA2. Panels c/f use the **factor(date)-adjusted** results (`PCoA{1,2}_plus_factordate.csv`) -- the deseasonalized association, not the permissive unadjusted one; f's subtitle flags PCoA2 as unstable across date controls (6/0/9 of 200 unadjusted/quadratic/factor(date)). Panels b/e read the shared `fungal_taxa_date_association.all_taxa.csv`. **REVISED 2026-09-02** to adopt a linear-term-takes-priority reading of peak timing (see "Peak-timing convention" below): panels a/b now color points by whether the taxon ALSO carries a significant quadratic term (not by linear significance, which the dashed q<0.10 line already shows) -- insect: 69/153 linear-significant, 15 of those also quadratic; fungal: 2,186/6,342 linear-significant, 403 of those also quadratic. Panels d/e now EXCLUDE any taxon with a significant linear term before plotting, since under the new priority a linear-significant taxon is called early/late-season regardless of its quadratic term -- insect: 23/84 significant among non-linear-trending taxa (of 153 total, 69 excluded); fungal: 304/4,156 (of 6,342 total, 2,186 excluded). Panels c/f unchanged. Only fig3 has been revised so far -- fig4/fig5/fig7/fig8 still use the OLD quadratic-priority `shape` reading and are pending the same revisit (see "Peak-timing convention" below). **REVISED 2026-09-08**: volcano x-axis switched from the raw permutation t-statistic to the fully-standardized partial slope (`beta_std`/`beta_std_quad`, `b*sd(x)/sd(y)`) -- a proper effect size for the continuous date/PCoA predictors rather than a two-group Cohen's d. Deterministic observed-fit quantity, retrofitted onto the existing screen CSVs (Lineage-A `general_workflow.r` + `fungal_community_seasonality.r` outputs via `compare_insect_fungi/add_standardized_slope_to_fig3_screens.R`; the two prevalence-filtered scripts re-run) with permutation p/q-values unchanged; `sign(beta_std)==sign(t_stat)` so the combo coloring / linear-exclusion in panels d/e are unaffected. Labeled hits are now ranked by `|beta_std|`. Lineage-A `presentation_items/fig3_volcano_plots.R` got the identical x-axis change the same day | -- (reads existing CSVs) | `figures/presentation_items_prevalence_filtered_insect/fig3_volcano_plots.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/figS1_date_shape_bin_examples.R` | **Supplemental.** Fungal-only, lineage-invariant (rebuilds the matched 69-sample dataset + prevalence filter + CLR transform used by `fungal_community_seasonality.r`), filed here because it exists to support the fig3 peak-timing redesign above. Classifies all 6,342 tested fungal taxa into the full 3x3 linear (later/earlier/n.s.) x quadratic (hump/dip/n.s.) crosstab (8 non-trivial bins), writes the bin-count table, and plots 2 example taxa per bin (CLR abundance vs. date, quadratic OLS trend line) so a reader can see what e.g. "later-season + hump" actually looks like. NOTE: this 8-bin crosstab is NOT the adopted peak-timing classification (which gives linear priority) -- it's the raw independent-terms view that motivated the redesign; kept for reference, not as competing labels | `data/presentation_items_prevalence_filtered_insect/figS1_date_shape_bin_counts.csv` | `figures/presentation_items_prevalence_filtered_insect/figS1_date_shape_bin_examples.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig4_association_taxonomic_breakdown.R` | Lineage-B counterpart of `presentation_items/fig4_association_taxonomic_breakdown.R`. **3 columns**: a insect~date and b fungal~date, by taxonomic group (subfamily / order+class respectively), each a full-height column; c (fungal~PCoA1) and d (fungal~PCoA2) stack in a third, narrower column, staying single-term (factor(date)-adjusted, unchanged design from fig3) since PCoA1/PCoA2 are different predictor axes. **REVISED 2026-09-02** (third pass) to the linear-priority "peak timing" reading (see "Peak-timing convention" below): panels a/b now show 4 mutually-exclusive dodged bars per group -- Peaks late season (blue) / Peaks early season (vermillion) / Peaks mid-season (orange, hump with no sig linear term) / Bimodal-early+late (sky-blue, dip with no sig linear term) -- instead of the old "up to 2 independent hit-instances per taxon" design; top-n fold ranking now uses plain non-overlapping taxon counts. Insect: 92 taxa (45 late + 24 early + 22 mid + 1 bimodal) across 55 subfamilies. Fungal: 2,490 taxa (1,272 late + 914 early + 215 mid + 89 bimodal) across 102 orders -- the 2,490 total is unchanged from the pre-revision "date-associated" count (linear OR quadratic significant), only how those taxa are split into bars changed. **Fourth pass, same day**: since the 4 categories are no longer a signed t-statistic, panels a/b's bars all now run the same direction (0 -> positive, no more diverging left/right layout), dodge order top-to-bottom within each group is early / bimodal / mid-season / late. Panels c/d unchanged throughout | `data/presentation_items_prevalence_filtered_insect/fig4_taxonomic_breakdown.*.csv` | `figures/presentation_items_prevalence_filtered_insect/fig4_association_taxonomic_breakdown.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig5_family_breakdown_by_order.R` | Lineage-B counterpart of `presentation_items/fig5_family_breakdown_by_order.R`. **3 panels**: a = family-within-order breakdown of fig4 panel b's fungal~date result; b/c = family-within-order breakdown of fig4 panel c/d's fungal~PCoA1/PCoA2 (single-term, unchanged). Revised 2026-09-01 alongside fig4, dropping from 4 panels to 3. **REVISED 2026-09-02** (third + fourth pass) alongside fig4 -- panel a uses the same 4-category peak-timing bars, non-overlapping order/family ranking, and same-direction/early-bimodal-mid-late bar layout as fig4 panel b (1,471 taxa across the top 14 fig4a/b orders, by family); panels b/c unchanged | `data/presentation_items_prevalence_filtered_insect/fig5_family_breakdown_by_order.*.csv` | `figures/presentation_items_prevalence_filtered_insect/fig5_family_breakdown_by_order.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig6_responsive_taxa_abundance_and_counts_by_site.R` | Lineage-B counterpart of `presentation_items/fig6_responsive_taxa_abundance_and_counts_by_site.R` (same 2x2 layout: a/b relative sequence abundance by site, c/d ASV count by site, date- vs. insect-associated). Rebuilt on the Lineage-B matched dataset (all-families, >=5-sample prevalence-filtered insect table, same 69-sample match as fig2). "Date-associated" and "insect-associated" are each a COMBINED (union) definition instead of Lineage A's single-test one: date = q<0.10 on linear OR quadratic date term (2490/15422 taxa); insect = q<0.10 on PCoA1 OR PCoA2, factor(date)-adjusted (26/15422 taxa, no overlap between the two hit sets) | -- (reads existing CSVs) | `figures/presentation_items_prevalence_filtered_insect/fig6_responsive_taxa_abundance_and_counts_by_site.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig7_association_trait_breakdown.R` | Lineage-B counterpart of `presentation_items/fig7_association_trait_breakdown.R` (FungalTraits/Polme et al. 2020 genus-join trophic-mode x growth-form breakdown). **2 columns**: a fungal trait vs. date, faceted by growth_form, full-height; b fungal trait vs. insect PCoA1 and c vs. PCoA2 (both single-term, factor(date)-adjusted) stack in a narrower column. "Unclassified genus"/"No FungalTraits match" kept as their own pseudo-facet bars (unchanged from Lineage A). **REVISED 2026-09-02** to the linear-priority "peak timing" reading, ported directly from fig4/fig5's finished design (4 mutually-exclusive same-direction dodged bars, early/bimodal/mid-season/late top-to-bottom order, plain non-overlapping taxon counts for the top-25-trait-group fold): 2,490 taxa (1,272 late + 914 early + 215 mid-season + 89 bimodal) across 42 trait groups -- same totals as fig4 panel b, different grouping | `data/presentation_items_prevalence_filtered_insect/fig7_trait_breakdown.*.csv` | `figures/presentation_items_prevalence_filtered_insect/fig7_association_trait_breakdown.{png,pdf}` |
| `presentation_items_prevalence_filtered_insect/fig8_association_trait_breakdown.R` | Lineage-B counterpart of `presentation_items/fig8_association_trait_breakdown.R` (same trait join, faceted by taxonomic Class instead of growth_form, `axis_id`-keyed since trait_group x Class is a genuine cross-tab, not a 1:1 nesting). **2 columns**, same layout as fig7: a fungal trait vs. date; b/c fungal trait vs. PCoA1/PCoA2 (single-term, factor(date)-adjusted) stack narrower. "Unclassified genus"/"No FungalTraits match" dropped entirely (unchanged from Lineage A). **REVISED 2026-09-02**, same peak-timing design as fig7: 1,747 taxa with a FungalTraits match and a significant peak-timing call (970 late + 576 early + 153 mid-season + 48 bimodal) across 20 lifestyles; the two `min_facet_hits_date` fold thresholds (drop a sparse Class facet; fold a sparse lifestyle-within-Class cell into "Other") switched from hit-instance sums to plain taxon counts | `data/presentation_items_prevalence_filtered_insect/fig8_trait_breakdown.*.csv` | `figures/presentation_items_prevalence_filtered_insect/fig8_association_trait_breakdown.{png,pdf}` |

**Peak-timing convention (2026-09-02, DONE for fig3/fig4/fig5/fig7/fig8):** after examining example taxa per
linear x quadratic combination (`figS1_date_shape_bin_examples.R`), a taxon's
non-linearity itself was judged not particularly noteworthy on its own -- what matters
for interpretation is WHEN it peaks. Adopted reading, **linear term takes priority
over quadratic** whenever the linear term is significant (the OPPOSITE priority from
the `shape` column `fungal_community_seasonality.r`/`insect_taxa_date_association.
prevalence_filtered_insect.r` write, which gives the quadratic term priority):
  - Significant negative linear trend -> **peaks early** (regardless of quadratic).
  - Significant positive linear trend -> **peaks late** (regardless of quadratic).
  - No significant linear term, significant negative quadratic (hump) -> **peaks
    mid-season**.
  - No significant linear term, significant positive quadratic (dip) -> **bimodal /
    early+late**.
Applied to `fig3_volcano_plots.R`, `fig4_association_taxonomic_breakdown.R`,
`fig5_family_breakdown_by_order.R`, `fig7_association_trait_breakdown.R`, and
`fig8_association_trait_breakdown.R` (all 2026-09-02) -- each taxon now gets exactly
ONE of the 4 categories above as a dodged bar (fig4/5/7/8) or a color/exclusion rule
(fig3), replacing the "up to 2 independent hit-instances per taxon" design; group/
family/trait-group ranking for every fold threshold now uses plain non-overlapping
taxon counts. fig4/5/7/8's dodged bars all run the same direction (0 -> positive,
not diverging left/right) and dodge top-to-bottom as early / bimodal / mid-season /
late. The underlying `shape` column in `fungal_taxa_date_association.all_taxa.csv`/
`insect_taxa_date_association.prevalence_filtered_insect.r`'s output is NOT updated
to this priority (still quadratic-priority) -- it isn't read by any of the revised
figures, but don't assume its `shape` label matches a taxon's peak-timing category.

**Bug caught and fixed during the fig4/fig5 bar-reorder (2026-09-02, see
`iterative_analysis_updates.md`'s "CORRECTION" for that date):** the peak-timing
`case_when()` blocks in fig4/fig5 originally referenced `cat_levels[1..4]` by
POSITION; reordering `cat_levels` for the bar layout above silently swapped the
"Peaks early season"/"Peaks mid-season" labels until fixed. All 4 figures now assign
literal category strings instead of indexing into `cat_levels`, so this can't recur.
If you add a 5th script using this pattern, follow that convention (literal strings,
never `cat_levels[i]`).

### Lineage C (full insect table, alpha-diversity-only by design -- see above)

| Script | What it does | Data output | Figure output |
|---|---|---|---|
| `compare_insect_fungi/insect_fungal_alpha_diversity.full_insect_table.r` | Same alpha-diversity comparison as Lineage A's version, but insect table = all families, **raw counts, no singleton filter as of 2026-09-09** (426 taxa, NOT the 153-taxon Lineage B table; was 277 under the `colSums > 1` filter through 2026-09-08 -- dropped for consistency with the asymptotic script below, see the Lineage C section and the 2026-09-09 log entry). Date effects tested linear+quadratic+cubic (`alpha_diversity_date_shape.csv`); also tests an insect x lure interaction on the insect~fungal correlation (`alpha_diversity_insect_fungal_interaction_by_lure.csv`) -- not significant for any metric | `data/compare_insects_fungi_alpha_diversity_full_insect_table/` | `figures/compare_insects_fungi_alpha_diversity_full_insect_table/` |
| `compare_insect_fungi/insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r` | Asymptotic (Chao-extrapolated) richness/Shannon(Hill q=1)/Simpson(Hill q=2) diversity via `iNEXT::ChaoRichness()`/`ChaoShannon()`/`ChaoSimpson()` (NOT `iNEXT()` itself -- its curve machinery didn't finish in 120s on the largest fungal sample; the point/SE/CI-only Chao*() functions reproduce its `$AsyEst` estimates and run in ~1-2s each). Insect table = all families, **RAW counts, NO singleton filter** (426 taxa) -- required here because Chao1 (`S_obs + f1^2/(2*f2)`) and the coverage-based `ChaoShannon`/`ChaoSimpson` all key on the singleton count `f1`, and a global singleton is by definition a within-sample singleton, so filtering these biases every estimate (measured up to ~2x on the worst-sampled insect samples; corrected 2026-09-09, see the 2026-09-09 entry in `iterative_analysis_updates.md`). The observed-diversity script above was switched to the same raw 426-taxon table the same day, so both Lineage C scripts now share one insect-table construction. Fungal side also raw -- Kingdom==Fungi counts, left UNRAREFIED (extrapolation itself corrects for uneven sequencing depth, so pre-rarefying would discard the rare-tail information the estimator needs). Per-sample granularity (user's choice). Mirrors the alpha-diversity script's cross-community correlation (both observed and asymptotic reported side by side), lure-interaction test, and site/lure/date covariate tests + date-shape classification, the latter two run on asymptotic values only. "simpson" here = Hill q=2 diversity (higher=more even), opposite direction from "simpson_dominance" in the alpha-diversity script -- don't compare directly. See `iterative_analysis_updates.md`'s 2026-09-08 and 2026-09-09 entries for results | `data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table/` | `figures/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table/` |
| `compare_insect_fungi/insect_abundance_fungal_alpha_diversity.full_insect_table.r` | Tests the non-targeted-dispersal-vector hypothesis: total insect ABUNDANCE (individuals/sample, `rowSums()` of the raw 426-taxon table) vs. fungal alpha diversity, per fungal metric, both raw-pooled Spearman and season-residualized (`log1p(abundance)` & fungal metric each residualized on `site + factor(date)`, Spearman on residuals, 9999-perm p); observed + asymptotic fungal diversity, read from the main script's CSVs. Flat null throughout (all \|rho\| <= 0.11, p_perm >= 0.38). Added 2026-09-10 | `data/compare_insects_fungi_alpha_diversity_full_insect_table/insect_abundance_{fungal_alpha_correlation,by_sample,context_correlations}.csv` | `figures/compare_insects_fungi_alpha_diversity_full_insect_table/insect_abundance_vs_fungal_alpha.{png,pdf}` |
| `presentation_items/fig9_alpha_diversity_by_lure.R` | Presentation figure: insect (panel a) + fungal (panel b) alpha diversity vs. collection date, faceted metric x lure, cubic OLS trend lines. Lives directly in `presentation_items/` (not a dedicated Lineage C folder) since alpha diversity has no Lineage B counterpart to disambiguate from | -- (reads existing CSVs) | `figures/presentation_items/fig9_alpha_diversity_by_lure.{png,pdf}` |
| `presentation_items/fig10_insect_fungal_alpha_diversity_correlation.R` | Presentation figure: per-sample insect vs. fungal alpha diversity, one panel per metric, colored by date, pooled Spearman correlation annotated. An earlier lure-faceted version was dropped after the interaction test above found no lure effect | -- (reads existing CSVs) | `figures/presentation_items/fig10_insect_fungal_alpha_diversity_correlation.{png,pdf}` |
| `presentation_items/fig11_asymptotic_diversity_by_lure.R` | Fig9's counterpart for ASYMPTOTIC (Chao-extrapolated, `insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r`) diversity -- same panel a/b (insect/fungal) x metric-row x lure-column layout, cubic OLS trend line. Row-strip metric labels shortened ("Shannon div. (q=1)" etc.) vs. fig9's, since the longer "Asymptotic Shannon diversity (Hill q=1)" wording clipped in the rotated switch="y" strip | -- (reads existing CSVs) | `figures/presentation_items/fig11_asymptotic_diversity_by_lure.{png,pdf}` |
| `presentation_items/fig12_insect_fungal_asymptotic_diversity_correlation.R` | Fig10's counterpart for asymptotic diversity -- same one-panel-per-metric layout, pooled (lure-blind) Spearman correlation annotated, reading the `value_type == "asymptotic"` rows of `asymptotic_diversity_cross_community_correlation.csv` | -- (reads existing CSVs) | `figures/presentation_items/fig12_insect_fungal_asymptotic_diversity_correlation.{png,pdf}` |

### Lineage D (bark-beetle x Ophiostomatales presence/absence co-occurrence -- see above)

All outputs under `data/insect_exploratory/{target_genera,scolytinae}/` and
`figures/insect_exploratory/{target_genera,scolytinae}/`. Full numbers +
interpretation in `lineage_D_insect_fungus_cooccurence_patterns.md`.

| Script | What it does | Data output | Figure output |
|---|---|---|---|
| `insect_exploratory/target_genus_phenology.R` | *Ips* + *Dendroctonus* individuals & species count vs. collection date, faceted by lure, pooled over sites (all insect trap data, not just the fungal-matched 69). Also the species x lure catch totals -- the lures partition the 3 target species (EA -> *D. valens* + *I. grandicollis*; Ips lure -> *I. pini* only) | `data/insect_exploratory/target_genera/target_{genus,species}_catch_by_date_lure.csv`, `target_species_totals_by_lure.csv` | `figures/insect_exploratory/target_genera/target_{genus,species}_phenology_by_lure.{png,pdf}`, `target_species_totals_by_lure.{png,pdf}` |
| `insect_exploratory/scolytinae_species_phenology.R` | All 45 Scolytinae species x date abundance heat map, faceted by lure (rows ordered by abundance-weighted mean date) + subfamily-level individuals/richness per date x lure | `data/insect_exploratory/scolytinae/scolytinae_species_abundance_by_date_lure.csv`, `scolytinae_subfamily_totals_by_date_lure.csv`, `scolytinae_species_phenology_order.csv` | `figures/insect_exploratory/scolytinae/scolytinae_species_abundance_heatmap.{png,pdf}`, `scolytinae_subfamily_totals_by_lure.{png,pdf}` |
| `insect_exploratory/target_species_ophiostomatales_cooccurrence.R` | Jaccard co-occurrence, 3 target species x Ophiostomatales (species/genus/ASV grids), `within_trap` null (1999 perms) + Fisher cross-check, BH within grid. Only FDR hit: *D. valens* x *Raffaelea* (genus, q_perm 0.09) | `data/insect_exploratory/target_genera/beetle_ophiostomatales_jaccard.{species,genus,asv}.csv`, `beetle_ophiostomatales_presence_absence.csv` | `figures/insect_exploratory/target_genera/beetle_ophiostomatales_jaccard.{species,genus,asv}.{png,pdf}` |
| `insect_exploratory/scolytinae_ophiostomatales_cooccurrence.R` | Same screen for all 34 Scolytinae species at >=3-sample prevalence, `within_trap` null. Nothing survives FDR over the larger grid; a block of early-season beetles (*Pityogenes hopkinsi* etc.) goes broadly positive -- the shared-phenology signature the residualized script below is built to test | `data/insect_exploratory/scolytinae/scolytinae_ophiostomatales_jaccard.{species,genus,asv}.csv`, `scolytinae_ophiostomatales_presence_absence.csv` | `figures/insect_exploratory/scolytinae/scolytinae_ophiostomatales_jaccard.{species,genus,asv}.{png,pdf}` |
| `insect_exploratory/scolytinae_ophiostomatales_cooccurrence.residualized.R` | Date-conditional re-test of an 11-beetle set (8-beetle top block by mean genus z + 3 reference rows incl. *Dryocoetes autographus* as a negative control): same observed Jaccard, three nulls (`within_trap` / `within_date` / `within_site_x_date`), genus + species + ASV grids. Season control halves the block; genus-level FDR survivors centre on *Pityogenes hopkinsi*, *Dendroctonus valens* x *Leptographium*, *Heteroborips seriatus* x *Grosmannia*/*Raffaelea*, *Ips grandicollis* x *Leptographium* | `data/insect_exploratory/scolytinae/scolytinae_ophiostomatales_jaccard.date_conditional_comparison.csv`, `..._date_conditional_long.csv` | `figures/insect_exploratory/scolytinae/scolytinae_ophiostomatales_jaccard.date_conditional.{genus,species,asv}.{png,pdf}`, `..._date_conditional.shrinkage.{png,pdf}` |
| `insect_exploratory/cooccurrence_lib.R` | Shared engine -- `jaccard_grid_test(beetle_pa, fungal_pa, perm_mat, grid_name)`: per-pair Jaccard + permutation null + Fisher + BH. Sourced by both `*_ophiostomatales_cooccurrence*.R` scripts | -- | -- |

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
  NOT prevalence-filter, and as of 2026-09-09 do NOT global-singleton-filter
  either: Lineage C's full-insect-table scripts (`insect_fungal_alpha_
  diversity.full_insect_table.r` and `insect_fungal_asymptotic_richness_
  iNEXT.full_insect_table.r`) both now use the **raw 426-taxon** insect
  table. These metrics are direct summaries of how many/how evenly
  distributed taxa are observed per sample; discarding taxa that occur in
  fewer than 5 samples across the whole dataset would systematically
  undercount richness for reasons unrelated to the sample being described.
  The milder `colSums > 1` global singleton filter (drop taxa with a single
  individual across the ENTIRE dataset) *was* kept through 2026-09-08 as
  noise-floor cleanup, but is now also dropped: it is **outright wrong for
  Chao-type asymptotic estimators** (`iNEXT` `ChaoRichness`/`ChaoShannon`/
  `ChaoSimpson`), which estimate undetected taxa from the observed
  singleton/doubleton counts (`f1`, `f2`) -- and a global singleton is
  necessarily a within-sample singleton, so dropping it corrupts the
  quantity the estimator is built on -- and was removed from the observed-
  diversity script too so "observed richness" is defined identically in
  both. (The Chao script had the filter through 2026-09-08; both were fixed
  2026-09-09 -- see `iterative_analysis_updates.md`. Lineage A's frozen
  family-filtered `insect_fungal_alpha_diversity.r` still applies it.)

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
  `t_stat_quad`/`p_perm_quad`/`q_value_quad`/`shape` (new). As of 2026-09-08
  also `beta_std`/`beta_std_quad` -- the fully-standardized partial slope
  `b * sd(x) / sd(y)` (QuantPsyc::lm.beta / effectsize::standardize_
  parameters(method="basic") convention) for the linear and quadratic date
  terms, from the OBSERVED fit only. Added as an effect-size companion to the
  t-statistic and used as the fig3 volcano x-axis (see the fig3 rows / the
  2026-09-08 entry in `iterative_analysis_updates.md`); it involves no
  permutation, so `p_perm`/`q_value` are untouched. `sign(beta_std) ==
  sign(t_stat)`, so `shape` and every direction/peak-timing reading are
  unaffected. The same two columns are emitted by the PCoA-axis per-taxon
  tests (`insect_fungal_coca_analysis.general_workflow.r`'s `test_one_taxon`
  and `insect_fungal_coca_factordate_check.prevalence_filtered_insect.r`) as
  a single `beta_std` for the axis under test. `shape` is
  "hump (peaks mid-season)" / "dip (troughs mid-season)" if the quadratic
  term is significant (quadratic takes priority over linear when both are),
  else "linear increase (late-season)" / "linear decrease (early-season)" if
  only the linear term is, else "no significant date pattern". **Gotcha
  already hit once**: `c(t_stat = obs["t_linear"], ...)` in R mangles the
  name to `"t_stat.t_linear"` when the right-hand side is already a named
  scalar -- always wrap in `unname()` before combining into a new named
  vector, or downstream `r["t_stat"]`-style extraction silently returns NA.

## Untracked large outputs (present locally, not in git)

As of 2026-09-09, four full-grid pairwise-association CSVs are `.gitignore`d
and NOT tracked in git (20-80MB each, over GitHub's soft large-file warning
threshold; they were briefly committed and pushed, then removed from all git
history via `git filter-repo`):

- `data/compare_insects_fungi_pairwise_taxa/fungal_insect_pairwise_full_grid.csv`
- `data/compare_insects_fungi_pairwise_taxa/fungal_insect_pairwise_full_grid.date_site_residualized.csv`
- `data/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect/fungal_insect_pairwise_full_grid.csv`
- `data/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect/fungal_insect_pairwise_full_grid.date_site_residualized.csv`

These are the exhaustive all-pairs output of the pairwise taxon-association
screens (`insect_fungal_pairwise_taxon_association.r` /
`.residualized.r` for Lineage A,
`insect_fungal_pairwise_taxon_association.prevalence_filtered_insect.r` /
`.residualized.prevalence_filtered_insect.r` for Lineage B) -- fully
regenerable by rerunning those scripts, so losing them from git costs
nothing but rerun time. Anyone who clones this repo fresh (including via
GitHub) will NOT have these four files and must rerun the corresponding
script to reproduce them locally. Every other output derived from them
(the BH-corrected summary CSVs, `original_vs_residualized_comparison.csv`,
the flagged-taxon figures) is small and remains tracked normally.

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
  shape classification), matching the same convention. §1-6 cover
  `insect_fungal_alpha_diversity.full_insect_table.r` and its two
  presentation figures (fig9, fig10); §7-12 (added 2026-09-08) extend this
  to the ASYMPTOTIC (Chao-extrapolated, `iNEXT`) diversity variant --
  `insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r` and fig11/
  fig12 -- with explicit observed-vs-asymptotic comparison throughout. Since
  Lineage C has no ordination/CoCA/Procrustes/pairwise-screen analog
  (alpha-diversity-only by design, see the Lineage C section above), this
  document's scope is these two scripts and their four figures -- update it
  if any of them change, rather than retaining outdated numbers.
- `lineage_D_insect_fungus_cooccurence_patterns.md` -- the Lineage D
  (bark-beetle x Ophiostomatales presence/absence co-occurrence) counterpart
  to the above: top-line numbers for the phenology figures, the Jaccard
  co-occurrence screens (3 target species + all Scolytinae, `within_trap`
  null), and the date-conditional re-test of the top block (`within_date` /
  `within_site_x_date` nulls, genus/species/ASV). Matching the same
  convention as the A/B/C interp docs -- numbered sections, tables, a bottom
  line, an underlying-files list. Update it when any of the six
  `insect_exploratory/` Lineage D scripts change. This line is still
  exploratory; the results in the doc are candidates (nothing survives FDR
  below genus resolution) pending an abundance-based or larger-sample
  follow-up.
