# Interpretation update (2026-09-01, extended 2026-09-08, §7-12 corrected 2026-09-09): top-line numbers for the full-insect-table alpha-diversity lineage (Lineage C)

This file is the Lineage C (full insect table, per `project_organization.md`) counterpart to `lineage_A_family_filtered_insect_top_level_interepretation.md` and `lineage_B_prevalence_filtered_insect_top_level_interpretation.md`. Same purpose: pull current, on-disk numbers into one place rather than requiring anyone to reconstruct them from the chat-style narrative in `iterative_analysis_updates.md`. Read `project_organization.md` first if you haven't -- it defines the three insect-table constructions this file assumes you already know about.

**Table construction reminder**: Lineage C = insect trap-catch table, ALL families (not just Curculionidae+Latridiidae), **raw counts** -- deliberately NOT the >=5-sample prevalence filter Lineage B uses, and (as of 2026-09-09) NOT the `colSums > 1` global singleton filter either. **426 taxa** across 98 families, 69 samples before matching to fungal data (68 after the fungal rarefaction step drops one low-depth sample; the asymptotic script keeps all 69 since it leaves the fungal side unrarefied). Both Lineage C scripts -- `insect_fungal_alpha_diversity.full_insect_table.r` (§1-6, observed) and `insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r` (§7-12, asymptotic) -- use this same raw table as of 2026-09-09; through 2026-09-08 both applied the singleton filter (277 taxa), which is outright wrong for the Chao estimators §7-12 use and was dropped from §1-6 too for consistency (see §1's note, project_organization.md's "When to prevalence-filter vs. not", and the 2026-09-09 log entry). **Lineage C is alpha-diversity-only by design** -- there is no Lineage C ordination/CoCA/Procrustes/pairwise-taxon-screen analog the way Lineage A and B both have; don't look for one.

§1-6 below come from `compare_insect_fungi/insect_fungal_alpha_diversity.full_insect_table.r` and the two presentation figures built on top of it (`presentation_items/fig9_alpha_diversity_by_lure.R`, `fig10_insect_fungal_alpha_diversity_correlation.R`), all using OBSERVED diversity (no extrapolation). §7-12 (added 2026-09-08; numbers corrected 2026-09-09 -- see below) extend this to ASYMPTOTIC (Chao-extrapolated, via the `iNEXT` package) diversity, from `compare_insect_fungi/insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r` and its two presentation figures (`presentation_items/fig11_asymptotic_diversity_by_lure.R`, `fig12_insect_fungal_asymptotic_diversity_correlation.R`) -- same permutation-testing machinery, but estimating the richness/diversity that WOULD be seen with infinite additional sampling effort rather than counting only what was observed, leaving the fungal side unrarefied since the asymptotic estimator itself is the correction for uneven sampling effort that rarefaction exists to provide elsewhere, and (corrected 2026-09-09) feeding the estimators a raw, unfiltered 426-taxon insect table -- the `colSums > 1` singleton filter both Lineage C scripts used through 2026-09-08 corrupts the very singleton/doubleton counts Chao estimators are built on, and was dropped from the §1-6 observed script the same day so both scripts share one table. All §7-12 numbers below are the post-correction values; see the 2026-09-09 entry in `iterative_analysis_updates.md` for the before/after and `project_organization.md`'s Lineage C section + "When to prevalence-filter" exception for the rationale.

---

## 1. Alpha-diversity value ranges (descriptive)

| Metric | Insect (full table, n=69) | Fungal (rarefied, n=68) |
|---|---|---|
| Shannon diversity | 0.96 - 3.45 (mean 2.65) | 2.07 - 5.61 (mean 3.90) |
| Simpson dominance | 0.045 - 0.648 (mean 0.155) | 0.009 - 0.466 (mean 0.115) |
| Richness (observed taxa) | 6 - 75 (mean 34.8) | 271 - 947 (mean 529) |

Insect and fungal diversity live on very different scales for every metric (most obviously richness: tens of insect taxa vs. hundreds of fungal ASVs per sample) -- this is why fig9/fig10 both need independent per-panel y-axes rather than a shared scale (see §4 below and each script's header comment for the two different ways that was solved).

**Insect-table change (2026-09-09):** these insect numbers are computed on the raw, unfiltered 426-taxon table -- the `colSums > 1` global singleton filter was dropped from `insect_fungal_alpha_diversity.full_insect_table.r` for consistency with the asymptotic script (§7), which cannot use it. Effect on this table: observed richness rose in 54 of 69 samples (median +2, max +8 taxa/sample: 4-68/mean 32.7 -> 6-75/mean 34.8); Shannon and Simpson dominance barely moved (global singletons carry negligible weight). See the 2026-09-09 entry in `iterative_analysis_updates.md`. Lineage A's frozen `insect_fungal_alpha_diversity.r` still applies the filter and is unchanged.

## 2. Insect vs. fungal alpha-diversity correlation (pooled, lure-blind Spearman)

The alpha-diversity analog of the whole-community Procrustes check used for community composition under Lineage A/B.

| Metric | rho | p-value | n |
|---|---|---|---|
| Shannon diversity | -0.268 | **0.027** (sig) | 68 |
| Simpson dominance | -0.131 | 0.285 (n.s.) | 68 |
| Richness | -0.118 | 0.336 (n.s.) | 68 |

Only Shannon diversity shows a significant relationship, and it's negative: samples with more diverse insect catches tend to have *less* diverse fungal communities, not more. This is a markedly weaker relationship than the whole-community COMPOSITION correlation (Procrustes r=0.604, p=0.001 under Lineage B) -- alpha diversity and beta diversity/composition are answering different questions here and needn't agree; a sample's insect and fungal communities can co-vary strongly in composition/timing while their raw diversity levels are only weakly (or not) linked.

## 3. Does lure modify the insect~fungal relationship?

Added specifically to check whether fig10 (below) should be split by lure. Model: `fungal_y ~ site + date + insect_x * lure` vs. the same model without the interaction, tested via the project's standard lure permutation (999 perms, lure shuffled among traps within site -- same scheme `test_term_lure` uses in §4 below and `insect_pcoa_lure_association.r` uses project-wide).

| Metric | F | p_perm |
|---|---|---|
| Shannon diversity | 0.152 | 0.900 |
| Simpson dominance | 0.321 | 0.812 |
| Richness | 0.056 | 0.829 |

**No metric shows any evidence of a lure-dependent insect~fungal relationship** (all F<0.35, all p_perm>=0.80). fig10 was built and then deliberately collapsed from a lure-faceted grid to a single pooled panel per metric on the strength of this result -- splitting the earlier draft by lure was cutting n roughly into thirds (~23/lure) without visualizing a real effect.

## 4. Within-community site/lure/date effects on each metric

Same trap-level permutation logic used throughout this project (date permuted within trap; lure permuted among traps within site; site permuted freely among all traps), applied here to a univariate diversity metric via `lm()` instead of a per-taxon abundance or distance matrix. Date is tested with linear + quadratic + cubic terms in one model (§5 covers the shape read); only the significant (p_perm<0.10) terms are shown below -- see `alpha_diversity_covariate_tests.csv` for the full table (30 rows: 2 communities x 3 metrics x 5 terms).

| Community | Metric | Term | Stat | p_perm |
|---|---|---|---|---|
| insect | Shannon diversity | lure | F=7.39 | **0.008** |
| insect | Shannon diversity | site | F=5.73 | **0.042** |
| insect | Shannon diversity | date (linear) | t=1.88 | 0.067 |
| insect | Simpson dominance | site | F=3.69 | 0.079 |
| insect | Richness | date (quadratic) | t=-5.41 | **0.001** |
| insect | Richness | date (cubic) | t=3.17 | **0.002** |
| insect | Richness | lure | F=13.16 | **0.008** |
| fungal | Shannon diversity | date (linear) | t=-4.40 | **0.001** |
| fungal | Shannon diversity | date (quadratic) | t=3.97 | **0.002** |
| fungal | Shannon diversity | date (cubic) | t=2.87 | 0.010 |
| fungal | Simpson dominance | date (linear) | t=3.56 | **0.001** |
| fungal | Simpson dominance | date (quadratic) | t=-3.09 | **0.004** |
| fungal | Simpson dominance | date (cubic) | t=-2.71 | 0.006 |
| fungal | Simpson dominance | site | F=3.85 | 0.046 |
| fungal | Richness | date (linear) | t=-4.79 | **0.001** |
| fungal | Richness | date (quadratic) | t=3.04 | 0.006 |
| fungal | Richness | date (cubic) | t=3.65 | **0.001** |
| fungal | Richness | site | F=2.85 | 0.025 |

Two clean patterns:
- **Lure is never a significant term on the fungal side** (p_perm 0.654-0.962 across all 3 fungal metrics, not shown above) but IS significant for insect Shannon diversity (p=0.008) and is the single largest effect size anywhere in this table for insect richness (F=13.16, p=0.008). Consistent with lure being an insect-trap attractant with a real effect on how diverse a trap's catch is, but no detectable knock-on effect on fungal diversity -- and consistent with §3's null interaction result.
- **Date dominates the fungal side completely** -- all 3 fungal metrics are significant on all 3 polynomial date terms (9/9 significant date tests). The insect side is patchier: richness has a strong quadratic+cubic date structure, insect Shannon diversity a weak late-season linear increase (t=1.88, p=0.067 -- this crossed the 0.10 line when the insect table was switched to raw counts on 2026-09-09, from p=0.103), and only insect Simpson dominance shows no date effect at any polynomial order (see §5).
- **Site is a real but always secondary effect** -- present at borderline-to-modest significance for 4 of the 6 community x metric combinations (insect Shannon p=0.042, insect Simpson p=0.079, fungal Simpson p=0.046, fungal richness p=0.025), never the largest term for any of them.

## 5. Date shape classification (linear vs. quadratic vs. cubic)

Extends this project's established linear+quadratic date-test convention (`fungal_community_seasonality.r`, `insect_fungal_coca_analysis.general_workflow.r`) one order further to cubic -- added specifically because fig9's trend lines showed more structure (e.g. an S-curve in insect Shannon diversity under the Ethanol lure) than a single hump/dip could capture. Shape classification cascades: a significant cubic term takes priority over quadratic, which takes priority over linear.

| Community | Metric | Shape | Linear t (p) | Quadratic t (p) | Cubic t (p) |
|---|---|---|---|---|---|
| insect | Shannon diversity | **linear increase (late-season)** | 1.88 (0.067) | -1.27 (0.208) | 0.39 (0.726) |
| insect | Simpson dominance | no significant date pattern | -1.42 (0.181) | 0.38 (0.734) | 0.38 (0.717) |
| insect | Richness | **complex (cubic significant)** | -1.03 (0.291) | -5.41 (0.001) | 3.17 (0.002) |
| fungal | Shannon diversity | **complex (cubic significant)** | -4.40 (0.001) | 3.97 (0.002) | 2.87 (0.010) |
| fungal | Simpson dominance | **complex (cubic significant)** | 3.56 (0.001) | -3.09 (0.004) | -2.71 (0.006) |
| fungal | Richness | **complex (cubic significant)** | -4.79 (0.001) | 3.04 (0.006) | 3.65 (0.001) |

**4 of the 6 community x metric date relationships are cubic-significant** (insect richness + all 3 fungal metrics); insect Shannon diversity shows a plain late-season linear increase (linear p=0.067, no higher-order term significant -- this crossed the 0.10 line with the 2026-09-09 raw-count switch, from p=0.103); only insect Simpson dominance shows no significant date structure at all, at any polynomial order. All three fungal metrics are cubic, meaning the simple "fungal diversity declines through the season" linear read (itself also significant, and in the expected direction: Shannon/richness falling, Simpson dominance rising) is an incomplete description -- there's a genuine dip-then-partial-recovery shape, matching fig9's visible mid-June trough for fungal Shannon diversity and richness across all three lures. Insect richness shows the same kind of higher-order structure (t=3.17 on the cubic term, p_perm=0.002) despite insect Simpson showing nothing and insect Shannon showing only a linear trend -- alpha-diversity metrics computed from the SAME table don't always agree on whether there's a date effect at all, let alone its shape.

## 6. Presentation figures

- **`fig9_alpha_diversity_by_lure.R`** -- insect (panel a) and fungal (panel b) alpha diversity vs. collection date, faceted metric (rows) x lure (columns), points colored by site, cubic OLS trend line per panel (`y ~ poly(date, 3)`) matching the permutation-tested cubic date model in §5. Output: `figures/presentation_items/fig9_alpha_diversity_by_lure.{png,pdf}`.
- **`fig10_insect_fungal_alpha_diversity_correlation.R`** -- per-sample insect vs. fungal alpha diversity, one panel per metric (NOT split by lure, per §3's null interaction result), colored by collection date, pooled Spearman rho/p from §2 annotated in each panel. Output: `figures/presentation_items/fig10_insect_fungal_alpha_diversity_correlation.{png,pdf}`.

## 7. Asymptotic (Chao-extrapolated) diversity: how much more is out there?

`insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r` estimates, per sample, the richness/Shannon(Hill q=1)/Simpson(Hill q=2) diversity expected with infinite additional sampling effort, using `iNEXT`'s Chao-type estimators (`ChaoRichness()`/`ChaoShannon()`/`ChaoSimpson()`). n=69 for both communities here (unlike §1-6's fungal n=68 -- the fungal side is unrarefied in this variant, so no rarefaction-depth sample drop). "Simpson" is Hill number q=2 DIVERSITY (higher = more even), the opposite direction from §1-6's Simpson DOMINANCE -- don't compare the two scripts' Simpson columns directly.

**Insect table note (corrected 2026-09-09):** this script feeds the estimators the **raw, unfiltered 426-taxon** insect table -- no singleton filter -- because Chao1 and the coverage-based `ChaoShannon`/`ChaoSimpson` all key on the observed singleton count `f1`, which the `colSums > 1` filter (used by both Lineage C scripts through 2026-09-08) corrupts (a global singleton is by definition a within-sample singleton). The §1-6 observed script was switched to the same raw table the same day, so this script's insect `richness_observed` (6-75, mean 34.8) once again matches §1's observed richness exactly. Both differ from the pre-2026-09-09 singleton-filtered values (4-68, mean 32.7); see the 2026-09-09 entry in `iterative_analysis_updates.md` for the before/after. Fungal-side values still do not match §1-6's fungal numbers (unrarefied here vs. rarefied-and-averaged there).

| Metric | Insect observed | Insect asymptotic | Fungal observed (unrarefied) | Fungal asymptotic |
|---|---|---|---|---|
| Richness | 6 - 75 (mean 34.8) | 13.7 - 236.8 (mean 81.6) | 364 - 4,794 (mean 1,723) | 490.7 - 5,546.7 (mean 2,053) |
| Shannon (Hill q=1) | 2.60 - 31.35 (mean 16.0) | 3.21 - 84.64 (mean 26.5) | 8.49 - 345.35 (mean 82.3) | 8.60 - 348.06 (mean 83.9) |
| Simpson (Hill q=2) | 1.54 - 22.26 (mean 9.33) | 1.55 - 30.92 (mean 11.1) | 2.15 - 113.75 (mean 20.4) | 2.15 - 113.82 (mean 20.4) |

Sample coverage (`iNEXT::DataInfo()`'s SC, the fraction of the true community estimated to be captured) tells the same story from a different angle: **insect coverage is 0.07 - 0.98 (mean 0.77)** -- trap catches (median ~100 individuals per sample) are genuinely incomplete samples of the local insect community, and the worst-sampled ones are very incomplete -- while **fungal coverage is 0.97 - 1.00 (mean 0.99)** -- amplicon sequencing depth (thousands to hundreds of thousands of reads per sample) is already close to exhaustive.

That asymmetry drives everything below: **median asymptotic/observed ratio is ~2.15x for insect richness but only ~1.20x for fungal richness** -- extrapolation roughly doubles the insect richness estimate but barely nudges the fungal one. Shannon (q=1, weighted toward common taxa) is far less sensitive to the correction than richness for either community (median ratio 1.44x insect, 1.01x fungal); Simpson (q=2, weighted even more toward dominant taxa) barely moves at all (1.08x insect, 1.00x fungal) -- expected, since higher-order Hill numbers are increasingly dominated by the already-well-sampled common taxa and increasingly insensitive to how many rare/undetected ones exist.

## 8. Does extrapolating change the insect~fungal correlation conclusion? No.

Same Spearman test as §2, computed on both observed and asymptotic values within the same script (a clean, same-metric-definition comparison unlike some of the §10 comparisons below):

| Metric | Observed rho (p) | Asymptotic rho (p) | Conclusion |
|---|---|---|---|
| Richness | 0.081 (0.511, n.s.) | 0.006 (0.962, n.s.) | Unchanged -- no relationship either way |
| Shannon (Hill q=1) | -0.258 (0.032, **sig**) | -0.258 (0.032, **sig**) | Unchanged -- same modest negative relationship |
| Simpson (Hill q=2) | -0.126 (0.300, n.s.) | -0.080 (0.513, n.s.) | Unchanged -- no relationship either way |

**Extrapolating to the asymptote does not change any of the three §2-style conclusions** -- correcting for undetected taxa (which, per §7, matters a lot more for insect richness than for anything else in this table) doesn't meaningfully move the correlation with fungal diversity. The weak-Shannon/null-richness-and-Simpson pattern from §2 is robust to whether "diversity" means observed counts or extrapolated estimates.

## 9. Does lure modify the insect~fungal relationship, under extrapolation? Still no for richness/Simpson -- but Shannon is now a nominally significant interaction

Same interaction test as §3 (`fungal ~ site+date+insect*lure`, 999-perm F-test on the interaction term), applied to asymptotic values:

| Metric | F | p_perm | vs. §3 (observed: F=0.02-0.34, p=0.80-0.95) |
|---|---|---|---|
| Richness | 0.008 | 0.994 | Consistent -- still clearly n.s. |
| Shannon (Hill q=1) | 0.694 | **0.013** | **Significant even at p<0.05** -- was p=0.913 (n.s.) in §3 |
| Simpson (Hill q=2) | 0.251 | 0.754 | Consistent -- still clearly n.s. |

Worth flagging plainly rather than glossing over: the asymptotic Shannon interaction (p_perm=0.013) is significant by any of this project's conventions, unlike every other lure-interaction test run anywhere in this document (both observed §3 and the other two asymptotic metrics here). It strengthened from p_perm=0.075 to ~0.01 when the insect table was corrected to raw counts on 2026-09-09 (§7 note; the corrected insect asymptotic Shannon values changed materially for several samples). Reasons to still hold it loosely: it is 1 of 3 uncorrected interaction tests, not a taxon screen; the F-statistic itself is small (0.694); and it rides entirely on the insect asymptotic Shannon estimate, which is built on poorly-sampled insect data (mean coverage 0.77, one q1>q0 Hill-ordering violation at NH.99 -- a mathematically impossible result that only over-extrapolation produces). So: a real nominal signal, stronger than before, but resting on the shakiest quantity in this analysis. fig12 (§12) still shows the pooled (non-lure-split) relationship; whether asymptotic Shannon specifically now warrants a lure-faceted panel is a live question for the next revision of this figure, and would be worth re-checking against a version of the interaction test that down-weights or excludes the lowest-coverage insect samples.

## 10. Within-community site/lure/date effects, asymptotic values

Same trap-level permutation logic as §4, run here on the asymptotic estimates only (observed-value covariate testing under these exact metric definitions -- Hill-number Shannon/Simpson, unrarefied fungal counts -- isn't separately re-derived; §4 already covers observed-value significance under its own, different, metric definitions). Only significant/borderline (p_perm<0.10) terms shown -- see `asymptotic_diversity_covariate_tests.csv` for the full 30-row table (2 communities x 3 metrics x 5 terms).

| Community | Metric | Term | Stat | p_perm |
|---|---|---|---|---|
| insect | Richness | date (quadratic) | t=-3.17 | **0.001** |
| insect | Shannon | site | F=9.28 | **0.015** |
| insect | Shannon | date (linear) | t=2.17 | **0.037** |
| insect | Simpson | date (quadratic) | t=2.06 | **0.050** |
| fungal | Richness | site | F=1.16 | **0.029** |
| fungal | Shannon | date (linear) | t=-4.28 | **0.001** |
| fungal | Shannon | date (quadratic) | t=3.23 | **0.003** |
| fungal | Shannon | date (cubic) | t=2.39 | **0.026** |
| fungal | Simpson | date (linear) | t=-3.00 | **0.002** |
| fungal | Simpson | date (quadratic) | t=2.78 | **0.005** |

(The `insect Shannon ~ lure` term was p_perm=0.070 here before the 2026-09-09 raw-count correction -- it dropped to F=2.60, p_perm~0.22 on the corrected insect table, i.e. it was an artifact of the singleton filter, and is no longer in this table. Two `~ site` terms sit just outside the table straddling the 0.10 line on permutation noise alone -- `fungal Simpson ~ site` (F=1.92, p_perm 0.09-0.10 across re-runs) and `fungal Shannon ~ site` (F=1.93, p_perm ~0.11); their F-statistics are deterministic and unchanged, only the permutation p wobbles. See the 2026-09-09 entry in `iterative_analysis_updates.md`.)

**The insect-richness lure/date comparison with §4**: as of the 2026-09-09 raw-count switch, both this script and §4's use the same raw 426-taxon insect table, so §4's observed insect richness and this script's `richness_observed` are once again the same numbers -- a clean before/after extrapolation comparison. §4 (observed richness) found a strong lure effect (F=13.16, **p=0.008**) and a significant cubic date term (t=3.17, **p=0.002**) for insect richness. Here (asymptotic richness), lure is nowhere near significant (F=0.30, p=0.807) and the cubic term is gone (t=1.22, p=0.219) -- only the quadratic (hump) term survives (t=-3.17, p=0.001). **Extrapolating to account for undetected insect taxa attenuates the lure effect on richness and simplifies its date shape from cubic to a plain mid-season hump.** A plausible reading: part of the observed lure-richness effect may reflect unequal SAMPLING COMPLETENESS across lures (some lures catch more individuals per trap, pushing observed richness up mechanically) rather than a true difference in the underlying taxon pool each lure attracts -- exactly the kind of confound asymptotic estimation exists to correct for. (Shannon/Simpson can't be compared this cleanly to §4, since their Hill-number definitions here differ from §4's raw entropy/dominance metrics, and the fungal side is unrarefied here vs. rarefied-and-averaged in §4 -- treat any resemblance or difference between those rows and §4's as coincidental, not evidence either way.)

## 11. Date shape classification, asymptotic values

Same linear/quadratic/cubic cascading classification as §5 (cubic takes priority over quadratic, which takes priority over linear), applied to asymptotic values:

| Community | Metric | Shape | Linear t (p) | Quadratic t (p) | Cubic t (p) |
|---|---|---|---|---|---|
| insect | Richness | **hump (peaks mid-season)** | 0.40 (0.695) | -3.17 (0.001) | 1.22 (0.236) |
| insect | Shannon | **linear increase (late-season)** | 2.17 (0.037) | 0.99 (0.339) | -0.25 (0.800) |
| insect | Simpson | **dip (troughs mid-season)** | 1.36 (0.183) | 2.06 (0.050) | -0.10 (0.923) |
| fungal | Richness | no significant date pattern | -0.86 (0.403) | -1.06 (0.319) | 1.06 (0.310) |
| fungal | Shannon | **complex (cubic significant)** | -4.28 (0.001) | 3.23 (0.003) | 2.39 (0.026) |
| fungal | Simpson | **dip (troughs mid-season)** | -3.00 (0.002) | 2.78 (0.005) | 1.29 (0.217) |

Three changes worth noting relative to §5's observed-value classification: **insect richness simplifies from "complex (cubic significant)" to a plain "hump"** (the same cubic-term loss already flagged in §10); **insect Shannon is a late-season linear increase in both** (§5's observed insect Shannon also classifies as "linear increase (late-season)" after the 2026-09-09 raw-count switch, so this is no longer a divergence -- both the observed and asymptotic insect Shannon carry a significant positive linear date term and nothing higher-order, asymptotic linear p=0.037); and **fungal richness loses its date signal entirely** (§5: "complex (cubic significant)," all 3 terms p<=0.03; here: no term reaches p<0.10) -- on the unrarefied, extrapolated fungal table, richness's date pattern doesn't hold up, even though fungal Shannon and Simpson's date patterns (still complex/dip respectively) do. Fungal Shannon remains the only "complex" (genuinely non-monotonic, all 3 terms significant) result in this table, same as §5.

## 12. Presentation figures (asymptotic)

- **`fig11_asymptotic_diversity_by_lure.R`** -- fig9 reproduced one-for-one on asymptotic values: insect (panel a) and fungal (panel b) asymptotic diversity vs. collection date, faceted metric (rows) x lure (columns), cubic OLS trend lines. Row-strip metric labels shortened ("Shannon div. (q=1)" etc.) vs. fig9's wording, which clipped in the rotated strip. Output: `figures/presentation_items/fig11_asymptotic_diversity_by_lure.{png,pdf}`.
- **`fig12_insect_fungal_asymptotic_diversity_correlation.R`** -- fig10 reproduced one-for-one on asymptotic values: per-sample insect vs. fungal asymptotic diversity, one panel per metric, pooled (not lure-split, per §9) Spearman rho/p from §8 annotated. Points colored by collection date AND shaped by lure (filled circle=Ethanol, triangle=Alpha-pinene_EtOH, square=Ips) -- descriptive only, doesn't affect the pooled trend line/correlation; this shape mapping was also back-ported to fig10 for consistency between the two figures. Output: `figures/presentation_items/fig12_insect_fungal_asymptotic_diversity_correlation.{png,pdf}`.

## Bottom line

1. **Fungal diversity shows strong, comprehensive date effects** (observed values, §1-6) -- all 3 metrics (Shannon, Simpson dominance, richness) are significant on all 3 polynomial date terms (linear, quadratic, AND cubic), and every one of them is a genuinely non-monotonic (cubic) shape, not a simple increase or decrease.
2. **Insect diversity's date signal is metric-specific and much weaker** -- richness shows a strong non-monotonic (cubic) date effect and Shannon diversity a weak late-season linear increase (p=0.067); only Simpson dominance shows no significant date effect at all under this full, all-families table. This contrasts with the strong insect~date signal already established at the COMMUNITY-COMPOSITION level under both Lineage A and B (PERMANOVA date R² ~13-15%) -- alpha diversity and community composition are different questions and don't have to move together.
3. **Lure drives insect diversity but not fungal diversity** -- significant for insect Shannon (p=0.008) and the largest effect size in the whole covariate table for insect richness (F=13.16, p=0.008), but never significant for any fungal metric (p=0.654-0.962) -- and §3 confirms lure doesn't even modify the insect~fungal relationship itself. Consistent with lure acting as an insect-trap attractant with no direct fungal-diversity consequence.
4. **Site is a real but always secondary effect** on both communities -- present at borderline-to-modest significance for roughly half the community x metric combinations, never the dominant term for any of them.
5. **Insect and fungal alpha diversity are only weakly linked, and only for one metric**: Shannon diversity shows a modest, significant negative correlation (rho=-0.27, p=0.03 -- more insect diversity co-occurs with LESS fungal diversity); Simpson dominance and richness show no significant correlation. This is a much weaker relationship than the strong whole-community Procrustes correlation for composition (r=0.60, p=0.001 under Lineage B) -- these two communities' raw diversity levels are far less tightly coupled than their compositions are.
6. **Lure does not modify the insect~fungal relationship for any metric under observed values** (all interaction p_perm>=0.80) -- fig10 shows the pooled, non-lure-faceted relationship rather than an artificially data-starved 3-way split.
7. **Extrapolating to asymptotic (Chao-estimated) diversity leaves the headline insect~fungal correlation conclusion unchanged** (§8) -- same weak-Shannon/null-richness-and-Simpson pattern whether using observed counts or diversity corrected for undetected taxa.
8. **How much "extra" diversity extrapolation adds is wildly uneven, and it's mostly about insect richness** (§7) -- insect sample coverage (mean 0.77) is far lower than fungal (mean 0.99), so insect richness roughly doubles on extrapolation (median ratio ~2.15x) while fungal richness moves only ~20%; insect Shannon moves moderately (~1.4x) and fungal Shannon barely (~1.0x); Simpson (q=2) barely moves for either community. Insect trap catches are genuinely incomplete community samples in a way fungal amplicon sequencing isn't.
9. **The one place extrapolation changes a covariate conclusion**: insect richness's strong lure effect (F=13.16, p=0.008) and cubic date shape under observed counts (§4) attenuate to no significant lure effect (F=0.30, p=0.807) and a plain quadratic hump (§10-11) once corrected for undetected taxa -- suggesting the observed lure effect on richness may partly reflect unequal sampling completeness across lures rather than a genuinely different taxon pool per lure. (A separate, spurious `insect Shannon ~ lure` term at p_perm=0.070 in the pre-2026-09-09 version of §10 turned out to be an artifact of the singleton filter and is gone on the corrected table.)
10. **One signal that is no longer borderline**: the lure x insect-diversity interaction on the fungal side, null everywhere else in this document (§3, and 2 of 3 asymptotic metrics in §9), is significant for asymptotic Shannon at p_perm=0.013 (it was p=0.075 before the 2026-09-09 raw-count correction). Still 1 of 3 uncorrected tests, still resting on the poorly-sampled (mean coverage 0.77, one Hill-ordering violation) insect asymptotic Shannon estimate, and not yet reflected in fig12's still-pooled design -- but a stronger nominal signal than previously recorded, and now worth actively deciding whether fig12 needs a lure-faceted Shannon panel (§9).

Underlying files: §1-6 -- `data/compare_insects_fungi_alpha_diversity_full_insect_table/insect_alpha_diversity.csv`, `fungal_alpha_diversity.csv`, `alpha_diversity_cross_community_correlation.csv`, `alpha_diversity_insect_fungal_interaction_by_lure.csv`, `alpha_diversity_covariate_tests.csv`, `alpha_diversity_date_shape.csv`; figures in `figures/compare_insects_fungi_alpha_diversity_full_insect_table/` (descriptive) and `figures/presentation_items/fig9_alpha_diversity_by_lure.{png,pdf}`, `fig10_insect_fungal_alpha_diversity_correlation.{png,pdf}`. §7-12 -- `data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table/insect_asymptotic_diversity.csv`, `fungal_asymptotic_diversity.csv`, `asymptotic_diversity_cross_community_correlation.csv`, `asymptotic_diversity_insect_fungal_interaction_by_lure.csv`, `asymptotic_diversity_covariate_tests.csv`, `asymptotic_diversity_date_shape.csv`; figures in `figures/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table/` (descriptive) and `figures/presentation_items/fig11_asymptotic_diversity_by_lure.{png,pdf}`, `fig12_insect_fungal_asymptotic_diversity_correlation.{png,pdf}`.
