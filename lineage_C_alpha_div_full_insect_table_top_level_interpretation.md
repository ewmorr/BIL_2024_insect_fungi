# Interpretation update (2026-09-01, extended 2026-09-08): top-line numbers for the full-insect-table alpha-diversity lineage (Lineage C)

This file is the Lineage C (full insect table, per `project_organization.md`) counterpart to `lineage_A_family_filtered_insect_top_level_interepretation.md` and `lineage_B_prevalence_filtered_insect_top_level_interpretation.md`. Same purpose: pull current, on-disk numbers into one place rather than requiring anyone to reconstruct them from the chat-style narrative in `iterative_analysis_updates.md`. Read `project_organization.md` first if you haven't -- it defines the three insect-table constructions this file assumes you already know about.

**Table construction reminder**: Lineage C = insect trap-catch table, ALL families (not just Curculionidae+Latridiidae), **singleton-taxon filter only** -- deliberately NOT the >=5-sample prevalence filter Lineage B uses. 277 taxa across 98 families, 69 samples before matching to fungal data (68 after the fungal rarefaction step drops one low-depth sample). This is the correct all-families construction *specifically for alpha-diversity metrics* -- prevalence-filtering before computing richness/Shannon/Simpson would discard real taxa and artificially deflate diversity (see project_organization.md's "When to prevalence-filter vs. not"). **Lineage C is alpha-diversity-only by design** -- there is no Lineage C ordination/CoCA/Procrustes/pairwise-taxon-screen analog the way Lineage A and B both have; don't look for one.

§1-6 below come from `compare_insect_fungi/insect_fungal_alpha_diversity.full_insect_table.r` and the two presentation figures built on top of it (`presentation_items/fig9_alpha_diversity_by_lure.R`, `fig10_insect_fungal_alpha_diversity_correlation.R`), all using OBSERVED diversity (no extrapolation). §7-12 (added 2026-09-08) extend this to ASYMPTOTIC (Chao-extrapolated, via the `iNEXT` package) diversity, from `compare_insect_fungi/insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r` and its two presentation figures (`presentation_items/fig11_asymptotic_diversity_by_lure.R`, `fig12_insect_fungal_asymptotic_diversity_correlation.R`) -- same full insect table, same permutation-testing machinery, but estimating the richness/diversity that WOULD be seen with infinite additional sampling effort rather than counting only what was observed, and (a deliberate departure from §1-6) leaving the fungal side unrarefied since the asymptotic estimator itself is the correction for uneven sampling effort that rarefaction exists to provide elsewhere. See that script's header comment and `project_organization.md`'s Lineage C section for the full methodological rationale.

---

## 1. Alpha-diversity value ranges (descriptive)

| Metric | Insect (full table, n=69) | Fungal (rarefied, n=68) |
|---|---|---|
| Shannon diversity | 0.87 - 3.38 (mean 2.59) | 2.07 - 5.61 (mean 3.90) |
| Simpson dominance | 0.047 - 0.673 (mean 0.161) | 0.009 - 0.466 (mean 0.115) |
| Richness (observed taxa) | 4 - 68 (mean 32.7) | 271 - 947 (mean 529) |

Insect and fungal diversity live on very different scales for every metric (most obviously richness: tens of insect taxa vs. hundreds of fungal ASVs per sample) -- this is why fig9/fig10 both need independent per-panel y-axes rather than a shared scale (see §4 below and each script's header comment for the two different ways that was solved).

## 2. Insect vs. fungal alpha-diversity correlation (pooled, lure-blind Spearman)

The alpha-diversity analog of the whole-community Procrustes check used for community composition under Lineage A/B.

| Metric | rho | p-value | n |
|---|---|---|---|
| Shannon diversity | -0.264 | **0.030** (sig) | 68 |
| Simpson dominance | -0.121 | 0.327 (n.s.) | 68 |
| Richness | -0.119 | 0.336 (n.s.) | 68 |

Only Shannon diversity shows a significant relationship, and it's negative: samples with more diverse insect catches tend to have *less* diverse fungal communities, not more. This is a markedly weaker relationship than the whole-community COMPOSITION correlation (Procrustes r=0.604, p=0.001 under Lineage B) -- alpha diversity and beta diversity/composition are answering different questions here and needn't agree; a sample's insect and fungal communities can co-vary strongly in composition/timing while their raw diversity levels are only weakly (or not) linked.

## 3. Does lure modify the insect~fungal relationship?

Added specifically to check whether fig10 (below) should be split by lure. Model: `fungal_y ~ site + date + insect_x * lure` vs. the same model without the interaction, tested via the project's standard lure permutation (999 perms, lure shuffled among traps within site -- same scheme `test_term_lure` uses in §4 below and `insect_pcoa_lure_association.r` uses project-wide).

| Metric | F | p_perm |
|---|---|---|
| Shannon diversity | 0.148 | 0.913 |
| Simpson dominance | 0.341 | 0.802 |
| Richness | 0.017 | 0.953 |

**No metric shows any evidence of a lure-dependent insect~fungal relationship** (all F<0.35, all p_perm>=0.80). fig10 was built and then deliberately collapsed from a lure-faceted grid to a single pooled panel per metric on the strength of this result -- splitting the earlier draft by lure was cutting n roughly into thirds (~23/lure) without visualizing a real effect.

## 4. Within-community site/lure/date effects on each metric

Same trap-level permutation logic used throughout this project (date permuted within trap; lure permuted among traps within site; site permuted freely among all traps), applied here to a univariate diversity metric via `lm()` instead of a per-taxon abundance or distance matrix. Date is tested with linear + quadratic + cubic terms in one model (§5 covers the shape read); only the significant (p_perm<0.10) terms are shown below -- see `alpha_diversity_covariate_tests.csv` for the full table (30 rows: 2 communities x 3 metrics x 5 terms).

| Community | Metric | Term | Stat | p_perm |
|---|---|---|---|---|
| insect | Shannon diversity | lure | F=8.07 | **0.013** |
| insect | Shannon diversity | site | F=4.49 | 0.071 |
| insect | Simpson dominance | site | F=3.27 | 0.096 |
| insect | Richness | date (quadratic) | t=-5.52 | **0.001** |
| insect | Richness | date (cubic) | t=3.25 | **0.001** |
| insect | Richness | lure | F=14.66 | **0.003** |
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
- **Lure is never a significant term on the fungal side** (p_perm 0.654-0.962 across all 3 fungal metrics, not shown above) but IS significant for insect Shannon diversity (p=0.013) and is the single largest effect size anywhere in this table for insect richness (F=14.66, p=0.003). Consistent with lure being an insect-trap attractant with a real effect on how diverse a trap's catch is, but no detectable knock-on effect on fungal diversity -- and consistent with §3's null interaction result.
- **Date dominates the fungal side completely** -- all 3 fungal metrics are significant on all 3 polynomial date terms (9/9 significant date tests). The insect side is much patchier: only richness shows any significant date effect (and only on the quadratic/cubic terms, not linear) -- insect Shannon diversity and Simpson dominance show NO significant date effect at all (see §5).
- **Site is a real but always secondary effect** -- present at borderline-to-modest significance for 4 of the 6 community x metric combinations, never the largest term for any of them.

## 5. Date shape classification (linear vs. quadratic vs. cubic)

Extends this project's established linear+quadratic date-test convention (`fungal_community_seasonality.r`, `insect_fungal_coca_analysis.general_workflow.r`) one order further to cubic -- added specifically because fig9's trend lines showed more structure (e.g. an S-curve in insect Shannon diversity under the Ethanol lure) than a single hump/dip could capture. Shape classification cascades: a significant cubic term takes priority over quadratic, which takes priority over linear.

| Community | Metric | Shape | Linear t (p) | Quadratic t (p) | Cubic t (p) |
|---|---|---|---|---|---|
| insect | Shannon diversity | no significant date pattern | 1.72 (0.103) | -1.33 (0.192) | 0.53 (0.623) |
| insect | Simpson dominance | no significant date pattern | -1.31 (0.205) | 0.49 (0.659) | 0.24 (0.823) |
| insect | Richness | **complex (cubic significant)** | -0.95 (0.346) | -5.52 (0.001) | 3.25 (0.001) |
| fungal | Shannon diversity | **complex (cubic significant)** | -4.40 (0.001) | 3.97 (0.002) | 2.87 (0.010) |
| fungal | Simpson dominance | **complex (cubic significant)** | 3.56 (0.001) | -3.09 (0.004) | -2.71 (0.006) |
| fungal | Richness | **complex (cubic significant)** | -4.79 (0.001) | 3.04 (0.006) | 3.65 (0.001) |

**4 of the 6 community x metric date relationships are cubic-significant** -- only insect Shannon diversity and insect Simpson dominance show no significant date structure at all, at any polynomial order. All three fungal metrics are cubic, meaning the simple "fungal diversity declines through the season" linear read (itself also significant, and in the expected direction: Shannon/richness falling, Simpson dominance rising) is an incomplete description -- there's a genuine dip-then-partial-recovery shape, matching fig9's visible mid-June trough for fungal Shannon diversity and richness across all three lures. Insect richness shows the same kind of higher-order structure (t=3.25 on the cubic term, tied for the single lowest p_perm anywhere in this analysis) despite insect Shannon/Simpson showing nothing -- alpha-diversity metrics computed from the SAME table don't always agree on whether there's a date effect at all, let alone its shape.

## 6. Presentation figures

- **`fig9_alpha_diversity_by_lure.R`** -- insect (panel a) and fungal (panel b) alpha diversity vs. collection date, faceted metric (rows) x lure (columns), points colored by site, cubic OLS trend line per panel (`y ~ poly(date, 3)`) matching the permutation-tested cubic date model in §5. Output: `figures/presentation_items/fig9_alpha_diversity_by_lure.{png,pdf}`.
- **`fig10_insect_fungal_alpha_diversity_correlation.R`** -- per-sample insect vs. fungal alpha diversity, one panel per metric (NOT split by lure, per §3's null interaction result), colored by collection date, pooled Spearman rho/p from §2 annotated in each panel. Output: `figures/presentation_items/fig10_insect_fungal_alpha_diversity_correlation.{png,pdf}`.

## 7. Asymptotic (Chao-extrapolated) diversity: how much more is out there?

`insect_fungal_asymptotic_richness_iNEXT.full_insect_table.r` estimates, per sample, the richness/Shannon(Hill q=1)/Simpson(Hill q=2) diversity expected with infinite additional sampling effort, using `iNEXT`'s Chao-type estimators (`ChaoRichness()`/`ChaoShannon()`/`ChaoSimpson()`). n=69 for both communities here (unlike §1-6's fungal n=68 -- the fungal side is unrarefied in this variant, so no rarefaction-depth sample drop). "Simpson" is Hill number q=2 DIVERSITY (higher = more even), the opposite direction from §1-6's Simpson DOMINANCE -- don't compare the two scripts' Simpson columns directly. Richness values match §1's exactly on the insect side (same raw table, "asymptotic Observed" = `specnumber()`); they do NOT match on the fungal side (unrarefied here vs. rarefied-and-averaged in §1-6).

| Metric | Insect observed | Insect asymptotic | Fungal observed (unrarefied) | Fungal asymptotic |
|---|---|---|---|---|
| Richness | 4 - 68 (mean 32.7) | 8.5 - 176.6 (mean 70.8) | 364 - 4,794 (mean 1,723) | 490.7 - 5,546.7 (mean 2,053) |
| Shannon (Hill q=1) | 2.38 - 29.32 (mean 15.2) | 2.79 - 62.11 (mean 24.0) | 8.49 - 345.35 (mean 82.3) | 8.60 - 348.06 (mean 83.9) |
| Simpson (Hill q=2) | 1.49 - 21.42 (mean 8.96) | 1.49 - 30.07 (mean 10.6) | 2.15 - 113.75 (mean 20.4) | 2.15 - 113.82 (mean 20.4) |

Sample coverage (`iNEXT::DataInfo()`'s SC, the fraction of the true community estimated to be captured) tells the same story from a different angle: **insect coverage is 0.18 - 0.98 (mean 0.79)** -- trap catches (median ~93 individuals per sample) are genuinely incomplete samples of the local insect community -- while **fungal coverage is 0.97 - 1.00 (mean 0.99)** -- amplicon sequencing depth (thousands to hundreds of thousands of reads per sample) is already close to exhaustive.

That asymmetry drives everything below: **median asymptotic/observed ratio is ~1.97x for insect richness but only ~1.20x for fungal richness** -- extrapolation roughly doubles the insect richness estimate but barely nudges the fungal one. Shannon (q=1, weighted toward common taxa) is far less sensitive to the correction than richness for either community (median ratio 1.40x insect, 1.01x fungal); Simpson (q=2, weighted even more toward dominant taxa) barely moves at all (1.08x insect, 1.00x fungal) -- expected, since higher-order Hill numbers are increasingly dominated by the already-well-sampled common taxa and increasingly insensitive to how many rare/undetected ones exist.

## 8. Does extrapolating change the insect~fungal correlation conclusion? No.

Same Spearman test as §2, computed on both observed and asymptotic values within the same script (a clean, same-metric-definition comparison unlike some of the §10 comparisons below):

| Metric | Observed rho (p) | Asymptotic rho (p) | Conclusion |
|---|---|---|---|
| Richness | 0.076 (0.534, n.s.) | 0.008 (0.950, n.s.) | Unchanged -- no relationship either way |
| Shannon (Hill q=1) | -0.253 (0.036, **sig**) | -0.251 (0.038, **sig**) | Unchanged -- same modest negative relationship |
| Simpson (Hill q=2) | -0.075 (0.539, n.s.) | -0.075 (0.539, n.s.) | Unchanged -- no relationship either way |

**Extrapolating to the asymptote does not change any of the three §2-style conclusions** -- correcting for undetected taxa (which, per §7, matters a lot more for insect richness than for anything else in this table) doesn't meaningfully move the correlation with fungal diversity. The weak-Shannon/null-richness-and-Simpson pattern from §2 is robust to whether "diversity" means observed counts or extrapolated estimates.

## 9. Does lure modify the insect~fungal relationship, under extrapolation? Mostly still no -- but Shannon crosses this project's significance line

Same interaction test as §3 (`fungal ~ site+date+insect*lure`, 999-perm F-test on the interaction term), applied to asymptotic values:

| Metric | F | p_perm | vs. §3 (observed: F=0.02-0.34, p=0.80-0.95) |
|---|---|---|---|
| Richness | 0.081 | 0.888 | Consistent -- still clearly n.s. |
| Shannon (Hill q=1) | 0.657 | **0.075** | **Crosses this project's p<0.10 convention** -- was p=0.913 (n.s.) in §3 |
| Simpson (Hill q=2) | 0.321 | 0.697 | Consistent -- still clearly n.s. |

Worth flagging plainly rather than glossing over: by this project's own uncorrected-p<0.10 convention for this kind of small, fixed covariate/interaction test set (see `project_organization.md`'s "Significance convention"), the asymptotic Shannon interaction is technically significant, unlike every other lure-interaction test run anywhere in this document (both observed §3 and the other two asymptotic metrics here). That said: it's a single test just under the 0.10 line (p=0.075), effect size is modest (F=0.657, smallest F of the three asymptotic tests is 0.081 but this is nowhere near the largest covariate F's seen in §4/§10), and it comes from 3 uncorrected tests, not a taxon screen -- treat as a "worth another look if the dataset grows," not a settled finding. fig12 (§12) still shows the pooled (non-lure-split) relationship, matching richness/Simpson's clean nulls and consistent with §3's overall conclusion; a future update could revisit whether Shannon specifically deserves a lure-faceted panel if this signal strengthens.

## 10. Within-community site/lure/date effects, asymptotic values

Same trap-level permutation logic as §4, run here on the asymptotic estimates only (observed-value covariate testing under these exact metric definitions -- Hill-number Shannon/Simpson, unrarefied fungal counts -- isn't separately re-derived; §4 already covers observed-value significance under its own, different, metric definitions). Only significant/borderline (p_perm<0.10) terms shown -- see `asymptotic_diversity_covariate_tests.csv` for the full 30-row table (2 communities x 3 metrics x 5 terms).

| Community | Metric | Term | Stat | p_perm |
|---|---|---|---|---|
| insect | Richness | date (quadratic) | t=-3.32 | **0.003** |
| insect | Shannon | site | F=8.04 | **0.022** |
| insect | Shannon | date (linear) | t=1.95 | 0.079 |
| insect | Shannon | lure | F=3.60 | 0.070 |
| insect | Simpson | date (quadratic) | t=1.97 | 0.056 |
| fungal | Richness | site | F=1.16 | **0.015** |
| fungal | Shannon | date (linear) | t=-4.28 | **0.001** |
| fungal | Shannon | date (quadratic) | t=3.23 | **0.001** |
| fungal | Shannon | date (cubic) | t=2.39 | **0.025** |
| fungal | Simpson | date (linear) | t=-3.00 | **0.003** |
| fungal | Simpson | date (quadratic) | t=2.78 | **0.008** |

**One directly comparable, genuinely interesting divergence**: insect richness is the one metric where "observed" here is IDENTICAL to §4's raw richness (same table, no transform), so §4's insect-richness row and this table's insect-richness row are a clean before/after extrapolation comparison. §4 (observed) found a strong lure effect (F=14.66, **p=0.003**) and a significant cubic date term (t=3.25, **p=0.001**) for insect richness. Here (asymptotic), lure is no longer significant (F=1.29, p=0.447) and the cubic term is gone (t=1.47, p=0.128) -- only the quadratic (hump) term survives (p=0.003). **Extrapolating to account for undetected insect taxa attenuates the lure effect on richness and simplifies its date shape from cubic to a plain mid-season hump.** A plausible reading: part of the observed lure-richness effect may reflect unequal SAMPLING COMPLETENESS across lures (some lures catch more individuals per trap, pushing observed richness up mechanically) rather than a true difference in the underlying taxon pool each lure attracts -- exactly the kind of confound asymptotic estimation exists to correct for. (Shannon/Simpson can't be compared this cleanly to §4, since their Hill-number definitions here differ from §4's raw entropy/dominance metrics, and the fungal side is unrarefied here vs. rarefied-and-averaged in §4 -- treat any resemblance or difference between those rows and §4's as coincidental, not evidence either way.)

## 11. Date shape classification, asymptotic values

Same linear/quadratic/cubic cascading classification as §5 (cubic takes priority over quadratic, which takes priority over linear), applied to asymptotic values:

| Community | Metric | Shape | Linear t (p) | Quadratic t (p) | Cubic t (p) |
|---|---|---|---|---|---|
| insect | Richness | **hump (peaks mid-season)** | 0.46 (0.621) | -3.32 (0.003) | 1.47 (0.128) |
| insect | Shannon | no significant date pattern | 1.95 (0.079)* | 0.97 (0.355) | 0.13 (0.904) |
| insect | Simpson | **dip (troughs mid-season)** | 1.28 (0.220) | 1.97 (0.056) | -0.05 (0.967) |
| fungal | Richness | no significant date pattern | -0.86 (0.423) | -1.06 (0.307) | 1.06 (0.330) |
| fungal | Shannon | **complex (cubic significant)** | -4.28 (0.001) | 3.23 (0.001) | 2.39 (0.025) |
| fungal | Simpson | **dip (troughs mid-season)** | -3.00 (0.003) | 2.78 (0.008) | 1.29 (0.223) |

*insect Shannon's linear term (p=0.079) doesn't reach significance under the shape classification's cascade (no quadratic/cubic term is significant either), so it's still called "no significant date pattern" here -- shown for completeness, not a contradiction.

Two changes worth noting relative to §5's observed-value classification: **insect richness simplifies from "complex (cubic significant)" to a plain "hump"** (the same cubic-term loss already flagged in §10), and **fungal richness loses its date signal entirely** (§5: "complex (cubic significant)," all 3 terms p<=0.03; here: no term reaches p<0.10) -- on the unrarefied, extrapolated fungal table, richness's date pattern doesn't hold up, even though fungal Shannon and Simpson's date patterns (still complex/dip respectively) do. Fungal Shannon remains the only "complex" (genuinely non-monotonic, all 3 terms significant) result in this table, same as §5.

## 12. Presentation figures (asymptotic)

- **`fig11_asymptotic_diversity_by_lure.R`** -- fig9 reproduced one-for-one on asymptotic values: insect (panel a) and fungal (panel b) asymptotic diversity vs. collection date, faceted metric (rows) x lure (columns), cubic OLS trend lines. Row-strip metric labels shortened ("Shannon div. (q=1)" etc.) vs. fig9's wording, which clipped in the rotated strip. Output: `figures/presentation_items/fig11_asymptotic_diversity_by_lure.{png,pdf}`.
- **`fig12_insect_fungal_asymptotic_diversity_correlation.R`** -- fig10 reproduced one-for-one on asymptotic values: per-sample insect vs. fungal asymptotic diversity, one panel per metric, pooled (not lure-split, per §9) Spearman rho/p from §8 annotated. Points colored by collection date AND shaped by lure (filled circle=Ethanol, triangle=Alpha-pinene_EtOH, square=Ips) -- descriptive only, doesn't affect the pooled trend line/correlation; this shape mapping was also back-ported to fig10 for consistency between the two figures. Output: `figures/presentation_items/fig12_insect_fungal_asymptotic_diversity_correlation.{png,pdf}`.

## Bottom line

1. **Fungal diversity shows strong, comprehensive date effects** (observed values, §1-6) -- all 3 metrics (Shannon, Simpson dominance, richness) are significant on all 3 polynomial date terms (linear, quadratic, AND cubic), and every one of them is a genuinely non-monotonic (cubic) shape, not a simple increase or decrease.
2. **Insect diversity's date signal is metric-specific and much weaker** -- only richness shows a significant (and non-monotonic, cubic) date effect; Shannon diversity and Simpson dominance show no significant date effect at all under this full, all-families table. This contrasts with the strong insect~date signal already established at the COMMUNITY-COMPOSITION level under both Lineage A and B (PERMANOVA date R² ~13-15%) -- alpha diversity and community composition are different questions and don't have to move together.
3. **Lure drives insect diversity but not fungal diversity** -- significant for insect Shannon (p=0.013) and the largest effect size in the whole covariate table for insect richness (F=14.66, p=0.003), but never significant for any fungal metric (p=0.654-0.962) -- and §3 confirms lure doesn't even modify the insect~fungal relationship itself. Consistent with lure acting as an insect-trap attractant with no direct fungal-diversity consequence.
4. **Site is a real but always secondary effect** on both communities -- present at borderline-to-modest significance for roughly half the community x metric combinations, never the dominant term for any of them.
5. **Insect and fungal alpha diversity are only weakly linked, and only for one metric**: Shannon diversity shows a modest, significant negative correlation (rho=-0.26, p=0.03 -- more insect diversity co-occurs with LESS fungal diversity); Simpson dominance and richness show no significant correlation. This is a much weaker relationship than the strong whole-community Procrustes correlation for composition (r=0.60, p=0.001 under Lineage B) -- these two communities' raw diversity levels are far less tightly coupled than their compositions are.
6. **Lure does not modify the insect~fungal relationship for any metric under observed values** (all interaction p_perm>=0.80) -- fig10 shows the pooled, non-lure-faceted relationship rather than an artificially data-starved 3-way split.
7. **Extrapolating to asymptotic (Chao-estimated) diversity leaves the headline insect~fungal correlation conclusion unchanged** (§8) -- same weak-Shannon/null-richness-and-Simpson pattern whether using observed counts or diversity corrected for undetected taxa.
8. **How much "extra" diversity extrapolation adds is wildly uneven, and it's mostly about insect richness** (§7) -- insect sample coverage (mean 0.79) is far lower than fungal (mean 0.99), so insect richness roughly doubles on extrapolation (median ratio ~2.0x) while fungal richness moves only ~20%; higher-order Hill numbers (Shannon, Simpson) barely move for either community. Insect trap catches are genuinely incomplete community samples in a way fungal amplicon sequencing isn't.
9. **The one place extrapolation changes a covariate conclusion**: insect richness's strong lure effect (F=14.66, p=0.003) and cubic date shape under observed counts (§4) attenuate to no significant lure effect (p=0.447) and a plain quadratic hump (§10-11) once corrected for undetected taxa -- suggesting the observed lure effect on richness may partly reflect unequal sampling completeness across lures rather than a genuinely different taxon pool per lure.
10. **One new borderline signal, flagged but not over-interpreted**: the lure x insect-diversity interaction on the fungal side, null everywhere else in this document (§3, and 2 of 3 asymptotic metrics in §9), crosses this project's p<0.10 convention for asymptotic Shannon specifically (p=0.075) -- a single test just under the line, not yet reflected in fig12's still-pooled design, worth re-checking as data accumulates.

Underlying files: §1-6 -- `data/compare_insects_fungi_alpha_diversity_full_insect_table/insect_alpha_diversity.csv`, `fungal_alpha_diversity.csv`, `alpha_diversity_cross_community_correlation.csv`, `alpha_diversity_insect_fungal_interaction_by_lure.csv`, `alpha_diversity_covariate_tests.csv`, `alpha_diversity_date_shape.csv`; figures in `figures/compare_insects_fungi_alpha_diversity_full_insect_table/` (descriptive) and `figures/presentation_items/fig9_alpha_diversity_by_lure.{png,pdf}`, `fig10_insect_fungal_alpha_diversity_correlation.{png,pdf}`. §7-12 -- `data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table/insect_asymptotic_diversity.csv`, `fungal_asymptotic_diversity.csv`, `asymptotic_diversity_cross_community_correlation.csv`, `asymptotic_diversity_insect_fungal_interaction_by_lure.csv`, `asymptotic_diversity_covariate_tests.csv`, `asymptotic_diversity_date_shape.csv`; figures in `figures/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table/` (descriptive) and `figures/presentation_items/fig11_asymptotic_diversity_by_lure.{png,pdf}`, `fig12_insect_fungal_asymptotic_diversity_correlation.{png,pdf}`.
