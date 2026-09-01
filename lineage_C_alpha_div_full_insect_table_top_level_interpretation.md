# Interpretation update (2026-09-01): top-line numbers for the full-insect-table alpha-diversity lineage (Lineage C)

This file is the Lineage C (full insect table, per `project_organization.md`) counterpart to `lineage_A_family_filtered_insect_top_level_interepretation.md` and `lineage_B_prevalence_filtered_insect_top_level_interpretation.md`. Same purpose: pull current, on-disk numbers into one place rather than requiring anyone to reconstruct them from the chat-style narrative in `iterative_analysis_updates.md`. Read `project_organization.md` first if you haven't -- it defines the three insect-table constructions this file assumes you already know about.

**Table construction reminder**: Lineage C = insect trap-catch table, ALL families (not just Curculionidae+Latridiidae), **singleton-taxon filter only** -- deliberately NOT the >=5-sample prevalence filter Lineage B uses. 277 taxa across 98 families, 69 samples before matching to fungal data (68 after the fungal rarefaction step drops one low-depth sample). This is the correct all-families construction *specifically for alpha-diversity metrics* -- prevalence-filtering before computing richness/Shannon/Simpson would discard real taxa and artificially deflate diversity (see project_organization.md's "When to prevalence-filter vs. not"). **Lineage C is alpha-diversity-only by design** -- there is no Lineage C ordination/CoCA/Procrustes/pairwise-taxon-screen analog the way Lineage A and B both have; don't look for one.

Everything below comes from `compare_insect_fungi/insect_fungal_alpha_diversity.full_insect_table.r` and the two presentation figures built on top of it (`presentation_items/fig9_alpha_diversity_by_lure.R`, `fig10_insect_fungal_alpha_diversity_correlation.R`).

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

## Bottom line

1. **Fungal diversity shows strong, comprehensive date effects** -- all 3 metrics (Shannon, Simpson dominance, richness) are significant on all 3 polynomial date terms (linear, quadratic, AND cubic), and every one of them is a genuinely non-monotonic (cubic) shape, not a simple increase or decrease.
2. **Insect diversity's date signal is metric-specific and much weaker** -- only richness shows a significant (and non-monotonic, cubic) date effect; Shannon diversity and Simpson dominance show no significant date effect at all under this full, all-families table. This contrasts with the strong insect~date signal already established at the COMMUNITY-COMPOSITION level under both Lineage A and B (PERMANOVA date R² ~13-15%) -- alpha diversity and community composition are different questions and don't have to move together.
3. **Lure drives insect diversity but not fungal diversity** -- significant for insect Shannon (p=0.013) and the largest effect size in the whole covariate table for insect richness (F=14.66, p=0.003), but never significant for any fungal metric (p=0.654-0.962) -- and §3 confirms lure doesn't even modify the insect~fungal relationship itself. Consistent with lure acting as an insect-trap attractant with no direct fungal-diversity consequence.
4. **Site is a real but always secondary effect** on both communities -- present at borderline-to-modest significance for roughly half the community x metric combinations, never the dominant term for any of them.
5. **Insect and fungal alpha diversity are only weakly linked, and only for one metric**: Shannon diversity shows a modest, significant negative correlation (rho=-0.26, p=0.03 -- more insect diversity co-occurs with LESS fungal diversity); Simpson dominance and richness show no significant correlation. This is a much weaker relationship than the strong whole-community Procrustes correlation for composition (r=0.60, p=0.001 under Lineage B) -- these two communities' raw diversity levels are far less tightly coupled than their compositions are.
6. **Lure does not modify the insect~fungal relationship for any metric** (all interaction p_perm>=0.80) -- fig10 shows the pooled, non-lure-faceted relationship rather than an artificially data-starved 3-way split.

Underlying files referenced above: `data/compare_insects_fungi_alpha_diversity_full_insect_table/insect_alpha_diversity.csv`, `fungal_alpha_diversity.csv`, `alpha_diversity_cross_community_correlation.csv`, `alpha_diversity_insect_fungal_interaction_by_lure.csv`, `alpha_diversity_covariate_tests.csv`, `alpha_diversity_date_shape.csv`; figures in `figures/compare_insects_fungi_alpha_diversity_full_insect_table/` (descriptive) and `figures/presentation_items/fig9_alpha_diversity_by_lure.{png,pdf}`, `fig10_insect_fungal_alpha_diversity_correlation.{png,pdf}`.
