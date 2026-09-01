# Interpretation update (2026-08-28): top-line numbers for the prevalence-filtered insect table lineage, vs. the family-filtered lineage

This file is the prevalence-filtered-lineage (Lineage B, per `project_organization.md`) counterpart to `interpretation_update_08172026.md`, which covers the family-filtered (Lineage A) results. Same purpose: pull current, on-disk numbers into one place rather than requiring anyone to reconstruct them from the chat-style narrative in `iterative_analysis_updates.md`. Read `project_organization.md` first if you haven't -- it defines the two (really three, counting the alpha-diversity-only full table) insect-table constructions this file assumes you already know about.

**Table construction reminder**: Lineage B = insect trap-catch table, ALL families (not just Curculionidae+Latridiidae), individual-taxon >=5-sample prevalence filter (same threshold used for fungal ASVs). 153 taxa, same 69-sample insect-fungal matched dataset as Lineage A. Everything below is this construction unless stated otherwise.

**This lineage is not fully built out yet.** Ported so far: insect NMDS/PERMANOVA, the full CoCA workflow (PCoA1-3 vs. fungal taxa), PCoA1-3 vs. date/site/lure (plus quadratic and factor(date) robustness checks), per-axis taxon drivers, and the direct insect~date screen. NOT yet ported: whole-community Procrustes, the all-against-all pairwise taxon screen, and all 8 presentation figures. See `project_organization.md`'s "Known gaps" section for the current checklist.

---

## 1. Insect community structure: NMDS + PERMANOVA

| | Lineage A (family-filtered, 52 taxa) | Lineage B (prevalence-filtered, 153 taxa) |
|---|---|---|
| NMDS distance/transform | raw Bray-Curtis, no autotransform | same |
| NMDS stress | ~0.14 | 0.162 (repeatable -- best solution hit in 5-18/20 tries across runs) |
| Site R² | 9.7-9.9% (p=0.001) | 9.5% (p=0.001) |
| Lure R² | 14.2-15.2% (p≤0.002) | 11.4% (p=0.001) |
| Date R² | 14.2-14.7% (p=0.001) | 12.7% (p=0.001) |
| Residual | 62.3-62.6% | 66.6% |

The three-way site/lure/date structure is essentially the same shape under both tables -- all three terms significant, lure and date the two largest, roughly comparable in size to each other, site a real but smaller third term. Lineage B's effect sizes are all a bit smaller and residual variance a bit larger, consistent with a broader taxon set diluting the signal with additional taxa that aren't as strongly lure- or season-responsive as the two originally-selected dominant families. Family composition of the 153 retained taxa: Curculionidae is still the largest single group (37 taxa) but no longer dominant -- Elateridae (15), Staphylinidae (12), Cerambycidae (6), Tenebrionidae (6), and dozens of smaller families are now represented, vs. only Curculionidae+Latridiidae under Lineage A.

## 2. Direct taxon~date screens (linear + quadratic, hump/dip shape)

Both lineages now use the linear+quadratic date test added 2026-08-28 (see `iterative_analysis_updates.md`). This is a within-lineage comparison, not something the prevalence filter itself changes the *method* of -- but the broader taxon set surfaces a substantially higher rate of both linear and non-monotonic seasonality:

| | Lineage A (52 insect taxa) | Lineage B (153 insect taxa) |
|---|---|---|
| Linear-significant (q<0.10) | 19 (36.5%) | 69 (45.1%) |
| Quadratic-significant (hump/dip, q<0.10) | 9 (17.3%) | 38 (24.8%) |
| -- hump (peaks mid-season) | 5 | 23 |
| -- dip (troughs mid-season) | 4 | 15 |
| -- linear increase (late-season) | 3 | 32 |
| -- linear decrease (early-season) | 11 | 22 |
| -- no significant pattern | 29 | 61 |

Lineage B finds proportionally MORE seasonally-structured insect taxa on every measure, not just more taxa in absolute terms -- consistent with the alpha-diversity finding (see `iterative_analysis_updates.md`'s full-insect-table entries) that the family restriction was suppressing a real seasonal signal carried by taxa outside Curculionidae+Latridiidae (click beetles, longhorn beetles, rove beetles, etc.).

The fungal side of this screen (`fungal_community_seasonality.r`'s unbiased 6,342-taxon test) does not depend on the insect table at all, so it's identical under both lineages: **2,186/6,342 (34.5%) linear-significant, 707/6,342 (11.1%) quadratic-significant** (377 hump, 330 dip). Not re-run or re-reported per lineage -- see `project_organization.md`.

## 3. insect_PCoA1-3 vs. date/site/lure (Lineage B specific -- no direct Lineage A equivalent for lure/site until `insect_pcoa_date_lure_association.r` is ported)

| Axis | %var | Date (r, partial R², p) | Lure (partial R², p) | Site (partial R², p, parametric) | Quadratic date (F, p) |
|---|---|---|---|---|---|
| PCoA1 | 19.3% | r=0.84, R²=0.80, p=0.001 | R²=0.36, p=0.029 | R²=0.10, p=0.078 | F=26.6, p=0.001 |
| PCoA2 | 9.5% | r=0.06, R²=0.008, p=0.465 (n.s.) | R²=0.27, p=0.005 | R²=0.16, p=0.012 | F=92.7, p=0.001 |
| PCoA3 | 7.3% | r=-0.13, R²=0.056, p=0.062 (n.s. after BH) | R²=0.69, p=0.003 | R²=0.028, p=0.617 | F=16.1, p=0.001 |

PCoA1 is the date axis (as under Lineage A), but now also carries a real secondary lure signal never tested for under Lineage A. PCoA3 is clearly the lure axis (comparable in kind, if not exact magnitude, to Lineage A's PCoA3/lure result, R²=0.81 in the original 10-axis screen). **PCoA2 is new** -- it did not show up as a distinct, characterizable axis under Lineage A's smaller taxon set. It has ~zero linear date association but a real lure association (R²=0.27) AND the single strongest quadratic date signal of the three axes (F=92.7) -- i.e. a real mid-season peak, not detectable by a linear test. All three axes have significant quadratic date components; this project's date-adjustment tests going forward should account for that (see §4).

## 4. Which insect taxa drive each axis (new for Lineage B -- not previously done as a formal script for either lineage)

- **PCoA1** drivers split cleanly by guild, same story as the original informal PCoA1-driver check under Lineage A: Scolytinae (*Xyleborinus attenuatus* r=-0.85, *Xylosandrus germanus* r=-0.58, *Cyclorhipidion pelliculosum* r=-0.53, *Pissodes strobi* r=-0.46) vs. fungivorous Latridiidae/Corticariinae (*Melanophthalma* sp. r=+0.83). Newly visible under Lineage B: several non-Curculionidae/Latridiidae taxa also load on PCoA1 (*Melanotus hyslopi* Elateridae r=+0.57, *Isorhipis obliqua* Eucnemidae r=+0.53, *Neoclytus acuminatus* Cerambycidae r=+0.50) -- taxa the family filter excluded entirely.
- **PCoA2** drivers are almost entirely WITHIN-Scolytinae, differentiating which bark/ambrosia beetle SPECIES dominate a sample: *Xylosandrus germanus* (r=0.68), *Ips grandicollis* (r=0.48), *Pissodes affinis*, *Anisandrus sayi*, *Dryocoetes affaber*, *Xylosandrus crassiusculus*, *Dendroctonus valens* (all positive) vs. *Xyleborinus attenuatus* again (r=-0.46, the same species anchoring PCoA1's negative end). Reads as a species-composition-within-guild axis, not a guild-vs-guild split like PCoA1.
- **PCoA3** drivers overlap partially with PCoA2's cast but with different, sometimes opposite signs (*Anisandrus sayi* r=+0.69, *Ips grandicollis* r=-0.62, *Hylobius pales* r=-0.58) -- consistent with PCoA2 and PCoA3 being two different (timing vs. lure-preference) ways of splitting the same set of bark-beetle species.
- **Reclassifying two headline taxa by shape, not just direction**: *Xyleborinus attenuatus* (PCoA1's strongest driver, described everywhere in this project as "early-season") is actually a **dip** shape (t_quad=9.3, q=0.014) -- high early, low mid-season, rising again late -- not a simple decline. *Melanophthalma* sp. (PCoA1's other headline taxon, "later-season") is a clean **hump** (t_quad=-3.4, q=0.022), peaking mid-season rather than simply increasing. Both reclassifications hold under Lineage A's smaller table too (see §2 above and `iterative_analysis_updates.md`).

## 5. CoCA insect-fungal association hits, and the Cytospora/Tympanis correction (Lineage B specific finding, most important item in this file)

Headline hit counts (top-200 CoCA candidates, with lure, same permutation scheme as Lineage A):

| Axis | no date | linear date | linear+quadratic date | factor(date) |
|---|---|---|---|---|
| PCoA1 | 59/200 | 9/200 | 19/200 | 18/200 |
| PCoA2 | 6/200 | 7/200 | **0/200** | **9/200** |
| PCoA3 | 0/200 (all four lure/date variants) | | | |

PCoA1's pattern is stable once any non-linear flexibility is allowed (9→19→18); PCoA3 remains completely null under Lineage B, same as Lineage A. **PCoA2 is unstable** -- 0 under quadratic, 9 under the maximally strict factor(date) control -- flagged as small-sample-sensitive, not a confirmed finding either way (see `iterative_analysis_updates.md` for the full discussion of why).

**Cytospora prunicola (ASV_1905), tracked across every date-adjustment strategy on PCoA1:**

| Model | t_stat | q_value |
|---|---|---|
| no date | -7.65 | 0.017 (sig) |
| linear date | -3.35 | 0.109 (fails) |
| linear+quadratic date | -2.41 | 0.160 (fails) |
| factor(date) | -2.01 | 0.257 (fails) |

**This is a correction to how Cytospora has been characterized throughout this project.** It fails every date-adjusted test under Lineage B, monotonically worse as the control gets stricter. It should NOT be cited as a date-robust insect-community-associated taxon going forward -- its standing rests entirely on (a) the unbiased fungal~date screen (§2, shared across lineages) and (b) its raw, no-date-adjustment PCoA1/PCoA2 loading, both of which are fundamentally seasonal observations, not evidence of a direct insect-mediated effect.

**Tympanis sp. (ASV_13242), same tracking, is the opposite pattern:**

| Model | t_stat | q_value |
|---|---|---|
| no date | 0.89 | 0.350 (n.s.) |
| linear date | 2.38 | 0.040 (sig) |
| linear+quadratic date | 2.42 | 0.088 (sig) |
| factor(date) | 2.55 | 0.094 (sig) |

Tympanis is null with no date control, then significant under every date-adjusted version -- the single most date-robust PCoA1 hit found under Lineage B so far. Combined with a second, different *Tympanis* ASV (ASV_5640) topping the PCoA2 factor(date)-adjusted list, and 4 different Valsaceae ASVs (2x *Cytospora*, 1x *Valsa*, 1x unresolved) also surviving that same PCoA2 test (though not ASV_1905 itself) -- **the recommended framing going forward is genus/family-level** ("Valsaceae/*Cytospora*-as-a-group" and "*Tympanis*-as-a-group" keep reappearing via different specific ASVs across axes and date-control strategies) **rather than any single ASV**, and specifically not *Cytospora prunicola* by name as a date-robust finding.

## Bottom line

1. **Community-level structure (PERMANOVA, NMDS) is qualitatively unchanged** between lineages -- same three significant terms, same rough ranking (lure/date > site), Lineage B just has smaller effect sizes and more residual variance from the broader, noisier taxon set.
2. **The broader taxon set reveals more seasonality, not less** -- both a higher linear hit rate (45% vs. 37%) and a higher hump/dip rate (25% vs. 17%) among insect taxa, consistent with the family filter having suppressed real signal in previously-excluded families.
3. **A new axis (PCoA2) emerged** that has no clean Lineage-A analog -- a within-Scolytinae, lure-and-mid-season-timing axis, distinct from PCoA1's guild-level date split and PCoA3's lure split. Its fungal-association hit count is unstable across date-control strategies and should be reported as a range (0-9/200), not a single number.
4. **The single most important correction from this line of work**: *Cytospora prunicola* (ASV_1905) -- previously the project's headline "most robust cross-cutting taxon" -- does not survive ANY date-adjusted test on any insect PCoA axis under Lineage B. *Tympanis* sp. is the new best candidate for a genuinely date-robust insect-community-linked fungal signal, and the Valsaceae/*Cytospora* and *Tympanis* genus-level groupings (not specific ASVs) are the more defensible framing overall.

Underlying files referenced above: `data/insect_exploratory_prevalence_filtered_insect/insect_prevalence_filter.permanova.csv`, `insect_prevalence_filter.family_composition.csv`; `data/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_taxa_date_association.csv`, `insect_pcoa_axes.date_site_lure_association.csv`, `insect_pcoa_axes.quadratic_date_test.csv`, `insect_pcoa_axes.taxon_drivers.csv`, `date_adjustment_comparison.all_strategies.csv`, `fungal_insect_association_results.PCoA{1,2,3}_*.csv`; `data/2024_fungi/fungal_taxa_date_association.all_taxa.csv` (shared, lineage-independent); `data/2024_insect_data/insect_taxa_date_association.csv` (Lineage A comparison).
