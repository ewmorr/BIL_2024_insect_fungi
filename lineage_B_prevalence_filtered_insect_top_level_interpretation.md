# Interpretation update (2026-08-28): top-line numbers for the prevalence-filtered insect table lineage, vs. the family-filtered lineage

This file is the prevalence-filtered-lineage (Lineage B, per `project_organization.md`) counterpart to `interpretation_update_08172026.md`, which covers the family-filtered (Lineage A) results. Same purpose: pull current, on-disk numbers into one place rather than requiring anyone to reconstruct them from the chat-style narrative in `iterative_analysis_updates.md`. Read `project_organization.md` first if you haven't -- it defines the two (really three, counting the alpha-diversity-only full table) insect-table constructions this file assumes you already know about.

**Table construction reminder**: Lineage B = insect trap-catch table, ALL families (not just Curculionidae+Latridiidae), individual-taxon >=5-sample prevalence filter (same threshold used for fungal ASVs). 153 taxa, same 69-sample insect-fungal matched dataset as Lineage A. Everything below is this construction unless stated otherwise.

**This lineage is now feature-complete with Lineage A.** Ported: insect NMDS/PERMANOVA, the full CoCA workflow (PCoA1-3 vs. fungal taxa), PCoA1-3 vs. date/site/lure (plus quadratic and factor(date) robustness checks), per-axis taxon drivers, the direct insect~date screen, the whole-community Procrustes (r=0.604, p=0.001, 68 matched samples -- see `iterative_analysis_updates.md`'s 2026-08-31 entry), the all-against-all pairwise taxon screen (§6), and all 8 presentation figures (fig1 lineage-invariant; fig2/fig3 2026-08-31; fig4-fig8 2026-09-01, see §7-§8 below). Remaining gaps (deliberately out of scope or minor variants) are in `project_organization.md`'s "Known gaps" section.

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

*Effect-size metric (2026-09-08):* fig3's volcano x-axis and all "top-ranked taxon" lists derived from these screens now use the fully-standardized partial slope β = b·sd(x)/sd(y) instead of the permutation t-statistic (the meaningful effect size for a continuous date/PCoA predictor -- not a two-group Cohen's d). β is an observed-fit quantity, so the counts and q-values above are untouched and `sign(β) == sign(t)`; only the *ranking within* the significant set changes, which shifts a few of fig3's labelled taxa (see the 2026-09-08 entry in `iterative_analysis_updates.md`).

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
- **Reclassifying two headline taxa by shape, not just direction**: *Xyleborinus attenuatus* (PCoA1's strongest driver, described everywhere in this project as "early-season") is actually a **dip** shape (standardized quadratic slope β_quad = +0.51, q=0.014; Lineage A +0.52, q=0.017) -- high early, low mid-season, rising again late -- not a simple decline. *Melanophthalma* sp. (PCoA1's other headline taxon, "later-season") is a clean **hump** (β_quad = -0.28, q=0.022; Lineage A -0.21, q=0.045), peaking mid-season rather than simply increasing. Both reclassifications hold under Lineage A's smaller table too (see §2 above and `iterative_analysis_updates.md`). (Effect sizes quoted as the standardized partial slope β = b·sd(x)/sd(y) as of 2026-09-08, replacing the raw t-statistic throughout the fig3 screens; sign and significance are unchanged.)

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

## 6. All-against-all pairwise fungal x insect taxon screen (2026-09-01 port)

Every prevalence-filtered fungal ASV tested against every prevalence-filtered
insect taxon individually, one model per pair
(`fungal_clr ~ site + lure + date + insect_taxon`), 999 trap-blocked
permutations, vectorized-OLS. Lineage B grid is **153 insect x 6,342 fungal =
970,326 pairs** (vs. Lineage A's 39 x 6,342 = 247,338). Fungal side is
identical between lineages; only the insect predictor set widens.

**The all-families table does not change the conclusion -- it just adds two
taxa that then fail the robustness check.**

| | Lineage A (39 insect taxa) | Lineage B (153 insect taxa) |
|---|---|---|
| Pairs q<0.10 within-insect-taxon | 1,918 | 2,402 |
| Pairs q<0.10 global BH | **0** | **0** (min global q = 0.37) |
| Insect taxa with >=1 significant partner | 3 | 5 |
| Which taxa | *Pityogenes hopkinsi* (1,593), *Xyleborinus attenuatus* (254), *Hylastes opacus* (71) -- all early-season Scolytinae | *P. hopkinsi* (1,592), *Asemum striatum* (266), *X. attenuatus* (213), *Orthoperus scutellaris* (198), *H. opacus* (133) |

- The **three Scolytinae are the same trio** Lineage A flagged, with nearly
  identical partner counts (*P. hopkinsi* 1,592 vs. 1,593 -- its Hellinger
  predictor column barely moves when the community basis widens, since it is
  an abundant taxon).
- The broader table adds **two non-Scolytine taxa**: *Asemum striatum*
  (Cerambycidae, a longhorn beetle) and *Orthoperus scutellaris*
  (Corylophidae, a minute fungus beetle) -- both families excluded by the
  Lineage-A Curculionidae+Latridiidae restriction.
- Signal is **extremely concentrated**: 5 of 153 insect taxa carry all 2,402
  hits; on the fungal side 1,800 of the 2,071 involved ASVs link to exactly
  one insect taxon (max 4). Top fungal genera among hits: *Aureobasidium*
  (117), *Dothiora* (83), *Cytospora* (66), *Taphrina* (42) -- a
  saprotroph/plant-surface-yeast cast, consistent with a shared seasonal
  bloom rather than beetle-specific symbiosis.
- **Not a single pair survives a global (all-970k-tests) BH correction**, same
  as Lineage A. The "within-insect-taxon" hits are best read as "this insect
  taxon's seasonal abundance curve happens to line up with a lot of fungal
  taxa's curves", not as 2,402 specific associations.

**Residualized follow-up** (`.residualized.prevalence_filtered_insect.r` --
each insect predictor stripped of `site + factor(date)`, the strictest
seasonal control the 6 sampling rounds allow; fungal-side model keeps
`site + lure + factor(date) + insect_taxon_resid`):

| Insect taxon | main-grid partners | residualized partners | survive |
|---|---|---|---|
| *Pityogenes hopkinsi* | 1,592 | **646** | 605 overlap (38%) |
| *Asemum striatum* | 266 | 0 | -- |
| *Xyleborinus attenuatus* | 213 | 0 | -- |
| *Orthoperus scutellaris* | 198 | 0 | -- |
| *Hylastes opacus* | 133 | 0 | -- |

Identical outcome to Lineage A (there: only *P. hopkinsi* survived, 717 of
1,593 = 45%). **Only *Pityogenes hopkinsi* retains any fungal partners once
non-linear season is removed from its predictor**; the other four collapse to
zero, i.e. their fungal "associations" were entirely a shared early-season
timing artifact. The two taxa the broader table added (*Asemum*,
*Orthoperus*) are therefore seasonal co-occurrence, not evidence of a
specific insect-fungal link -- the all-families table adds no new robust
pairwise association over Lineage A.

*P. hopkinsi*'s surviving 646 pairs remain uncorrected-globally and are not
pursued further here; per the parent script's note, targeted with/without
partialling checks on that one taxon's candidates would be the follow-up if
this thread is picked up.

Files: `data/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect/`
-- `fungal_insect_pairwise_full_grid.csv`,
`fungal_insect_pairwise_significant_hits.csv`,
`fungal_partners_per_insect_taxon.csv`,
`insect_partners_per_fungal_taxon.csv`, `fungal_family_tally.csv`, and the
`*.date_site_residualized.csv` + `original_vs_residualized_comparison.csv`
counterparts; figures in the matching `figures/` directory.

## 7. Higher-level summaries: Order/Class and trophic mode (fig4/fig5/fig7/fig8, 2026-09-01 port; all four revised 2026-09-02)

Lineage-B's fig4/fig5/fig7/fig8 differ structurally from Lineage A's in one way beyond
just the broader taxon set: the direct insect~date and fungal~date screens now carry
a quadratic term (§2). insect_PCoA1 and insect_PCoA2 stay as separate single-term
panels, factor(date)-adjusted (§5's strictest test) rather than Lineage A's permissive
no-lure/no-date test. See `project_organization.md` and `iterative_analysis_updates.md`'s
2026-09-01 entries for the full design rationale.

**2026-09-02 revision (all four figures):** fig4/fig5/fig7/fig8 now read the linear
and quadratic terms with the linear term taking PRIORITY whenever it's significant,
giving each taxon exactly ONE "peak timing" category (late / early / mid-season /
bimodal), not up to 2 independent hit-instances as in the first port -- see `project_
organization.md`'s "Peak-timing convention". All numbers below (both the Order/Class
section and the trophic-mode section) now use this reading; the two sections are
directly comparable to each other again.

### Order/Class level, date-associated hits (fig4/fig5, n=2,490)

2,490/6,342 fungal taxa have a significant peak-timing call -- same file as Lineage A
(`fungal_taxa_date_association.all_taxa.csv` is lineage-invariant); this total is
unchanged from the pre-revision "linear OR quadratic significant" count, only how
those 2,490 taxa split into categories changed. Split: **1,272 peak late, 914 peak
early, 215 peak mid-season, 89 bimodal (early+late)**.

Class-level split (fig4 facets; mutually-exclusive peak-timing taxon counts):

| Class | Peaks late | Peaks early | Peaks mid-season | Bimodal |
|---|---|---|---|---|
| Dothideomycetes (Mycosphaerellales/Dothideales/Capnodiales) | 436 | 127 | 62 | 34 |
| Agaricomycetes (Polyporales/Russulales) | 222 | 111 | 6 | 5 |
| Sordariomycetes (Diaporthales/Xylariales/Hypocreales) | 70 | 203 | 8 | 6 |
| Tremellomycetes (Tremellales) | 132 | 5 | 70 | 0 |
| Leotiomycetes (Helotiales) | 46 | 110 | 20 | 5 |
| Eurotiomycetes (Chaetothyriales) | 22 | 88 | 0 | 12 |
| Taphrinomycetes (Taphrinales) | 37 | 2 | 14 | 0 |
| Lecanoromycetes (Lecanorales) | 23 | 61 | 1 | 7 |

The same two-class opposition Lineage A found is unchanged: **Dothideomycetes skews
later-season, Sordariomycetes skews earlier-season**. Under the linear-priority
reading, a meaningful minority of each class's signal is now cleanly separated out as
NOT actually early/late-trending at all, but genuinely peaking mid-season with no
significant linear direction: **Tremellomycetes** is the clearest case -- 70 of its
207 significant taxa (34%) are "peaks mid-season," nearly matching its 132 "peaks
late" -- a taxon-by-taxon split the old hit-instance table (74 hump vs. 123 later,
allowed to overlap) couldn't show directly. **Dothideomycetes** (62/659, 9%) and
**Taphrinomycetes** (14/53, 26%) show the same pattern at smaller scale. These
mid-season taxa should be described as peaking mid-season, not lumped in with the
"later-season" majority of their class.

### Order/Class level, insect-axis-associated hits (fig4/fig5, factor(date)-adjusted)

**PCoA1: 17/200**, spread thin across 10 orders/9 classes -- no single real order
carries more than 3 hits (Chaetothyriales: 2 earlier + 1 later; the 4-hit
"Unclassified" bucket is unresolved-order taxa, not a real order). Direction split
**11 PCoA1+ / 6 PCoA1-**, the reverse tilt from Lineage A's 77-hit, 62-negative/
15-positive no-date-adjustment result (§5 already covers why direct comparison is
risky here: the permissive PCoA1 test collapses toward 0 once ANY date control is
added, and only a non-linear one recovers a comparable count -- 9 linear -> 19
quadratic -> 17-18 factor(date) -- so this is a much smaller, differently-selected
hit set, not a scaled-up version of Lineage A's).

**PCoA2: 9/200** across 4 orders/4 classes, still flagged unstable (§5). Direction:
**8 PCoA2-** (Diaporthales x4, Helotiales x3, Pleosporales x1) vs. **1 PCoA2+**
(Microbotryomycetes_ord_Incertae_sedis). The Diaporthales+Helotiales concentration on
the PCoA2- side is the same Valsaceae/*Cytospora* + Tympanidaceae/*Tympanis* cluster
identified by ASV in §5.

### Trophic mode / FungalTraits (fig7/fig8; "Unclassified genus"/"No FungalTraits match" excluded)

Pooled across classes (fig7), the date signal is again dominated by **plant pathogen**
(303 peak late + 221 peak early + 45 peak mid-season + 12 bimodal) and **wood
saprotroph** (188 peak late + 157 peak early + 29 peak mid-season + 4 bimodal) --
both close to even on the late/early split alone, same as Lineage A's linear-only
read (247/205 and 176/152 respectively). Splitting "plant pathogen" by Class (fig8)
reproduces Lineage A's two-distinct-signals finding almost exactly:

| Class | Lifestyle | Peaks late | Peaks early | Peaks mid-season | Bimodal |
|---|---|---|---|---|---|
| Dothideomycetes | plant pathogen (*Ramularia*-dominated) | 175 | 29 | 12 | 5 |
| Sordariomycetes | plant pathogen (*Cytospora*-dominated) | 36 | 124 | 5 | 1 |

(Lineage A's linear-only read: 173/28 and 34/123 -- essentially the same late/early
numbers, confirming the broader taxon set adds volume but doesn't change this
pattern; the mid-season/bimodal columns are new information Lineage A's linear-only
test can't produce.) "Wood saprotroph" splits similarly by Class:

| Class | Peaks late | Peaks early | Peaks mid-season | Bimodal |
|---|---|---|---|---|
| Agaricomycetes | 129 | 79 | 4 | 2 |
| Sordariomycetes | 10 | 27 | 2 | 0 |
| Dothideomycetes | 34 | 27 | 22 | 2 |

Agaricomycetes skews late, Sordariomycetes skews early, and **Dothideomycetes' wood-
saprotroph signal is substantially mid-season** (22 of 85 significant taxa, 26% --
nearly as many as its "peaks late" bar) -- the same mid-season-peaking pattern §7's
Order/Class section already flags for Dothideomycetes overall, now visible within
this one lifestyle specifically.

For the PCoA-axis trait breakdowns (fig7/fig8 panels b/c), both are too small (17 and
9 total hits, 9 and 7 with a FungalTraits match) to show a clear trait-level trend the
way Lineage A's 77-hit PCoA1 result did -- PCoA1's 9 FungalTraits-matched hits spread
across 5 distinct lifestyle x growth-form combinations (4 lifestyles once growth-form
is dropped, matching fig8 panel b) with no more than 2 hits in any one bucket. PCoA2's
7 FungalTraits-matched hits are all PCoA2-: 6 "plant pathogen : filamentous mycelium"
(Sordariomycetes x4, Leotiomycetes x2 per fig8 panel c) plus 1 "wood saprotroph :
filamentous mycelium" (Dothideomycetes) -- the trait-level restatement of the same
Diaporthales/Helotiales-dominated PCoA2 pattern noted above.

## 8. Responsive-taxa abundance and richness by site (fig6, 2026-09-01 port)

Same matched 69-sample dataset as fig2 (15,422 Kingdom==Fungi ASVs), split into
"date-associated" (q<0.10 linear or quadratic, 2,490/15,422 = 16.2%) vs.
"insect-associated" (q<0.10 PCoA1 or PCoA2, factor(date)-adjusted, 26/15,422 = 0.17%,
no overlap between the two axes' hit sets) categories. Exact per-site numbers are in
`iterative_analysis_updates.md`'s 2026-09-01 fig6 entry (added there in response to a
request for citable exact figures); headline pattern:

- **Date-associated taxa carry 81-88% of relative sequence abundance at every site**
  (lowest at Manchester Airport, 81.4%) despite being only 16% of taxa -- the same
  disproportionate-abundance pattern the informal `responsive_taxa_abundance_by_site.R`
  script first surfaced, now with an exact per-site breakdown and richness counts
  (2,126-2,401 date-associated ASVs per site) alongside the abundance figures.
- **Insect-associated taxa are a much smaller, near-invisible slice** -- 0.005-1.03%
  of relative abundance and 5-26 ASVs per site. Pease stands out numerically (1.03%
  abundance, 26 ASVs, both site maxima), but all four sites carry the full
  Ethanol/Alpha-pinene_EtOH/Ips lure set, so this isn't a lure-design artifact --
  not investigated further.

## Bottom line

1. **Community-level structure (PERMANOVA, NMDS) is qualitatively unchanged** between lineages -- same three significant terms, same rough ranking (lure/date > site), Lineage B just has smaller effect sizes and more residual variance from the broader, noisier taxon set.
2. **The broader taxon set reveals more seasonality, not less** -- both a higher linear hit rate (45% vs. 37%) and a higher hump/dip rate (25% vs. 17%) among insect taxa, consistent with the family filter having suppressed real signal in previously-excluded families.
3. **A new axis (PCoA2) emerged** that has no clean Lineage-A analog -- a within-Scolytinae, lure-and-mid-season-timing axis, distinct from PCoA1's guild-level date split and PCoA3's lure split. Its fungal-association hit count is unstable across date-control strategies and should be reported as a range (0-9/200), not a single number.
4. **The single most important correction from this line of work**: *Cytospora prunicola* (ASV_1905) -- previously the project's headline "most robust cross-cutting taxon" -- does not survive ANY date-adjusted test on any insect PCoA axis under Lineage B. *Tympanis* sp. is the new best candidate for a genuinely date-robust insect-community-linked fungal signal, and the Valsaceae/*Cytospora* and *Tympanis* genus-level groupings (not specific ASVs) are the more defensible framing overall.
5. **The pairwise (individual taxon x taxon) screen adds nothing over Lineage A**: same picture at the broader taxon set -- a handful of early-season Scolytinae (above all *Pityogenes hopkinsi*) account for essentially all q<0.10 pairs, zero pairs survive a global correction, and only *P. hopkinsi* keeps any partners after a strict `factor(date)` control on its predictor. The two extra taxa the all-families table surfaces (*Asemum striatum*, *Orthoperus scutellaris*) both collapse entirely under that control, so they are seasonal co-occurrence, not specific associations.
6. **Order/Class and trophic-mode date patterns are essentially unchanged from Lineage A** -- same Dothideomycetes-later/Sordariomycetes-earlier class split, same *Ramularia*-vs-*Cytospora* two-signal plant-pathogen story (fig8 numbers match Lineage A's linear-only read almost exactly: 175/29 vs. 173/28, 36/124 vs. 34/123). Under fig4/fig5/fig7/fig8's linear-priority peak-timing reading (§7, all four figures revised 2026-09-02), a genuinely new observation stands out: **some of what reads as "later-season" is really mid-season peaking with no significant linear trend at all** -- 34% of Tremellomycetes' significant taxa (70/207), 26% of Taphrinomycetes' (14/53), and 9% of Dothideomycetes' (62/659) are "peaks mid-season," not early or late (§7); the same pattern shows up within a single lifestyle too -- 26% of Dothideomycetes' significant wood-saprotroph taxa (22/85) peak mid-season rather than late (§7 trophic-mode subsection).
7. **The insect-axis (PCoA1/PCoA2) trait/taxonomic breakdowns are far thinner under the strict factor(date)-adjusted hit sets** (17 and 9 hits) than Lineage A's permissive 77-hit PCoA1 result -- no single order or lifestyle dominates PCoA1's spread-thin hits, while PCoA2's hits concentrate on "plant pathogen : filamentous mycelium" (Diaporthales/Helotiales, PCoA2- direction), the same Valsaceae/*Cytospora*+Tympanidaceae/*Tympanis* cluster §5 already flags as the most defensible finding on that axis.
8. **Date-associated fungal taxa carry the large majority of both sequence abundance (81-88%) and richness at every site**, while insect-associated taxa are a near-invisible slice (<=1.03% abundance, <=26 ASVs) -- consistent with, and a direct site-by-site quantification of, the much stronger seasonal vs. insect-community signal established throughout this file (fig6, §8).

Underlying files referenced above: `data/insect_exploratory_prevalence_filtered_insect/insect_prevalence_filter.permanova.csv`, `insect_prevalence_filter.family_composition.csv`; `data/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_taxa_date_association.csv`, `insect_pcoa_axes.date_site_lure_association.csv`, `insect_pcoa_axes.quadratic_date_test.csv`, `insect_pcoa_axes.taxon_drivers.csv`, `date_adjustment_comparison.all_strategies.csv`, `fungal_insect_association_results.PCoA{1,2,3}_*.csv`; `data/2024_fungi/fungal_taxa_date_association.all_taxa.csv` (shared, lineage-independent); `data/2024_insect_data/insect_taxa_date_association.csv` (Lineage A comparison); `data/presentation_items_prevalence_filtered_insect/fig4_taxonomic_breakdown.*.csv`, `fig5_family_breakdown_by_order.*.csv`, `fig7_trait_breakdown.*.csv`, `fig8_trait_breakdown.*.csv` (fig4/5/7/8 source counts, §7); `compare_insects_fungi_procrustes_prevalence_filtered_insect/` (Procrustes r/p).
