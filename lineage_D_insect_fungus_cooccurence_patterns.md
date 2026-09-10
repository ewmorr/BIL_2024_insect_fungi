# Interpretation update (2026-09-10): Lineage D -- bark-beetle x Ophiostomatales presence/absence co-occurrence

This file introduces **Lineage D** and is the co-occurrence counterpart to
`lineage_A_family_filtered_insect_top_level_interepretation.md`,
`lineage_B_prevalence_filtered_insect_top_level_interpretation.md`, and
`lineage_C_alpha_div_full_insect_table_top_level_interpretation.md`. Same
purpose as those: pull the current, on-disk numbers into one place. Read
`project_organization.md` first for the Lineage A/B/C insect-table
constructions -- Lineage D reuses one of them (see §1) rather than defining a
fourth.

**Status of the "Lineage D" label**: as of this writing it lives only in this
document. `project_organization.md`'s script/output index and
`iterative_analysis_updates.md`'s chronological log do **not** yet describe
this line of work; the scripts sit in `insect_exploratory/` and their outputs
in `data/insect_exploratory/{target_genera,scolytinae}/` and
`figures/insect_exploratory/{target_genera,scolytinae}/`. Wire Lineage D into
those two top-level docs when this line is taken past exploratory.

---

## 1. What Lineage D is

Lineage D is a **targeted, hypothesis-driven co-occurrence analysis** of the
beetle taxa the EDRR traps are actually baited for -- the bark/ambrosia
beetles (subfamily **Scolytinae**) -- against the fungal order most tied to
their natural history, **Ophiostomatales** (the ophiostomatoid "blue-stain"
fungi). It is not a whole-community method (Lineage A/B ordination/CoCA) and
not a diversity method (Lineage C); it asks, one beetle taxon against one
fungal taxon, whether the two are found in the same samples more (or less)
than expected.

- **Metric**: Jaccard similarity `J = a / (a + b + c)` on presence/absence,
  where `a` = samples with both present, `b` = beetle only, `c` = fungus
  only. Shared absence (`d`) is deliberately excluded -- co-presence is the
  question and joint absence is uninformative when most samples lack most
  rare taxa. This is what the user asked for over an abundance correlation:
  "focus on shared presence, downweight shared absence."
- **Comparisons are beetle x fungal-taxon only.** No beetle-vs-beetle and no
  fungus-vs-fungus pairs are computed.
- **Insect table**: the full all-families trap-catch table restricted to
  Subfamily == Scolytinae, presence/absence. This is a *slice of the Lineage
  C raw table*, not the Lineage A family-filter or the Lineage B >=5-sample
  prevalence filter -- Scolytinae is the biological unit of interest, so the
  filter is taxonomic, not statistical. A per-taxon >=3-sample prevalence
  floor (in the 69 matched samples) is applied only to decide which beetles
  carry enough presences to test.
- **Fungal table**: raw (un-rarefied) ASV table, Kingdom-invariant filter to
  Order == `o__Ophiostomatales`, presence/absence. Raw rather than rarefied
  so the rare-tail Ophiostomatales detections -- the whole point here -- are
  kept; `num_reads` spans ~100x across samples, so this read-depth caveat
  bites hardest at ASV resolution (§9). Rolled up three ways, each run as its
  own grid: **species** (named Ophiostomatales species, prevalence >= 3),
  **genus** (the 8 genera + an "unclassified genus" bin + an "Any
  Ophiostomatales" union column), **ASV** (individual ASVs, prevalence >= 3).
  A genus/species/order unit is "present" in a sample if any constituent ASV
  is.
- **Matched dataset**: the same 69 samples with both insect and fungal data
  used by `compare_insect_fungi/insect_fungal_pairwise_taxon_association.r`.
  68 Ophiostomatales ASVs occur in these 69 samples; "Any Ophiostomatales" is
  present in 44 of 69.

Note: the taxonomy contains `Graphilbum ipis-grandicollis` (ASV_3518) --
literally named for *Ips grandicollis* -- but that ASV was dropped during ASV
table construction and is absent from `ASV_tab.csv`, so it cannot be tested
here.

Scripts: `insect_exploratory/target_genus_phenology.R`,
`target_species_ophiostomatales_cooccurrence.R`,
`scolytinae_species_phenology.R`,
`scolytinae_ophiostomatales_cooccurrence.R`,
`scolytinae_ophiostomatales_cooccurrence.residualized.R`, and the shared
engine `cooccurrence_lib.R` (`jaccard_grid_test()`).

## 2. Target-genus phenology (Ips, Dendroctonus)

`target_genus_phenology.R`, all insect trap data (not restricted to the 69
fungal-matched samples): 6 collection dates 1-May to 10-Jul 2024, 4 sites, 3
lures (E = Ethanol, EA = Alpha-pinene + Ethanol, I = Ips lure), pooled across
sites.

Only **three** target species exist in the data: *Ips grandicollis* (92
individuals), *Ips pini* (8), *Dendroctonus valens* (199). Season totals by
lure:

| Species | Ethanol | Alpha-pinene + Ethanol | Ips lure |
|---|--:|--:|--:|
| *Dendroctonus valens* | 2 | **197** | 0 |
| *Ips grandicollis* | 0 | **92** | 0 |
| *Ips pini* | 1 | 0 | **7** |

**The lures partition the three species almost perfectly.** Alpha-pinene +
ethanol takes essentially all of *D. valens* and *every* *I. grandicollis*
individual; the Ips lure takes *I. pini* and nothing else (0 *I.
grandicollis*, 0 *D. valens*); ethanol catches almost nothing. Both
EA-caught species peak sharply on 15-May and decline through the season. So
the earlier read that "the Ips lure barely works" was wrong at the genus
level -- it is **species-specific to *I. pini***, consistent with its
ipsdienol-type pheromone blend vs. *I. grandicollis*'s ipsenol-based
aggregation pheromone (not in this lure), so *I. grandicollis* is picked up
only as bycatch on the host-volatile lure. Caveat: *I. pini* (n=8) and the
ethanol column are too thin for more than presence/absence.

## 3. Scolytinae-wide phenology

`scolytinae_species_phenology.R`, all insect trap data. **45 Scolytinae
species, 28 genera, 5,280 individuals.** Dominants: *Xylosandrus germanus*
(2,104), *Anisandrus sayi* (522), *Xyleborinus attenuatus* (582),
*Xylosandrus crassiusculus* (347), *Xyleborinus saxesenii* (285),
*Cyclorhipidion pelliculosum* (243), *Dendroctonus valens* (199),
*Heteroborips seriatus* (154), *Orthotomicus caelatus* (129), *Dryocoetes
autographus* (113).

**Scolytinae as a group are strongly spring-weighted**: subfamily catch and
species richness peak 15-May to 29-May on the Ethanol and Alpha-pinene +
Ethanol lures (~24-25 species per collection on EA at the peak), then
collapse after 12-Jun. The Ips lure catches far less throughout (peak ~17
species). A handful of species are July-only (*Corthylus columbianus*,
*Lymantor decipiens*). This spring concentration is the source of the
seasonality confound that dominates §7.

Outputs: `scolytinae_species_abundance_heatmap.{png,pdf}` (species x date,
faceted by lure, rows ordered by abundance-weighted mean date),
`scolytinae_subfamily_totals_by_lure.{png,pdf}`.

## 4. Method: Jaccard co-occurrence + restricted permutation

For each (beetle, fungal-unit) pair: observed Jaccard, then a null built by
**permuting the beetle presence/absence vector only** (the fungal vector is
held observed) and recomputing Jaccard, 1,999 times.

- `p_pos` = add-one `Pr(J_perm >= J_obs)` (co-occurrence), `p_neg` = add-one
  `Pr(J_perm <= J_obs)` (avoidance), `p_perm = 2 * min(p_pos, p_neg)` capped
  at 1.
- `z_score = (J_obs - mean(J_perm)) / sd(J_perm)` -- a standardized effect
  size, analogous to this project's `beta_std` convention for the regression
  screens; it is what the figures colour on, so the eye tracks
  excess-over-null rather than raw prevalence.
- BH-FDR (`q_perm`) within each grid; project convention is **q < 0.10**.
- Fisher's exact test on the 2x2 is reported alongside as a structure-blind
  cross-check (`q_fisher`), but it ignores the site/trap/season structure and
  overstates significance -- treat `p_perm`/`q_perm` as primary.

Permutation schemes (via `permute::how`):

| scheme | blocks | what it holds fixed | used in |
|---|---|---|---|
| `within_trap` | site, plots = trap (type none), within = free | site + trap + lure; **NOT season** | §5, §6 |
| `within_date` | collection date | **season** (each date's beetle prevalence) | §8, §9 |
| `within_site_x_date` | site x date | season + site (permutes beetle among the 3 lure-traps of a site-date; strict, low power) | §8, §9 |

`within_trap` and `within_date` cannot be combined: each trap has exactly one
sample per date, so "within trap and within date" leaves nothing to permute.

## 5. Three target species x Ophiostomatales (`target_species_ophiostomatales_cooccurrence.R`, within_trap null)

Beetle presence in the 69 matched samples: *I. grandicollis* 19, *I. pini* 5,
*D. valens* 14. Grids: 20 species, 10 genus units, 39 ASVs.

**Only one pair survives FDR: *Dendroctonus valens* x *Raffaelea* (genus),
q_perm = 0.09** (J = 0.263 vs null 0.085, z = 3.77; both present in 5 of *D.
valens*'s 14 samples; Fisher p = 0.024). At species level the same signal
appears as *Raffaelea arxii* (p_perm = 0.003) but does not clear BH over 60
tests.

*D. valens* co-occurs positively with **most** Ophiostomatales genera (the
entire *D. valens* column is z > 0), strongest after *Raffaelea*:
*Grosmannia* (z = 3.0), *Ophiostoma* (z = 2.2), Any Ophiostomatales
(z = 2.2). **Neither *Ips* species shows any co-occurrence beyond chance** --
*I. grandicollis* has high *raw* Jaccard with common genera (*Leptographium*
0.39, *Sporothrix* 0.29) but z near 0 or negative once trap structure is
accounted for; *I. pini* (n = 5) is uninformative. The broad, non-specific
positivity of the *D. valens* row is already a hint of a shared driver rather
than a specific symbiosis (see §7).

Outputs: `beetle_ophiostomatales_jaccard.{species,genus,asv}.csv`,
`beetle_ophiostomatales_jaccard.{species,genus,asv}.{png,pdf}`.

## 6. All Scolytinae x Ophiostomatales (`scolytinae_ophiostomatales_cooccurrence.R`, within_trap null)

Beetle set: 34 Scolytinae species at >= 3-sample prevalence in the 69 matched
samples. Grids: 20 species, 10 genus units, 39 ASVs.

**Nothing survives FDR** over the larger grid (34 beetles x 10/20/39 fungal
units; the *D. valens* x *Raffaelea* pair from §5 is now q_perm ~ 0.14, same
z, larger correction). Signal is concentrated in one block of beetles.
Ranked by mean z across the 9 Ophiostomatales genera:

| Beetle | mean z | genera at p_pos < 0.05 |
|---|--:|--:|
| *Pityogenes hopkinsi* | 3.11 | 8 / 9 |
| *Xyleborinus attenuatus* | 2.25 | 4 |
| *Cyclorhipidion pelliculosum* | 2.08 | 4 |
| ***Dendroctonus valens*** | 2.00 | 4 |
| *Heteroborips seriatus* | 1.99 | 6 |
| *Xyloterinus politus* | 1.88 | 5 |
| *Hylesinus aculeatus* | 1.61 | 2 |
| *Hylastes opacus* | 1.59 | 3 |
| ... | | |
| *Dryocoetes autographus* | **-1.71** | 0 (avoidance) |

*Pityogenes hopkinsi* stands out -- broad positive co-occurrence with
Ophiostomatales, several genera at q_perm ~ 0.11 (*Graphilbum*,
*Leptographium*, *Raffaelea*, Any Ophiostomatales). This is consistent with
*P. hopkinsi* dominating the project's CLR-abundance residualized pairwise
screen (`insect_fungal_pairwise_taxon_association.residualized.r`, see
`project_organization.md`). The fungal taxon most often flagged as a strong
partner is *Grosmannia francke-grosmanniae* (species-grid J = 0.62 with
*Heteroborips seriatus*, 0.48 with *Xyleborinus attenuatus*, 0.44 with
*Cyclorhipidion pelliculosum*). Both *Ips* species sit near zero, matching
§5.

Outputs: `scolytinae_ophiostomatales_jaccard.{species,genus,asv}.csv` and
`.{png,pdf}`.

## 7. The seasonality problem with the within_trap null

The within_trap null permutes only the beetle vector; the fungal vector keeps
its full seasonal signal. Its null hypothesis is *"beetle presence is
date-uniform within a trap."* Rejecting it flags any date-structured
association -- **including two partners that simply both peak in spring**
(§3), which is not a direct association. The uniformly-positive block in §6 is
exactly that signature: a shared seasonal driver, not 8 separate symbioses.
Controlling a confounder by permutation requires preserving its relationship
with *both* variables, so the fix is to permute the beetle vector *within
collection-date strata* (§4, `within_date` / `within_site_x_date`) -- the
presence/absence analogue of the project's `.residualized.r` scripts that
residualize on `site + factor(date)`.

## 8. Date-conditional re-test of the top block, genus level (`scolytinae_ophiostomatales_cooccurrence.residualized.R`)

Re-tests an 11-beetle set with the **same observed Jaccard** under all three
nulls. The set: the 8-beetle **block** (mean genus z >= 1.5 from §6) plus 3
**reference** rows carried through as contrast, not part of the block's FDR
set:

- *Ips grandicollis*, *Ips pini* -- flat anchors (near-zero throughout §5-6).
- ***Dryocoetes autographus*** -- a **negative control**. It was the single
  most negative beetle in §6 (mean z = -1.71, apparent avoidance); it is a
  mid-season EA-lure species peaking after the Ophiostomatales decline, so
  the "avoidance" looked like the same seasonal confound running in reverse.
  Watching it collapse under date control is a two-sided check on the null.

Mean z per beetle across the 9 Ophiostomatales genera:

| Beetle | in block | z within_trap | z within_date | z site x date | verdict |
|---|:--:|--:|--:|--:|---|
| *Pityogenes hopkinsi* | Y | 3.11 | **2.67** | **1.81** | date-robust |
| *Dendroctonus valens* | Y | 2.00 | **1.88** | **1.65** | date-robust |
| *Heteroborips seriatus* | Y | 1.99 | 1.62 | 1.19 | partly retained |
| *Xyleborinus attenuatus* | Y | 2.25 | 1.55 | 1.51 | partly retained |
| *Cyclorhipidion pelliculosum* | Y | 2.08 | 0.88 | 0.82 | mostly seasonal |
| *Xyloterinus politus* | Y | 1.88 | **-0.17** | -0.52 | entirely seasonal |
| *Hylesinus aculeatus* | Y | 1.61 | 0.49 | 0.25 | seasonal |
| *Hylastes opacus* | Y | 1.59 | 0.21 | -0.21 | seasonal |
| *Ips pini* | ref | 0.89 | 0.27 | -0.08 | flat |
| *Ips grandicollis* | ref | 0.17 | **1.56** | **1.56** | signal emerges under date control |
| *Dryocoetes autographus* | ref | **-1.71** | 0.27 | 0.48 | "avoidance" was seasonal |

Controlling season roughly halves the block. *Xyloterinus politus*, *Hylastes
opacus*, *Hylesinus aculeatus* and most of *Cyclorhipidion pelliculosum* were
seasonal coincidence; the *Dryocoetes autographus* avoidance was too. And a
date-conditional signal **appears** for *Ips grandicollis* that the
within_trap null had been absorbing.

**Pairs surviving BH (q_perm < 0.10) under `within_date` -- all genus level (8 of 110 pairs):**

| Beetle | Fungal genus | J | z | q_perm |
|---|---|--:|--:|--:|
| *Pityogenes hopkinsi* | *Graphilbum* | 0.56 | 5.57 | 0.028 |
| *Pityogenes hopkinsi* | *Leptographium* | 0.44 | 4.10 | 0.028 |
| *Pityogenes hopkinsi* | Any Ophiostomatales | 0.43 | 3.76 | 0.028 |
| *Dendroctonus valens* | *Leptographium* | 0.48 | 4.76 | 0.028 |
| *Pityogenes hopkinsi* | *Ceratocystiopsis* | 0.43 | 4.49 | 0.037 |
| *Heteroborips seriatus* | *Grosmannia* | 0.50 | 3.94 | 0.037 |
| *Heteroborips seriatus* | *Raffaelea* | 0.41 | 3.41 | 0.083 |
| *Ips grandicollis* | *Leptographium* | 0.39 | 3.87 | 0.083 |

Under the strict `within_site_x_date` null, only two clear BH:
**\*Dendroctonus valens\* x \*Leptographium\* (q_perm = 0.055)** and **\*Ips
grandicollis\* x \*Leptographium\* (q_perm = 0.055)**; *Pityogenes* x
*Graphilbum*/*Leptographium*/*Ophiostoma*, *Heteroborips* x
*Grosmannia*/*Ophiostoma*/Any, and *Xyleborinus attenuatus* x *Leptographium*
persist at uncorrected p_perm < 0.05 but not FDR (that null has little
power).

Outputs: `scolytinae_ophiostomatales_jaccard.date_conditional_comparison.csv`
(wide, one row per pair, all three nulls), `..._date_conditional_long.csv`,
`scolytinae_ophiostomatales_jaccard.date_conditional.genus.{png,pdf}`,
`..._date_conditional.shrinkage.{png,pdf}`.

## 9. Date-conditional at species and ASV resolution

FDR power drops sharply with resolution -- of the 11-beetle grid, `within_date` hits:

| grid | pairs | p_perm < 0.05 | q_perm < 0.10 |
|---|--:|--:|--:|
| genus | 110 | 14 | **8** |
| species | 220 | 18 | **0** |
| ASV | 429 | 30 | **0** |

So genus is the only level where anything clears BH; species/ASV are read on
uncorrected p_perm **and on consistency across resolutions and nulls**.

**Species level, `within_date`, p_perm < 0.05 (recurring partners):**

| Fungal species | Beetles (z) | survives site x date? |
|---|---|---|
| *Leptographium gracile* | *D. valens* (4.0), *Ips grandicollis* (3.3), *Xyleborinus attenuatus* (3.1) | yes, all three |
| *Grosmannia francke-grosmanniae* | *Heteroborips seriatus* (**6.0**, J = 0.62 -- strongest single pair anywhere), *Xyleborinus attenuatus* (2.8) | yes |
| *Raffaelea arxii* | *Heteroborips seriatus* (4.1) | yes |
| *Ceratocystiopsis manitobensis* | *Pityogenes hopkinsi* (4.2) | -- |
| *Ophiostoma pallidulum* / *O. denticulatum* | *D. valens*, *Ips grandicollis* | partly |
| *Sporothrix eucastaneae* | *Hylesinus aculeatus* (3.5), *Ips pini* (3.4) | -- (both thin, shared_n 5 / 3) |
| *Sporothrix rossii* | *Pityogenes hopkinsi* (3.1) | -- |

**ASV level -- each species-level signal resolves to essentially one ASV**,
not spread across variants:

| Fungal species (ASV) | Beetles |
|---|---|
| *Leptographium gracile* (**ASV_6117**) | *D. valens*, *Ips grandicollis*, *Xyleborinus attenuatus* (all also survive site x date) |
| *Grosmannia francke-grosmanniae* (**ASV_660**) | *Heteroborips seriatus*, *Xyleborinus attenuatus* |
| *Raffaelea arxii* (**ASV_1122**) | *Heteroborips seriatus* |
| *Ceratocystiopsis manitobensis* (**ASV_8895**) | *Pityogenes hopkinsi* |
| *Sporothrix rossii* (**ASV_8274**) | *Pityogenes hopkinsi*, *D. valens*, *Ips grandicollis* |
| *Ophiostoma pallidulum* (**ASV_15258**) / *O. denticulatum* (**ASV_13249**) | *D. valens*, *Ips grandicollis* |
| *Sporothrix* sp. (**ASV_5162**) | *Cyclorhipidion pelliculosum* -- the one *Cyclorhipidion* signal that survives date control, even though its genus-level mean z collapsed |

That the signals are carried by single ASVs (rather than diluted across
co-occurring variants of a species) is itself a result; the strongest of them
-- *L. gracile* / ASV_6117 with the pine-associated beetles -- is firm at
genus, species, and ASV resolution and under the strict null.

Outputs: `scolytinae_ophiostomatales_jaccard.date_conditional.{species,asv}.{png,pdf}`;
species/ASV rows are in the same
`..._date_conditional_comparison.csv` / `..._long.csv` as §8, with
`asv_genus` / `asv_species` annotation.

## Bottom line

1. **The traps sample three target species, and the lures cleanly partition
   them** (§2): alpha-pinene + ethanol -> all of *Dendroctonus valens* and
   *Ips grandicollis*; Ips lure -> *Ips pini* only; ethanol -> almost
   nothing. Any co-occurrence analysis of "Ips + Dendroctonus" is really an
   analysis of these three species, dominated by *D. valens* (n = 199) and
   *I. grandicollis* (n = 92).
2. **Scolytinae are strongly spring-weighted** (§3): subfamily catch and
   richness peak mid-to-late May, collapse after mid-June. This is the
   confound that governs everything downstream.
3. **The within-trap permutation null does not control season** (§7). Under
   it, a whole block of early-season Scolytinae -- *Pityogenes hopkinsi*,
   *Xyleborinus attenuatus*, *Cyclorhipidion pelliculosum*, *Dendroctonus
   valens*, *Heteroborips seriatus*, *Xyloterinus politus*, *Hylesinus
   aculeatus*, *Hylastes opacus* -- turns uniformly positive against nearly
   every Ophiostomatales genus (§6), which is the signature of shared
   phenology, not specific association. Nothing survives FDR at the
   subfamily scale under this null.
4. **Under the within-trap null, the only FDR-surviving pair in the 3-species
   analysis is *Dendroctonus valens* x *Raffaelea* (genus, q_perm = 0.09)**
   (§5). Neither *Ips* species shows anything.
5. **Controlling season (permute within collection date) cuts the block
   roughly in half** (§8). Confirmed seasonal artifacts: *Xyloterinus
   politus* (z 1.88 -> -0.17), *Hylastes opacus*, *Hylesinus aculeatus*,
   most of *Cyclorhipidion pelliculosum*, and *Dryocoetes autographus*'s
   apparent avoidance.
6. **What survives season control, FDR, and (where noted) the strict site x
   date null is a short, interpretable list**:
   - ***Pityogenes hopkinsi*** -- broad Ophiostomatales carrier (*Graphilbum*,
     *Leptographium*, *Ceratocystiopsis*, Any Ophiostomatales all q_perm <
     0.10), sharpening at ASV level to *Ceratocystiopsis manitobensis*
     (ASV_8895) and *Sporothrix rossii* (ASV_8274). Consistent with its
     status in the project's CLR pairwise screen.
   - ***Dendroctonus valens* x *Leptographium* (gracile / ASV_6117)** --
     genus q_perm = 0.028, and one of only two pairs to clear BH under the
     strict site x date null (q_perm = 0.055).
   - ***Heteroborips seriatus*** with ***Grosmannia francke-grosmanniae***
     (ASV_660; J = 0.62, z = 6.0 -- the strongest single pair) and
     ***Raffaelea arxii*** (ASV_1122).
   - ***Ips grandicollis* x *Leptographium* (gracile / ASV_6117)** -- a
     date-conditional signal that the within-trap null had masked (§8); the
     other pair to clear BH under the strict null (q_perm = 0.055).
7. **The firmest single theme is *Leptographium gracile* (ASV_6117) with the
   pine bark/ambrosia beetles** -- *D. valens*, *Ips grandicollis*,
   *Xyleborinus attenuatus* -- robust at genus, species, and ASV resolution
   and under the strict null.
8. **Caveats**: nothing survives FDR below genus resolution, so species/ASV
   entries in §9 are candidates, not confirmed. Read-depth bias (100x range
   in `num_reads`) is worst at ASV level. Several flagged cells are thin
   (shared_n 2-3 with a high z, e.g. *Xyleborinus attenuatus* x ASV_6117).
   The strict `within_site_x_date` null has limited power (many degenerate
   3-sample strata). And this is presence/absence only on 69 samples -- an
   abundance-based or larger-sample follow-up would be the way to firm up any
   of the §6-9 candidates.

Underlying files:
- Phenology: `data/insect_exploratory/target_genera/target_species_totals_by_lure.csv`,
  `target_{genus,species}_catch_by_date_lure.csv`;
  `data/insect_exploratory/scolytinae/scolytinae_species_abundance_by_date_lure.csv`,
  `scolytinae_subfamily_totals_by_date_lure.csv`,
  `scolytinae_species_phenology_order.csv`; figures in
  `figures/insect_exploratory/{target_genera,scolytinae}/`.
- Co-occurrence (within_trap): `data/insect_exploratory/target_genera/beetle_ophiostomatales_jaccard.{species,genus,asv}.csv`,
  `data/insect_exploratory/scolytinae/scolytinae_ophiostomatales_jaccard.{species,genus,asv}.csv`,
  `*_presence_absence.csv`.
- Date-conditional: `data/insect_exploratory/scolytinae/scolytinae_ophiostomatales_jaccard.date_conditional_comparison.csv`,
  `..._date_conditional_long.csv`; figures
  `figures/insect_exploratory/scolytinae/scolytinae_ophiostomatales_jaccard.date_conditional.{genus,species,asv}.{png,pdf}`,
  `..._date_conditional.shrinkage.{png,pdf}`.
