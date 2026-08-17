# Interpretation update (2026-08-17): taxonomic summaries and ecological interpretation on the current, Kingdom==Fungi-filtered analyses

`interpretation.md`'s first ~120 lines (initial CoCA taxonomic breakdown, "what's driving insect_PCoA1", and the first PCoA1-vs-date pass) and the PCoA2 section (lines 85-119) were all written **before** `Kingdom == "k__Fungi"` filtering was added to the fungal ASV table (see `interpretation.md:227-258` for the mechanics of that fix) and before several other iterative corrections landed on top of them (the permutation-scheme fix that moved PCoA1-no-date hits from 96/200 to 81/200, the 2-axis -> 3-axis CoCA change, the lure/no-lure comparison). None of those early numbers, taxon lists, or ASV IDs should be quoted going forward. This file re-derives the same taxonomic/ecological summaries directly from the current output files (all post-filter, current as of this writing) so there's one place with numbers that match what's on disk right now. It doesn't restate the mechanistic/methods sections (distance choice, lure partialling, Procrustes, etc.) from `interpretation.md`, which remain valid and unaffected by the filter.

Core current-state numbers (community level, all post-filter, from `interpretation.md:240-249` — reproduced here for reference):

| Result | Value |
|---|---|
| Fungal ASVs retained after Kingdom==Fungi + prevalence filter | 6,342 (of 15,422 confirmed Fungi; 1,717 non-fungal/unidentified ASVs dropped from 17,139 total) |
| Whole-community PERMANOVA (fungal) | site R²=15.0% (p=0.001), lure R²=2.8% (p=0.14), date R²=16.5% (p=0.001), residual 65.7% |
| Whole-community PERMANOVA (insect) | site R²=9.7%, lure R²=14.2%, date R²=14.2% (all p≤0.002) |
| Procrustes (insect NMDS vs. fungal NMDS) | r = 0.56, p = 0.001 |

---

## 1. Insect taxa driving PCoA1 (the seasonal axis)

Insect PCoA1 is still overwhelmingly a date/seasonal-turnover axis (r ≈ 0.78 with collection date, by far the largest term in the insect PERMANOVA). The clean, current-data way to ask "which taxa drive it" is the direct per-taxon date test (`data/2024_insect_data/insect_taxa_date_association.csv`), which supersedes the earlier informal PCoA1-correlation check — this is unaffected by the fungal Kingdom filter (insect side wasn't touched), but the earlier top-level summary understated it slightly (20/52 was reported informally; the current file gives **18/52 significant at q<0.10**).

By subfamily, the direction split is clean and matches the fig4 panel-a breakdown exactly:

| Subfamily | Earlier-season (t<0) | Later-season (t>0) |
|---|---|---|
| Scolytinae (bark/ambrosia beetles) | 11 | 2 |
| Cossoninae | 0 | 2 |
| Corticariinae | 0 | 1 |
| Molytinae | 1 | 1 |

Top hits by \|t\|:

| Taxon | Subfamily | t | Direction |
|---|---|---|---|
| Melanophthalma sp. | Corticariinae | +9.15 | later season |
| Stenoscelis brevis | Cossoninae | +8.79 | later season |
| Xyleborinus attenuatus | Scolytinae | -8.50 | earlier season |
| Xyloterinus politus | Scolytinae | -5.08 | earlier season |
| Pissodes strobi | Molytinae | -4.68 | earlier season |
| Cyclorhipidion pelliculosum | Scolytinae | -4.49 | earlier season |
| Corthylus columbianus | Scolytinae | +4.07 | later season |
| Xylosandrus germanus | Scolytinae | -3.43 | earlier season |
| Himatium errans | Cossoninae | +3.22 | later season |

**Interpretation (unchanged from before the filter, since it doesn't depend on the fungal side):** PCoA1 captures a spring ambrosia/bark-beetle-dominated catch (Scolytinae: *Xyleborinus*, *Xylosandrus*, *Cyclorhipidion*, *Pissodes*) giving way to a later-season catch dominated by fungivorous/detritivorous beetles (*Melanophthalma*, Corticariinae "minute brown scavenger beetles"; *Stenoscelis*, Cossoninae).

---

## 2. Fungal taxa that respond to insect PCoA1 (current 3-axis CoCA screen)

Source: `data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv` — the top-200 CoCA-loading candidates (selected jointly over insect PCoA1-3), tested against PCoA1 alone, lure dropped from the model (established to have no detectable independent effect on fungi; dropping it changes almost nothing: 63/200 sig with lure vs. 77/200 without). This is the variant fig4/fig5/fig7/fig8 downstream all use.

**77 of 200 significant at q<0.10** (up from the pre-filter range of 46-60/200 reported in `interpretation.md:188`). Direction split: **62 taxa negative** (associated with the bark-beetle-dominated/early-season end of PCoA1) vs. **15 positive** (associated with the fungivorous-beetle/later-season end) — still a strongly one-sided pattern, consistent with the original story, just with different specific ASV membership post-filter.

Family tally among the 77 hits (all but 2 of the "+" hits are singleton families):

| Family | Direction | n |
|---|---|---|
| Valsaceae | PCoA1- | 7 |
| Phaffomycetaceae | PCoA1- | 8 |
| Phleogenaceae | PCoA1- | 5 |
| Ophiostomataceae | PCoA1- | 3 |
| Rhytismataceae | PCoA1- | 3 |
| Diatrypaceae | PCoA1- | 2 |
| Spirographaceae | PCoA1- | 2 |
| Tympanidaceae | PCoA1- | 2 |
| Dermateaceae | PCoA1- | 2 |
| Unclassified | PCoA1- / PCoA1+ | 11 / 2 |

Genus-level, with an internal-consistency check (same logic as the original pre-filter writeup — independent ASVs from the same genus/species moving the same direction is a good sign the signal isn't noise):

- **Cytospora** (Valsaceae, canker fungus) — 6 ASVs, all PCoA1- (bark-beetle end). Species where resolved: *C. prunicola*, *C. predappioensis*, *C. phialidica*. Plus *Valsa* sp. (1 ASV), same family, same direction — 7/7 Valsaceae hits agree in sign. This is the same genus flagged pre-filter, though the ASV membership has changed (was 7 ASVs pre-filter, now 6 — not the same set).
- **Wickerhamomyces alni** (Phaffomycetaceae, yeast) — 5 ASVs, all PCoA1-, all the same species. Not mentioned at all in the pre-filter writeup — this is a new top-loading taxon under the current filtered data. *W. alni* is a bark/wood-associated yeast; several *Wickerhamomyces* species are documented from bark-beetle galleries and frass, so this is a biologically plausible companion to the Cytospora signal rather than an isolated artifact.
- **Phleogena faginea** (Phleogenaceae, Atractiellales) — 5 ASVs, all PCoA1-, all the same species. Also new relative to the pre-filter list. *P. faginea* is a bark-surface basidiomycete on dead hardwood — again consistent with the bark-beetle-associated end of the axis.
- **Tympanis** (Tympanidaceae, canker/twig fungus) — 2 ASVs, both PCoA1- (down from 3 pre-filter).
- **Spirographa** (Spirographaceae) — 2 ASVs, both PCoA1- (down from 5 pre-filter; species-level resolution collapsed post-filter, no longer all confidently *S. fusisporella*).
- **Cyberlindnera**, **Pyrenopeziza**, **Lophodermium** — 2 ASVs each, all PCoA1-.
- 15 of 77 hits lack genus-level assignment (all but one on the PCoA1- side).

**Same caveat as before still holds under current data:** these are associations with PCoA1 alone. Adding date as a covariate to the identical 200-candidate set collapses the result to **0/200 significant** (`fungal_insect_association_results.PCoA1_plus_date.csv` / `.no_lure.csv`), exactly as it did pre-filter. The Cytospora/bark-beetle association is the one piece of this that gets a second, independent line of support (see §3) — a direct fungal~date test, run without ever conditioning on the insect axes, points to the same genus and the same seasonal end.

---

## 3. Taxa that respond to date directly (top-loading ASVs, both communities)

### Fungal (unbiased screen, all 6,342 prevalence-filtered Kingdom==Fungi taxa, no pre-selection)

Source: `data/2024_fungi/fungal_taxa_date_association.all_taxa.csv`. **2,116 / 6,342 significant at q<0.10 (33.4%)** — up slightly from 32.6% (2,248/6,887) pre-filter, i.e. filtering out non-fungal ASVs modestly *increased* the seasonal fraction, consistent with those excluded ASVs being mostly low-abundance noise rather than a source of real seasonal signal. Direction split: 1,236 later-season (+) vs. 880 earlier-season (-).

Top 20 hits by \|t\| split almost entirely into the same two genera identified pre-filter, now with corrected ASV IDs and t-statistics:

| Direction | Genus / species | Family / Order / Class | n in top 20 |
|---|---|---|---|
| Later season (+) | *Ramularia* spp. | Mycosphaerellaceae / Mycosphaerellales / Dothideomycetes | 11 |
| Earlier season (-) | *Cytospora* spp. (incl. *C. prunicola*, *C. beilinensis*) | Valsaceae / Diaporthales / Sordariomycetes | 7 |
| Earlier season (-) | *Piptoporus betulinus* (birch polypore) | Fomitopsidaceae / Polyporales / Agaricomycetes | 1 |
| Earlier season (-) | *Ganoderma* sp. (later in top 20) | Ganodermataceae / Polyporales / Agaricomycetes | 1 |

*Ramularia* (foliar plant-pathogen genus, leaf-spot disease) building up over the season as leaf area increases is textbook phenology. *Cytospora* (canker fungus) being higher earlier in the season is now supported two independent ways: this direct date test, and the PCoA1-association test in §2 — both point to the same genus and the same seasonal/insect-community end, without one being derived from the other. Wood-decay taxa (*Piptoporus*, *Ganoderma*) tracking the same early-season/bark-associated direction as *Cytospora* is a third, weaker line pointing the same way.

### Insect (all 52 taxa, direct date test)

See §1 for the full table — 18/52 significant, split cleanly by subfamily (Scolytinae early, Cossoninae/Corticariinae late).

### Where both communities agree

The insect and fungal date signals point at the same underlying gradient (spring bark/ambrosia-beetle activity coinciding with elevated Cytospora canker-fungus abundance; both decline into summer as *Ramularia* and fungivorous beetles increase) — but this is a **shared seasonal trend**, not evidence the two are directly interacting beyond what's already captured in §2's PCoA1-plus-date null result. The Procrustes whole-community concordance (r=0.56) is consistent with — and largely explained by — this same shared site:date structure in both community PERMANOVAs, not proof of a taxon-to-taxon relationship.

---

## 4. Higher-level summaries: Order/Class and trophic mode (current, post-filter data)

These come straight from the current fig4/fig5 (taxonomic) and fig7/fig8 (FungalTraits trophic-mode) output files, all of which post-date the Kingdom filter and were already computed on the current 6,342-taxon table — reproduced/verified here directly from the CSVs rather than from `interpretation.md`'s prose, since a couple of numbers there conflate the class-resolved (fig8) and class-pooled (fig7) versions of the same trait tally.

### Order/Class level, date-associated hits (fig4/fig5, n=2,116)

Cleanly one-directional orders:
- **Mycosphaerellales** (*Ramularia*-dominated): 169 later vs. 12 earlier
- **Tremellales**: 113 later vs. 4 earlier
- **Dothideales**: 97 later vs. 8 earlier
- **Diaporthales** (*Cytospora*-dominated): 87 earlier vs. ~6 later
- **Helotiales**: 76 earlier vs. relatively few later

Class-level split (fig4 facets): **Dothideomycetes** (Mycosphaerellales, Dothideales, Capnodiales) skews later-season; **Sordariomycetes** (Diaporthales, Xylariales, Hypocreales) skews earlier-season — the two dominant fungal classes point in opposite seasonal directions.

### Order/Class level, PCoA1-associated hits (fig4/fig5, n=77 — matches §2's candidate set)

Almost the entire signal sits on the PCoA1-/bark-beetle side: **Saccharomycetales** (10), **Diaporthales** (7, = the Cytospora/Valsa hits), **Atractiellales** (6, = the *Phleogena* hits), **Helotiales** (5), **Ophiostomatales** (3, ambrosia-beetle-associated fungi — *Grosmannia*, *Sporothrix*, *Raffaelea*), **Ostropales** (3, = *Spirographa*), **Rhytismatales** (3). Only a handful of singleton hits (Polyporales, Tremellales, Capnodiales, a few others — 2 each at most) sit on the PCoA1+/later-season side. Class-level: Sordariomycetes (Diaporthales, Ophiostomatales, Xylariales) and Saccharomycetes/Atractiellomycetes are essentially entirely on the PCoA1- side.

### Trophic mode / FungalTraits (fig7/fig8; "Unclassified genus" and "No FungalTraits match" excluded, ~36% of date hits and ~27% of PCoA1 hits)

Pooled across all classes (fig7), the date signal is dominated by two trait categories in both directions: **plant pathogen** (247 later vs. 205 earlier) and **wood saprotroph** (176 later vs. 152 earlier). Pooling masks a real split, though — splitting "plant pathogen" by Class (fig8) shows it is **not one signal**, it's two taxonomically distinct, oppositely-directed ones:

| Class | Lifestyle | Later (+) | Earlier (-) |
|---|---|---|---|
| Dothideomycetes | plant pathogen (*Ramularia*-dominated) | 173 | 28 |
| Sordariomycetes | plant pathogen (*Cytospora*-dominated) | 34 | 123 |

This is the trait-level restatement of the Order/Class pattern above — foliar leaf-spot pathogens build up over the growing season; canker pathogens track the early-season bark-beetle-associated end. Pooling across Class (as the fig7 view does) makes "plant pathogen" look only mildly seasonal (247 vs. 205, close to even) when in fact it's two strong, opposite, taxonomically separable trends that largely cancel out.

For PCoA1 hits (fig7/fig8, n=77), trait categories skew heavily toward the earlier-season/bark-beetle end almost everywhere: plant pathogen 19 earlier vs. 1 later; wood saprotroph 10 vs. 2; litter/nectar-tap saprotroph yeasts (5 vs. 0, matching the *Wickerhamomyces*/§2 finding) all earlier. Sordariomycetes plant pathogens alone account for 11 of the 19 earlier-season "plant pathogen" PCoA1 hits — i.e. the PCoA1 trait pattern is substantially the same *Cytospora*-driven signal seen in §2 and §3, viewed through a different lens.

---

## Bottom line (current, post-filter state)

The core conclusions from before the filter still hold, and in most cases are now on slightly stronger or cleaner footing:

1. **Date is the dominant driver of both communities** (insect PCoA1 r=0.78 with date; fungal community PERMANOVA date R²=16.5%, the largest fungal term).
2. **No fungal taxon's association with insect-community structure (PCoA1 or PCoA3) survives adjusting for date** — 77/200 -> 0/200 for PCoA1; PCoA3 null throughout, with or without date.
3. **The one relationship that shows up two independent ways — a direct fungal~date test (§3) and a fungal~PCoA1 test never conditioned on date (§2) — is *Cytospora* (Valsaceae canker fungi) tracking the early-season, bark/ambrosia-beetle-dominated end of the community.** This remains the most defensible candidate for a real (not purely date-confounded) ecological link, though it still isn't proof of direct interaction.
4. **New under the current filtered data**: *Wickerhamomyces alni* and *Phleogena faginea*, both bark/wood-associated fungi absent from the pre-filter top-hit lists, now sit alongside *Cytospora* among the strongest PCoA1-/bark-beetle-associated taxa (5 ASVs each, single-species, same-sign internal consistency) — worth folding into the ecological narrative alongside Cytospora going forward.
5. **Plant-pathogen seasonality is two distinct signals, not one** — Dothideomycetes/*Ramularia* (foliar, later-season) vs. Sordariomycetes/*Cytospora* (canker, earlier-season/bark-beetle-associated) — visible only once trait breakdowns are split by Class (fig8) rather than pooled (fig7).

Underlying files referenced above (all current, post-filter): `data/2024_insect_data/insect_taxa_date_association.csv`; `data/2024_fungi/fungal_taxa_date_association.all_taxa.csv`; `data/compare_insects_fungi_top3axes/fungal_insect_association_results.PCoA1_only.no_lure.csv` and `.PCoA1_plus_date.no_lure.csv`; `data/presentation_items/fig4_taxonomic_breakdown.*.csv`, `fig5_family_breakdown_by_order.*.csv`, `fig7_trait_breakdown.*.csv`, `fig8_trait_breakdown.*.csv`; `data/2024_fungi/ASVs_taxonomy.tsv` (taxonomy join).
