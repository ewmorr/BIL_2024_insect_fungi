### ###############################################################
### Initial Co-correspondence analysis (CoCA) plus per-taxon lms
### against insect_PCoA1 to identify fungal taxa that are 
### correlated with insect community variance
### ###############################################################

Family/genus composition of the 96 significant CoCA hits:
- Valsaceae is the top family (10 hits), driven mostly by Cytospora (7 ASVs) — a canker-causing fungal genus (Diaporthales) commonly associated with bark beetle/weevil vectoring in trees. Given the insect data here is Curculionidae (weevils, including Scolytinae bark beetles) and Latridiidae, this is a biologically very plausible hit, not just statistical noise.
- Spirographa (Spirographaceae, 5 ASVs) — all annotated as S. fusisporella, likely the same organism split across ASVs.
- Tympanis (Tympanidaceae, 3 ASVs) — another canker/twig-associated genus.
- Cladosporium, Jattaea, Valsa, Pyrenopeziza, Bjerkandera each appear among the top 20 by q-value.
- 23 of the 96 significant hits lack family-level taxonomy (unclassified/low-confidence assignments).

Internal consistency check (a good sign the signal is real rather than noise): all 7 Cytospora ASVs, all 5 Spirographa ASVs, and all 3 Tympanis ASVs share the same sign on their t-statistic — independent ASVs from the same genus are moving together with insect community composition, which you would not expect by chance if these were spurious hits.

Overall 75 of 96 significant taxa share one sign and 21 the other — worth checking later whether that split lines up with a specific direction on the insect ordination (e.g., something distinguishing certain traps/lures) rather than a generic gradient.

### ###############################################################
### Identify what is driving insect_PCoA1 
### ###############################################################

PCoA1 is overwhelmingly a seasonal/date gradient (r = 0.78, F = 168, p < 2e-16 — by far the largest term), with weaker but real Lure (F=22) and Site (F=3.9) effects. Let me find which insect taxa actually drive that axis.


This clarifies the whole story. Insect PCoA1 is a seasonal turnover axis (r=0.78 with collection date, explaining far more variance than Site or Lure), and it's driven by a clear compositional shift:

┌──────────────────┬─────────────────────────────────────────────────────────────────────────────────────────────────┬────────────────────────────────────┐
│    Direction     │                                          Driving taxa                                           │              Pattern               │
├──────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────┼────────────────────────────────────┤
│ Positive (later  │ Melanophthalma sp. (cor = 0.96, n=1780)                                                         │ Latridiidae — fungivorous "minute  │
│ season)          │                                                                                                 │ brown scavenger beetles"           │
├──────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────┼────────────────────────────────────┤
│ Negative         │ Xyleborinus attenuatus (-0.82), Xylosandrus germanus (-0.56, n=2104), Cyclorhipidion            │ Scolytinae — ambrosia/bark         │
│ (earlier season) │ pelliculosum, Pissodes strobi, Dendroctonus valens, Xyloterinus politus,                        │ beetles, typical spring-flight     │
│                  │ Pityophthorus/Pityogenes spp.                                                                   │ taxa                               │
└──────────────────┴─────────────────────────────────────────────────────────────────────────────────────────────────┴────────────────────────────────────┘

So PCoA1 essentially captures a spring bark/ambrosia-beetle-dominated catch giving way to a later-season fungivorous-beetle-dominated catch.

This ties directly back to the CoCA fungal hits: 75 of 96 significant fungal taxa (negative t-statistic) — including all the Cytospora, Spirographa, and Tympanis ASVs — increase in abundance toward the negative (bark-beetle-dominated, early-season) end of PCoA1. Cytospora is a canker fungus frequently associated with bark/ambrosia beetle activity on stressed trees, so this isn't just a seasonal read-count coincidence — it lines up with a real, plausible insect-vectoring or shared-substrate relationship. The 21 positive-direction hits (e.g. Sporobolomyces, one Didymellaceae ASV) would be worth checking for association with the Melanophthalma/later-season end instead.

One caution: because PCoA1 is so strongly collinear with date, and date isn't in the per-taxon test's fixed-effect blocking (only site and lure are), what looks like an "insect community" signal could partly just be a shared seasonal trend in both insect and fungal abundance rather than a direct biotic interaction. Want me to add date as a covariate in the per-taxon model to see how many hits survive that stricter test?

### ###############################################################
###  after adding date as a covariate in the per-taxon models
### ###############################################################

Summary: 0 of 200 significant at q < 0.10 (down from 96/200), only 4 at q < 0.20

The four closest (q = 0.15, all still short of significance):

┌──────────┬──────────────────────────────────────────────────┬────────┐
│  Taxon   │             Family / Genus / Species             │ t_stat │
├──────────┼──────────────────────────────────────────────────┼────────┤
│ ASV_199  │ Phaffomycetaceae, Cyberlindnera misumaiensis     │ -3.56  │
├──────────┼──────────────────────────────────────────────────┼────────┤
│ ASV_2007 │ Saccharomycetales inc. sedis, Candida ponderosae │ -3.09  │
├──────────┼──────────────────────────────────────────────────┼────────┤
│ ASV_1905 │ Valsaceae, Cytospora prunicola                   │ -3.53  │
├──────────┼──────────────────────────────────────────────────┼────────┤
│ ASV_534  │ Valsaceae, Valsa sp.                             │ -3.65  │
└──────────┴──────────────────────────────────────────────────┴────────┘

What this means: the earlier "48% hit rate" was almost entirely an artifact of insect_PCoA1 and fungal seasonal abundance both tracking collection date independently — not evidence of a direct fungal/insect-community relationship. Once date is properly partialled out, the signal collapses. This is exactly the caution flagged earlier, now confirmed quantitatively.

It's worth noting Cytospora prunicola (ASV_1905) and the Valsaceae group still sit at the top even after this much stricter test (q=0.15, close-ish to conventional significance) — not proof of anything, but a mild "still worth a closer look with more data" signal rather than one that should be dropped outright. Everything else in the earlier list (Spirographa, Tympanis, the broader Cytospora set) does not survive.

### ###############################################################
### Explain Euclidean distance usage (it's actually Hellinger)
### ###############################################################

 Why did you use euclidean distance in the insect community PCoA instead of Bray-Curtis distance?

That's a comment already baked into the CoCA script (insect_dist <- vegdist(insect_hel, method = "euclidean")  # Hellinger + Euclidean ~ chord distance, compare_insect_fungi/insect_fungal_coca_analysis.general_workflow.r:151) — it's not Euclidean on raw counts, it's Euclidean on Hellinger-transformed data, which is a specific, deliberate combination:

- Hellinger transform = square-root of relative abundances per sample (decostand(x, "hellinger")). This down-weights dominant taxa and converts the data to proportions before any distance is computed.
- Euclidean distance on Hellinger-transformed data is mathematically equivalent to the Hellinger distance, one of the transformations Legendre & Gallagher (2001) specifically recommend as a Bray-Curtis substitute for methods that require true Euclidean geometry — PCA, RDA, and classical (metric) PCoA via cmdscale(), which is what's used here.

The reason it matters here specifically: Bray-Curtis is a semi-metric, not a proper Euclidean-embeddable distance. Feeding it into cmdscale() (classical/metric PCoA) commonly produces negative eigenvalues and a non-Euclidean embedding that needs a correction step (e.g. Cailliez) to interpret cleanly. Hellinger+Euclidean avoids that — it's guaranteed to embed properly — while still behaving like a community-ecology-appropriate distance (down-weighting dominants, handling double-zeros reasonably), unlike raw Euclidean on counts.

This is a different choice than insect_ords.R, which uses Bray-Curtis directly with metaMDS (non-metric MDS) — Bray-Curtis is the right tool there because NMDS doesn't need a Euclidean-embeddable distance, it just ranks dissimilarities. The CoCA script needed the axis as a numeric covariate in a cmdscale-based PCoA, so it uses the Hellinger/Euclidean combination instead.

If you'd rather it match the Bray-Curtis convention used elsewhere in the project (e.g., substituting a PCoA on Bray-Curtis with a Cailliez correction, or non-metric MDS scores instead of classical PCoA), I can swap that in — but the current choice wasn't arbitrary, it's the standard fix for this exact use case.

### ###############################################################
### Adding a per-taxon test against PCoA2 not just PCoA1
### ###############################################################

Adding insect_PCoA2 as a second tested axis, and adding a direct per-taxon test of association with collection date (independent of the insect-community axes), for both insect and fungal taxa.

While generalizing the per-taxon permutation test to cover PCoA2, fixed a subtle issue in the original permutation scheme: previously insect_PCoA1 and insect_PCoA2 were shuffled together (same permutation index) even when only testing PCoA1's coefficient. Each axis is now permuted independently while the other stays fixed at its observed value -- the correct way to isolate a single term's marginal effect. This is why the PCoA1-no-date hit count changed from the earlier 96/200 to 81/200 here (stricter, more correct null, not a regression).

PCoA1 x date:
- PCoA1 without date: 81 / 200 significant (q<0.10)
- PCoA1 with date: 0 / 200 significant (best q=0.27)

PCoA2 x date -- answering "did anything come up significant for PCoA2?":
- PCoA2 without date: 0 / 200 significant (best q=0.20)
- PCoA2 with date: 0 / 200 significant (best q=0.20)
No, nothing comes up for PCoA2 either way -- it never had a significant fungal association to begin with, consistent with it being a much weaker axis (13.5% variance vs. PCoA1's 31%) not tied to a single strong driver the way PCoA1 is tied to date.

Insect taxa vs. date (direct test, all 52 taxa): 20/52 significant at q<0.10, confirming and formalizing the earlier informal correlation check:
- Later-season: Melanophthalma sp. (t=+9.2), Stenoscelis brevis (+8.8), Corthylus columbianus (+4.1), Dryocoetes autographus (+2.7), Himatium errans (+3.2), Hylobius pales (+2.4)
- Earlier-season: Xyleborinus attenuatus (t=-8.5), Xyloterinus politus (-5.1), Pissodes strobi (-4.7), Cyclorhipidion pelliculosum (-4.5), Xylosandrus germanus (-3.4), Heteroborips seriatus, Pityogenes hopkinsi, Pityophthorus cariniceps, Dendroctonus valens

Fungal taxa vs. date: 200/200 tested significant at q<0.10 -- caveat: those 200 were pre-screened as the taxa with the largest raw correlation with date out of 6,887, so near-universal significance among them is expected by construction, not a surprise finding. It does confirm those specific ASVs have a real (not chance) seasonal signal under the trap-blocked permutation null, with very strong t-statistics (~11-14) -- but "200/200" shouldn't be read as "the whole fungal community is seasonal."

### [for later] We should add a test for how much of the fungal community is seasonal, not just the pre-selected taxa.

Fungal taxa vs. date, with taxonomy -- a very clean pattern. The top 20 hits split almost entirely into two genera pulling in opposite directions:
- Mycosphaerellaceae / Ramularia -- increasing over the season (later): 11 of top 20
- Valsaceae / Cytospora -- decreasing over the season (earlier): 6 of top 20

Ramularia is a foliar plant-pathogen genus (leaf-spot fungi) -- building up over the growing season as leaf area increases is textbook phenology for a foliar pathogen. Cytospora (canker fungus) being higher earlier in the season lines up with the earlier CoCA finding that tied Cytospora to the bark-beetle-dominated (early-season) end of the insect community axis -- this is now a second, independent line of evidence (direct date test, not mediated through the insect axis) for the same seasonal pattern. Also notable in the early-season group: Piptoporus betulinus (birch polypore, wood decay) and Ganoderma sp. (wood decay), both consistent with bark/wood-associated fungi peaking alongside bark beetle activity.

Overall takeaway: no confirmed direct fungal-taxon / insect-community-composition relationship once season is accounted for (either axis). But there is a strong, independent seasonal signal in both communities, with a plausible shared ecological driver (Cytospora canker fungus tracking early-season bark/ambrosia beetle activity, evidenced two separate ways) worth following up with a larger dataset or a study design that better decouples date from insect community turnover.

Output files (data/2024_insect_data/): fungal_insect_association_results.PCoA1_only.csv, fungal_insect_association_results.PCoA1_plus_date.csv, fungal_insect_association_results.PCoA2_only.csv, fungal_insect_association_results.PCoA2_plus_date.csv, insect_taxa_date_association.csv, fungal_taxa_date_association.csv, fungal_taxa_date_association.taxonomy.csv.

### ###############################################################
### Test of the effect of partialling lure in the CoCA+ analyses
### ###############################################################

Lure had a strong effect on the insect community (Lure F=21.9 in insect_PCoA1 ~ site+lure+date), as expected since it's the trap attractant. But there's no biological reason to expect fungi respond to lure independently (all lures contain ethanol), and a companion dataset showed minimal fungal~lure effect there too. Quick check: does including lure in the fungal~insect_PCoA1 model partial out real signal?

Reran the PCoA1 permutation test (same 200 CoCA-ranked candidate taxa) with and without lure in the formula, no-date and date-adjusted:
- PCoA1, with lure, no date: 81/200 significant (q<0.10)
- PCoA1, without lure, no date: 82/200 significant
- PCoA1, with lure + date: 0/200 significant
- PCoA1, without lure + date: 0/200 significant

Dropping lure barely changes anything (81 vs 82 taxa, no-date; 0 vs 0, date-adjusted) -- if lure were partialing out real signal, removing it should have recovered noticeably more hits, and it didn't. Correlation of insect_PCoA1 t-statistics with vs without lure (no-date models): r=0.96, i.e. essentially identical effect sizes either way.

Direct check -- do the candidate fungal taxa respond to lure at all, independent of the insect axes? lure is constant within trap (trap = site x lure), so the within-trap permutation scheme used elsewhere can't test it meaningfully; used a parametric, site-controlled F-test instead (descriptive only, not permutation-corrected): 0 of 200 candidate fungal taxa show even a nominal lure effect (q<0.10).

Conclusion: lure is a genuine, strong driver of the insect community but has no detectable independent effect on the fungal community here, and including it in the models neither suppresses nor inflates the fungal~insect_PCoA1 relationship. The date-driven conclusion above stands on its own, unaffected by whether lure is in the model.

This check is now folded into the main script as two extra PCoA1 test variants (no-lure, no-date/date-adjusted) plus a direct fungal~lure parametric check, rather than a one-off side analysis. Output files added: fungal_insect_association_results.PCoA1_only.no_lure.csv, fungal_insect_association_results.PCoA1_plus_date.no_lure.csv, fungal_taxa_lure_direct_check.csv.

### ###############################################################
### Test of date against whole fungal community
### ###############################################################

The 200/200 fungal taxa found significant vs. date earlier were pre-selected for strong |correlation with date|, so that number can't be read as "how seasonal is the fungal community" -- it's true by construction. Added a new script, compare_insect_fungi/fungal_community_seasonality.r, to answer that question two unbiased ways.

1) Per-taxon date test on every prevalence-filtered fungal taxon (no pre-selection), same trap-blocked permutation scheme as before (fungal_clr ~ site + lure + date, permuting date within trap): 2116 of 6342 taxa (32.6%) significant at q<0.10. That's a much more modest, and more trustworthy, estimate of the seasonal fraction of the fungal community than "200/200" implied -- roughly a third of taxa carry a detectable seasonal signal, not "basically everything."

2) Whole-community view, mirroring insect_ords.R's approach for the insect side: the raw ASV table (sequencing depth varies ~150x across the 69 samples) was repeatedly rarefied to depth 5000 (100 iterations, all samples retained) via multiple_subsamples(), the resulting Bray-Curtis distance matrices averaged via avg_matrix_list(), then PERMANOVA (~ site + site:lure + site:date) and NMDS run on the averaged matrix, same term structure as the insect community PERMANOVA in insect_ords.R.

PERMANOVA (fungal community): site R2=14.7% (F=4.96, p=0.001), site:lure R2=10.1% (F=1.27, p=0.05), site:date R2=22.7% (F=5.72, p=0.001), residual 52.5%. Site:date is the single largest term -- directly parallel to the insect community result (Site R2=9.9%, Site:Lure R2=20.1%, Site:mdy R2=19.7%), except here date/season dominates over lure, whereas lure dominates on the insect side. That's a coherent picture: lure structures which insects get caught, while season structures which fungi are present, and the per-taxon result above (32.6%) is the taxon-level face of the same site:date term.

NMDS (stress ~0.17) plotted with date as a color gradient and site as point shape: figures/fungal_community_NMDS.date_and_site.pdf.

Output files added: data/2024_insect_data/fungal_taxa_date_association.all_taxa.csv, data/2024_insect_data/fungal_community_permanova.csv, figures/fungal_community_NMDS.date_and_site.pdf.

### ###############################################################
### We observed a lure effect in the insect NMDS, but it does not
### appear in the above CoCA+ analyses, which is the lure axis?
### ###############################################################

Open question from the fungal PERMANOVA above: site:lure explains 10.1% of fungal community variance (p=0.05), but no lure-associated taxa turned up among the CoCA hits, which were only ever tested against insect_PCoA1/PCoA2. Hypothesis: the lure signal in the insect community lives on a lower-variance axis (3rd, 4th, ...) that a 2-axis PCoA screen wouldn't catch.

Built a new script, compare_insect_fungi/insect_pcoa_lure_association.r, to check this directly -- same Hellinger+Euclidean PCoA strategy as the CoCA script, but keeping all (not just the top 2) PCoA axes, then testing each axis's association with lure. Lure is a trap-level constant (exactly one trap per site x lure combination, sampled repeatedly over dates), so the within-trap date-permutation used elsewhere can't test it -- instead permuted lure labels among the 3 traps within each site (matches how lure actually varies in this design: across traps within a site, not across samples within a trap).

Result, ranking insect PCoA axes by lure association (F-test, site+date-adjusted, trap-within-site permutation, n=999):

| Axis | % variance | partial R2 (lure) | q-value |
|---|---|---|---|
| PCoA3 | 10.2% | 0.81 | 0.013 |
| PCoA2 | 13.5% | 0.36 | 0.013 |
| PCoA4 | 6.5% | 0.34 | 0.013 |
| PCoA1 | 31.1% | 0.41 | 0.088 |
| PCoA5-10 | -- | -- | n.s. |

Hypothesis confirmed: PCoA3, never tested in the original CoCA per-taxon screen, is by far the strongest lure-associated axis (R2=0.81) -- much stronger than PCoA2, which *was* tested, and stronger than PCoA1's borderline association (q=0.088, likely secondary/incidental given PCoA1 is fundamentally the date axis). PCoA4 is also significant but left out of the follow-up per scope (one axis at a time).

Output: data/2024_insect_data/insect_pcoa_axes.lure_association.csv, figures/insect_pcoa_scree.lure_highlighted.png, figures/insect_pcoa_top_lure_axis.by_lure.png.

### ###############################################################
### This section describes the final CoCA+ analysis which includes
### 3 insect axes in the taxon selection and modeling steps
### ###############################################################

Extended the CoCA workflow to test insect_PCoA3 as a predictor of fungal abundance, following up on the PCoA3/lure finding above. Also bumped the CoCA algorithm itself (Step 3, the coca() call that ranks candidate taxa by loading strength) from n_axes_use=2 to 3, so the top-200 candidate list is now selected using CoCA axes 1-3 rather than 1-2.

insect_PCoA3 was added as a third covariate alongside PCoA1/PCoA2 in every per-taxon model (all three are mutually controlled for, same logic as the original PCoA1/PCoA2 design: the PCoA axes are orthogonal by construction, so this isn't correcting collinearity between axes, it's making each axis's permutation test specific to "does this taxon track axis k given what's already known about the other axes" rather than a less precise univariate test). Ran the same with/without-date x with/without-lure spec matrix already used for PCoA1 -- the without-lure comparison matters more for PCoA3 than it did for PCoA1, since PCoA3 is itself ~81% explained by lure, so the with-lure models are only testing the ~19% of PCoA3 that isn't lure.

Result: PCoA3 comes up completely null. 0 of 200 CoCA-ranked candidate taxa significant (q<0.10) in all four variants -- with lure/no date, with lure/with date, without lure/no date, without lure/with date. Contrast with PCoA1 over the same candidates: 46-60/200 significant without date (46 with lure, 60 without), collapsing to 0/200 once date is added, as found previously.

So despite insect_PCoA3 being the single strongest lure-associated axis in the whole insect ordination, no fungal taxon among the top-200 CoCA-loading candidates tracks it -- whether or not lure is partialled out. Combined with the direct parametric fungal~lure check (also 0/200, both before and after this update), the small lure effect PERMANOVA detects in the fungal community (10.1%, p=0.05) does not appear to be mediated through this specific insect community axis. Caveat: the candidate list is still selected by CoCA loading strength (now over 3 axes), not by a PCoA3- or lure-specific screen, so this rules out "PCoA3 explains the CoCA hits" but doesn't fully rule out a lure-associated fungal signal outside that candidate set -- though the null parametric check makes that less likely.

Output files added: fungal_insect_association_results.PCoA3_only.csv, fungal_insect_association_results.PCoA3_plus_date.csv, fungal_insect_association_results.PCoA3_only.no_lure.csv, fungal_insect_association_results.PCoA3_plus_date.no_lure.csv, figures/fungal_insect_volcano.PCoA3.with_lure.png, figures/fungal_insect_volcano.PCoA3.no_lure.png.

### #################################################################################
### Testing the effect of selecting fungal taxa by direct comparison to date/lure
### #################################################################################

The 3-axis CoCA loading-strength screen (results moved to data/compare_insects_fungi_top3axes/) selects candidate fungal taxa using CoCA axes 1-3 jointly. Set up two alternative single-axis candidate-selection strategies to compare against it: PCoA1 (the date-associated axis) and PCoA3 (the lure-associated axis), each isolating one specific driver instead of the joint CoCA structure.

cocorresp::coca() requires non-negative "community-style" data for both blocks, so the signed insect_PCoA1/PCoA3 axis scores can't be plugged in directly as a single-column CoCA predictor. Candidate selection instead reuses/extends the unbiased, all-taxa per-taxon screens already used elsewhere in this project: fungal_taxa_date_association.all_taxa.csv (already run, fungal_community_seasonality.r) for the PCoA1 strategy, and a new parallel script, fungal_community_lure_association.r, for PCoA3 -- same design (F-test of the lure term, site+date vs. site+lure+date, trap-within-site permutation, same scheme insect_pcoa_lure_association.r used to screen insect PCoA axes for lure), run on all 6887 prevalence-filtered fungal taxa with no pre-selection. Result: 0 of 6887 fungal taxa individually significant for lure (q<0.10) -- consistent with the earlier finding that fungal PERMANOVA detects a real but small, diffuse lure effect (10.1%, p=0.05) not concentrated in specific taxa. The ranking is still usable as a candidate screen even with no individually-significant taxa, same as how the original CoCA loading-strength ranking never required its top-200 to be individually significant either.

Two new scripts (insect_fungal_coca_analysis.PCoA1_date_axis.r, insect_fungal_coca_analysis.PCoA3_lure_axis.r) take the top 200 taxa from each screen and run the same per-taxon permutation-test machinery as the original CoCA script (site + lure? + date? + insect_PCoA1/2/3, permuting only the target axis within trap), restricted to their respective target axis (4 variants: with/without date x with/without lure). Results: data/compare_insects_fungi_PCoA1_date/, data/compare_insects_fungi_PCoA3_lure/; figures in the parallel figures/ subfolders.

PCoA1/date-screen results: 199/200 candidates significant vs. insect_PCoA1 without date (expected -- these taxa were selected FOR date association, and PCoA1 is largely the date axis), dropping to 21/200 once date is added as a covariate (50/200 if lure is also dropped). Contrast with the original 3-axis-CoCA candidate set over the same PCoA1-plus-date test: 0/200. So a candidate list built specifically around date association retains some signal beyond a linear date covariate that the CoCA-selected list does not -- but this is likely at least partly circular (candidates were pre-selected for date association, and insect_PCoA1's collinearity with date is r=0.78, not 1.0, so imperfect collinearity alone could produce partial signal survival) rather than clean evidence of a direct insect-community effect. Should not be over-read as new biology without a stricter check.

PCoA3/lure-screen results: 0/200 candidates significant vs. insect_PCoA3 in all 4 variants -- same null result as the original 3-axis-CoCA candidate set. Even a candidate list specifically built to surface lure-associated fungal taxa shows no detectable individual-taxon relationship with the lure-associated insect axis. One difference: the direct parametric fungal~lure check (site-controlled, uncorrected) finds 10/200 candidates nominally lure-associated in this screen, vs. 0/200 for the CoCA-selected candidates -- consistent with the community-level lure signal being diffuse (spread thinly across many taxa, each below the multiple-testing threshold) rather than concentrated, and more findable by a direct univariate screen than by CoCA loading strength.

Candidate-set overlap comparison (compare_selection_strategies.r, data/compare_selection_strategies/, figures/compare_selection_strategies/candidate_overlap_regions.png) is the more striking result: the 3-axis CoCA list barely overlaps with either single-axis screen (4/200 shared with PCoA1_date, 4/200 shared with PCoA3_lure, Jaccard ~0.01 both times), while the PCoA1_date and PCoA3_lure screens overlap substantially with EACH OTHER (76/200 shared, Jaccard 0.235) despite nominally targeting different drivers (date vs. lure). Only 1 taxon is common to all three lists. Two read on this: (1) joint CoCA loading strength (which reflects taxa best explaining insect-fungal covariation via a multivariate PLS-type decomposition) is picking up on a genuinely different signal than "correlates with one univariate covariate," so the original 3-axis screen and these single-axis screens are not interchangeable; (2) the substantial date/lure screen overlap suggests both univariate per-taxon tests may be partly surfacing generically high-effect-size/high-variance taxa (strong site or trap-driven differences) rather than drivers cleanly specific to date vs. lure -- worth treating the "top-loading" label for any of these three strategies as method-dependent, not a fixed property of a taxon.

Overlap among each strategy's own SIGNIFICANT hits (not just candidate membership) is also low even where both strategies found hits: for PCoA1 without date, top3axes found 46/200 significant and PCoA1_date-screen found 199/200, but only 4 taxa are significant in both (Jaccard 0.017) -- the two methods are largely flagging different specific ASVs as "significant," not just different-sized supersets of the same core signal.

Output files added: data/2024_fungi/fungal_taxa_lure_association.all_taxa.csv; data/compare_insects_fungi_PCoA1_date/ (fungal_insect_association_results.PCoA1_*.csv, lure_partialling_check.summary.csv, fungal_taxa_lure_direct_check.csv); data/compare_insects_fungi_PCoA3_lure/ (fungal_insect_association_results.PCoA3_*.csv, lure_partialling_check.summary.csv, fungal_taxa_lure_direct_check.csv); data/compare_selection_strategies/ (candidate_taxa_membership.csv, candidate_set_overlap_summary.csv, significant_hits_overlap.csv); figures/compare_insects_fungi_PCoA1_date/fungal_insect_volcano.PCoA1.png; figures/compare_insects_fungi_PCoA3_lure/fungal_insect_volcano.PCoA3.{with,no}_lure.png; figures/compare_selection_strategies/candidate_overlap_regions.png.

### #################################################################################
### Procrustes (whole community)

All the analysis above (CoCA, per-taxon date/lure tests) asks a taxon-level question. Added compare_insect_fungi/insect_fungal_procrustes.r to ask the whole-community version directly with Procrustes analysis: do samples that are close together in insect-ordination space also tend to be close together in fungal-ordination space, taken as full communities rather than one taxon at a time?

Same 69-sample matched dataset (Curculionidae + Latridiidae insect trap catch x fungal ITS2 ASVs) used throughout this comparison. Insect NMDS: raw Bray-Curtis, no autotransform (insect_ords.R's setting for this table, stress ~0.14). Fungal NMDS: rarefaction-averaged Bray-Curtis (same depth-5000 x 100-iteration rarefy-and-average approach as fungal_community_seasonality.r, stress ~0.17) -- deliberately built independently of the insect data, so the comparison isn't circular.

Result: Procrustes correlation = 0.56, p = 0.001 (999 permutations, symmetric rotation) -- a real, statistically significant whole-community concordance. This is consistent with, and largely expected given, the two lines of evidence already established: both community PERMANOVAs are dominated by a site:date term (insect Site:mdy R2=19.7%; fungal Site:date R2=22.7%), so a shared seasonal/site structure alone should pull the two ordinations into partial alignment even without any direct taxon-to-taxon relationship. The Procrustes result confirms that shared structure is strong enough to show up as significant whole-community concordance -- it does NOT on its own imply a direct biotic interaction, and shouldn't be read as contradicting the per-taxon finding that no individual fungal taxon's association with the insect community axes survives adjusting for date (PCoA1: 81/200 -> 0/200 once date is added; PCoA3: null throughout). A moderate, imperfect correlation (0.56, not close to 1) is itself consistent with that picture -- real shared structure, but far from a 1:1 mapping.

Worst-matched (highest-residual) samples span several sites and lures rather than clustering in one condition (top of the list: NH.34, Durham, Alpha-pinene_EtOH, 2024-05-15; NH.99 and NH.46, Manchester Cedar Swamp; NH.113, Pease) -- no obvious single site/lure/date driving the discordance, though this hasn't been formally tested.

Output files: data/2024_insect_data/insect_fungal_procrustes_summary.csv (correlation, sum-of-squares, p-value), data/2024_insect_data/insect_fungal_procrustes_residuals.csv (per-sample residuals with metadata); figures/insect_fungal_procrustes.pdf (standard Procrustes vector diagram + residuals-by-sample plot), figures/insect_fungal_procrustes.NMDS_overlay.png (ggplot overlay of both NMDS configurations post-rotation, colored by community and shaped by site).

### #################################################################################
### Non-fungal ASV filtering: the fungal community table was never restricted
### to Kingdom == Fungi
### #################################################################################

Discovered that ASV_tab.csv (and therefore every core script's fungal community table) had never been filtered to confirmed fungal taxa -- ASVs_taxonomy.tsv includes ASVs identified as Rhizaria, Viridiplantae (plant), and Metazoa (animal), plus a large unidentified-Kingdom (NA) bucket, all left in alongside real fungi.

Quick check on the 69-sample matched dataset's prevalence-filtered taxa (>=5 samples, the filter shared by all 5 core scripts): of 6,887 retained ASVs, 545 (7.9%) were not confirmed Fungi -- 470 completely unidentified (NA Kingdom), 48 Viridiplantae, 20 Rhizaria, 7 Metazoa. By read count these 545 ASVs were a much smaller share, 2.09% of total reads (210,303 / 10,069,212) -- mostly low-abundance ASVs, not dominating the community signal. Spot-checking the top-200 candidates actually used in the date-association permutation test found 8 with no Kingdom assignment at all, two of which (ASV_1079, ASV_4243) were among the top "significant" seasonal hits at q=0.001 -- i.e. taxonomically unidentified sequences were sitting inside a headline result, not just padding out the tail.

Fix: added a `Kingdom == "k__Fungi"` filter (via ASVs_taxonomy.tsv) to all 5 core scripts immediately after loading and sample-matching the ASV table, before the existing prevalence filter -- insect_fungal_coca_analysis.general_workflow.r, .PCoA1_date_axis.r, .PCoA3_lure_axis.r, fungal_community_seasonality.r, fungal_community_lure_association.r. Consistent across all 5 (same 69-sample matched dataset): 15,422 of 17,139 ASVs confirmed Kingdom == k__Fungi (1,717 dropped), of which 6,342 pass the prevalence filter (vs. 6,887 before).

Before vs. after comparison, rerunning all 5 scripts:

| Result | Before | After |
|---|---|---|
| Whole-community PERMANOVA -- site R2 | 14.8%, p=0.001 | 15.0%, p=0.001 |
| -- lure R2 | 2.81%, p=0.129 | 2.80%, p=0.143 |
| -- date R2 | 16.3%, p=0.001 | 16.5%, p=0.001 |
| Date-association, unbiased screen (all prevalence-filtered taxa) | 2,248/6,887 sig. (32.6%) | 2,116/6,342 sig. (33.4%) |
| Lure-association, unbiased screen | 0/6,887 sig. | 0/6,342 sig. |
| insect_PCoA1 only (top-200 CoCA candidates) | 199/200 sig. | 199/200 sig. |
| insect_PCoA1 + date | 21/200 sig. | 27/200 sig. |
| insect_PCoA3 (with/without date, with/without lure) | 0/200 sig. | 0/200 sig. |

Every prior conclusion holds -- nothing flips significance direction, effect sizes shift by low single-digit percentage points at most. Concretely: ASV_1079 and ASV_4243, the two unidentified-Kingdom ASVs flagged above, are correctly excluded from the date-association results now. One knock-on effect worth flagging: the whole-community rarefied PERMANOVA/NMDS lost a sample (69 -> 68) -- one sample's total read count fell below the common rarefaction depth once its non-fungal reads were stripped out, since a chunk of what had been counted toward its sequencing depth wasn't fungal to begin with.

While rerunning, also fixed two pre-existing, unrelated bugs surfaced by the rerun rather than caused by it:
- fungal_community_seasonality.r was writing fungal_taxa_date_association.all_taxa.csv to data/2024_insect_data/, while every script that reads it (insect_fungal_coca_analysis.PCoA1_date_axis.r, compare_selection_strategies.r, presentation_items/fig3_volcano_plots.R) expected it at data/2024_fungi/ -- a stale manually-relocated copy at the expected path was silently going stale on every rerun. Fixed the write path to match.
- insect_fungal_coca_analysis.general_workflow.r's per-axis-model comparison outputs (the 10 fungal_insect_association_results.*.csv files, fungal_taxa_lure_direct_check.csv, fungal_taxa_date_association.csv, and the CoCA biplot / volcano figures) were writing to data/2024_insect_data/ and figures/ directly, unlike the two single-axis scripts which each write to their own data/compare_insects_fungi_PCoA1_date(or PCoA3_lure)/ and figures/compare_insects_fungi_PCoA1_date(or PCoA3_lure)/ subfolders. Redirected these to data/compare_insects_fungi_top3axes/ and figures/compare_insects_fungi_top3axes/ for structural parity across all three model variants. insect_taxa_date_association.csv stays in data/2024_insect_data/ since it isn't axis-model-specific.

Output files updated in place (same filenames, taxonomy-filtered content): data/2024_fungi/fungal_community_permanova.csv, data/2024_fungi/fungal_taxa_date_association.all_taxa.csv, data/2024_fungi/fungal_taxa_lure_association.all_taxa.csv, data/compare_insects_fungi_top3axes/*, data/compare_insects_fungi_PCoA1_date/*, data/compare_insects_fungi_PCoA3_lure/*, figures/compare_insects_fungi_top3axes/*, figures/fungal_community_NMDS.date_and_site.pdf.

### #################################################################################
### Taxonomic breakdown of the significant hits behind fig3's three volcano panels
### #################################################################################

fig3 shows every tested taxon as one point; added presentation_items/fig4_association_taxonomic_breakdown.R to ask a follow-up question -- of just the significant (q<0.10) hits, which taxonomic groups do they belong to, and which direction (later-season / earlier-season, or PCoA1+/PCoA1-) do they point? Diverging bar charts, one per fig3 panel: insect~date (panel a), fungal~date unbiased screen (panel b), fungal~insect_PCoA1 no-lure/no-date CoCA-candidate screen (panel c, matches fig3 panel c exactly).

Taxonomic level was chosen per taxon group to balance a manageable number of bars against ecological resolution: insects use Subfamily (Family collapses to 2 groups since Curculionidae is ~98% of the catch; Subfamily gives 4 groups and is the level already used above to distinguish Scolytinae from Corticariinae/Cossoninae). Fungi use Family, capped at the top 15 by hit count + "Other" (2116 date hits span 259 families and 77 PCoA1 hits span 38 -- Family is too fine to show unfolded, but coarser levels like Order/Class would lose the Valsaceae/Mycosphaerellaceae-level resolution this project already relies on).

Result is a clean confirmation of the patterns described above, now quantified across ALL significant hits rather than just the top 20: Scolytinae is 11/13 earlier-season insect hits (vs. Cossoninae/Corticariinae/Molytinae split evenly or later); Valsaceae is 72/73 earlier-season fungal date hits and 7/7 tracking the PCoA1- (bark-beetle-dominated) end; Mycosphaerellaceae is 165/170 later-season fungal date hits. A large "Unclassified" (family-level) bucket shows up top-ranked in both fungal panels (406/2116 date hits, 13/77 PCoA1 hits) -- these ASVs are confirmed Kingdom==Fungi (post the filter fix above) but unresolved below Order/Class, not the excluded non-fungal ASVs from that earlier fix.

Output: figures/presentation_items/fig4_association_taxonomic_breakdown.png/.pdf.

### #################################################################################
### fig4's fungal panels moved from Family to Order level; new fig5 supplement
### #################################################################################

Two follow-up exploratory scripts (compare_insect_fungi/fungi_date_order_breakdown.R, fungi_insect_order_breakdown.R) looked at the same significant hits at Order level instead of Family, to see whether whole orders track cleanly in one direction. They do, more cleanly than Family-level suggested: e.g. Mycosphaerellales (169 later-season vs. 12 earlier), Tremellales (113 vs. 4), and Dothideales (97 vs. 8) are each almost entirely one-directional for date; Diaporthales, Atractiellales, Helotiales, and Saccharomycetales are each almost entirely on the PCoA1-/bark-beetle-associated side. Family-level obscured this because families like Phleogenaceae/Atheliaceae (Atractiellales) or Rhytismataceae/Dermateaceae/Tympanidaceae (Helotiales) split the same order's signal across several small bars.

Based on this, fig4_association_taxonomic_breakdown.R panels b/c were changed from Family to Order (panel b: top 15 of 103 orders + "Other"; panel c: all 31 orders present, no folding needed -- see the script's updated header for the full rationale). A new supplemental figure, fig5_family_breakdown_by_order.R, now carries the family-level detail that fig4 used to show: family-level diverging bars, grouped into facets by their subsuming order (facet strip = order name) and ordered to match fig4's order ranking, so it reads as a direct drill-down of fig4 panels b/c. Panel b covers the same 15 major orders as fig4b, with each order's own family list folded to its top 6 + "Other" where needed (Pleosporales alone spans 24 families) -- 88 family-level bars total. Panel c covers all 31 orders from fig4c unfolded (no order there spans more than 3 families) -- 41 bars.

Output: figures/presentation_items/fig4_association_taxonomic_breakdown.png/.pdf (revised), fig5_family_breakdown_by_order.png/.pdf (new); data/presentation_items/fig4_taxonomic_breakdown.fungal_date.by_order.csv, fig4_taxonomic_breakdown.fungal_insectPCoA1.by_order.csv, fig5_family_breakdown_by_order.fungal_date.csv, fig5_family_breakdown_by_order.fungal_insectPCoA1.csv; exploratory figures figures/fungi_date_order_breakdown.png, figures/fungi_insect_order_breakdown.png.

Also added the same facet_grid technique to fig4 itself: panels b/c's order-level bars are now grouped into Class facets (one level up from fig5's order facets around family bars). This makes the class-level pattern immediately visible: in panel b, Dothideomycetes (Mycosphaerellales, Dothideales, Capnodiales -- all later-season) and Sordariomycetes (Diaporthales, Xylariales, Hypocreales -- all earlier-season) each cleanly split by direction at the class level, not just the order level. In panel c, Sordariomycetes (Diaporthales, Ophiostomatales, Xylariales) and Saccharomycetes/Atractiellomycetes are entirely on the PCoA1-/bark-beetle side.

#### Consider adding a chi-sqaured test for the insect-associated classes/orders wherein we compare insect associated/not versus date associated/not (could make this multi-level by comapring date-assocaited/not to whole community

### #################################################################################
### Responsive taxa carry disproportionate sequence abundance
### #################################################################################

The 2133 date- or insect_PCoA1-associated ("responsive") fungal ASVs (q<0.10) are only ~14% of Kingdom==Fungi ASVs in the 69-sample matched dataset, but make up 60-70% of average per-sample relative sequence abundance at every site (Durham 64%, Pease Airport 61%, Manchester Cedar Swamp 70%, Manchester Airport 68%) -- responsiveness concentrates in the more abundant taxa rather than being spread evenly across the community, and that concentration is fairly consistent across sites (compare_insect_fungi/responsive_taxa_abundance_by_site.R, figures/responsive_taxa_abundance_by_site.png).

Splitting that combined "responsive" set into its two component screens (presentation_items/fig6_responsive_taxa_abundance_and_counts_by_site.R, figures/presentation_items/fig6_responsive_taxa_abundance_and_counts_by_site.png) shows the disproportion is almost entirely a date effect: date-associated ASVs are ~20-25% of taxa present per site but ~60-70% of relative abundance, the same small-count/large-abundance pattern as above. Insect-associated taxa alone (only 77 total, 60 of which also overlap with the date-associated set) are a small share on both axes at every site -- the dramatic split in the combined figure is carried almost entirely by the date signal, not the insect signal.

### #################################################################################
### Trophic mode (FungalTraits) breakdown of the significant hits -- fig7/fig8
### #################################################################################

Same question as fig4/fig5 (taxonomic breakdown of the date- and insect_PCoA1-significant hits), but by ecological trait instead of taxonomy: ASV genus joined to FungalTraits (Polme et al. 2020) primary_lifestyle, same join fig1 panel d uses. Two companion figures, same underlying data: fig7 (presentation_items/fig7_association_trait_breakdown.R) facets by growth_form with lifestyle as the bar label; fig8 (presentation_items/fig8_association_trait_breakdown.R) facets by taxonomic Class instead, with growth_form dropped from the label since it adds little once Class is already shown. "Unclassified genus" and "No FungalTraits match" ASVs are dropped from both (they carry no trait information to plot), which removes ~36% of the 2116 date hits and ~27% of the 77 PCoA1 hits from the FungalTraits-based figures -- a substantial blind spot worth remembering when comparing these panel totals (1491/56 matched taxa) against fig4/fig5's full counts.

The headline pattern, visible in both figures: plant pathogen and wood saprotroph are by far the two largest trait categories in both directions of the date signal, and PCoA1 hits skew heavily toward the earlier-season/bark-beetle end across nearly every trait category (plant pathogen 19 earlier vs. 1 later; wood saprotroph 10 vs. 2), consistent with the CoCA/Cytospora/bark-beetle story already established above.

fig8's Class-faceted cross-tab adds a layer the growth-form view (fig7) and fig1 panel d's site-pooled bars can't show: "plant pathogen" is not one uniform seasonal signal, it's at least two taxonomically distinct ones pointing in opposite directions. Dothideomycetes plant pathogens (dominated by Mycosphaerellaceae/Ramularia-type leaf-spot fungi, per the Order/Family analysis above) skew later-season (247 later vs. 205 earlier), consistent with foliar disease building up as leaf area increases over the growing season. Sordariomycetes plant pathogens (dominated by Valsaceae/Cytospora canker fungi) skew strongly earlier-season (123 earlier vs. 34 later), consistent with the canker-fungus/bark-beetle association already established via the CoCA and direct-date analyses above. Pooling across Class (as fig7, fig1 panel d, and any Family/Order-blind trait summary would) mostly cancels these two signals out and would suggest "plant pathogen" is only mildly seasonal, when in fact it contains two strong, opposite, taxonomically separable trends -- Class is doing real explanatory work here, not just adding bars.

A second finding came from checking FungalTraits' Secondary_lifestyle field: taxa whose primary lifestyle is litter or wood saprotroph but whose secondary lifestyle is also plant_pathogen (133 of 1491 matched date-hit taxa, ~9%) are overwhelmingly later-season, and disproportionately so relative to the same primary lifestyle's un-tagged taxa -- e.g. in Dothideomycetes, plain "litter saprotroph" is a small, roughly balanced bar (5 earlier/5 later) once the secondary-plant-pathogen-tagged taxa are split out, while "litter saprotroph (+plant pathogen)" is 1 earlier/70 later. Digging into the single largest contributor to that tagged bucket (litter saprotroph : thallus photosynthetic, 26 taxa, 26/26 later-season) found it collapses almost entirely to one genus, Sphaerulina (18/26 ASVs resolved to species are S. pelargonii) -- a reminder that a clean-looking trait-level bar can be carried by one prolific genus detected as many ASV variants, not a taxonomically diverse guild.

Output: figures/presentation_items/fig7_association_trait_breakdown.png/.pdf, fig8_association_trait_breakdown.png/.pdf; data/presentation_items/fig7_trait_breakdown.fungal_date.csv, fig7_trait_breakdown.fungal_insectPCoA1.csv, fig8_trait_breakdown.fungal_date.csv, fig8_trait_breakdown.fungal_insectPCoA1.csv.

### #################################################################################
### All-against-all individual fungal taxon x individual insect taxon screen --
### negative result, does not support broad specific-symbiosis signal
### #################################################################################

Every analysis above tests fungal taxa against OVERALL insect community variance (insect_PCoA1-3). That's the right test for community-level covariation, but it could miss a fungal taxon that forms a tight association with ONE insect species whose signal doesn't dominate any top PCoA axis. compare_insect_fungi/insect_fungal_pairwise_taxon_association.r instead tests every prevalence-filtered fungal taxon against every prevalence-filtered insect taxon individually: 6,342 fungal ASVs x 39 insect taxa (Curculionidae + Latridiidae, >=5-sample prevalence filter, same thresholds used throughout) = 247,338 pairs, one model per pair (fungal_taxon_clr ~ site + lure + date + insect_taxon), trap-blocked 999-permutation test -- first pass skips the with/without lure/date partialling suite used for the PCoA-axis tests (both included together throughout) per scope decision going in.

Computational note: a direct port of the existing per-pair lm()+apply()-over-permutations pattern (used for the PCoA-axis tests and the all-6,342-taxa date screen) would have taken multiple hours on this much larger grid. Since every fungal taxon shares the same design matrix for a given insect-taxon predictor (and a given permutation of it), all 6,342 fungal responses can be fit at once via vectorized OLS matrix algebra -- (t(X) %*% X)^-1 %*% t(X) %*% Y with Y as the full multi-column response matrix -- instead of one lm() call per fungal taxon. Verified to reproduce lm()'s t-statistics exactly; cut the full rigorous 999-permutation x 247,338-pair grid to ~3 minutes single-core, so no subsampling or pre-screening of fungal candidates was needed.

Headline result: 1,918 of 247,338 pairs significant at q<0.10 (FDR corrected WITHIN each insect taxon's 6,342 tests, the same convention used for the PCoA-axis tests) -- but 0 survive a global correction across the whole grid. More strikingly, only 3 of the 39 insect taxa have ANY significant fungal partner, and all three are early-season Scolytinae (bark/ambrosia beetles): Pityogenes hopkinsi (1,593 hits, 25% of every fungal taxon tested), Xyleborinus attenuatus (254), Hylastes opacus (71). The other 36 insect taxa: zero.

That concentration -- one taxon accounting for 83% of all hits -- prompted several checks before trusting it:
- Not a single-sample leverage artifact: P. hopkinsi's raw counts have one outlier sample (15 individuals vs. a max of 4 elsewhere among 19 nonzero samples), but removing that sample from the model strengthens rather than weakens the association.
- Only partially explained by the fungal community's own direct date-association signal (fungal_taxa_date_association.all_taxa.csv, 2,116/6,342 taxa): significantly enriched (Fisher OR=1.4, p=2e-8) but this overlap only accounts for 39% of P. hopkinsi's hits.
- Taxonomically incoherent for a specific symbiosis: alongside the expected Cytospora/Valsaceae bark-beetle fungi (30 hits) are lichens (Physcia), foliar pathogens (Ramularia), wood-decay polypores (Bjerkandera, Xylodon), and epiphytic yeasts (Aureobasidium, Vishniacozyma) -- guilds with no shared biology.

Follow-up (compare_insect_fungi/insect_fungal_pairwise_taxon_association.residualized.r): does the signal survive a stricter seasonal control on the insect-predictor side? With LINEAR date in both the insect predictor and the fungal model, the fungal-side coefficient on insect_taxon is mathematically identical whether insect_taxon enters raw or pre-residualized against site+date (Frisch-Waugh-Lovell) -- so that version would just reproduce the original result. Since this dataset has only 6 distinct collection dates, residualizing each insect taxon against site + factor(date) (a 6-level factor, not a line) removes ANY date effect the data can resolve, linear or not -- a genuinely stricter test, not a redundant one.

Result: 2 of the 3 flagged taxa collapse entirely to zero -- Xyleborinus attenuatus and Hylastes opacus were pure seasonal/site artifacts, not real associations. Pityogenes hopkinsi drops from 1,593 to 717 hits (45% survive, 42% are literally the same fungal taxa as before) -- so roughly half of its original signal was also seasonal, but a substantial chunk survives even this strictest available deseasonalizing. That surviving chunk, though, still doesn't look like 717 independent symbioses: still concentrated in exactly one insect taxon (0/38 others), still 0 significant under global correction, still taxonomically broad (yeasts, wood-decay fungi, plant pathogens, a shrunk-but-present Cytospora/Valsaceae core of ~30 hits), and a PCA of just these 717 taxa's CLR profiles shows PC1 alone explains 31% of their variance (24 PCs needed for 80%) -- i.e. there's a real dominant shared axis behind a lot of this, meaning it behaves more like a handful of correlated community-level gradients riding along with P. hopkinsi than 717 independent point-to-point relationships. BH-FDR assumes roughly independent tests; a shared latent axis like this inflates the apparent hit count past what the q-values suggest at face value.

**Conclusion:** this screen does not support a broad layer of specific insect-fungus symbioses across the insect community sampled here. No insect taxon shows a clean, robust, taxon-specific fungal partner set that survives scrutiny. The one partial exception, Pityogenes hopkinsi, retains a real but modest association after the strictest deseasonalizing available, but it reads as community/environment-mediated (part of, and substantially overlapping with, the seasonal bark-beetle/Cytospora pattern already established via insect_PCoA1 above) rather than a validated point symbiosis. Kept here for the record since this was a real analysis path explored, not because it's expected to appear in the primary paper.

Output: data/compare_insects_fungi_pairwise_taxa/fungal_insect_pairwise_full_grid.csv (all 247,338 pairs), fungal_insect_pairwise_significant_hits.csv, fungal_partners_per_insect_taxon.csv, insect_partners_per_fungal_taxon.csv, fungal_family_tally.csv, fungal_insect_pairwise_full_grid.date_site_residualized.csv, fungal_insect_pairwise_significant_hits.date_site_residualized.csv, original_vs_residualized_comparison.csv; figures/compare_insects_fungi_pairwise_taxa/fungal_partners_per_insect_taxon.png, family_by_insect_taxon_heatmap.png, fungal_specialist_generalist_distribution.png, original_vs_residualized_comparison.png.

### #################################################################################
### Alpha diversity comparison: insect vs. fungal (Shannon, Simpson dominance, richness)
### #################################################################################

Every prior comparison asks a beta-diversity question (composition). Added compare_insect_fungi/insect_fungal_alpha_diversity.r to ask the alpha-diversity version: do insect and fungal per-sample diversity covary, and does each community's own diversity respond to site/lure/date the way its composition does?

Community tables match the ordinations exactly: insect uses the raw Curculionidae+Latridiidae trap-catch table (same as insect_ords.R/fig2, no rarefaction); fungal uses the same rarefy-to-5000-x-100-iterations approach as fungal_community_seasonality.r/fig2, but averaging per-sample alpha diversity across iterations instead of averaging Bray-Curtis distance matrices (1 sample, NH.54, dropped for falling below depth 5000, same as the beta-diversity version). Shannon and richness (observed taxa, no extrapolation) via vegan; Simpson **dominance** (D = sum(p_i^2), not vegan's 1-D "simpson" diversity index) computed directly since dominance was what was requested.

Two tests, per user request:

1. **Per-sample insect-vs-fungal correlation** (Spearman, alpha-diversity analog of the whole-community Procrustes check, unconstrained permutation/correlation like that check): all three metrics point the same weakly-positive direction but only two clear q<0.10-equivalent hits -- Simpson dominance (rho=0.22, p=0.068) and richness (rho=0.21, p=0.080) are marginal; Shannon (rho=0.19, p=0.127) is not significant. Consistent with, but much weaker than, the whole-community Procrustes result (r=0.56) -- alpha diversity is a coarser summary than full composition, so a weaker cross-community signal here is expected, not a contradiction.

2. **Site/lure/date effects within each community** (same permutation logic used throughout this project -- date permuted within trap, lure permuted among traps within site; site tested the same way as lure, one level up, by permuting the site label freely among all traps since site, like lure, is trap-constant): mirrors the compositional PERMANOVA story closely.
   - **Insect**: lure dominates all three metrics (Shannon F=32.6 R2=0.51 p=0.005; Simpson dominance F=21.9 R2=0.42 p=0.003; richness F=37.5 R2=0.55 p=0.005) -- the Ips lure drives markedly lower diversity/richness and higher dominance than Alpha-pinene_EtOH or Ethanol, visible directly in the boxplots. Richness also shows a real site effect (F=3.73, R2=0.17, p=0.007); date is only marginal (p=0.025-0.06).
   - **Fungal**: date dominates all three metrics (Shannon t=-4.02 R2=0.21 p=0.001; Simpson dominance t=2.48 R2=0.09 p=0.02; richness t=-3.25 R2=0.15 p=0.004 -- diversity/richness fall and dominance rises later in the season), lure has no detectable effect on any metric (p=0.66-0.96), and site is a modest but real hit for dominance and richness (p=0.02-0.05) but not Shannon (p=0.11).

This is a clean alpha-diversity restatement of the whole-community pattern already established: **lure structures the insect community (both composition and diversity), season structures the fungal community (both composition and diversity)**, with lure having essentially zero measurable effect on fungal diversity, matching the null lure-taxon results from the per-taxon and PERMANOVA work above.

Output: data/compare_insects_fungi_alpha_diversity/{insect,fungal}_alpha_diversity.csv (per-sample values), alpha_diversity_cross_community_correlation.csv, alpha_diversity_covariate_tests.csv; figures/compare_insects_fungi_alpha_diversity/alpha_diversity_by_community.png, alpha_diversity_by_site_lure.png, alpha_diversity_by_date.png, alpha_diversity_cross_community_correlation.png.

### #################################################################################
### Alpha diversity variant: full insect table (all families, not just Curculionidae+Latridiidae)
### #################################################################################

The alpha-diversity script above restricts insect data to Curculionidae+Latridiidae, matching every other insect/fungal comparison in this project. Added a variant, compare_insect_fungi/insect_fungal_alpha_diversity.full_insect_table.r, that reruns the identical analysis (same fungal side, same permutation schemes, same three metrics) but with NO family restriction on the insect table -- all 426 taxa across 98 families in insect_species_tab.csv, singleton taxa and resulting empty samples dropped the same way (277 taxa retained, same 69-sample matched dataset as the family-restricted version). Kept as a separate script/output set alongside (not replacing) the Curculionidae+Latridiidae version for reference.

**Result changes substantially and in a biologically informative way once the full insect community is used:**

- **Insect diversity now tracks season strongly, not just lure.** Date becomes the dominant term for all three insect metrics (Shannon t=6.19 R2=0.38 p=0.001; Simpson dominance t=-3.10 R2=0.13 p=0.005; richness t=4.46 R2=0.24 p=0.001) -- and the direction is the OPPOSITE of the Curculionidae+Latridiidae-only result: insect diversity/richness *rises* over the season here (full community), whereas the family-restricted analysis showed a (weaker, marginal) decline. Lure remains significant but much weaker than before (Shannon R2=0.20 p=0.013 vs. 0.51 p=0.005 in the family-restricted version; richness R2=0.32 vs. 0.55); site drops to marginal/non-significant everywhere (p=0.07-0.21, vs. a real richness hit at p=0.007 before).
- **Fungal side is unchanged** (same table, same results as the reference analysis): date dominates all three fungal metrics, lure has no effect, site is a modest hit for dominance/richness only.
- **Cross-community correlation flips sign and becomes significant for Shannon**: insect vs. fungal Shannon diversity rho=-0.26, p=0.030 (vs. +0.19, p=0.13 in the family-restricted version) -- insect Shannon rises across the season while fungal Shannon falls (visible directly in the by-date figure), so the two communities' diversity are seasonally anti-correlated once the full insect catch is used. Simpson dominance (rho=-0.12, p=0.33) and richness (rho=-0.12, p=0.34) are weak/non-significant in the same direction.

Interpretation: restricting to Curculionidae+Latridiidae (done everywhere else in this project because those are the taxonomically well-resolved, ecologically targeted groups) evidently also filters out most of the insect community's own seasonal diversity signal -- the two bark/fungus-associated families examined elsewhere apparently have a flatter or even opposite seasonal diversity trend than the broader trap catch. The full-table result (insects becoming more diverse as fungal diversity declines through the season) is a cleaner, more conventional phenological story (community turnover/diversification over the season) but is not directly comparable to the beta-diversity/CoCA/PERMANOVA results elsewhere in this project, which are all built on the Curculionidae+Latridiidae-only insect table. Treat this as a complementary alpha-diversity-only check, not a replacement for the family-restricted analyses used throughout the rest of the workflow.

Output: data/compare_insects_fungi_alpha_diversity_full_insect_table/{insect,fungal}_alpha_diversity.csv, alpha_diversity_cross_community_correlation.csv, alpha_diversity_covariate_tests.csv; figures/compare_insects_fungi_alpha_diversity_full_insect_table/alpha_diversity_by_community.png, alpha_diversity_by_site_lure.png, alpha_diversity_by_date.png, alpha_diversity_cross_community_correlation.png.

### #################################################################################
### Insect NMDS variant: individual-taxon prevalence filter instead of family filter
### #################################################################################

insect_ords.R restricted the insect trap-catch table to Curculionidae + Latridiidae because the full, unfiltered table gave a high-stress (>0.23), non-repeatable NMDS -- too many taxa caught in only 1-2 samples. Added a new, clean script (insect_ords.R is left as-is, per user request, since it's exploratory/messy) -- insect_exploratory/insect_ords.prevalence_filter.R -- that asks whether an individual-taxon prevalence filter (present in >=5 samples, the SAME threshold already used for fungal ASVs throughout compare_insect_fungi/) fixes the same problem without discarding taxa outside those two families.

Result: yes, cleanly. 153 of 425 named taxa (72 raw samples, one all-zero/no-ID row dropped as in insect_ords.R) pass the >=5-sample filter; 71 samples survive after dropping resulting empty rows. Family composition of the retained taxa is much broader than the two-family table -- Curculionidae still the largest single group (37 taxa) but Elateridae (15), Staphylinidae (12), Cerambycidae (6), Tenebrionidae (6), Nitidulidae (5), and others now included. Confirmed the same methodological choice insect_ords.R made still holds here: raw Bray-Curtis with autotransform=FALSE gives stress=0.162 (repeatable, best solution hit 5/20 tries) vs. stress=0.233 (barely repeatable, 1/21) under vegan's default sqrt+Wisconsin autotransform.

PERMANOVA (same site-blocked/trap-stratified permutation scheme as insect_ords.R's date/site test): site R2=9.5% p=0.001, lure R2=11.4% p=0.001, date R2=12.7% p=0.001, residual 66.6% -- all three terms significant, same three-way structure as the family-filtered table (there: site 9.7%, lure 14.2%, date 14.2%, residual 62.3%), just with somewhat smaller lure/date effect sizes and more residual variance now that more (noisier, less lure/season-responsive) taxa are included. The NMDS plot shows a visually clean date gradient along NMDS1, matching that PERMANOVA result.

Output: data/2024_insect_data/insect_prevalence_filter.family_composition.csv, insect_prevalence_filter.NMDS.RDS, insect_prevalence_filter.permanova.csv; figures/insect_exploratory/insect_ords.prevalence_filter.NMDS.png.

### #################################################################################
### CoCA general workflow rerun with the prevalence-filtered insect table
### #################################################################################

Added compare_insect_fungi/insect_fungal_coca_analysis.prevalence_filtered_insect.r -- the full CoCA -> per-taxon insect_PCoA1-3 association -> date-screen workflow from insect_fungal_coca_analysis.general_workflow.r, rerun on the >=5-sample individual-taxon prevalence-filtered insect table (153 taxa, all families -- insect_exploratory/insect_ords.prevalence_filter.R) instead of the Curculionidae+Latridiidae family restriction. Same 69-sample matched dataset (one fewer than the 71-sample prevalence-filtered insect table alone, since 2 more samples lack fungal data), same fungal filtering/CLR transform, same n_axes_use=3, same per-taxon permutation scheme and full with/without-date x with/without-lure test matrix. Output goes to its own directory (data/compare_insects_fungi_top3axes_prevalence_filtered_insect/, figures/...) so the original family-filtered results are untouched.

**One numerical snag, handled non-fatally:** the LOO cross-validation diagnostic (crossval(), used only for informal inspection -- n_axes_use is hardcoded to 3 either way, matching the general workflow script) hit a LAPACK SVD convergence failure (`La.svd` error code 1) on one particular leave-one-out fold with this broader insect table. Wrapped that call in tryCatch to skip gracefully; it doesn't affect anything downstream since n_axes_use was never derived from it.

**Results are broadly similar in shape to the family-filtered version, with a stronger PCoA1 signal and a new PCoA2 signal:**

- **insect_PCoA1 (no date, with lure): 62/200 significant** (vs. 46/200 for the original Curculionidae+Latridiidae 3-axis CoCA screen over the same candidate-count) -- and **without lure: 74/200** (vs. 60/200 originally). Dropping to **9/200 with date** (17/200 without lure) -- still collapses sharply once date is added, same story as before, just retaining slightly more (9 vs. 0) under the strictest with-lure+date variant. Top hit is again ASV_1905 (*Cytospora prunicola*, Valsaceae) at t=-7.65 (no date) / t=-3.35 (with date) -- the same taxon flagged as the strongest, most date-independent signal throughout this project's prior work, now recovered independently under a completely different insect-table construction.
- **insect_PCoA2 now shows a real signal for the first time**: 6/200 (with lure) and 7/200 (with lure+date) significant -- previously null under the family-filtered table (0/200 throughout iterative_analysis_updates.md's PCoA2 checks). Top hits (ASV_13242, ASV_12553, ASV_616, ASV_3449) are consistent across the no-date/with-date variants (same taxa, same sign), suggesting this isn't a date-collinearity artifact. Worth a taxonomic follow-up look, since this is new relative to everything documented above.
- **insect_PCoA3 remains completely null** (0/200 in all four with/without-lure x with/without-date variants), consistent with the family-filtered table's result -- PCoA3 still doesn't pick up any fungal signal even with the broader insect community behind it.
- **Lure-partialling check**: dropping lure barely changes the PCoA1 hit count (62->74 no-date, 9->17 with date; t-statistic correlation with vs. without lure = 0.968) -- same conclusion as before, lure isn't suppressing or inflating the PCoA1 signal. PCoA3's with/without-lure t-statistic correlation is much weaker (0.388), but the hit count is 0 either way so this doesn't change any conclusion.
- **Direct insect~date screen, all 153 taxa: 71/153 (46%) significant at q<0.10** -- much higher than the family-filtered table's 18/52 (35%), consistent with the broader table's fig described in the alpha-diversity full-insect-table entry above (the wider insect community carries a stronger, cleaner seasonal signal than Curculionidae+Latridiidae alone). Top hits mix familiar Scolytinae/Latridiidae names (*Xyleborinus attenuatus*, *Melanophthalma* sp., *Pissodes strobi*) with several taxa outside the two-family table entirely (*Melanotus hyslopi* click beetle, *Isorhipis obliqua*, *Neoclytus acuminatus* longhorn beetle, *Paromalus teres*, *Ampedus melanotoides*) -- the family restriction was excluding real, strongly date-associated taxa, not just noise.
- **Fungal~date screen (same 200-candidate approach): 200/200 significant**, same as the general workflow script (expected by construction, candidates are pre-selected for date correlation) -- unchanged since it doesn't depend on the insect table.

**Bottom line:** the prevalence-filtered insect table doesn't overturn any core conclusion (Cytospora/PCoA1/date story intact, PCoA3/lure still null, date-adjustment still collapses most of the PCoA1 signal) but it (a) surfaces a new, non-date-collinear PCoA2 signal absent from the family-filtered analysis, and (b) shows the family restriction was leaving real seasonal insect signal on the table (71/153 vs. 18/52 date-associated taxa). Worth treating the PCoA2 result as a genuine follow-up candidate rather than a construction artifact, given it isn't just riding date.

Output: data/compare_insects_fungi_top3axes_prevalence_filtered_insect/ (10 fungal_insect_association_results.*.csv files, fungal_taxa_lure_direct_check.csv, insect_taxa_date_association.csv, fungal_taxa_date_association.csv); figures/compare_insects_fungi_top3axes_prevalence_filtered_insect/ (coca_biplot.pdf, fungal_insect_volcano.png, fungal_insect_volcano.PCoA3.{with,no}_lure.png).

### #################################################################################
### What IS insect_PCoA2 under the prevalence-filtered table, and what fungal taxa track it?
### #################################################################################

Added compare_insect_fungi/insect_pcoa_date_lure_association.prevalence_filtered_insect.r -- a prevalence-filtered-table rebuild of insect_pcoa_date_lure_association.r, to follow up on the new PCoA2 signal found above. Two questions: (1) does insect_PCoA2 track date, site, or lure under this table, and (2) what are the fungal taxa driving the PCoA2 CoCA hit?

**(1) insect_PCoA2 is a lure axis, not a date axis -- and PCoA3 is still (even more clearly) a lure axis too:**

| Axis | %var | date (r, partial R2, p) | lure (partial R2, p) | site (partial R2, p, parametric) |
|---|---|---|---|---|
| PCoA1 | 19.3% | r=0.84, R2=0.80, **p=0.001** | R2=0.36, **p=0.029** | R2=0.10, p=0.078 |
| PCoA2 | 9.5% | r=0.06, R2=0.008, p=0.465 (n.s.) | R2=0.27, **p=0.005** | R2=0.16, **p=0.012** |
| PCoA3 | 7.3% | r=-0.13, R2=0.056, p=0.062 (n.s. after BH) | R2=0.69, **p=0.003** | R2=0.028, p=0.617 |

So under the broader (prevalence-filtered, all-families) insect table, PCoA1 is still clearly the date axis (as in the family-filtered table) but now ALSO carries a real secondary lure signal (R2=0.36, p=0.029 -- new; PCoA1 was not tested for lure in earlier scripts). PCoA2 has essentially zero date association (r=0.06, matching the CoCA finding that its fungal hits don't change under a date covariate) but a real lure association (R2=0.27) plus a modest site effect (R2=0.16, parametric/uncorrected). PCoA3 remains the strongest lure axis by far (R2=0.69, even a bit stronger than the family-filtered table's PCoA3, which had R2=0.81 in the original screen restricted to the top 10 axes -- comparable magnitude, different specific table). The by-lure boxplot (figures/.../insect_pcoa_top3axes.by_lure.png) shows this cleanly: PCoA2 separates Alpha-pinene_EtOH (high) from Ips (low) with Ethanol intermediate, and PCoA3 separates Ethanol (high) from Alpha-pinene_EtOH (low).

**Interpretation:** the new PCoA2 <-> fungal signal found in the prevalence-filtered CoCA workflow is not a date-confounded artifact (confirmed twice now -- unchanged by adding date as a covariate, and PCoA2 itself doesn't correlate with date) but it IS a lure-correlated axis. This doesn't invalidate the fungal association, though: every PCoA2 per-taxon test already includes lure as a separate covariate in the same model (site + lure + insect_PCoA1/2/3), and the permutation scheme permutes only insect_PCoA2 within trap while holding lure fixed -- so the reported association is specifically with the part of PCoA2 that ISN'T explained by lure, the same logic already used for PCoA3 throughout this project's PCoA3/lure work. Given PCoA2's more moderate lure correlation (R2=0.27, vs. PCoA3's 0.69-0.81), this is a cleaner separation than PCoA3 ever had, and PCoA3's fungal association remains completely null regardless (0/200 throughout) -- so a real fungal effect showing up specifically on PCoA2 rather than PCoA3 isn't simply "the lure axis, whichever one it happens to be."

**(2) Taxonomy of the fungal taxa behind the PCoA2 hit** (6/200 without date, 7/200 with date -- 6 taxa shared between the two, all same sign/negative t-statistic, i.e. all track the same end of PCoA2):

| Taxon | Family / Order | Genus / species | Notes |
|---|---|---|---|
| ASV_616 | Valsaceae / Diaporthales | *Cytospora prunicola* | The same species/genus flagged as the single most robust cross-cutting taxon throughout this entire project (top PCoA1 hit, top direct-date hit) -- now also the top PCoA2 hit. Independent third line of evidence for this taxon. |
| ASV_13242, ASV_12553 | Tympanidaceae / Helotiales | *Tympanis* sp. | Canker/twig fungus, flagged in the very first (pre-filter) CoCA screen at the top of this log as a top PCoA1-associated genus (3 ASVs then) -- reappearing here on a completely different axis. |
| ASV_9178 | Orbiliaceae / Orbiliales | *Orbilia eucalypti* | New to this project's hit lists -- nematode-trapping/saprobic genus. |
| ASV_12734 | unclassified family / Agaricomycetes | *Mycobernardia incrustans* (Corticiaceae, per the no-date model's join) | New. |
| ASV_15958 | Trichomeriaceae / Chaetothyriales | *Trichomerium* sp. | Only significant in the date-adjusted model; sooty-mold-associated family, new to this project. |

All 6-7 hits share the same negative-t-statistic direction (the internal-consistency check used throughout this project as a sign of real, non-noise signal) -- i.e. all track the same end of PCoA2 (the Ips-lure-associated / low-PCoA2 end per the boxplot). Re-finding *Cytospora* and *Tympanis* -- both flagged independently on PCoA1/date in the original, unfiltered, very first CoCA pass at the top of this log -- via a completely different axis and insect-table construction is a meaningful third line of support for those two taxa specifically being real, recurring signals rather than construction-dependent artifacts.

Output: data/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_pcoa_axes.date_site_lure_association.csv, PCoA2_significant_taxa_taxonomy.csv; figures/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_pcoa_top3axes.by_date.png, insect_pcoa_top3axes.by_lure.png.

### #################################################################################
### Correction: insect_PCoA1-3's date relationship is non-linear (hump-shaped), not absent
### #################################################################################

User caught, from the insect_pcoa_top3axes.by_date.png scatterplot, that insect_PCoA2's relationship with date looks like a hump (rising from early May to a peak around late May/early June, then falling through July) rather than a flat/null relationship or a straight line. **This matters because the linear date test above (r=0.06, p=0.465, reported as "no date association") has essentially zero power to detect a symmetric hump -- a hump's linear correlation with date is close to zero by construction regardless of how strong the real relationship is.** That earlier "PCoA2 has no date association" conclusion needs a correction, not just a caveat.

Extended insect_pcoa_date_lure_association.prevalence_filtered_insect.r (Part 6) with a proper test: add a quadratic date term (centered date, `+ I(date_c^2)`) to the same site+lure-adjusted model and ask whether it explains additional variance beyond the linear term, using the same trap-blocked date permutation as the linear test. Result -- **all three axes have a highly significant quadratic date component (q=0.001 each)**, and it is by far strongest for PCoA2 (F=92.7, vs. F=26.6 for PCoA1 and F=16.1 for PCoA3) -- i.e., PCoA2's relationship with date is not just present, it's the single strongest non-linear date signal of the three axes, despite having ~zero linear date correlation. A loess-fit version of the by-date scatterplot (figures/.../insect_pcoa_top3axes.by_date.loess.png) confirms this visually: PCoA2 has a clean hump peaking in late May/early June; PCoA3 shows a weaker version of the same hump shape; PCoA1's already-known near-monotonic rise is essentially unchanged by the loess fit.

**Revised interpretation of insect_PCoA2**: it's a seasonal-timing axis after all, just capturing a mid-season peak/decline rather than a directional trend -- on top of the already-established lure association (Part 1 above). This is worth keeping in mind for the fungal CoCA hits on PCoA2 (Cytospora, Tympanis, etc., see above): those tests used a LINEAR date covariate (matching the general workflow script's convention throughout this project), so they've only been checked against linear date confounding, not this newly-discovered non-linear component. The fungal PCoA2 hits should be treated as not yet fully deconfounded from date until a quadratic-date-adjusted version of that per-taxon test is run.

Output: data/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_pcoa_axes.quadratic_date_test.csv; figures/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_pcoa_top3axes.by_date.loess.png.

### #################################################################################
### Which insect taxa drive insect_PCoA2 (and PCoA1/PCoA3, for completeness)?
### #################################################################################

Also added Part 7 to the same script: Pearson correlation of each of the 153 prevalence-filtered insect taxa's Hellinger-transformed abundance (the exact data the PCoA was built from) against each axis's scores -- the formalized version of the informal PCoA1-driver check from earlier in this project ("Melanophthalma sp., cor=0.96"), now covering all three axes and the full all-families taxon set.

**insect_PCoA2's top drivers are almost entirely Scolytinae (bark/ambrosia beetles), but a DIFFERENT set of species than PCoA1**, split into two opposing groups:
- Positive end: *Xylosandrus germanus* (r=0.68, also a top PCoA1 negative/early-season driver), *Ips grandicollis* (r=0.48), *Pissodes affinis* (r=0.47), *Anisandrus sayi* (r=0.45), *Dryocoetes affaber* (r=0.42), *Xylosandrus crassiusculus* (r=0.38), *Dendroctonus valens* (r=0.38) -- plus a few non-Scolytinae (*Asemum striatum* Cerambycidae r=0.47, *Isarthrus calceatus* Eucnemidae r=0.43, *Cis fuscipes* Ciidae r=0.42, *Tetropium schwarzianum* Cerambycidae r=0.41).
- Negative end: *Xyleborinus attenuatus* (r=-0.46) -- the single strongest NEGATIVE driver of PCoA1 as well (r=-0.85 there), i.e. this species anchors the early/bark-beetle end of both axes.

So PCoA2 reads as a within-Scolytinae/ambrosia-beetle axis -- differentiating WHICH bark/ambrosia beetle species dominate a sample (several different Xylosandrus/Ips/Dryocoetes/Anisandrus species vs. Xyleborinus attenuatus specifically) -- consistent with both of Part 1's findings above: a real lure association (different bark beetle species have different lure preferences -- e.g. Ips grandicollis and the Ips-branded lure) and now a real mid-season timing peak (different bark beetle species have staggered emergence/flight phenology within the same May-July window). PCoA1, by contrast, differentiates Scolytinae as a whole from the fungivorous Latridiidae/Corticariinae group (Melanophthalma, r=+0.83) -- a broader guild-level split. PCoA3's top drivers (*Anisandrus sayi* r=0.69, *Ips grandicollis* r=-0.62, *Hylobius pales* r=-0.58, *Xyleborinus saxesenii* r=0.58) overlap partially with PCoA2's cast of characters but with different, sometimes opposite signs -- e.g. *Anisandrus sayi* and *Ips grandicollis* have the SAME sign relationship on PCoA2 (both positive) but OPPOSITE signs on PCoA3, meaning the two axes separate these species differently (plausibly the lure-preference axis, PCoA3, vs. the timing axis, PCoA2).

Output: data/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_pcoa_axes.taxon_drivers.csv (all 153 taxa x 3 axes); figures/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_pcoa2_top_taxon_drivers.png.

### #################################################################################
### CORRECTION: the insect_PCoA2 fungal signal (Cytospora, Tympanis, etc.) does NOT
### survive proper (linear+quadratic) date adjustment -- it was a date artifact
### #################################################################################

Direct follow-up to the quadratic-date finding above. Added compare_insect_fungi/insect_fungal_coca_pcoa2_quadratic_date_check.prevalence_filtered_insect.r: reran the exact same per-taxon permutation test insect_fungal_coca_analysis.prevalence_filtered_insect.r used for insect_PCoA2 (same top-200 CoCA candidate set, same trap-blocked permutation of insect_PCoA2, same site+lure+insect_PCoA1/2/3 covariates), but replacing the linear-only date covariate with a centered linear + quadratic term (`date_c + I(date_c^2)`) -- the proper adjustment given insect_PCoA2's date relationship is a hump, not a line (previous entry).

**Result: 0 of 200 candidates significant at q<0.10 -- down from 6/200 (no date) and 7/200 (linear date only).** None of the 7 taxa significant under linear-date adjustment (ASV_13242/*Tympanis*, ASV_12553/*Tympanis*, ASV_616/*Cytospora prunicola*, ASV_3449, ASV_9178/*Orbilia eucalypti*, ASV_12734/*Mycobernardia*, ASV_15958/*Trichomerium*) survive the quadratic adjustment -- 0/7. The best surviving q-value in the quadratic-adjusted model is 0.20 (several taxa, including ASV_616/*Cytospora* itself at q=0.23), nowhere near the q<0.10 threshold used throughout this project.

**This reverses the previous entry's interpretation.** The earlier writeup treated *Cytospora*'s appearance on insect_PCoA2 as "a third, independent line of evidence" (beyond PCoA1 and the direct date screen) for that taxon being a real, non-date-confounded signal. That conclusion does not hold: insect_PCoA2's association with *Cytospora* (and every other PCoA2 candidate) collapses exactly the same way insect_PCoA1's date-adjusted signal always has throughout this project (81/200 -> 0/200, 62/200 -> 9/200, etc.) once date is adjusted for correctly. The apparent "independence" from date in the earlier (linear-only) check was an artifact of testing the wrong functional form of date, not evidence of a real insect-community-mediated effect. **PCoA1's already-established date-adjusted survivors (the handful of taxa surviving 81/200->0/200-style linear adjustment, none of which include Cytospora at q<0.10 either) remain the only defensible fungal/insect-axis associations in this project, and even those are already flagged elsewhere as weak/marginal.** *Cytospora*'s status as "the most robust cross-cutting taxon in this project" now rests on exactly two lines of evidence, not three: the direct fungal~date screen, and its (linear-date-collapsing) PCoA1 loading -- both of which are fundamentally date-driven observations, not independent confirmations of an insect-community relationship.

**General lesson for this project going forward**: any future PCoA axis found to correlate with a fungal taxon should be checked for a non-linear date relationship (as done here) before its date-adjusted survival is treated as meaningful -- a linear-only date covariate is not a sufficient deconfounding check when the underlying seasonal relationship isn't linear, which insect_PCoA1's already-known near-linear date trend obscured as the "obvious" functional form to test against.

Output: data/compare_insects_fungi_top3axes_prevalence_filtered_insect/fungal_insect_association_results.PCoA2_plus_quadratic_date.csv, PCoA2_date_adjustment_comparison.csv.

### #################################################################################
### PCoA1's quadratic-date-adjusted survivors: Cytospora still does NOT survive
### #################################################################################

Follow-up question: insect_PCoA1 also has a significant quadratic date component (F=26.6, p=0.001, from the earlier quadratic-date-test entry) on top of its already-known strong linear trend (r=0.84) -- does PCoA1's *Cytospora prunicola* (ASV_1905) hit survive the same linear+quadratic date adjustment that erased the entire PCoA2 signal? Added compare_insect_fungi/insect_fungal_coca_pcoa1_quadratic_date_check.prevalence_filtered_insect.r, identical logic to the PCoA2 version, applied to PCoA1's top-200 candidate set.

**Answer: no. ASV_1905 (*Cytospora prunicola*) does not survive under any date adjustment tested so far, including the original linear one:**

| Model | t_stat | q_value |
|---|---|---|
| PCoA1_only (no date) | -7.65 | 0.017 |
| PCoA1_plus_date (linear) | -3.35 | **0.109** (already fails q<0.10, by a hair) |
| PCoA1_plus_quadratic_date | -2.41 | **0.160** (fails more) |

This is a correction to how the linear-date PCoA1 result was characterized in the very first prevalence-filtered-table CoCA entry above, which reported "9/200 with lure+date" and listed ASV_1905 in the top-15 preview table without flagging that its q=0.109 already misses the q<0.10 cutoff used everywhere else in this project. *Cytospora* was never actually among PCoA1's significant date-adjusted survivors under the prevalence-filtered table, linear or quadratic.

**Unlike PCoA2 (whose signal collapsed entirely, 7->0), PCoA1's hit count actually INCREASES under quadratic adjustment: 9 -> 19 significant, and all 9 linear-adjusted hits remain significant (a strict superset, nothing lost).** This makes sense given PCoA1 has both a strong linear AND a real quadratic date component (unlike PCoA2, which is *only* quadratic) -- modeling the true (curved) date relationship removes more residual date-driven noise than a straight-line fit does, revealing additional taxa whose association was previously masked by date-model misspecification rather than being spurious.

Two taxa worth flagging among the newly-revealed quadratic-only hits:
- **ASV_192 (Valsaceae, genus unresolved)** -- q=0.248 under linear date, q=0.095 under quadratic (newly crosses the threshold). Same family as *Cytospora*/*Valsa* (the genus most associated with the bark-beetle end of this project's story throughout), though not resolved enough to confirm the same genus.
- **ASV_13242 (Tympanidaceae, *Tympanis* sp.)** -- significant under BOTH linear (q=0.04) and quadratic (q=0.0875) date adjustment on PCoA1. This is the same ASV that was the top PCoA2 hit before that entire signal collapsed under quadratic adjustment -- so *Tympanis* (ASV_13242) has a real, quadratic-date-robust association with PCoA1 specifically, not PCoA2.

**Bottom line:** *Cytospora prunicola*'s cross-cutting-taxon status in this project rests entirely on the direct fungal~date screen and its raw (no-date) PCoA1/PCoA2 loadings -- it has never actually survived a proper date-adjusted test on any insect-community axis under the prevalence-filtered table. *Tympanis* sp. (ASV_13242), by contrast, is the one taxon so far with a genuinely date-robust (both linear- and quadratic-adjusted) association with an insect PCoA axis (PCoA1) -- worth promoting ahead of *Cytospora* as the strongest remaining candidate for a real (not purely date-confounded) insect-community-linked fungal signal, at least under this insect-table construction.

Output: data/compare_insects_fungi_top3axes_prevalence_filtered_insect/fungal_insect_association_results.PCoA1_plus_quadratic_date.csv, PCoA1_date_adjustment_comparison.csv.

### #################################################################################
### factor(date) (maximally strict deseasonalizing) check -- PCoA2's signal is
### NOT stably zero; it depends on the specific functional form of the date control
### #################################################################################

Added compare_insect_fungi/insect_fungal_coca_factordate_check.prevalence_filtered_insect.r: the strictest possible deseasonalizing check, replacing the date covariate with `factor(date)` (6 levels, absorbing ANY functional form of date, not just a line or a parabola) -- the same logic already used elsewhere in this project (insect_fungal_pairwise_taxon_association.residualized.r's `factor(date)` residualization). Reran both PCoA1's and PCoA2's top-200 candidate sets.

**Full comparison across all four date-adjustment strategies:**

| Axis | no date | linear date | linear+quadratic date | factor(date) |
|---|---|---|---|---|
| PCoA1 | 59/200 | 9/200 | 19/200 | 18/200 |
| PCoA2 | 6/200 | 7/200 | **0/200** | **9/200** |

**PCoA1 is stable across the stricter controls** (9 -> 19 -> 18, roughly consistent once *any* non-linear flexibility is allowed) -- quadratic and factor(date) agree closely, both a large step up from linear-only, both well below the untouched 59. This is the reassuring, expected pattern: once date is modeled with enough flexibility to capture its real shape, the result stabilizes.

**PCoA2 is NOT stable -- it swings from 0 (quadratic) back up to 9 (factor(date)), a bigger count than either the no-date or linear-date versions.** This is a real result, not a bug (spot-checked candidate lists and formulas), but it means PCoA2's fungal association is sensitive to the exact functional form used to control for date, in a way PCoA1's is not. The likely explanation: factor(date) absorbs not just the smooth hump quadratic captures but also any idiosyncratic, single-date-specific irregularity in insect_PCoA2 (weather, a trap anomaly, etc. on one particular collection date) that a 2-parameter quadratic curve cannot represent. Depending on how a given fungal taxon's abundance relates to that idiosyncratic component vs. the smooth component, removing it can either suppress or (as apparently happened here) sharpen the taxon's estimated association with the AXIS's remaining within-date (largely lure-driven, per the earlier date/lure/site test) variation. With only 6 dates and 69 samples, this is exactly the kind of small-sample sensitivity that should be flagged, not overinterpreted as a confirmed finding either way.

**Specific taxa, tracked across all four PCoA1 strategies:**

| Taxon | no date | linear | quadratic | factor(date) |
|---|---|---|---|---|
| ASV_1905 (*Cytospora prunicola*) | q=0.017 (sig) | q=0.109 | q=0.160 | q=0.257 |
| ASV_13242 (*Tympanis* sp.) | q=0.350 (n.s.) | q=0.040 (sig) | q=0.0875 (sig) | q=0.094 (sig) |

*Cytospora prunicola* (ASV_1905) fails every date-adjusted test, consistently, and its q-value gets monotonically worse as the date control gets stricter -- the clearest, most consistent "this is a date artifact" pattern in the whole comparison. *Tympanis* sp. (ASV_13242) is the opposite: null with no date control, then significant under all three date-adjusted versions (linear, quadratic, factor) -- the single most date-robust PCoA1 hit found so far.

**But the genus/family-level picture is more robust than any single ASV, and re-emerges on PCoA2 under factor(date):** the 9 PCoA2 factor(date) survivors include **4 separate Valsaceae ASVs** -- ASV_192 (unresolved genus), ASV_9007 (*Cytospora*), ASV_6273 (*Cytospora*), ASV_11019 (*Valsa*) -- plus **ASV_5640, a SECOND, different *Tympanis* ASV** (Tympanidaceae) as the single strongest hit (t=-4.25, q=0.067). Notably, ASV_616 (the specific *C. prunicola* OTU) is NOT among these survivors (q=0.20) -- it's other members of the same genus/family that show up. So across this whole investigation, no single specific ASV survives every test, but **Valsaceae (*Cytospora*/*Valsa*) and *Tympanis* as genus/family-level groups keep reappearing across different axes, different ASVs within the group, and different date-control strategies** -- a pattern more consistent with "this fungal family has a real, if diffuse, relationship with the insect community that isn't fully reducible to any single ASV or any single functional form of date" than with either "it's all just Cytospora prunicola" (the original framing) or "it's nothing, purely date" (the quadratic-only framing from the previous entry).

**Bottom line / recommended framing going forward:** don't cite *Cytospora prunicola* (ASV_1905) specifically as a date-robust insect-community signal -- it isn't, under any adjustment tested. Do treat Valsaceae/*Cytospora*-as-a-group and *Tympanis*-as-a-group as the more defensible candidates, since they reappear (via different specific ASVs) across PCoA1 and PCoA2 and across every date-control strategy from quadratic onward -- but flag that PCoA2's specific hit count is unstable across date-control choices and shouldn't be quoted as a single fixed number (e.g. "9/200") without noting that alternative reasonable specifications gave anywhere from 0 to 9.

Output: data/compare_insects_fungi_top3axes_prevalence_filtered_insect/fungal_insect_association_results.PCoA1_plus_factordate.csv, fungal_insect_association_results.PCoA2_plus_factordate.csv, date_adjustment_comparison.all_strategies.csv.

### #################################################################################
### Direct taxon~date screens updated with a quadratic term: which taxa are linear
### (early-/late-season) vs. hump/dip-shaped?
### #################################################################################

Extended the two main direct taxon~date screens -- fungal_community_seasonality.r (all 6,342 prevalence-filtered fungal taxa, unbiased) and insect_fungal_coca_analysis.general_workflow.r (all 52 Curculionidae+Latridiidae insect taxa, plus its own top-200-by-correlation fungal date screen) -- to test a quadratic date term alongside the existing linear one, so each taxon can be classified as a linear early-/late-season trend or a real hump/dip (non-monotonic) seasonal pattern, per user request. Also added a new standalone script, compare_insect_fungi/insect_taxa_date_association.prevalence_filtered_insect.r, running the identical test on the >=5-sample prevalence-filtered (all-families, 153-taxon) insect table -- broken out on its own rather than rerunning the full prevalence-filtered CoCA workflow (which also fits CoCA/LOO-CV/the 10-variant PCoA grid, all irrelevant to this specific question), so it can be rerun quickly.

**Backward compatibility, by design:** the linear-term columns (`t_stat`, `p_perm`, `q_value`) are UNCHANGED in both meaning and values from before -- both output files (`fungal_taxa_date_association.all_taxa.csv`, `insect_taxa_date_association.csv`) feed a dozen+ downstream scripts (fig3-fig8, compare_selection_strategies.r, insect_fungal_coca_analysis.PCoA1_date_axis.r, responsive_taxa_abundance_by_site.R, etc.) that all filter on `q_value < 0.10`, so the quadratic test was added as NEW columns (`t_stat_quad`, `p_perm_quad`, `q_value_quad`, `shape`) rather than replacing anything. Spot-checked: 34.5% fungal linear hit rate (2186/6342) matches the previously-reported ~33.4% within normal run-to-run noise; nothing downstream needs to change.

One implementation bug caught and fixed before running: `c(t_stat = obs["t_linear"], ...)` -- combining a named scalar into a new named vector element -- mangles the name to `"t_stat.t_linear"` in R, which would have made every `r["t_stat"]` extraction silently return `NA`. Fixed with `unname()` on each `obs[...]` element before combining; verified with a standalone test before running the full 6,342-taxon screen.

**Results, all three tables (q<0.10 for both linear and quadratic terms, independently BH-corrected):**

| Screen | n taxa | Linear-sig | Quadratic-sig (hump/dip) | Hump | Dip | Linear late | Linear early | n.s. |
|---|---|---|---|---|---|---|---|---|
| Fungal, unbiased (fungal_community_seasonality.r) | 6,342 | 2,186 (34.5%) | 707 (11.1%) | 377 | 330 | 1,035 | 748 | 3,852 |
| Insect, family-filtered (general_workflow.r) | 52 | 19 (36.5%) | 9 (17.3%) | 5 | 4 | 3 | 11 | 29 |
| Insect, prevalence-filtered, all families (new script) | 153 | 69 (45.1%) | 38 (24.8%) | 23 | 15 | 32 | 22 | 61 |
| Fungal, top-200 date-correlation screen (general_workflow.r, pre-selected, not unbiased) | 200 | 200 (100%, by construction) | 107 (53.5%) | 17 | 90 | 58 | 35 | 0 |

Non-trivial fractions of every taxon set show a real hump/dip shape rather than a monotonic trend -- roughly a third to a half of each screen's "date-associated" taxa (by whichever measure) turn out non-monotonic once tested for it, confirming this was a real gap in every prior "X% of taxa are seasonal" statement in this project (all of which only ever tested linear date).

**Notable patterns:**
- The fungal unbiased screen's hump/dip split is roughly even (377 hump vs. 330 dip), but the pre-selected top-200 (candidates chosen for the STRONGEST raw |correlation| with date, which favors monotonic taxa almost by construction) skews heavily toward dip (90) over hump (17) among its quadratic-significant subset -- a selection-method artifact worth remembering when comparing hit-rate percentages across these two fungal screens.
- *Xyleborinus attenuatus* (the single strongest linear PCoA1/date driver throughout this project) is a clean "dip" taxon in BOTH insect screens (family-filtered: t_quad=10.0, q=0.017; prevalence-filtered: t_quad=9.3, q=0.014) -- i.e. its actual trajectory is high-early, low-mid-season, rising again late, not a simple monotonic decline the way it's been described everywhere in this project so far ("earlier-season" was the linear characterization). *Melanophthalma* sp. (the other headline PCoA1 taxon) is a clean "hump," peaking mid-season rather than simply "later-season."
- *Anisandrus sayi* is a good example of why the quadratic term matters: its LINEAR term is completely null (p=0.218, family-filtered table) but its quadratic term is highly significant (q=0.017, hump-shaped) -- a taxon the original linear-only screen would have called "not date-associated" at all.

Output: data/2024_fungi/fungal_taxa_date_association.all_taxa.csv (updated in place, new columns appended), data/2024_insect_data/insect_taxa_date_association.csv (updated in place), data/compare_insects_fungi_top3axes/fungal_taxa_date_association.csv (updated in place), data/compare_insects_fungi_top3axes_prevalence_filtered_insect/insect_taxa_date_association.csv (new, with Family/Subfamily/Order context columns).

### #################################################################################
### 2026-08-31 -- Lineage B port: whole-community Procrustes + fig2 (NMDS/Procrustes)
### #################################################################################

Started porting the "Known gaps" list in project_organization.md to Lineage B
(prevalence-filtered, all-families, >=5-sample insect table). Two scripts this round.

**1. compare_insect_fungi/insect_fungal_procrustes.prevalence_filtered_insect.r** --
Lineage-B rebuild of insect_fungal_procrustes.r. Same whole-community Procrustes
concordance test (insect NMDS vs. an independently-built rarefied fungal NMDS), same
69-sample matched dataset, only the insect table changes: 153-taxon >=5-sample
prevalence filter (all families) instead of the Curculionidae+Latridiidae restriction.
Two deliberate deviations from the Lineage-A standalone script, both bringing it in
line with presentation_items/fig2_nmds_procrustes.R (which already does both):
  - Added the `Kingdom == "k__Fungi"` ASV filter. The Lineage-A standalone procrustes
    script was never one of the "5 core scripts" that filter was retrofitted into, so
    it still runs on the unfiltered ASV table; this port uses the filter (15,422 of
    17,139 ASVs confirmed fungal).
  - Handles the rarefaction sample drop (NH.54 falls below depth 5000 once non-fungal
    reads are stripped -- the same 69->68 drop documented in the Kingdom-filter entry
    above) by restricting the Procrustes fit to the 68 survivors, instead of the
    Lineage-A script's `stopifnot` that assumed no drop.
Also computes and writes the matched-69-sample insect community PERMANOVA
(site + lure + date, marginal, permute::how() with the within-trap/site-blocked
scheme from insect_ords.prevalence_filter.R) so the fig2 port can annotate its insect
panel without re-running an ordination -- same role fungal_community_seasonality.r
plays as the (lineage-invariant) source of the fungal PERMANOVA.

Output dir: data/compare_insects_fungi_procrustes_prevalence_filtered_insect/,
figures/compare_insects_fungi_procrustes_prevalence_filtered_insect/ (dedicated
Lineage-B dir -- the Lineage-A script wrote CSVs into data/2024_insect_data/ and
figures/ root, the same anti-pattern the Kingdom-filter cleanup fixed for
general_workflow.r; not replicated here).

**Results -- the Procrustes conclusion is unchanged under Lineage B:**
  - Insect NMDS stress 0.163 (repeatable; matches insect_ords.prevalence_filter.R's
    0.162 for this table). Fungal NMDS stress 0.166.
  - **Procrustes correlation = 0.604, p = 0.001** (999 perms, symmetric) -- vs.
    Lineage A's 0.56. Essentially the same moderate, significant whole-community
    concordance: real shared structure, far from a 1:1 mapping. Consistent, as in
    Lineage A, with both community PERMANOVAs being dominated by a shared seasonal/
    site term (see below + fungal Site:date R2=22.7% established earlier), and with
    the per-taxon finding that no individual fungal taxon's association with an insect
    PCoA axis survives proper date adjustment -- a ~0.6 whole-community correlation is
    what shared season+site structure alone should produce, not evidence of direct
    taxon-to-taxon interaction.
  - Matched-dataset insect PERMANOVA: site R2=9.5% (F=3.00, p=0.001), lure R2=11.0%
    (F=5.18, p=0.001), date R2=13.5% (F=12.69, p=0.001), residual 65.8%. All three
    terms significant, same three-way structure as (a) the Lineage-A fig2 insect
    PERMANOVA (site 9.7%, lure 14.2%, date 14.2%) and (b) insect_ords.prevalence_
    filter.R's 71-sample version (site 9.5%, lure 11.4%, date 12.7%) -- the 69-sample
    match trims lure/date effect sizes slightly but changes nothing qualitatively.
  - Worst-matched (highest-residual) samples: NH.34 (Durham, Alpha-pinene_EtOH,
    2024-05-15), NH.99 & NH.46 (Manchester Cedar Swamp), NH.113 (Pease) -- all
    early-season (May), and nearly the same set the Lineage-A script flagged (NH.34,
    NH.99, NH.46, NH.113 all appeared there too). Still spans several sites/lures, no
    single obvious driver.

**2. presentation_items_prevalence_filtered_insect/fig2_nmds_procrustes.R** --
Lineage-B counterpart of presentation_items/fig2_nmds_procrustes.R, same filename
(dedicated-folder convention), same 3-panel layout / aesthetics / fig-number mapping.
Panel a insect NMDS + PERMANOVA annotation from script 1's output; panel b fungal
NMDS annotated from the shared data/2024_fungi/fungal_community_permanova.csv (fungal
side is lineage-invariant, so this is the exact file Lineage-A fig2 reads); panel c
Procrustes overlay (r = 0.60, p = 0.001, 68 samples). One insect NMDS + one fungal
NMDS run feed all three panels; the fungal rarefaction/NMDS/Procrustes fit is cached
to data/presentation_items_prevalence_filtered_insect/fig2_ordination_cache.rds
(force_recompute <- T by default, as in Lineage A).

Output: figures/presentation_items_prevalence_filtered_insect/fig2_nmds_procrustes.
{png,pdf}; data/compare_insects_fungi_procrustes_prevalence_filtered_insect/
(insect_fungal_procrustes_summary.csv, insect_fungal_procrustes_residuals.csv,
insect_community_permanova.csv); figures/compare_insects_fungi_procrustes_prevalence_
filtered_insect/ (insect_fungal_procrustes.pdf, insect_fungal_procrustes.NMDS_
overlay.png).

### #################################################################################
### 2026-08-31 -- fig1 panel c: top insect genera OVERALL, not just within the
### two focal families
### #################################################################################

presentation_items/fig1_taxonomic_breakdown.R panel c previously showed the 10 most
abundant genera *within Curculionidae + Latridiidae*, renormalized inside that
two-family subset. Changed it to the 10 most abundant insect genera across ALL
families, renormalized against each sample's total catch (the same denominator as
panel a), so panel c now reads as a finer-grained panel a. Panels a, b, d and the
rest of the script are unchanged.

First checked whether this actually changes what's shown -- i.e. whether the top
genera overall are just the focal-family genera anyway. They are not:

| Rank | Genus | Individuals | % of catch (n=10,962) | Family | In focal 2? |
|---|---|---|---|---|---|
| 1 | *Xylosandrus* | 2451 | 22.4% | Curculionidae | yes |
| 2 | *Melanophthalma* | 1829 | 16.7% | Latridiidae | yes |
| 3 | *Xyleborinus* | 867 | 7.9% | Curculionidae | yes |
| 4 | *Anisandrus* | 522 | 4.8% | Curculionidae | yes |
| 5 | *Enoclerus* | 313 | 2.9% | **Cleridae** | **no** |
| 6 | (genus-less bucket) | 291 | 2.7% | >=13 families | -- |
| 7 | *Cyclorhipidion* | 258 | 2.4% | Curculionidae | yes |
| 8 | *Dendroctonus* | 199 | 1.8% | Curculionidae | yes |
| 9 | *Pissodes* | 182 | 1.7% | Curculionidae | yes |
| 10 | *Heteroborips* | 154 | 1.4% | Curculionidae | yes |
| 11 | *Asemum* | 144 | 1.3% | **Cerambycidae** | **no** |

So the focal-family restriction was hiding *Enoclerus* (Cleridae -- checkered
beetles, common bark-beetle predators) at rank 5. The rank-6 slot is a genus-level-
unclassified bucket that pools Staphylinidae/Anthocoridae/Chrysomelidae/Cicadellidae/
etc. -- not a real genus.

Handling (per user decision): the genus-less rows are **excluded from the ranking**
(they aren't a genus and span many families) but still count toward each sample's
total and fall into "Other". With that exclusion the shown set is 10 real genera --
the 8 focal-family genera above plus *Enoclerus* (Cleridae, rank 5) and *Asemum*
(Cerambycidae, promoted into rank 10). Switching to the overall ranking also drops
*Orthotomicus* and *Dryocoetes* out of the shown set relative to the old panel c.
Because the denominator is now the whole catch, panel c's "Other" bar is large
(~44-58% per site): it absorbs every non-top genus across all ~98 families plus the
genus-less individuals.

Note fig1 is entirely lineage-invariant -- every panel is built from the raw trap
catch / raw ASV table, none of it touches the 52-taxon (Lineage A) or 153-taxon
(Lineage B) community table -- so this edit is not lineage-specific and there is no
separate Lineage-B fig1 to keep in sync.

Output: figures/presentation_items/fig1_taxonomic_breakdown.{png,pdf} (regenerated).

### #################################################################################
### 2026-08-31 -- Lineage B port: fig3 volcano plots, expanded to carry the
### quadratic date term
### #################################################################################

presentation_items_prevalence_filtered_insect/fig3_volcano_plots.R. The Lineage-A
fig3 was a 1x3 row of volcanoes on the LINEAR date term (a insect~date, b
fungal~date) plus one fungal~insect_PCoA1 panel (c). Under Lineage B the direct
date screens now also carry a QUADRATIC date term and a hump/dip `shape` column
(the 2026-08-31 taxon~date entries above), and there is no clean way to show two
t-statistics per taxon on one volcano. Resolved by expanding to a **2x3 grid,
rows = functional form of the seasonal test**:

| | insect ~ date | fungal ~ date | fungal ~ insect axis |
|---|---|---|---|
| **Row 1 -- linear** | a  linear-term volcano | b  linear-term volcano | c  PCoA1, factor(date)-adjusted |
| **Row 2 -- quadratic** | d  quadratic-term volcano | e  quadratic-term volcano | f  PCoA2, factor(date)-adjusted |

- Volcano y = -log10(q) throughout. Row 1 a/b x = linear t (wings = earlier / later
  season). Row 2 d/e x = quadratic t, with the sign mapped to shape in the axis
  label: `t < 0: hump` (concave-down, mid-season peak) / `t > 0: dip`. So d/e read
  with the same left/right grammar as a/b.
- Panel f = fungal ~ insect_PCoA2 was the user's call: it makes row 2 cohere as
  "the non-linear seasonal story" -- d/e are the quadratic date terms, and
  insect_PCoA2 is the axis whose fungal hits are what first surfaced the quadratic
  date trend (the "what IS insect_PCoA2" entry above). c and f are a matched pair:
  PCoA1 (the ~linear date axis) in row 1, PCoA2 (the hump-shaped date axis) in row 2.
- **Panels c and f use the factor(date)-adjusted results**
  (`fungal_insect_association_results.PCoA{1,2}_plus_factordate.csv`), not the
  permissive unadjusted ones, per user request and per
  `interpretation_update_prevalence_filtered_insect_08282026.md` sec. 5. factor(date)
  is a 6-level term absorbing any functional form of season, so these panels show
  the DESEASONALIZED insect-community association -- the defensible version -- rather
  than the date-collinear one. Both are the with-lure variant (no `.no_lure`
  factordate file exists), so c and f are also matched on lure handling.
  - Panel c (PCoA1): 17/200 sig. PCoA1's hit count is stable once any non-linear
    date flexibility is allowed (9 linear -> 19 quadratic -> 17 factor(date)), so no
    instability caveat needed. Note *Cytospora prunicola* (ASV_1905) is NOT among the
    labeled hits here -- it fails factor(date) adjustment (q=0.257), consistent with
    the CORRECTION entries. The labeled hits are a mix of both directions (Exophiala,
    Diatrype stigma on the negative side; Sarcotrochila, a Verrucariaceae taxon,
    several unclassified on the positive side).
  - Panel f (PCoA2): 9/200 sig, headline count, BUT the subtitle flags PCoA2 as
    unstable across date controls -- 6/200 unadjusted, 0/200 under a quadratic date
    term, 9/200 under factor(date) -- and small-sample-sensitive (see text). Labeled
    hits include *Tympanis* sp. and two *Cytospora* sp. ASVs -- the genus/family-level
    Valsaceae + *Tympanis* reappearance that sec. 5 of the interpretation doc
    identifies as the defensible reading, as opposed to *C. prunicola* by name.
- `make_volcano()` generalized with `t_col`/`q_col` args so the one builder draws
  both the linear (t_stat/q_value) and quadratic (t_stat_quad/q_value_quad)
  volcanoes; significance and label-ranking use whichever term is plotted. c and f
  use the default t_stat/q_value columns (the PCoA-axis coefficient from the
  factor(date)-adjusted model).

Dynamic panel counts on this run: a 69/153 sig (linear), b 2186/6342, c 17/200
(factor(date)-adjusted), d 38/153 (quadratic q<0.10), e 707/6342, f 9/200
(factor(date)-adjusted; 6/200 unadjusted, 0/200 quadratic). Sanity checks against
the log: *Xyleborinus attenuatus* lands on the DIP side of panel d (t_quad>0) and
*Anisandrus sayi* on the HUMP side with a near-null linear term -- both as described
in the 2026-08-31 taxon~date entry.

Panels b and e read the shared, lineage-invariant
`data/2024_fungi/fungal_taxa_date_association.all_taxa.csv`; if Lineage-A fig3 is
ever given the same quadratic expansion, its b/e panels would be identical to these.

The 2x3 is assembled with `gtable` `rbind`/`cbind` (size="max"), matching fig1/fig2
rather than the Lineage-A fig3's `gridExtra::arrangeGrob` -- each column is an
rbind of its two panels (equalizes y-axis/panel widths so the boxes and x-tick
labels line up top-to-bottom), then the three columns are cbind'd (equalizes
row heights left-to-right, so panel f's 2-line subtitle bumps the whole bottom
subtitle band uniformly instead of pushing f's panel box out of line). The
freeze-to-PNG / replay-to-PDF step for ggrepel label-position consistency is
unchanged.

Output: figures/presentation_items_prevalence_filtered_insect/fig3_volcano_plots.
{png,pdf}.

### #################################################################################
### 2026-09-01 -- Lineage B port: all-against-all pairwise fungal x insect taxon
### screen (main + site/factor(date)-residualized)
### #################################################################################

Ported the last outstanding non-figure analysis to Lineage B:
compare_insect_fungi/insect_fungal_pairwise_taxon_association.prevalence_filtered_insect.r
and .residualized.prevalence_filtered_insect.r, mirroring the two Lineage-A scripts
of the same name. Method, model (fungal_clr ~ site + lure + date + insect_taxon),
vectorized-OLS solve, 999 trap-blocked permutations, and the two-tier FDR
(within-insect-taxon + global BH) are all unchanged -- the only swap is the insect
table: all families + >=5-sample prevalence filter (153 taxa) instead of
Curculionidae+Latridiidae + singleton filter (39 taxa after the parent's post-match
re-filter). Insect predictor columns are Hellinger values computed on the full
153-taxon community then subset, matching how insect_hel is built in the Lineage-B
CoCA script; the post-match >=5 re-filter is kept (Lineage-A parent does the same)
but changes nothing here -- all 153 pass, because 69 of 71 insect samples match a
fungal sample. Fungal side is lineage-invariant: same 6,342 prevalence-filtered
k__Fungi ASVs. Grid = 153 x 6,342 = 970,326 pairs, ~11 min single-core (main),
~16 min (residualized). Outputs ->
data/  and figures/compare_insects_fungi_pairwise_taxa_prevalence_filtered_insect/.

RESULT -- the broader table reproduces the Lineage-A story and adds nothing robust:

- Main grid: 2,402 pairs at q<0.10 within-insect-taxon (vs. 1,918 for Lineage A),
  but 0 pairs at q<0.10 global BH (min global q = 0.37), same as Lineage A.
- Only 5 of 153 insect taxa carry ANY significant partner:
  Pityogenes hopkinsi (1,592), Asemum striatum (266), Xyleborinus attenuatus (213),
  Orthoperus scutellaris (198), Hylastes opacus (133). The 3 Scolytinae are exactly
  the trio Lineage A flagged, with near-identical counts (P. hopkinsi 1,592 vs.
  1,593 -- an abundant taxon's Hellinger column barely shifts when the community
  basis widens). Asemum striatum (Cerambycidae) and Orthoperus scutellaris
  (Corylophidae) are new, both from families the Lineage-A filter excluded.
- Fungal side spread thin: 1,800 of 2,071 involved ASVs link to exactly 1 insect
  taxon (max 4). Top fungal genera among hits are a saprotroph / plant-surface-yeast
  cast (Aureobasidium 117, Dothiora 83, Cytospora 66, Taphrina 42) -- shared
  seasonal bloom, not beetle-specific symbiosis.

- Residualized follow-up (each insect predictor stripped of site + factor(date),
  the strictest seasonal control 6 sampling rounds allow; fungal-side model keeps
  site + lure + factor(date) + insect_taxon_resid): 646 pairs at q<0.10 within, and
  ONLY Pityogenes hopkinsi retains partners (646, overlap 605 = 38% of its main-grid
  hits survive). Asemum, Xyleborinus, Orthoperus, Hylastes all collapse to 0 --
  their fungal "associations" were purely a shared early-season timing artifact.
  Lineage A had the identical outcome (only P. hopkinsi survived, 717/1,593 = 45%).
  The residualized script's flagged-taxon list is data-driven here (top 5 by
  main-grid hit count) rather than the hard-coded 3-Scolytinae list the Lineage-A
  parent used.

Bottom line: the all-families insect table does not change the pairwise-association
picture -- a few early-season bark beetles, above all P. hopkinsi, drive everything;
nothing survives global correction; only P. hopkinsi's signal partially survives a
strict seasonal control; and the two extra taxa the broader table surfaces are
seasonal co-occurrence, not specific links. See
lineage_B_prevalence_filtered_insect_top_level_interpretation.md sec. 6.

Output: data/ and figures/compare_insects_fungi_pairwise_taxa_prevalence_filtered_
insect/ (full grid + significant hits + per-insect / per-fungal / family summaries
for the main run; *.date_site_residualized.csv + original_vs_residualized_
comparison.{csv,png} for the follow-up).

### #################################################################################
### 2026-09-01 -- Lineage B port: fig4/fig5 taxonomic breakdowns, with linear +
### quadratic date terms combined onto one set of axes as dodged bars
### #################################################################################

presentation_items_prevalence_filtered_insect/fig4_association_taxonomic_breakdown.R
and fig5_family_breakdown_by_order.R. Lineage-A's fig4/fig5 broke down the LINEAR-
term hits only (fig4: insect~date by subfamily, fungal~date by order/class,
fungal~insect_PCoA1 by order/class; fig5: family-within-order detail for the two
fungal panels). Ported forward with the quadratic date term (added to the direct
taxon~date screens on 2026-08-28, same as fig3) included alongside the linear term,
per user request, on the SAME taxonomic-group axes rather than in separate panels --
each group gets up to 4 dodged (side-by-side) bars: existing blue (#0072B2, later
season/linear t>0) and vermillion (#D55E00, earlier season/linear t<0), plus two new
Okabe-Ito colors for the quadratic term, sky blue (#56B4E9, dip/t>0) and orange
(#E69F00, hump/t<0). This only applies where linear and quadratic are two terms of
the SAME test (insect~date, fungal~date) -- fig4 panels c/d (fungal~insect_PCoA1,
fungal~insect_PCoA2) are different predictor axes, not a linear/quadratic pair of one
test (no quadratic variant of a PCoA-axis coefficient exists in this project), so
they stay single-term, using the same factor(date)-adjusted CoCA files and direction
framing fig3 established (panel c keeps PCoA1's guild-level Scolytinae-vs-
fungivorous-Latridiidae/Corticariinae language; panel d uses generic PCoA2+/- labels
since PCoA2 is a within-Scolytinae species-composition axis, not a guild split, and
its subtitle carries the same instability caveat as fig3 -- 9/200 factor(date)-
adjusted, 6/200 unadjusted, 0/200 quadratic-date).

- fig4 is **3 columns**: a insect~date (linear+quadratic combined, dodged, by
  subfamily) and b fungal~date (linear+quadratic combined, dodged, by order/class)
  each get a full-height column; c fungal~PCoA1 and d fungal~PCoA2 (both single-term,
  factor(date)-adjusted) stack in a third, narrower column, the same layout idea
  fig5 uses for its busy vs. narrow panels.
- The top-n fold cutoff for a/b (top 15 + "Other") ranks groups by the COMBINED
  linear+quadratic hit count, per user request -- a taxon significant in both terms
  counts twice toward its group's total. Panel a needed this folding for the first
  time in this lineage: Lineage A's insect~date panel never needed it (only 4
  subfamilies, family-filtered table), but Lineage B spans 53 families/74
  subfamilies. Subfamily fallback-to-Family logic was generalized to match:
  insect_taxa_date_association.csv carries baked-in Family/Subfamily/Order columns
  (no join needed, unlike Lineage A), but Subfamily uses two upstream conventions --
  literal "NA" (unresolved, Family also "NA") and "NoSubfamily_<x>" (Family resolved,
  no subfamily rank) -- both collapse to Family, extending Lineage A's original
  is.na(Subfamily)-only fallback. Panels a/b stay unfaceted at the Order level for
  insects specifically (b IS faceted by fungal Class) -- Lineage B insect taxa are
  92% Coleoptera (140/153), so an order facet would produce one dominant facet and
  several near-empty ones, unlike the fungal panels' more even Class split.
- fig5 is **3 panels**: a = family-within-order breakdown of fig4 panel a/b's
  fungal~date COMBINED result (same dodge design and combined-count fold rule as
  fig4), b = fig4 panel c's fungal~PCoA1 breakdown (single-term, unchanged design),
  c = fig4 panel d's fungal~PCoA2 breakdown (single-term). Panel c (PCoA2, only 9
  taxa across 4 orders) is small but genuinely informative: Helotiales alone splits
  into 3 different families (incl. Tympanidaceae) that fig4 panel d's order-level bar
  collapses into one, and Diaporthales resolves to Valsaceae specifically -- the
  "Tympanis-as-a-group" / "Valsaceae/Cytospora-as-a-group" framing the interpretation
  doc already recommends for this axis (sec. 6). Panel a's layout uses a 2-column
  arrangeGrob with panel a (the busiest, combined-term) wide on the left and b/c
  stacked narrower on the right.
- Every kept group gets a complete 4-category block via `tidyr::complete(...,
  fill = list(n = 0))` even where a category has zero hits, so `position_dodge2()`
  reserves the same 4 dodge slots for every group/facet row -- without this, groups
  with fewer non-zero categories would dodge to different widths/offsets than groups
  with all 4, breaking cross-row alignment. Zero-count bars render at zero height
  with a blank (not "0") text label.
- Both fig4/fig5 hit the same gtable-alignment limitation fig5-Lineage-A already
  found: the gtable rbind/cbind technique fig3 uses (column = rbind of its row
  panels, then cbind the columns) requires every column to have the same total row
  count, and here panels are facet_grid'd with DIFFERING numbers of Class/Order
  facets -- so cbind fails. Used gridExtra::arrangeGrob instead (Lineage A's
  original approach for these two figures), which places each panel as an opaque
  grob and sidesteps the row-count requirement.
- **Bug found and fixed**: the combined-panel builder initially computed `group =
  factor(group, levels = ...)` BEFORE the `left_join(group_class, ...)` used for
  Class faceting. `dplyr::left_join` silently coerces a factor join column back to
  character when the other side is character, which drops the factor's levels -- so
  `levels(counts$group)` came back NULL, the label lookup vector built from it was
  empty, and every faceted row's axis label rendered as literal "NA" (correct bar
  lengths/positions, wrong text -- isolated to the facet_by_class path since the
  unfaceted insect panel was unaffected). Fixed by moving the join before the factor
  conversion, matching the order the working single-term builder already used, and
  applied consistently in fig5's combined builder (join Order/Family back onto
  `fam_id` while still character, convert to factor after).
- A couple of subtitle-clipping issues (text running past a narrow column's width,
  e.g. panel d's instability caveat) were fixed with explicit "\n" line breaks at
  natural clause boundaries, verified by re-rendering.

Dynamic counts on this run -- fig4: a 92/153 taxa sig in >=1 term (69 linear + 38
quadratic = 107 hit-instances, 55 subfamilies), b 2490/6342 (2186 linear + 707
quadratic = 2893 hit-instances, 102 orders), c 17/200 (PCoA1 factor(date), 10
orders), d 9/200 (PCoA2 factor(date), 4 orders, unstable across date controls --
6/200 unadjusted, 0/200 quadratic-date). fig5: a 1467 taxa/14 orders (1280 linear +
467 quadratic = 1747 hit-instances), b 13 taxa/9 orders, c 9 taxa/4 orders.

Both scripts run standalone end-to-end (verified with Rscript, output inspected
visually, spot-checked e.g. fig4 panel a's Scolytinae row -- 12 linear hits split 10
earlier/2 later, 6 quadratic hits split 3 hump/3 dip -- against the raw CSV) and
write only into figures/ and data/presentation_items_prevalence_filtered_insect/ --
no Lineage-A file touched.

Output: figures/presentation_items_prevalence_filtered_insect/fig4_association_
taxonomic_breakdown.{png,pdf}, fig5_family_breakdown_by_order.{png,pdf}; per-group
counts as data/presentation_items_prevalence_filtered_insect/fig4_taxonomic_
breakdown.{insect_date_combined.by_subfamily,fungal_date_combined.by_order,
fungal_insectPCoA1_factordate.by_order,fungal_insectPCoA2_factordate.by_order}.csv
and fig5_family_breakdown_by_order.{fungal_date_combined,fungal_insectPCoA1_
factordate,fungal_insectPCoA2_factordate}.csv.

### #################################################################################
### 2026-09-01 -- Lineage B port: fig6 (responsive-taxa abundance/count by site),
### using a combined linear-or-quadratic date definition and combined PCoA1-or-PCoA2
### insect definition
### #################################################################################

presentation_items_prevalence_filtered_insect/fig6_responsive_taxa_abundance_and_
counts_by_site.R. Same 2x2 layout as Lineage A's fig6 (a/b mean relative sequence
abundance by site, c/d ASV count/richness by site; columns = date-associated vs.
insect-associated), rebuilt on the Lineage-B matched dataset (all-families, >=5-
sample prevalence-filtered insect table, same construction as fig2's Lineage-B port
-- 69 matched samples, 15,422 Kingdom==Fungi ASVs) and, per user request, on a
COMBINED significance definition for each association type rather than Lineage A's
single-test one:

- "Date-associated" = q<0.10 on EITHER the linear or quadratic date term
  (fungal_taxa_date_association.all_taxa.csv, lineage-invariant -- same file Lineage
  A's fig6 reads, but Lineage A only used the linear term since the quadratic term
  didn't exist yet when that figure was built). 2490/15422 taxa (2186 linear + 707
  quadratic, matching fig4 panel b's numbers exactly, since it's the same union).
- "Insect-associated" = q<0.10 on EITHER insect_PCoA1 or insect_PCoA2, factor(date)-
  adjusted (fungal_insect_association_results.PCoA{1,2}_plus_factordate.csv -- same
  deseasonalized files fig3/fig4/fig5 use, replacing Lineage A's single permissive
  no-lure PCoA1 test). Both files test the identical top-200 CoCA-loading candidate
  set, and PCoA1's 17 hits and PCoA2's 9 hits don't overlap at all, so the union is a
  clean 26/15422 taxa, no double-counting to worry about.
- Site labels use the raw `Site` metadata column (Durham/Pease/Manchester Cedar
  Swamp/Manchester Airport), matching fig2's Lineage-B convention, rather than
  Lineage A's fig6 "Pease" -> "Pease Airport" display rename.

Panel c and panel d's total bar heights are identical per site by construction (both
partition the SAME 15,422-taxon presence data into "associated"/"not", just under a
different definition of "associated") -- confirmed visually and expected, not a
scaling artifact: Durham 2401+7443 = 20+9824 = 9844; Pease 2363+6871 = 26+9208 = 9234;
Manchester Cedar Swamp 2356+6779 = 11+9124 = 9135; Manchester Airport 2126+5474 =
5+7595 = 7600.

Headline pattern matches Lineage A qualitatively: a small taxon share (2490/15422 =
16.2% date-associated; 26/15422 = 0.17% insect-associated) carries disproportionate
sequence abundance and richness. Exact per-site numbers underlying each panel
(site order matches the figure's site_levels; abundance = mean relative sequence
abundance across that site's samples):

**Panel a -- relative abundance, date-associated:**

| Site | Date-associated | Not date-associated |
|---|---|---|
| Durham | 86.9108% | 13.0892% |
| Pease | 88.3081% | 11.6919% |
| Manchester Cedar Swamp | 88.2896% | 11.7104% |
| Manchester Airport | 81.4339% | 18.5661% |

**Panel b -- relative abundance, insect-associated:**

| Site | Insect-associated | Not insect-associated |
|---|---|---|
| Durham | 0.183432% | 99.8381% |
| Pease | 1.03079% | 99.2117% |
| Manchester Cedar Swamp | 0.00944035% | 99.9958% |
| Manchester Airport | 0.00509525% | 99.9988% |

**Panel c -- ASV count, date-associated:**

| Site | Date-associated | Not date-associated | Total |
|---|---|---|---|
| Durham | 2401 | 7443 | 9844 |
| Pease | 2363 | 6871 | 9234 |
| Manchester Cedar Swamp | 2356 | 6779 | 9135 |
| Manchester Airport | 2126 | 5474 | 7600 |

**Panel d -- ASV count, insect-associated:**

| Site | Insect-associated | Not insect-associated | Total |
|---|---|---|---|
| Durham | 20 | 9824 | 9844 |
| Pease | 26 | 9208 | 9234 |
| Manchester Cedar Swamp | 11 | 9124 | 9135 |
| Manchester Airport | 5 | 7595 | 7600 |

The insect-associated slice is small enough (26 taxa total) to be barely visible in
panel b and only a thin sliver in panel d. Pease stands out numerically as carrying
the most of both its abundance (1.03%, ~5.6x Durham's 0.18%, the next highest) and
count (26 ASVs, also the site max) in the insect-associated category, but all four
sites carry the full Ethanol/Alpha-pinene_EtOH/Ips lure set (no site is lure-
restricted), so this is a site-level pattern, not a lure-design artifact -- not
further investigated here. Overall the insect-community-association signal is far
weaker and narrower than the seasonal one under Lineage B, matching the sec. 5/6
finding.

Output: figures/presentation_items_prevalence_filtered_insect/fig6_responsive_taxa_
abundance_and_counts_by_site.{png,pdf}.

### #################################################################################
### 2026-09-01 -- Lineage B port: fig7/fig8 (trait breakdowns), same combined-term
### dodged-bar design as fig4/fig5 -- Lineage B presentation-figure set now complete
### #################################################################################

presentation_items_prevalence_filtered_insect/fig7_association_trait_breakdown.R and
fig8_association_trait_breakdown.R. Same FungalTraits/Polme et al. 2020 genus-join
trophic-mode breakdown as Lineage A (fig7 faceted by growth_form, fig8 by taxonomic
Class), rebuilt per user request with the same design already established for fig4/
fig5: the linear and quadratic date terms share one set of trait-group axes as 4
dodged bars (blue/vermillion linear later/earlier season, sky-blue/orange quadratic
dip/hump), and insect_PCoA1/PCoA2 stay as separate single-term panels (factor(date)-
adjusted, replacing Lineage A's permissive no-lure PCoA1-only test) rather than
merged, since they're different predictor axes with different direction semantics
(PCoA1 guild-level, PCoA2 within-Scolytinae composition -- same reasoning as fig4c/d).

- Both figures are now **2 columns**: panel a (fungal trait vs. date, combined dodge)
  full-height on the left; panels b (PCoA1) and c (PCoA2), both single-term, stacked
  narrower on the right -- the same a-wide/b+c-stacked layout fig5 uses, since fig7/
  fig8 have only one "combined-term" panel (trait vs. date) rather than fig4's two
  (insect~date AND fungal~date both combined), unlike fig4's 3-column layout.
- fig7 (growth_form facets): panel a folds to the top 25 trait groups + "Other",
  ranked by COMBINED linear+quadratic hit count (2490 taxa, 2186 linear + 707
  quadratic = 2893 hit-instances, 42 trait groups, matching fig4 panel b's date
  numbers exactly since it's the same taxon set under a different grouping).
  "Unclassified genus" / "No FungalTraits match" kept as their own pseudo-facet bars,
  unchanged from Lineage A.
- fig8 (Class facets, trait_group x Class is a genuine cross-tab, not a 1:1 nesting
  like fig4's Order-within-Class): panel a's two Lineage-A thresholds -- drop a Class
  facet with too few total hits, fold a (lifestyle, Class) cell with too few hits
  into that facet's "Other" -- now use the COMBINED hit count via
  `min_facet_hits_date = 6` (unchanged value from Lineage A). This run: 10 of 25
  classes dropped (30 hit-instances, all <6), 29 lifestyle-within-class cells folded
  (65 hit-instances). "Unclassified genus" / "No FungalTraits match" are DROPPED
  entirely here (not their own bars) -- unchanged from Lineage A, opposite of fig7's
  treatment, same rationale as before (thin Class-by-Class splitting vs. one legible
  bucket).
- Implementation: both combined-panel builders reuse fig4's join-before-factor
  ordering (facet/axis_id joined while columns are still character, converted to
  factor only in the final mutate) to avoid the "NA" axis-label bug fig4 hit and
  fixed on its first pass -- both scripts rendered correctly on the first run this
  time, no follow-up fix needed. fig8's builder also carries forward fig4/fig7's
  zero-fill-via-complete()/crossing() trick so every (group,facet) row gets all 4
  dodge slots reserved even when a category has zero hits, keeping cross-row/cross-
  facet alignment consistent.

Dynamic counts on this run -- fig7: a 2490/6342 (2186 linear + 707 quadratic, 42
trait groups), b 17/200 (PCoA1 factor(date), 7 trait groups), c 9/200 (PCoA2
factor(date), 3 trait groups). fig8 (FungalTraits-matched subset only, "Unclassified
genus"/"No FungalTraits match" excluded): a 1747 taxa (1546 linear + 478 quadratic =
2024 hit-instances) across 20 surviving lifestyles/15 surviving classes, b 9 taxa
across 4 lifestyles/5 classes, c 7 taxa across 2 lifestyles/3 classes (2 of PCoA2's 9
hits fall in the excluded Unclassified-genus/No-FungalTraits-match buckets, hence 7
not 9).

This completes the Lineage-B presentation-figure set (fig1 lineage-invariant, no
port needed; fig2-fig8 all built) -- see project_organization.md's "Known gaps"
section, now updated to reflect Lineage B as feature-complete with Lineage A other
than the PCoA1_date_axis/PCoA3_lure_axis/compare_selection_strategies variants and
alpha diversity (both deliberately out of scope, see that section).

Output: figures/presentation_items_prevalence_filtered_insect/fig7_association_
trait_breakdown.{png,pdf}, fig8_association_trait_breakdown.{png,pdf}; per-group
counts as data/presentation_items_prevalence_filtered_insect/fig{7,8}_trait_
breakdown.{fungal_date_combined,fungal_insectPCoA1_factordate,
fungal_insectPCoA2_factordate}.csv.

### #################################################################################
### 2026-09-01 -- Lineage C presentation: fig9 (alpha diversity by lure) + fig10
### (insect~fungal correlation), plus a cubic date term and a lure-interaction test
### #################################################################################

New presentation_items/fig9_alpha_diversity_by_lure.R and fig10_insect_fungal_
alpha_diversity_correlation.R -- the first two presentation figures built for
Lineage C (full insect table, alpha-diversity-only; see project_organization.md).
Both live directly in presentation_items/ alongside the Lineage A fig1-fig8 (not a
dedicated Lineage C folder) since alpha diversity has no Lineage B counterpart to
disambiguate against.

**fig9**: insect (panel a) + fungal (panel b) alpha diversity vs. collection date,
two side-by-side facet_grid blocks (rows = Shannon/Simpson dominance/richness,
columns = lure), points colored by site. Row (metric) strip moved to panel a's left
side via facet_grid(switch="y") with its grey background blanked, since panel b's
own row strip would otherwise duplicate the label. Trend lines iterated loess (span
0.75, then widened to span 1 for a smoother curve) before being swapped to a cubic
OLS fit (y ~ poly(date, 3)) once the underlying stats (below) confirmed a matching
permutation-tested cubic date term, so the visual trend line and the formal test now
agree on functional form.

**Stats (compare_insect_fungi/insect_fungal_alpha_diversity.full_insect_table.r)**:
the date-effect test (sec. 6) was extended from linear-only to linear+quadratic+
cubic in ONE model (date centered, same `date_c`/`I(date_c^2)` pattern already
established in fungal_community_seasonality.r and insect_fungal_coca_analysis.
general_workflow.r, extended one order further here specifically because fig9's
trend lines showed more structure -- e.g. an S-curve in insect Shannon diversity
under Ethanol -- than a single hump/dip could capture; this cubic extension is a
targeted addition for this alpha-diversity check, not (yet) a project-wide
convention). `alpha_diversity_covariate_tests.csv` now carries 3 date rows
(date_linear/date_quadratic/date_cubic) per community x metric instead of 1; new
`alpha_diversity_date_shape.csv` classifies each community x metric's shape (cubic
significant takes priority over quadratic, which takes priority over linear, same
cascading rule fungal_community_seasonality.r uses one level down). Result: 4 of the
6 community x metric date tests show a significant cubic term (insect richness;
fungal Shannon diversity, Simpson dominance, richness) -- only insect Shannon
diversity and Simpson dominance show no significant higher-order structure.

**fig10, first pass**: built faceted by lure, mirroring fig9's metric x lure grid,
one point per sample as insect value (x) vs. fungal value (y), colored by date. Hit
and fixed a real ggplot2 semantics bug along the way: plain
`facet_grid(metric ~ lure, scales="free")` ties `free_x` to COLUMN position and
`free_y` to ROW position, NOT to an arbitrary facet variable -- verified with a
small synthetic facet_grid test before fixing. With metric as the row variable, all
three metrics sharing a lure column were forced onto the SAME x-range, squashing
Shannon/Simpson diversity into a sliver next to richness's much wider scale. First
fix was a manual 3-block `gtable::rbind` stack (one block per metric, each
facet_grid'd by lure only); replaced afterward with `ggh4x::facet_grid2(scales=
"free", independent="all")` (ggh4x installed for this) once available, since it
gives every one of the 9 panels a genuinely independent x AND y range while keeping
the normal facet_grid row/column strip layout -- shorter and no gtable workaround.

**fig10, collapsed to pooled**: per user's suspicion that splitting by lure wasn't
well-motivated, added a new sec. 5b to insect_fungal_alpha_diversity.full_insect_
table.r testing whether lure modifies the insect~fungal relationship: `fungal ~
site+date+insect*lure` vs. the same model without the interaction, permutation
F-test on the interaction term (999 perms, lure shuffled among traps within site --
same scheme test_term_lure already uses in sec. 6). NO metric showed a significant
interaction (Shannon diversity p_perm=0.913, F=0.148; Simpson dominance p_perm=
0.802, F=0.341; richness p_perm=0.953, F=0.017) -- so fig10 was collapsed to a
single, non-lure-faceted row of 3 panels (one per metric, facet_wrap with per-panel
free scales, no ggh4x needed once lure was dropped as a second facet dimension)
showing the pooled Spearman correlation already computed in sec. 5 (`alpha_
diversity_cross_community_correlation.csv`): Shannon diversity rho=-0.26, p=0.03
(the only significant one); Simpson dominance rho=-0.12, p=0.33; richness rho=-0.12,
p=0.34. The stale by-lure correlation CSV that the since-removed intermediate
version of fig10 had written was deleted.

Output: figures/presentation_items/fig9_alpha_diversity_by_lure.{png,pdf},
fig10_insect_fungal_alpha_diversity_correlation.{png,pdf}; data/compare_insects_
fungi_alpha_diversity_full_insect_table/alpha_diversity_covariate_tests.csv
(updated), alpha_diversity_date_shape.csv (new), alpha_diversity_insect_fungal_
interaction_by_lure.csv (new).

### #################################################################################
### 2026-09-02 -- "Peak timing" convention adopted for fig3, then ported to
### fig4/fig5/fig7/fig8: linear term takes priority over quadratic
### #################################################################################

After examining example taxa from each linear x quadratic significance combination
(new supplemental script `figS1_date_shape_bin_examples.R`, prompted by a question
about whether fig4's panel-b hit-instance counts were double-countable -- they are:
a taxon significant in both the linear and quadratic date terms contributed to 2
separate bars under the original combined-dodge design), the non-linearity itself
wasn't judged particularly interesting on its own. What matters for interpretation
is WHEN a taxon peaks. Adopted reading (see `project_organization.md`'s "Peak-timing
convention"): **the linear term takes priority over the quadratic term** whenever
the linear term is significant -- the OPPOSITE priority from the `shape` column
`fungal_community_seasonality.r`/`insect_taxa_date_association.prevalence_filtered_
insect.r` write (which gives the quadratic term priority). A taxon now gets exactly
ONE of 4 mutually-exclusive categories: peaks early (significant negative linear,
regardless of quadratic), peaks late (significant positive linear, regardless of
quadratic), peaks mid-season (no significant linear, significant quadratic hump),
bimodal/early+late (no significant linear, significant quadratic dip).

- **fig3** (`fig3_volcano_plots.R`, first to adopt this): panels a/b (linear-term
  volcanoes) now color points by whether the SAME taxon also carries a significant
  quadratic term (4-way combo: grey neither / blue linear-only / orange quadratic-
  only / purple both), replacing the old linear-significance-only coloring. Panels
  d/e (quadratic-term volcanoes) now EXCLUDE any taxon with a significant linear
  term first, since under the priority rule those taxa are already classified by
  their linear direction. Point z-order fixed at the same time so the sparser
  colored categories always draw on top of the grey majority (they were being
  buried by row order before).
- **fig4/fig5**: each taxon's up-to-2 hit-instances collapsed to exactly 1 dodged
  bar (Peaks late / Peaks early / Peaks mid-season / Bimodal), group/family ranking
  for the top-n fold switched from hit-instance sums to plain taxon counts. Per user
  request, the 4 bars were then changed to all run the same direction (0 -> positive
  -- a diverging left/right layout doesn't map onto 4 unordered categories the way
  it did onto a signed t-statistic), ordered top-to-bottom within each group as
  early / bimodal / mid-season / late.
- **fig7/fig8**: same design ported over, including the same-direction/early-
  bimodal-mid-late bar layout from the start (no separate diverging-bars pass).
  fig8's two `min_facet_hits_date` fold thresholds (drop a sparse Class facet; fold
  a sparse lifestyle-within-Class cell into "Other") switched from hit-instance sums
  to plain taxon counts, same as fig4/fig5's top-n fold.

**CORRECTION -- a real bug, caught before either affected figure was reported as
final**: implementing fig4's bar-reorder (early/bimodal/mid-season/late) required
reordering `cat_levels`, and the `case_when()` blocks that assign each taxon's
category referenced `cat_levels[1]`/`[2]`/`[3]`/`[4]` by POSITION rather than by
name. Reordering `cat_levels` silently swapped the "Peaks early season" and "Peaks
mid-season" labels in fig4 and fig5's generated figures/CSVs (Diaporthales, e.g.,
briefly showed 89 taxa as "mid-season" that are actually early-season -- caught
because that reading didn't match Diaporthales/*Cytospora*'s well-established early-
season profile from earlier in this file). Fixed by switching every `case_when()` in
fig4/fig5/fig7/fig8 to assign literal category strings ("Peaks late season", etc.)
instead of indexing into `cat_levels` -- makes the category assignment immune to
`cat_levels`' display-order changes by construction, rather than relying on the two
staying manually in sync. Caught and fixed within the same session before any
downstream interpretation-doc numbers were drawn from the buggy output (§7 of
`lineage_B_prevalence_filtered_insect_top_level_interpretation.md` was computed
independently from a correct one-off script, not from the buggy plot code, so it
did not need correcting).

Output: figures/presentation_items_prevalence_filtered_insect/fig{3,4,5,7,8}_*.
{png,pdf}, figS1_date_shape_bin_examples.{png,pdf}; data/presentation_items_
prevalence_filtered_insect/fig{4,5,7,8}_*.csv (all regenerated), figS1_date_shape_
bin_counts.csv (new).

### #################################################################################
### 2026-09-08 -- New analysis: asymptotic (Chao-extrapolated) richness via iNEXT,
### on the Lineage C full insect table; fed into the same cross-community
### correlation/covariate pipeline as insect_fungal_alpha_diversity.full_insect_table.r
### #################################################################################

New script `compare_insect_fungi/insect_fungal_asymptotic_richness_iNEXT.
full_insect_table.r`, per user request to estimate extrapolated species richness
with the `iNEXT` package, using the same full (all-families, singleton-filtered
only, no prevalence filter) insect table as Lineage C -- asymptotic estimation is
itself an alpha-diversity method (estimating how many taxa exist, same quantity
richness/Shannon/Simpson summarize), so the existing "don't prevalence-filter for
alpha-diversity" rule applies.

Per-sample granularity (user's choice, offered against site/site-x-lure/whole-
dataset pooling alternatives): iNEXT's Chao-type estimators run on each of the 69
matched insect/fungal samples individually, same unit as the existing observed-
richness pipeline.

Two methodological departures from insect_fungal_alpha_diversity.full_insect_
table.r, both explained in the new script's header:
- **Fungal table left unrarefied** (raw Kingdom==Fungi-filtered ASV counts, no
  rarefy-and-average). Rarefying first would throw away exactly the deep-sample
  information the Chao estimator needs to extrapolate the rare tail, and would
  artificially cap "observed" richness at the rarefaction depth -- asymptotic
  estimation already IS the principled correction for uneven sequencing depth that
  rarefaction exists to provide elsewhere in this project.
- **`ChaoRichness()`/`ChaoShannon()`/`ChaoSimpson()` used instead of `iNEXT()`**.
  `iNEXT()`'s full sample-size-based rarefaction/extrapolation curve computation
  (default 40 knots + bootstrap) did not finish in 120s on even a single large
  fungal sample (~800k reads); the lighter Chao*() point/SE/CI estimators (same
  underlying formulas, verified to reproduce `iNEXT()`'s `$AsyEst` point estimates
  exactly on a test subset) run in ~1-2s each, making all ~140 samples (69 insect +
  69 fungal) tractable (~80s total runtime). Two `iNEXT` package quirks worked
  around: `ChaoShannon()`/`ChaoSimpson()` need `transform = TRUE` to return Hill-
  number ("effective species count") form rather than raw entropy/Gini-Simpson
  index; `ChaoShannon()`'s SE column is named `"Est_s.e"` (no period) vs. the other
  two functions' `"Est_s.e."` (with period) -- extracted by column position, not
  name, to avoid a silent NA.

Naming note: "simpson" in this script's output is Hill number q=2 (inverse Simpson
DIVERSITY, higher = more even), the OPPOSITE direction from "simpson_dominance"
(higher = more dominated) used in `insect_fungal_alpha_diversity*.r` elsewhere in
this project -- don't compare the two scripts' "simpson" columns directly.

Per user follow-up request ("feed these estimates into the same type of
statistical analysis pipeline as the raw estimates to look for correlations
between insect and fungal richness"), the script mirrors insect_fungal_alpha_
diversity.full_insect_table.r's sections 5/5b/6: cross-community Spearman
correlation (reported for BOTH observed and asymptotic values side by side, so
the effect of extrapolation on the conclusion is visible directly), the lure-
modifies-the-relationship interaction test, and the site/lure/date permutation
covariate tests + linear/quadratic/cubic date-shape classification -- the latter
two run on the asymptotic values only (the new content this script adds; observed-
value covariate significance already exists under different metric definitions in
the original script, so wasn't re-derived a second time here).

Results (69 matched samples): median sample coverage (fraction of the true
community captured, per `iNEXT::DataInfo()`) was 0.809 for insect trap catches vs.
0.995 for fungal ASVs -- insect sampling is far from its asymptote (trap catches
median ~93 individuals across up to dozens of taxa per sample) while fungal
amplicon sequencing is close to saturated per sample, consistent with the very
different sampling mechanisms. Cross-community correlation on ASYMPTOTIC values:
richness rho=0.008, p=0.95 (no relationship, same conclusion as observed:
rho=0.076, p=0.53); Shannon rho=-0.251, p=0.038 (weak negative, same conclusion as
observed: rho=-0.253, p=0.036); Simpson rho=-0.075, p=0.54 (n.s., same as observed:
rho=-0.116, p=0.34) -- extrapolating to the asymptote did not change any of the
three conclusions. Lure does not modify the insect~fungal relationship for any
asymptotic metric (all p_perm >= 0.075). Within-community covariate tests on
asymptotic values: insect richness shows a significant quadratic (hump, mid-season
peak) date term (p=0.003); insect Shannon and fungal richness each show a
significant site effect (p=0.022, p=0.015); fungal Shannon and Simpson both show
significant linear (early-season) + quadratic + (Shannon only) cubic date terms,
classified "complex" under the cubic-takes-priority shape rule.

Output: data/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table/
insect_asymptotic_diversity.csv, fungal_asymptotic_diversity.csv,
asymptotic_diversity_cross_community_correlation.csv,
asymptotic_diversity_insect_fungal_interaction_by_lure.csv,
asymptotic_diversity_covariate_tests.csv, asymptotic_diversity_date_shape.csv (all
new); figures/compare_insects_fungi_asymptotic_richness_iNEXT_full_insect_table/
asymptotic_diversity_by_community.png, asymptotic_richness_observed_vs_
extrapolated.png, asymptotic_diversity_by_site_lure.png, asymptotic_diversity_by_
date.png, asymptotic_diversity_cross_community_correlation.png (all new).

### #################################################################################
### 2026-09-08 (same day) -- fig9/fig10 reproduced as fig11/fig12 on the
### asymptotic (iNEXT) diversity estimates
### #################################################################################

Per user request, `presentation_items/fig9_alpha_diversity_by_lure.R` and
`fig10_insect_fungal_alpha_diversity_correlation.R` reproduced one-for-one as
`fig11_asymptotic_diversity_by_lure.R` and `fig12_insect_fungal_asymptotic_
diversity_correlation.R`, reading the asymptotic-diversity CSVs from the new
iNEXT script above instead of the observed-diversity CSVs fig9/fig10 read. Same
panel layout, trend lines, and color scheme in both cases -- no new design, just
the metric swapped (and metric labels updated: "Asymptotic species richness",
"Shannon/Simpson diversity (Hill q=1/q=2)" in place of "Richness (observed
taxa)"/"Shannon diversity"/"Simpson dominance"; fig11/fig12's "Simpson" is Hill
q=2 diversity, opposite direction from fig9/fig10's Simpson dominance).

One fix needed: fig11's row-strip metric labels (switch="y", rotated text along
the left edge of panel a) clipped with the full "Asymptotic Shannon diversity
(Hill q=1)"-style wording fig12 uses without issue (fig12's strip labels are
horizontal top strips, unconstrained) -- shortened to "Shannon div. (q=1)",
"Simpson div. (q=2)", "Asymptotic richness" for fig11 only; fig12 kept the full
wording.

fig12's in-panel Spearman labels match the iNEXT script's asymptotic-value
correlation numbers already reported above (richness rho=0.01 p=0.95, Shannon
rho=-0.25 p=0.038, Simpson rho=-0.08 p=0.54) -- read directly from `asymptotic_
diversity_cross_community_correlation.csv`'s `value_type == "asymptotic"` rows,
not recomputed.

Output: figures/presentation_items/fig11_asymptotic_diversity_by_lure.{png,pdf},
fig12_insect_fungal_asymptotic_diversity_correlation.{png,pdf} (both new).

**Follow-up, same day**: per user request, fig12 now also shapes each point by
lure (filled circle = Ethanol, filled triangle = Alpha-pinene_EtOH, filled square
= Ips -- shapes 21/24/22, chosen so fill continues to carry the date gradient
alongside shape) in addition to the existing date-fill coloring; fig12 still pools
across lure for the trend line/Spearman annotation (unchanged rationale -- no
lure interaction, see above), the shape is descriptive only. fig9/fig11 already
facet by lure so were not touched.

Output: figures/presentation_items/fig12_insect_fungal_asymptotic_diversity_
correlation.{png,pdf} (updated).

**Follow-up, same day**: same shape-by-lure change ported back to fig10
(`fig10_insect_fungal_alpha_diversity_correlation.R`) for consistency between
the observed- and asymptotic-diversity versions of this figure -- identical
shapes/legend (circle=Ethanol, triangle=Alpha-pinene_EtOH, square=Ips), same
pooled-trend-line/Spearman-annotation rationale unchanged.

Output: figures/presentation_items/fig10_insect_fungal_alpha_diversity_
correlation.{png,pdf} (updated).

### #################################################################################
### 2026-09-08 (same day) -- Asymptotic-diversity summary added to lineage_C_
### alpha_div_full_insect_table_top_level_interpretation.md (§7-12)
### #################################################################################

Per user request, extended `lineage_C_alpha_div_full_insect_table_top_level_
interpretation.md` -- previously scoped to only `insect_fungal_alpha_diversity.
full_insect_table.r` and fig9/fig10 (observed diversity) -- with 6 new sections
(§7-12) covering the asymptotic-diversity work above: value ranges + sample
coverage (§7), observed-vs-asymptotic cross-community correlation (§8, unchanged
conclusion for all 3 metrics), the lure-interaction test (§9), within-community
covariate tests (§10), date-shape classification (§11), and fig11/fig12 (§12).
Bottom-line list extended from 6 to 10 points.

Two things worth surfacing here since they're genuine findings, not just
restated numbers already reported above:
- **§9/bottom-line point 10**: the asymptotic Shannon lure-interaction test
  (F=0.657, p_perm=0.075) crosses this project's uncorrected p<0.10 convention
  for small fixed covariate/interaction test sets -- unlike every other lure-
  interaction test in the document (observed-value §3: all p>=0.80; asymptotic
  richness/Simpson: p=0.888/0.697). Flagged explicitly as a single borderline
  result worth watching, not treated as a settled finding -- fig12 still shows
  the pooled (non-lure-split) relationship.
- **§10/bottom-line point 9**: insect richness is the one metric where "observed"
  is identical between the two scripts (same raw table, no transform), giving a
  clean before/after-extrapolation comparison. The strong lure effect (F=14.66,
  p=0.003) and cubic date shape (p=0.001) found on observed richness (§4)
  attenuate to no significant lure effect (p=0.447) and a plain quadratic hump
  once extrapolated (§10-11) -- suggesting part of the observed lure-richness
  effect may reflect unequal sampling completeness across lures rather than a
  genuinely different taxon pool per lure. Shannon/Simpson comparisons between
  the two scripts are explicitly flagged as NOT clean (different metric
  definitions -- Hill number vs. raw entropy/dominance, unrarefied vs. rarefied
  fungal table) to avoid over-reading superficially similar/different p-values.

Also updated project_organization.md's `lineage_C_alpha_div_full_insect_table_
top_level_interpretation.md` bullet to describe the doc's now-6-script/4-figure
scope (was 1 script/2 figures).

Output: lineage_C_alpha_div_full_insect_table_top_level_interpretation.md
(extended, not regenerated -- no new data/figures, all numbers pulled from
already-existing CSVs); project_organization.md (updated).
