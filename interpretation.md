Family/genus composition of the 96 significant CoCA hits:
- Valsaceae is the top family (10 hits), driven mostly by Cytospora (7 ASVs) — a canker-causing fungal genus (Diaporthales) commonly associated with bark beetle/weevil vectoring in trees. Given the insect data here is Curculionidae (weevils, including Scolytinae bark beetles) and Latridiidae, this is a biologically very plausible hit, not just statistical noise.
- Spirographa (Spirographaceae, 5 ASVs) — all annotated as S. fusisporella, likely the same organism split across ASVs.
- Tympanis (Tympanidaceae, 3 ASVs) — another canker/twig-associated genus.
- Cladosporium, Jattaea, Valsa, Pyrenopeziza, Bjerkandera each appear among the top 20 by q-value.
- 23 of the 96 significant hits lack family-level taxonomy (unclassified/low-confidence assignments).

Internal consistency check (a good sign the signal is real rather than noise): all 7 Cytospora ASVs, all 5 Spirographa ASVs, and all 3 Tympanis ASVs share the same sign on their t-statistic — independent ASVs from the same genus are moving together with insect community composition, which you would not expect by chance if these were spurious hits.

Overall 75 of 96 significant taxa share one sign and 21 the other — worth checking later whether that split lines up with a specific direction on the insect ordination (e.g., something distinguishing certain traps/lures) rather than a generic gradient.

#
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

# after adding date
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

############
 Why did you use euclidean distance in the insect community PCoA instead of Bray-Curtis distance?

That's a comment already baked into the CoCA script (insect_dist <- vegdist(insect_hel, method = "euclidean")  # Hellinger + Euclidean ~ chord distance, compare_insect_fungi/insect_fungal_coca_analysis.general_workflow.r:151) — it's not Euclidean on raw counts, it's Euclidean on Hellinger-transformed data, which is a specific, deliberate combination:

- Hellinger transform = square-root of relative abundances per sample (decostand(x, "hellinger")). This down-weights dominant taxa and converts the data to proportions before any distance is computed.
- Euclidean distance on Hellinger-transformed data is mathematically equivalent to the Hellinger distance, one of the transformations Legendre & Gallagher (2001) specifically recommend as a Bray-Curtis substitute for methods that require true Euclidean geometry — PCA, RDA, and classical (metric) PCoA via cmdscale(), which is what's used here.

The reason it matters here specifically: Bray-Curtis is a semi-metric, not a proper Euclidean-embeddable distance. Feeding it into cmdscale() (classical/metric PCoA) commonly produces negative eigenvalues and a non-Euclidean embedding that needs a correction step (e.g. Cailliez) to interpret cleanly. Hellinger+Euclidean avoids that — it's guaranteed to embed properly — while still behaving like a community-ecology-appropriate distance (down-weighting dominants, handling double-zeros reasonably), unlike raw Euclidean on counts.

This is a different choice than insect_ords.R, which uses Bray-Curtis directly with metaMDS (non-metric MDS) — Bray-Curtis is the right tool there because NMDS doesn't need a Euclidean-embeddable distance, it just ranks dissimilarities. The CoCA script needed the axis as a numeric covariate in a cmdscale-based PCoA, so it uses the Hellinger/Euclidean combination instead.

If you'd rather it match the Bray-Curtis convention used elsewhere in the project (e.g., substituting a PCoA on Bray-Curtis with a Cailliez correction, or non-metric MDS scores instead of classical PCoA), I can swap that in — but the current choice wasn't arbitrary, it's the standard fix for this exact use case.

############
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

# We should add a test for how much of the fungal community is seasonal.

Fungal taxa vs. date, with taxonomy -- a very clean pattern. The top 20 hits split almost entirely into two genera pulling in opposite directions:
- Mycosphaerellaceae / Ramularia -- increasing over the season (later): 11 of top 20
- Valsaceae / Cytospora -- decreasing over the season (earlier): 6 of top 20

Ramularia is a foliar plant-pathogen genus (leaf-spot fungi) -- building up over the growing season as leaf area increases is textbook phenology for a foliar pathogen. Cytospora (canker fungus) being higher earlier in the season lines up with the earlier CoCA finding that tied Cytospora to the bark-beetle-dominated (early-season) end of the insect community axis -- this is now a second, independent line of evidence (direct date test, not mediated through the insect axis) for the same seasonal pattern. Also notable in the early-season group: Piptoporus betulinus (birch polypore, wood decay) and Ganoderma sp. (wood decay), both consistent with bark/wood-associated fungi peaking alongside bark beetle activity.

Overall takeaway: no confirmed direct fungal-taxon / insect-community-composition relationship once season is accounted for (either axis). But there is a strong, independent seasonal signal in both communities, with a plausible shared ecological driver (Cytospora canker fungus tracking early-season bark/ambrosia beetle activity, evidenced two separate ways) worth following up with a larger dataset or a study design that better decouples date from insect community turnover.

Output files (data/2024_insect_data/): fungal_insect_association_results.PCoA1_only.csv, fungal_insect_association_results.PCoA1_plus_date.csv, fungal_insect_association_results.PCoA2_only.csv, fungal_insect_association_results.PCoA2_plus_date.csv, insect_taxa_date_association.csv, fungal_taxa_date_association.csv, fungal_taxa_date_association.taxonomy.csv.

############
Lure had a strong effect on the insect community (Lure F=21.9 in insect_PCoA1 ~ site+lure+date), as expected since it's the trap attractant. But there's no biological reason to expect fungi respond to lure independently, and a companion dataset showed minimal fungal~lure effect there too. Quick check: does including lure in the fungal~insect_PCoA1 model partial out real signal?

Reran the PCoA1 permutation test (same 200 CoCA-ranked candidate taxa) with and without lure in the formula, no-date and date-adjusted:
- PCoA1, with lure, no date: 81/200 significant (q<0.10)
- PCoA1, without lure, no date: 82/200 significant
- PCoA1, with lure + date: 0/200 significant
- PCoA1, without lure + date: 0/200 significant

Dropping lure barely changes anything (81 vs 82 taxa, no-date; 0 vs 0, date-adjusted) -- if lure were partialing out real signal, removing it should have recovered noticeably more hits, and it didn't. Correlation of insect_PCoA1 t-statistics with vs without lure (no-date models): r=0.96, i.e. essentially identical effect sizes either way.

Direct check -- do the candidate fungal taxa respond to lure at all, independent of the insect axes? lure is constant within trap (trap = site x lure), so the within-trap permutation scheme used elsewhere can't test it meaningfully; used a parametric, site-controlled F-test instead (descriptive only, not permutation-corrected): 0 of 200 candidate fungal taxa show even a nominal lure effect (q<0.10).

Conclusion: lure is a genuine, strong driver of the insect community but has no detectable independent effect on the fungal community here, and including it in the models neither suppresses nor inflates the fungal~insect_PCoA1 relationship. The date-driven conclusion above stands on its own, unaffected by whether lure is in the model.

This check is now folded into the main script as two extra PCoA1 test variants (no-lure, no-date/date-adjusted) plus a direct fungal~lure parametric check, rather than a one-off side analysis. Output files added: fungal_insect_association_results.PCoA1_only.no_lure.csv, fungal_insect_association_results.PCoA1_plus_date.no_lure.csv, fungal_taxa_lure_direct_check.csv.

############
The 200/200 fungal taxa found significant vs. date earlier were pre-selected for strong |correlation with date|, so that number can't be read as "how seasonal is the fungal community" -- it's true by construction. Added a new script, compare_insect_fungi/fungal_community_seasonality.r, to answer that question two unbiased ways.

1) Per-taxon date test on every prevalence-filtered fungal taxon (no pre-selection), same trap-blocked permutation scheme as before (fungal_clr ~ site + lure + date, permuting date within trap): 2248 of 6887 taxa (32.6%) significant at q<0.10. That's a much more modest, and more trustworthy, estimate of the seasonal fraction of the fungal community than "200/200" implied -- roughly a third of taxa carry a detectable seasonal signal, not "basically everything."

2) Whole-community view, mirroring insect_ords.R's approach for the insect side: the raw ASV table (sequencing depth varies ~150x across the 69 samples) was repeatedly rarefied to depth 5000 (100 iterations, all samples retained) via multiple_subsamples(), the resulting Bray-Curtis distance matrices averaged via avg_matrix_list(), then PERMANOVA (~ site + site:lure + site:date) and NMDS run on the averaged matrix, same term structure as the insect community PERMANOVA in insect_ords.R.

PERMANOVA (fungal community): site R2=14.7% (F=4.96, p=0.001), site:lure R2=10.1% (F=1.27, p=0.05), site:date R2=22.7% (F=5.72, p=0.001), residual 52.5%. Site:date is the single largest term -- directly parallel to the insect community result (Site R2=9.9%, Site:Lure R2=20.1%, Site:mdy R2=19.7%), except here date/season dominates over lure, whereas lure dominates on the insect side. That's a coherent picture: lure structures which insects get caught, while season structures which fungi are present, and the per-taxon result above (32.6%) is the taxon-level face of the same site:date term.

NMDS (stress ~0.17) plotted with date as a color gradient and site as point shape: figures/fungal_community_NMDS.date_and_site.pdf.

Output files added: data/2024_insect_data/fungal_taxa_date_association.all_taxa.csv, data/2024_insect_data/fungal_community_permanova.csv, figures/fungal_community_NMDS.date_and_site.pdf.

############
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

############
Extended the CoCA workflow to test insect_PCoA3 as a predictor of fungal abundance, following up on the PCoA3/lure finding above. Also bumped the CoCA algorithm itself (Step 3, the coca() call that ranks candidate taxa by loading strength) from n_axes_use=2 to 3, so the top-200 candidate list is now selected using CoCA axes 1-3 rather than 1-2.

insect_PCoA3 was added as a third covariate alongside PCoA1/PCoA2 in every per-taxon model (all three are mutually controlled for, same logic as the original PCoA1/PCoA2 design: the PCoA axes are orthogonal by construction, so this isn't correcting collinearity between axes, it's making each axis's permutation test specific to "does this taxon track axis k given what's already known about the other axes" rather than a less precise univariate test). Ran the same with/without-date x with/without-lure spec matrix already used for PCoA1 -- the without-lure comparison matters more for PCoA3 than it did for PCoA1, since PCoA3 is itself ~81% explained by lure, so the with-lure models are only testing the ~19% of PCoA3 that isn't lure.

Result: PCoA3 comes up completely null. 0 of 200 CoCA-ranked candidate taxa significant (q<0.10) in all four variants -- with lure/no date, with lure/with date, without lure/no date, without lure/with date. Contrast with PCoA1 over the same candidates: 46-60/200 significant without date (46 with lure, 60 without), collapsing to 0/200 once date is added, as found previously.

So despite insect_PCoA3 being the single strongest lure-associated axis in the whole insect ordination, no fungal taxon among the top-200 CoCA-loading candidates tracks it -- whether or not lure is partialled out. Combined with the direct parametric fungal~lure check (also 0/200, both before and after this update), the small lure effect PERMANOVA detects in the fungal community (10.1%, p=0.05) does not appear to be mediated through this specific insect community axis. Caveat: the candidate list is still selected by CoCA loading strength (now over 3 axes), not by a PCoA3- or lure-specific screen, so this rules out "PCoA3 explains the CoCA hits" but doesn't fully rule out a lure-associated fungal signal outside that candidate set -- though the null parametric check makes that less likely.

Output files added: fungal_insect_association_results.PCoA3_only.csv, fungal_insect_association_results.PCoA3_plus_date.csv, fungal_insect_association_results.PCoA3_only.no_lure.csv, fungal_insect_association_results.PCoA3_plus_date.no_lure.csv, figures/fungal_insect_volcano.PCoA3.with_lure.png, figures/fungal_insect_volcano.PCoA3.no_lure.png.

############
The 3-axis CoCA loading-strength screen (results moved to data/compare_insects_fungi_top3axes/) selects candidate fungal taxa using CoCA axes 1-3 jointly. Set up two alternative single-axis candidate-selection strategies to compare against it: PCoA1 (the date-associated axis) and PCoA3 (the lure-associated axis), each isolating one specific driver instead of the joint CoCA structure.

cocorresp::coca() requires non-negative "community-style" data for both blocks, so the signed insect_PCoA1/PCoA3 axis scores can't be plugged in directly as a single-column CoCA predictor. Candidate selection instead reuses/extends the unbiased, all-taxa per-taxon screens already used elsewhere in this project: fungal_taxa_date_association.all_taxa.csv (already run, fungal_community_seasonality.r) for the PCoA1 strategy, and a new parallel script, fungal_community_lure_association.r, for PCoA3 -- same design (F-test of the lure term, site+date vs. site+lure+date, trap-within-site permutation, same scheme insect_pcoa_lure_association.r used to screen insect PCoA axes for lure), run on all 6887 prevalence-filtered fungal taxa with no pre-selection. Result: 0 of 6887 fungal taxa individually significant for lure (q<0.10) -- consistent with the earlier finding that fungal PERMANOVA detects a real but small, diffuse lure effect (10.1%, p=0.05) not concentrated in specific taxa. The ranking is still usable as a candidate screen even with no individually-significant taxa, same as how the original CoCA loading-strength ranking never required its top-200 to be individually significant either.

Two new scripts (insect_fungal_coca_analysis.PCoA1_date_axis.r, insect_fungal_coca_analysis.PCoA3_lure_axis.r) take the top 200 taxa from each screen and run the same per-taxon permutation-test machinery as the original CoCA script (site + lure? + date? + insect_PCoA1/2/3, permuting only the target axis within trap), restricted to their respective target axis (4 variants: with/without date x with/without lure). Results: data/compare_insects_fungi_PCoA1_date/, data/compare_insects_fungi_PCoA3_lure/; figures in the parallel figures/ subfolders.

PCoA1/date-screen results: 199/200 candidates significant vs. insect_PCoA1 without date (expected -- these taxa were selected FOR date association, and PCoA1 is largely the date axis), dropping to 21/200 once date is added as a covariate (50/200 if lure is also dropped). Contrast with the original 3-axis-CoCA candidate set over the same PCoA1-plus-date test: 0/200. So a candidate list built specifically around date association retains some signal beyond a linear date covariate that the CoCA-selected list does not -- but this is likely at least partly circular (candidates were pre-selected for date association, and insect_PCoA1's collinearity with date is r=0.78, not 1.0, so imperfect collinearity alone could produce partial signal survival) rather than clean evidence of a direct insect-community effect. Should not be over-read as new biology without a stricter check.

PCoA3/lure-screen results: 0/200 candidates significant vs. insect_PCoA3 in all 4 variants -- same null result as the original 3-axis-CoCA candidate set. Even a candidate list specifically built to surface lure-associated fungal taxa shows no detectable individual-taxon relationship with the lure-associated insect axis. One difference: the direct parametric fungal~lure check (site-controlled, uncorrected) finds 10/200 candidates nominally lure-associated in this screen, vs. 0/200 for the CoCA-selected candidates -- consistent with the community-level lure signal being diffuse (spread thinly across many taxa, each below the multiple-testing threshold) rather than concentrated, and more findable by a direct univariate screen than by CoCA loading strength.

Candidate-set overlap comparison (compare_selection_strategies.r, data/compare_selection_strategies/, figures/compare_selection_strategies/candidate_overlap_regions.png) is the more striking result: the 3-axis CoCA list barely overlaps with either single-axis screen (4/200 shared with PCoA1_date, 4/200 shared with PCoA3_lure, Jaccard ~0.01 both times), while the PCoA1_date and PCoA3_lure screens overlap substantially with EACH OTHER (76/200 shared, Jaccard 0.235) despite nominally targeting different drivers (date vs. lure). Only 1 taxon is common to all three lists. Two read on this: (1) joint CoCA loading strength (which reflects taxa best explaining insect-fungal covariation via a multivariate PLS-type decomposition) is picking up on a genuinely different signal than "correlates with one univariate covariate," so the original 3-axis screen and these single-axis screens are not interchangeable; (2) the substantial date/lure screen overlap suggests both univariate per-taxon tests may be partly surfacing generically high-effect-size/high-variance taxa (strong site or trap-driven differences) rather than drivers cleanly specific to date vs. lure -- worth treating the "top-loading" label for any of these three strategies as method-dependent, not a fixed property of a taxon.

Overlap among each strategy's own SIGNIFICANT hits (not just candidate membership) is also low even where both strategies found hits: for PCoA1 without date, top3axes found 46/200 significant and PCoA1_date-screen found 199/200, but only 4 taxa are significant in both (Jaccard 0.017) -- the two methods are largely flagging different specific ASVs as "significant," not just different-sized supersets of the same core signal.

Output files added: data/2024_fungi/fungal_taxa_lure_association.all_taxa.csv; data/compare_insects_fungi_PCoA1_date/ (fungal_insect_association_results.PCoA1_*.csv, lure_partialling_check.summary.csv, fungal_taxa_lure_direct_check.csv); data/compare_insects_fungi_PCoA3_lure/ (fungal_insect_association_results.PCoA3_*.csv, lure_partialling_check.summary.csv, fungal_taxa_lure_direct_check.csv); data/compare_selection_strategies/ (candidate_taxa_membership.csv, candidate_set_overlap_summary.csv, significant_hits_overlap.csv); figures/compare_insects_fungi_PCoA1_date/fungal_insect_volcano.PCoA1.png; figures/compare_insects_fungi_PCoA3_lure/fungal_insect_volcano.PCoA3.{with,no}_lure.png; figures/compare_selection_strategies/candidate_overlap_regions.png.
