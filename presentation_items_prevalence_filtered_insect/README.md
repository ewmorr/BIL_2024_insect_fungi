# presentation_items_prevalence_filtered_insect/

Parallel to `presentation_items/`, rebuilt on the >=5-sample individual-taxon
prevalence-filtered insect table (all families) instead of the
Curculionidae+Latridiidae family restriction. See `project_organization.md`
at the repo root for the full rationale, table-construction details, and the
current status of each figure in this lineage.

Convention: scripts here keep the SAME names as their `presentation_items/`
counterparts (`fig1_taxonomic_breakdown.R`, `fig2_nmds_procrustes.R`, ...) --
the containing directory is what marks them as the prevalence-filtered
variant, so the filename doesn't need to repeat that. Outputs go to
`data/presentation_items_prevalence_filtered_insect/` and
`figures/presentation_items_prevalence_filtered_insect/`.

Status as of 2026-09-01: all 8 figures are done (fig1 lineage-invariant, no
port needed; fig2/fig3 2026-08-31; fig4-fig8 2026-09-01). Nothing in
`presentation_items/` has been touched or overwritten.

**2026-09-02**: fig3, fig4, fig5, fig7, and fig8 all revised to a linear-
term-takes-priority reading of peak timing (early/late/mid-season/bimodal),
replacing the old quadratic-priority `shape` reading and (for fig4/5/7/8)
the "up to 2 independent hit-instances per taxon" combined-dodge design --
each taxon now gets exactly ONE mutually-exclusive category. fig4/5/7/8's
bars all run the same direction (0 -> positive) in early/bimodal/mid-
season/late top-to-bottom order. Added `figS1_date_shape_bin_examples.R`, a
supplemental figure/table showing example taxa and counts for the full
linear x quadratic significance crosstab that motivated the change. A
positional-indexing bug introduced partway through this revision (fig4/fig5
briefly had "early"/"mid-season" swapped) was caught and fixed the same
day -- see `iterative_analysis_updates.md`'s 2026-09-02 CORRECTION entry.
See `project_organization.md`'s "Peak-timing convention" note for the full
rationale and current status.
