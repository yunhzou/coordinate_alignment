# Minimum-event alternatives on the 140 XYZ/WBO cases

[Open the offline bond-edit viewer](viewer.html), or download it and open locally.
`comparison.pdf`, `comparison.svg` and `comparison.png` show the aggregate result;
`per_case.csv` and `per_case.json` retain every case and its membership status.

## AAM versus SLAP with the single-edge sweep

At the **138 tied benchmark scores**, the verified AAM catalogue contains
**204** distinct bond-edit patterns; the complete saved
SLAP-sweep minimum-score label-family catalogue contains **207**.
There are **204 verified shared patterns**, **0
AAM-only patterns**, **3 patterns proved present only
in SLAP's saved families**, and **0 unresolved cross-method memberships**.
Counts for an incompletely enumerated AAM catalogue are lower bounds.

Both methods have more than one verified minimum-score pattern in
**40** tied cases. Their observed pattern sets agree in
135 cases; complete enumeration certifies equal
catalogues in 91 cases. A difference is certified
in 3 tied cases: [11, 64, 101].

| Case | Name | Common minimum score | Verified AAM patterns | SLAP-sweep patterns | Shared | SLAP-only proved |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 11 | pr10.carborane_ts3 | 4 | 1 | 2 | 1 | 1 |
| 25 | pr14.Pd_hydroamination_JOC2025_TS11_step1_reductive_coupling | 5 | ≥4 | 4 | 4 | 0 |
| 64 | pr16.carbocation_ts5 | 5 | 2 | 3 | 2 | 1 |
| 101 | pr7.V.dodh_ts19 | 4 | 1 | 2 | 1 | 1 |
| 129 | TS_04 | 8 | ≥16 | 16 | 16 | 0 |
| 132 | TS_08 | 4 | 5 | 5 | 5 | 0 |

The original one-representative comparison hid alternatives in both methods.
In case 129, for example, the terminal/optimized representatives gave 12 AAM
patterns and seven SLAP-sweep patterns. Family enumeration and membership
queries establish at least 16 shared patterns in the saved outputs. A smaller representative
list does not establish that the algorithm missed those alternatives.

Against **native SLAP without the sweep**, 134 scores tie. The
catalogues contain 165 verified shared patterns,
34 AAM-only patterns, 2
SLAP-only patterns and 0 unresolved memberships on
those tied cases. The native and swept protocols remain separate.

## What constitutes an alternative

An alternative is the full set of **broken, formed, strengthened and weakened
bonds**, pulled back to original reactant atom indices. Order changes use the
same 0.5 WBO threshold and bond presence uses the same 0.2 floor as the benchmark.
Explicit H is included. Multiple atom maps with the same edits count once.
All edits are jointly canonicalized on the full reactant graph, with exact signed
event-response colors. Symmetry-equivalent atom relabelings are merged; different
affected bonds or event types remain distinct. Geometry and stereochemistry do
not distinguish patterns in this analysis. These are candidate bond-edit
mechanisms, not experimentally established reaction pathways.

## Completeness and checks

The analysis scanned **317,913 full AAM terminal
witnesses** from the frozen cap-1000, one-seed, tolerance-1.0 archives. It then
queried the recorded correlated path families and enumerated event patterns at
the **previously recorded benchmark minimum score**. It does not replace that
score with an independently optimized global minimum over all implicit mappings.

SLAP enumeration is complete for all 140 saved label-family catalogues. It blocks
whole event sets and removes only safe source/target twin-atom permutations,
preserving every event class. AAM enumeration is complete for
96/140 cases. Incomplete tied cases are
[1, 4, 6, 18, 19, 25, 26, 28, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 47, 48, 49, 53, 54, 55, 56, 58, 68, 70, 71, 76, 77, 86, 87, 90, 91, 97, 105, 106, 109, 121, 122, 123, 129, 134]; their AAM counts carry a lower-bound marker in the
viewer. An incomplete enumeration is never used to assert absence. The separate
cross-family membership checks can still prove exclusion of an individual pattern
by scanning every applicable saved family.

Every published mapping was independently rescored with the scalar bond-event
routine (595 method/witness checks).
Direct enumeration of 13,872 complete assignments in ten small real SLAP families
exactly matches the symbolic event-pattern sets. `finite_checks.json` records
these tests and event-type/symmetry checks. All analysis is offline; none of the
mapping searches or reported mapping timings was rerun.

Full analysis and frozen sources: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_minimum_event_patterns_20260910`.
Compressed snapshots, membership witnesses, enumeration records, input/source
hashes, and job submissions accompany this report. The fixed manuscript baseline
is separate from this cap-1000 follow-up.

## Reproduction

The array submission records retain the exact environment, commands and budgets.
`analysis_sources.tar.gz` includes the frozen worker versions, published analysis
scripts and targeted-reproduction checks. Mapping inputs and large AAM archives
are identified by the preceding benchmark reports and their hash manifests.

The targeted case-25 SLAP retry and case-129 enrichment can be reproduced with
`bench/complete_minimum_event_checks.py retry_slap --index 25` and
`bench/complete_minimum_event_checks.py enrich_aam --index 129`, supplying
`--run` and a fresh `--destination`. Both wrappers were checked against the saved
results. The default retry uses a 30-second solver-query limit.

To rebuild this report, run `bench/publish_minimum_event_analysis.py --run`
with the source directory above and a fresh `--destination`, using the frozen
source paths in the submission records and the manuscript plotting dependencies.
`slurm_accounting.tsv` records analysis workers, including cancelled allocations
that never started. These postprocessing costs do not alter the mapping-speed
comparison in the preceding benchmark report.
