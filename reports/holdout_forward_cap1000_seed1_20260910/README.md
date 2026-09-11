# Forward-only comparison on 140 elementary-step cases

**Reactant → product only.** AAM uses one seed order per cut, branch cap 1000,
tolerance 1.0 and the uncut/single-edge sweep. SLAP is shown as the native XYZ
call and as the uncut plus every individual reactant-edge deletion, including H.
The frozen engine is `98b01b175eeed31f70d13e7cbf178b80bf07c9e0`.

All three protocols provide full mappings for all 140 cases. AAM has
6 lower scores and
134 ties against native SLAP;
against SLAP with the sweep it has
**4 lower scores,
136 ties and
0 higher scores**.
No reference mappings are available, so these are minimum-change scores rather
than mapping-accuracy percentages.

## Mapping time

| Method | Total CPU s | Mean CPU s/case | Median CPU s/case | Mean observed wall s/case |
| --- | ---: | ---: | ---: | ---: |
| AAM + sweep | 942.576 | 6.733 | 2.249 | 1.760 |
| Native SLAP | 12.145 | 0.087 | 0.076 | 0.108 |
| SLAP + sweep | 693.756 | 4.955 | 2.016 | 4.958 |

AAM uses **1.359×** the mapping CPU of SLAP+sweep.
These are actual saved forward-direction measurements, never halved totals.
Each comparison uses the same host per case. AAM CPU sums the parent and eight
workers and excludes measured checkpoint persistence/loading. SLAP sweep CPU
includes XYZ preparation, graph construction, native mapping and export; it
excludes frame encoding, persistence and subsequent candidate aggregation.
The forward sweep contains 8,370 native calls and 10,765 native outputs.

**Wall time is an observed execution comparison with different core counts:**
AAM used eight workers; SLAP used one. AAM wall includes archive I/O; SLAP-sweep
wall is the sum of measured preparation/construction/mapping/export durations.
It is not an equal-core speedup or full end-to-end latency comparison.

Ranking/scoring is separate from mapping: the saved forward AAM scorer used
117.533 CPU s and native SLAP scoring used
11.083 CPU s. The saved sweep
scoring timer covers the deduplicated bidirectional candidate union. It cannot
be assigned to the forward subset, so forward SLAP-sweep scoring CPU is null.
No total mapping-plus-scoring speedup is claimed.

## Alternative minimum-event patterns

At the 136 AAM / SLAP-sweep tied scores, the catalogues contain
196 verified AAM patterns and 195 SLAP-sweep
patterns: **192 shared, 4 AAM-only,
3 SLAP-only proved, and 0
unresolved memberships**. Both methods have multiple verified alternatives in
39 tied cases. Observed sets agree in
132 cases; equality is certified for
112 complete catalogues.

| Case | AAM / sweep score | Verified AAM patterns | Sweep patterns | Shared | Proved AAM-only / sweep-only |
| --- | ---: | ---: | ---: | ---: | ---: |
| 11 | 4 / 4 | 1 | 2 | 1 | 0 / 1 |
| 64 | 5 / 5 | 2 | 1 | 0 | 2 / 1 |
| 101 | 4 / 4 | 1 | 2 | 1 | 0 / 1 |
| 129 | 8 / 8 | 10 | 8 | 8 | 2 / 0 |

SLAP catalogues are complete for all 140 saved forward output families. AAM
catalogues are complete in **120/140 cases**;
remaining counts are lower bounds. A complete bidirectional catalogue can bound
the forward subset at the same score: every class is then tested for forward
membership. Reverse witnesses are never accepted as forward detections. Other
cases use bounded enumeration of the forward AAM families. An incomplete scan
never proves exclusion. Comparisons of alternatives are inapplicable when the
two recorded minimum scores differ.

Against native SLAP, 134 scores tie, with 165 shared
patterns, 28 AAM-only patterns and
2 SLAP-only patterns proved.

Patterns distinguish broken, formed, strengthened and weakened bonds, including
explicit H, modulo joint symmetries of the reactant event-response graph. Bond
presence uses WBO > 0.2 and order changes use magnitude > 0.5. These counts describe
candidate bond edits; they do not validate a physical reaction pathway.

## Effect of dropping the reverse direction

AAM retains its bidirectional best score on every case. SLAP sweep loses the
better scores contributed by reverse searches in cases 31 and 96. Consequently,
the earlier two AAM-lower scores become four in the forward-only comparison.

| Case | Name | AAM | Native SLAP | SLAP sweep |
| --- | --- | ---: | ---: | ---: |
| 31 | pr14.Pd_hydroamination_JOC2025_TS3_step2_alkene_inserion | 7 | 9 | 9 |
| 59 | pr16.carbocation_ts11 | 2 | 4 | 2 |
| 64 | pr16.carbocation_ts5 | 5 | 9 | 5 |
| 69 | pr17.carbene.ins_ts8 | 5 | 6 | 6 |
| 96 | pr7.V.dodh_ts13 | 1 | 3 | 3 |
| 104 | pr7.V.dodh_ts32 | 5 | 7 | 7 |

## Verification and artifacts

The audit scans 96,481 forward full-terminal AAM
witnesses and independently rescored 579
published method/pattern witnesses with the scalar bond-event routine. The
event-enumeration engine was also checked against 13,872 exhaustive assignments
in ten small real SLAP families in the preceding audit. Candidate origins and
all forward call journals are checked before accepting scores or times.
Direct forward-family enumeration independently confirms the subset-catalogue
results for cases 11, 64 and 101 (`directional_checks.json`).

[Open the offline bond-edit viewer](viewer.html). `comparison.pdf`, `.png` and
`.svg` show the scores, mapping CPU and alternative overlap. `per_case.csv` and
`per_case.json` contain all 140 records. Compressed case outputs, source hashes,
frozen analysis scripts, job records and the preceding finite checks accompany
the report. Full source run: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_forward_cap1000_seed1_20260910`.

`summary.json:new_analysis_timing` records case-worker analysis separately
(6374.888 CPU s; publication,
validation and process startup are excluded). It reuses
prior catalogues and exclusion proofs, so it is not a fresh alternative-extraction
speed benchmark. Mapping searches were not rerun, and the earlier bidirectional
report is preserved as a distinct protocol.

Presentation cleanup: the bundled HTML template now embeds the shared reaction
style. `presentation-migration.json` inside the archive records its old and new
hashes; original provenance hashes describe the historical publication. All other
archived source members, including the frozen numerical analysis, are unchanged.
