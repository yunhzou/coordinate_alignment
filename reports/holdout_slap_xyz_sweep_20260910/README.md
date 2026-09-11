# 140 cases: one-seed cap-1000 AAM versus native SLAP and SLAP with cuts

All three configurations return full mappings on all 140 XYZ/WBO cases.
Against native XYZ SLAP, AAM has six lower minimum-change scores, 134 ties,
and no higher scores. Against **SLAP plus the single-edge sweep**, AAM has
**2 lower scores, 138 ties,
and 0 higher scores**.

The sweep improves on native SLAP in 4 cases.
Compared with the bidirectional uncut SLAP control, cuts improve
4 cases.

| Configuration | Total mapping CPU s | Mean CPU s/case | Median CPU s/case | Separate scoring CPU s |
| --- | ---: | ---: | ---: | ---: |
| AAM: one seed, cap 1000 | 2238.525 | 15.989 | 5.841 | 243.277 |
| Native XYZ SLAP | 12.145 | 0.087 | 0.076 | 11.083 |
| SLAP + single-edge sweep | 1321.806 | 9.441 | 4.036 | 47.046 |

AAM uses 1.694 times the recorded mapping CPU
of SLAP plus the sweep. The sweep uses 108.832
times the native SLAP mapping CPU. Per-case timings are in `per_case.csv`.

## Protocol and timing scope

AAM is the unchanged stable engine `98b01b175eeed31f70d13e7cbf178b80bf07c9e0`,
with one deterministic seed order per cut, cap 1000, tolerance 1.0, both directions,
and the existing uncut/single-edge sweep. Its results and native SLAP timings are
from the fresh paired experiment in `../holdout_cap1000_seed1_20260910/`.

The added SLAP sweep keeps native XYZ graph perception (`geoms2lgp`, original
component files, bond scale 1.2) and `binary=True`, with heavy-atom symmetry
breaking. It runs the unmodified mapper on the uncut graph and every individual
source-edge deletion, including H bonds, in both directions. Initial WL labels
are rebuilt after each cut. Native XYZ edges are binary; this is a cut ablation
on the existing XYZ comparator, separate from the Golden SMILES binary/weighted
union. Every output is scored on the original full-H WBO matrices, with graph
threshold 0.2 and bond-order-change threshold 0.5. No reference mappings or
previous scores enter the search. These are best saved representative scores,
not reference accuracy or a proof of global optimality.

Each sweep case ran on the same host as the paired AAM/native-SLAP case. AAM
uses eight workers, while SLAP uses one; the table reports additive CPU, not
equal-core wall latency. AAM CPU excludes measured checkpoint persistence and
loading. Sweep CPU includes initial XYZ graph preparation, per-cut graph
construction, native mapping and label export; compressed archival I/O and
offline scoring are separate. Scoring workloads differ: AAM also canonicalizes
heavy mapping classes, while SLAP optimizes H assignments within its native
label groups. Do not interpret the mapping-only ratio as full application cost.

## Audit

All 16,748 planned SLAP mapping calls completed without
errors or timeouts. All 21,325 native outputs
are accounted for through 1,660 unique
label partitions. Identical partitions are deduplicated across cuts/directions;
all discovering cut identities are retained. Every scored representative is a
full element-preserving bijection and was independently rescored with the scalar
bond-event implementation. Every H-label optimization certified its final score.
All 140 uncut R-to-P controls reproduce the exact native label partitions and
minimum scores of the preceding comparison.

`audited_witnesses.json.gz` contains every scored mapping and its discovering
cuts. `native_call_records.json.gz` records all committed cuts, timings and
native-frame hashes. Full native LAP archives remain at `/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_slap_xyz_sweep_20260910`;
`archive_hashes.json` identifies them. Frozen sources and input hashes are in
`manifest.json`. There are no unresolved calls or score optimizations.
Case 17 had one allocation stall before starting; it was canceled and submitted
again on the same host, with no mapping work repeated. Both allocation records
are retained in `submissions.json` and the Slurm accounting export.

## Cases with differing scores

Indices are zero-based, matching the earlier reports.

| Case | Name | AAM | Native SLAP | Uncut SLAP both directions | SLAP sweep |
| --- | --- | ---: | ---: | ---: | ---: |
| 31 | pr14.Pd_hydroamination_JOC2025_TS3_step2_alkene_inserion | 7 | 9 | 9 | 7 |
| 59 | pr16.carbocation_ts11 | 2 | 4 | 4 | 2 |
| 64 | pr16.carbocation_ts5 | 5 | 9 | 6 | 5 |
| 69 | pr17.carbene.ins_ts8 | 5 | 6 | 6 | 6 |
| 96 | pr7.V.dodh_ts13 | 1 | 3 | 3 | 1 |
| 104 | pr7.V.dodh_ts32 | 5 | 7 | 7 | 7 |
