# One seed, branch cap 1000: 140-case comparison with native XYZ SLAP

Fresh experiment on the same 140 XYZ/WBO inputs, using frozen stable AAM engine
`98b01b175eeed31f70d13e7cbf178b80bf07c9e0`. AAM uses one deterministic seed order per cut,
branch cap 1000, tolerance 1.0, both directions, and the original uncut/single-edge
cut sweep. All other search settings match the stable ten-seed/cap-100 baseline.
SLAP is rerun through its unmodified native XYZ interface with original components,
binary adjacency, `break_sym='heavy'`, and default bond scale 1.2.
This comparator is separate from the Golden-set SLAP single-cut sweep ablation.

| Minimum full-H WBO change score | Cases |
| --- | ---: |
| AAM lower | 6 |
| Equal | 134 |
| SLAP lower | 0 |

Both methods return full mappings for all 140 cases. All 840 execution/scoring
phases succeed, and all SLAP H-label score optimizations finish with certified
bounds. These cases have no reference atom maps: this is an objective-score
comparison, not a chemical-accuracy estimate or proof of a global minimum.
Case 123 scores 3 for AAM and 3 for SLAP.

| Recorded time | AAM | Native XYZ SLAP |
| --- | ---: | ---: |
| Total search CPU seconds | 2238.525 | 12.145 |
| Mean search CPU seconds/case | 15.989 | 0.087 |
| Median search CPU seconds/case | 5.841 | 0.076 |
| Separately measured scoring CPU seconds | 243.277 | 11.083 |

AAM uses **184.31 times** the recorded mapping CPU of
native XYZ SLAP. Each paired case runs on the same host. Case 0 is the retained
eight-CPU local pilot; the other 139 pairs run in eight-CPU Slurm allocations.
AAM uses eight workers and SLAP one, so wall times do not establish equal-core
latency. AAM CPU includes parent and child search, process setup, IPC and profiling,
excluding measured checkpoint persistence/loading. SLAP CPU covers `map_3d`,
including native XYZ processing. Imports and scheduling are excluded.
Scoring is separate: AAM scoring also canonicalizes saved heavy mapping classes;
SLAP scoring resolves native H-label assignments. These are different output
processing workloads, so the main ratio uses the recorded search timers.

The prior cap-100/ten-seed comparison gave six AAM-lower cases, one SLAP-lower
case and 133 ties. Fresh SLAP scores differ from the archived comparison in
0 cases.

| Case index (zero-based) | Name | AAM score | SLAP score |
| --- | --- | ---: | ---: |
| 31 | pr14.Pd_hydroamination_JOC2025_TS3_step2_alkene_inserion | 7 | 9 |
| 59 | pr16.carbocation_ts11 | 2 | 4 |
| 64 | pr16.carbocation_ts5 | 5 | 9 |
| 69 | pr17.carbene.ins_ts8 | 5 | 6 |
| 96 | pr7.V.dodh_ts13 | 1 | 3 |
| 104 | pr7.V.dodh_ts32 | 5 | 7 |

`per_case.csv` contains all 140 comparisons and timings. `case_metrics.json`
retains directional completion, cap flags, phase status and hardware. The
compressed witnesses were independently rescored with the scalar bond-event
implementation; the published scores use the original vectorized evaluator.
`archive_hashes.json` identifies full local search archives at `/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_cap1000_seed1_20260910`.
The search driver is `bench/holdout_cap_seed_benchmark.py`; this report is built
by `bench/publish_holdout_cap_seed.py`. Existing manuscript baseline panels are
unchanged by this separate experiment.
