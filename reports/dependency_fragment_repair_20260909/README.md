# Dependency repair + native AAM scheduling

**Result: 1.24x overall compute-CPU speedup, 1.43x in search, and lower peak
worker memory, with identical tested compressed results. This is useful but not
the requested order-of-magnitude improvement. Production defaults are unchanged.**

Branch: `experiment_dependency_fragment_repair`, based on `3a4d459`.
[Design and exact reuse conditions](design.md).

## Final paired comparison

Four directional cases, ten parallel seed workers per case, every sweep cut,
branch cap 100, matching tolerance 1.0, event tolerance 0.5, bond floor 0.2,
explicit hydrogen. **Existing cut-specific random seed orders are unchanged.**
Each fresh/optimized pair ran sequentially in the same ten-CPU allocation on
bosque7 (Xeon Gold 6248R); pair order was counterbalanced. All four allocations
finished successfully in under 146 seconds, including both methods and I/O.
Each method had a 300-second watchdog; allocations had ten-minute limits.

| Phase | Current Python scheduler + native growth | Native scheduler + dependency repair |
|---|---:|---:|
| Search, CPU seconds | 827.71 | 580.48 |
| Symmetry finalization, CPU seconds | 501.61 | 492.13 |
| Setup + terminal-witness scoring, CPU seconds | 6.19 | 6.36 |
| **Total compute CPU seconds** | **1335.51** | **1078.96** |
| Maximum single-worker high-water RSS, MiB | 1576.00 | 921.44 |

This is **19.2% less total compute time**, **29.9% less search time** and **41.5%
less maximum worker RSS**. RSS is not aggregate ten-worker job memory.

Per-case totals (seconds):

| Case / direction | Fresh compute CPU | Optimized compute CPU | Fresh ten-worker elapsed including I/O | Optimized elapsed including I/O |
|---|---:|---:|---:|---:|
| 25 / R-to-P | 539.19 | 373.14 | 74.91 | 55.58 |
| 76 / P-to-R | 530.18 | 470.06 | 75.29 | 67.33 |
| 77 / P-to-R | 264.54 | 232.94 | 40.50 | 35.12 |
| 114 / R-to-P | 1.61 | 2.83 | 0.87 | 1.66 |

**The small vanadium case was slower in this final paired observation.** Earlier
paired observations differed substantially at that scale. We do not claim a
universal speedup or a repeated-trial timing distribution. The aggregate is
dominated by the three long tails, all of which improved in this comparison.

CPU time sums actual worker process CPU, not reserved cores multiplied by wall
time. Queue, imports and input loading are excluded. Audit JSON encoding and
archive persistence are separately recorded, excluded from compute totals.
Parallel elapsed time is shown only with its I/O-inclusive label; it is not
obtained by subtracting summed worker I/O from job duration. This benchmark does
not include the full downstream publication pattern-collection workflow.

## What was tested and learned

- Dependency repair reused **590,199 / 928,772 conditional growth calls** (63.5%).
  Of these, 589,350 crossed different mapping/island histories; 111,409 rebased
  unrelated inherited boundary lists. These are whole compressed results, not
  representative-only matches. Recorded dependency checks still recompute calls
  when relevant inputs change.
- It reused 1,690,712 of 3,327,124 logical extensions and 862,680 of 1,938,984
  canonicalization calls. Operations have unequal costs; these ratios are not
  claimed as CPU speedups.
- The first ten-seed, same-node four-way ablation measured: repair alone
  **1.179x**, the initial native scheduler alone **1.014x**, and both **1.190x**.
  The initial port copied long frontier keys and repeated no-change admission.
  Cached hashes, immutable-key references and exact no-work-step skipping led
  to the final **1.238x** paired result above. These runs are separate observations.
- Original cProfile diagnostics put roughly 42-60% of long-tail search time
  inside the native growth binding, with substantial Python branch/conversion
  work. Profiling inflates Python call overhead and is **not speed evidence**.
  A final native diagnostic on case 77 / seed 0 measured 8.49 CPU seconds in
  native scheduling/export and 3.00 in Python graph construction, under profiling.
- Symmetry finalization is now **45.6%** of the measured remaining compute CPU.
  We did not add the earlier large shared-symmetry workspace to this comparison.

## Correctness and scope

- Final paired comparison: **414 cut digests**, each aggregating ten individual
  seed graphs: **4140 seed/cut graph comparisons, zero mismatches**.
- Earlier four-way ten-seed ablation: **12,420** seed/cut comparisons, zero
  mismatches. The first one-seed, both-direction pilot checked **1648** graphs,
  also without mismatches.
- Comparisons include every state, transition, fragment, representative mapping,
  compressed symmetry record, finalized group, terminal and cap/stop record.
  They are not merely equal minimum scores.
- Both previously saved vanadium event patterns remain present in every final
  policy, verified using saved mappings only. `mechanism_check_114.json` records
  the exact event signatures, witnesses and source hashes.
- **450 tests passed in 93.58 seconds.** Tests include dependent/independent
  mapping changes, newly available placements, boundary rebasing, cap stages,
  anchors, requested-core stopping, nonidentity atom labels, incomplete endpoint
  composition, symmetry, and native/Python full-DAG equivalence.
- A 21-call kernel probe against the saved pre-edit binary had identical outputs;
  diagnostic CPU totals were 1.819 versus 1.784 seconds. This small probe is not
  evidence of a universal default-path speedup.

These are native XYZ/WBO elementary-step examples, not a retro bank scan or a new
Golden publication run. There are no independent curated reference atom mappings
for this holdout. Historical vanadium patterns are not chemical ground truth.

**Important limitation:** dependency repair is implemented as a reusable
conditional-fragment primitive inside the existing sequential search. We have
not replaced seed-times-cut exploration with a new dynamic-backtracking search
policy. Therefore this experiment does not remove that outer workload multiplier,
prove exhaustive discovery, or justify a 100x forecast. The native scheduler
returns the existing Python `AAMSearchGraph`; the Python abstractions remain.

## Artifacts and reproducibility

`evidence.tar.gz` contains compact per-case/per-seed summaries, complete profiling
tables, selected vanadium mappings, input data, frozen experimental sources,
manifests, native hashes, build/test logs and Slurm accounting. Full compressed
graphs remain in these cluster folders:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/fragment_dependency_20260909_TVriiX
/project/yunhengzou/coordinate_alignment/aam_benchmarks/fragment_repair_pilot_20260909_bHAccT
/project/yunhengzou/coordinate_alignment/aam_benchmarks/fragment_native_paired_20260909_tqGrV4
/project/yunhengzou/coordinate_alignment/aam_benchmarks/fragment_native_final_20260909_JGOrLU
```

Full ten-seed graphs:
`results/<index>/<direction>/<policy>/seed_XX/cut_NNNN.pkl.gz`.
No search rerun is needed for inspection or further analysis. `pairs.json`
defines scheduled tasks; prepared-but-unscheduled task specifications are not
missing results. Frozen drivers and native hashes identify each implementation.

The first folder preserves the pre-edit engine and profile attempts. Initial
profiling failed before search because the driver passed a JSON list instead of
a NumPy WBO matrix; corrected attempts are in `profile_retry`. Early native test
logs also retain a stale-build signature failure and a test-fixture error; both
were corrected before accepted benchmark runs. No failed attempts are included
in the timing tables above.

Scripts: `bench/cut_replay_pilot.py` prepares, submits, compares and archives the
ablations; `bench/profile_fragment_pipeline.py` records phase/function profiles.
Always use a new run directory for fresh computation. For saved-data analysis:

```bash
PYTHONPATH=bench .venv/bin/python bench/cut_replay_pilot.py report --run /project/yunhengzou/coordinate_alignment/aam_benchmarks/fragment_native_final_20260909_JGOrLU
```

Full report path:
`/h/399/yunhengzou/coordinate_alignment/reports/dependency_fragment_repair_20260909/README.md`
