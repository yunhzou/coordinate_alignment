# Exact AAM reuse: measured improvement and rejected optimization

**1.58x compute-CPU speedup (36.5% less compute)** over the original scheduler,
with **zero mismatches in 24,840 complete cut-graph comparisons**. This is a
four-case, full ten-seed, two-repeat stress experiment, not a new Golden accuracy
benchmark. It does not achieve the requested order-of-magnitude improvement.

## Full ten-seed experiment

Each of 80 paired tasks runs one of the original ten seed orders across every
cut. Modes within a pair use the same node, one CPU, fresh processes and
counterbalanced ordering. Seeds, cuts, branch cap 100, explicit H and tolerance
1.0 are unchanged. Comparisons include every state, transition, ordered symmetry
generator and cap/stop record; they do not compare only sampled bijections.

| Mode | Compute CPU seconds | Search | Symmetry |
|---|---:|---:|---:|
| Frozen original scheduler/native growth | 2552.72 | 1609.36 | 942.73 |
| Prior native scheduler + whole-fragment reuse | 2020.36 | 1117.32 | 901.47 |
| Prior reuse + shared Python symmetry workspace | 2316.41 | 1315.84 | 999.01 |
| Prior reuse + native conditioned symmetry | 1620.40 | 1213.70 | 404.51 |

The additional gain over the **prior reuse optimization** is **1.25x**, not
1.58x. Symmetry processing itself is 2.33x faster than the original. Sharing the
existing Python cache alone did not help relative to prior reuse. Both repeats
agree: original/new totals were 1289.05/821.09 and 1263.67/799.32 CPU seconds.

Per-case compute, averaged across two repeats of the complete ten-seed policy:

| Case / direction | Original CPU s | Optimized CPU s | Speedup |
|---|---:|---:|---:|
| 25: Pd hydroamination TS11 / R to P | 514.23 | 264.79 | 1.94x |
| 76: Suzuki Ni TS11 / P to R | 502.32 | 366.45 | 1.37x |
| 77: Suzuki Ni TS14 / P to R | 257.18 | 177.08 | 1.45x |
| 114: V dodh TS910 / R to P | 2.62 | 1.89 | 1.39x |

These are **summed compute CPU times**, excluding queue, imports, input loading,
audit encoding and persistence. They include setup, search and symmetry but not
the public API's cross-cut merge or downstream pattern collection. They are not
measured parallel wall times. Maximum worker RSS was 1501 MiB original versus
986 MiB optimized (prior reuse: 870 MiB); these are workload maxima, not total
cluster memory or cache sizes. Shared-Python-mode RSS is unavailable.

## What changed

- **AAM execution:** the existing native sequential scheduler uses
  dependency-validated whole-fragment results across cuts/seeds within a worker.
- **Symmetry conditioning:** one native target graph reuses exact role-color
  strings and ordered partitions, returning the same correlated generators.
  Python fragments, branches and the final search-graph representation remain.
- **Individual growth-extension cache:** implemented and tested, but rejected
  for this backend because it costs more than it saves. It is disabled.

Profiling 309 group requests on one saved long-tail cut found only about 13% of
instrumented symmetry time inside pynauty autgrp; most was preparation of role
labels/partitions. The optimization targets that preparation. It is **not** a
general permutation-group stabilizer algorithm, nor enumeration of automorphisms.
Native conditioning has a 64 MiB cache budget; this is not a process-memory cap.

The extension experiment caches raw children before the existing frontier
deduplication/cap decisions, including negative answers and all relevant
constraints. Source-edge reads remain dependency-tracked. Even with refined
immutable keys, repair alone took 98.91 CPU s versus 105.69 with extension caching;
conditioned symmetry alone took 77.34 versus 82.46 with both. This refined pilot
has another 2,484 exact comparisons with no mismatches. Its unchanged default
engine took 122.92 versus frozen baseline 121.16 CPU s (+1.45% in this small run);
this does not establish a systematic regression or prove zero overhead.

## Public API integration

Opt in explicitly; the default and publication protocol have not changed:

```python
result = search_aam(
    problem, config, workers=10,
    execution="reused_native",
    intermediate_dir=output_directory,
    archive_format="checkpoint",
)
```

This requires the built native engine and never switches algorithms on failure.
The returned object is still `AAMResult`; cross-backend checkpoint resume is
tested. The individual-extension cache remains off.

A separate full public-API check uses 10 CPUs, all ten seeds and all cuts, saving
every cut plus the final combined result. These two cases are integration checks,
not a broad end-to-end speed claim:

| Case | Original/new compute CPU s | Original/new elapsed s |
|---|---:|---:|
| 77: Suzuki Ni TS14 | 259.55 / 194.83 | 54.00 / 45.59 |
| 114: V dodh TS910 | 1.09 / 0.92 | 0.252 / 0.245 |

Compute includes parent/worker CPU, restoration, IPC and merge, minus measured
checkpoint-write CPU. Actual elapsed **includes saving**; summed worker-write
wall times are not subtracted from elapsed. Both exclude subsequent audit hashing.
Raw timing records retain both inclusive and write-excluded CPU counts.

Both final combined graph hashes agree exactly: case 77 contains 501,736 states
and 54,555 terminal records; case 114 contains 1,955 states and 463 terminal
records. Case 77 is branch-capped in **both** engines, identically; equivalence
does not mean exhaustive search. Audit hashing ran after mapping and persistence
and is not included in the timing table.

## Reproducibility and evidence

Final complete suite: **470 passed**, 23 multiprocessing/fork deprecation
warnings, 105.23 s. Targeted tests cover both worker counts, tolerances 0.5/1.0,
tiny and normal caches, correlated ring symmetry, partial composition, anchors,
different cut/graph floors, cap records and cross-backend checkpoint resume.
Earlier test attempts exposed an inherited Slurm environment assumption in a
local-worker fixture (fixed), then missing files in the frozen test bundle
(bundle corrected). Their logs are retained; these are not mapper failures.

Baseline core: `f158864`, frozen from the real-TS run. Experimental stages:
`446d82a`, `92a8035`, followed by the public-API integration committed with this
report. Frozen source/binaries, input hashes, settings, manifests, original and
resumed job logs, timings, full graphs and compact analyses are retained.

One benchmark statistics-export bug occurred **after** shared-Python-mode graphs
and phase timers had been saved. The reporter was repaired, summaries recovered
from those checkpoints, and only previously unstarted modes ran on the same
original node. No completed matching was rerun. Missing recovered RSS/cache
counters are explicitly null. Consequently campaign job elapsed time includes
interruptions and is not used as compute time. All mapping work units have a
300-second watchdog; allocations have a ten-minute limit.

Full artifacts on the cluster (never rerun mapping just to inspect them):

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/conditioned_reuse_pilot_20260910_KciqQq
/project/yunhengzou/coordinate_alignment/aam_benchmarks/conditioned_reuse_refined_20260910_QdViEx
/project/yunhengzou/coordinate_alignment/aam_benchmarks/conditioned_reuse_full_20260910_0o09px
/project/yunhengzou/coordinate_alignment/aam_benchmarks/public_reuse_20260910_kgCQdQ
```

Drivers: `bench/conditioned_reuse_pilot.py` and `bench/public_reuse_pilot.py`.
Compact reports and provenance accompany this file; large raw graphs and frozen
engines remain at the paths above. The outer seed-by-cut multiplier is still
present; these measurements do not support a 10–100x claim from caching alone.
