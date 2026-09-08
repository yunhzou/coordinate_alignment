# AAM memory-fix regression testing

## Outcome

The broad comparison found a real JSON checkpoint speed regression, which is fixed in this change. The corrected implementation preserves the exact search graphs in all 83 tested configurations. Core search aggregate time is effectively unchanged (+0.19%); aggregate pipeline cold-search time is 0.45% lower. These are regression checks, not a claim of universal speedup or a new full-bank timing.

## What was tested

- Baseline: pre-memory-fix `c2e7efe`. Initial comparison: `e3e8d2f`. Final comparison: `e3e8d2f` plus the JSON/serial-transport correction committed with this report. The historical manifest's `after_commit` is the base commit, not the exact final source; `source_sha256.json` fingerprints the frozen Python and native source/binaries actually tested.
- 59 saved full-AAM augmented inputs, from one previous BIAN experiment (not 59 independent reactions). Original explicit hydrogens, anchors, tolerance 0.5, cap 100 and seed settings retained.
- Four pipeline problems: TEMPO (57 atoms), tetraphenyl (45), tetratbu (53), and unequal-composition `CCBr.O` to `CCO`, with explicit hydrogens. Each tested with one/four workers and in-memory/JSON/compact checkpoints: 24 configurations. Default three seeds, cap 100, tolerance 0.5 and full cut sweep.
- Three fresh-process repetitions per version/configuration, alternating before/after order on the same Slurm allocation. 498 completed measured processes in the final campaign. Three noisy configurations received seven additional paired repeats each on isolated nodes (42 processes).
- Exact canonical full-graph hashes, state/transition/terminal counts and cap flags; cross-worker and cross-format graph equivalence; warm resume with zero repeated search; JSON and compact write/read round trips; stage timings, own-process peak RSS, and sampled process-tree RSS.
- Per-process watchdog 180 seconds, job budget below ten minutes. Node configuration stalls with no running benchmark were rescheduled; scheduler waiting is not included in the timings.

“Cold search” means no saved search intermediates, not a flushed filesystem cache. Process-tree RSS sums shared pages across processes and is sampled every 0.1 seconds, so it is not an exact physical-memory measurement.

## Regression found and corrected

Streaming `json.dump` made raw checkpoint writing slower than the previous C-backed `json.dumps`. Serial searches also wrote a graph and unnecessarily read it back. In the initial paired campaign, serial JSON cold search increased about 16% on TEMPO and tetratbu, and 70% on the tiny unequal case (only about 58 ms absolute).

The correction restores fast JSON encoding and retains the already-owned graph for true serial execution. Multiprocess workers still persist their output and return references, avoiding graph buffering behind slow sibling tasks. Compact checkpoints still preserve sharing. No seed, symmetry, branch cap, cut sweep, search or matching semantics were changed.

JSON encoding temporarily materializes a string: use compact checkpoints for large memory-sensitive runs. This is an explicit format tradeoff, not a fallback.

## Final measured timings

Seconds, sums of per-configuration medians; not concurrent job wall time:

| Stage | Before | After | Change |
|---|---:|---:|---:|
| Core search, 59 inputs | 87.622 | 87.789 | +0.19% |
| Cold pipeline search, 24 configurations | 149.017 | 148.352 | -0.45% |
| Warm resume, 16 persistent configurations | 23.429 | 22.761 | -2.85% |
| JSON write | 10.194 | 10.179 | -0.15% |
| JSON read | 9.826 | 9.877 | +0.52% |
| Compact write | 2.530 | 2.335 | -7.71% |
| Compact read | 1.773 | 1.601 | -9.72% |

Cold pipeline timings (before → after, seconds):

| Problem | Workers | In memory | JSON checkpoints | Compact checkpoints |
|---|---:|---:|---:|---:|
| TEMPO | 1 | 4.449 → 4.419 | 5.076 → 5.024 | 4.933 → 4.765 |
| TEMPO | 4 | 1.707 → 1.684 | 2.097 → 2.081 | 1.854 → 1.959 |
| Tetraphenyl | 1 | 2.592 → 2.594 | 2.828 → 2.817 | 2.937 → 3.122 |
| Tetraphenyl | 4 | 0.890 → 0.891 | 1.032 → 1.053 | 0.960 → 1.116 |
| Tetratbu | 1 | 27.601 → 27.457 | 30.026 → 29.804 | 29.035 → 28.760 |
| Tetratbu | 4 | 9.321 → 9.286 | 11.154 → 11.161 | 9.986 → 9.827 |
| Unequal | 1 | 0.062 → 0.060 | 0.084 → 0.085 | 0.080 → 0.082 |
| Unequal | 4 | 0.088 → 0.089 | 0.115 → 0.116 | 0.109 → 0.102 |

Do not hide the noisy first measurements: compact TEMPO/four workers and tetraphenyl/one and four workers initially looked 6%, 6%, and 16% slower. Seven additional paired repeats on exclusive nodes gave respectively **1.673 → 1.674**, **2.756 → 2.707**, and **0.973 → 0.969** seconds. Those cold-search regressions did not reproduce.

Small residual stage outliers remain in the main campaign: TEMPO/four-worker standalone JSON write was 0.215 → 0.252 s; tiny unequal/serial JSON resume was 0.034 → 0.041 s. These have not received separate isolated confirmation. The results do not prove every individual IO operation is faster.

## Correctness and memory

Full native-enabled test suite: **374 passed in 91.34 seconds**, using `PYTHONPATH=src:bench RXN_CORE_NATIVE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m pytest -q`. This includes checkpoint transport, mixed-format resume, graph sharing and existing AAM/retro regression tests. `git diff --check` also passes.

All 83 configurations have identical before/after full-graph hashes across all repetitions. Within each pipeline problem, hashes also agree across worker counts and archive modes. All tested warm resumes and JSON/compact round trips preserve those graphs. This verifies preservation of the returned compressed alternatives and histories, not merely one sampled mapping. It does not establish ground-truth chemical accuracy beyond the baseline.

Memory is workload-dependent. Median own-process peak RSS across the 59 core inputs is 75.65 → 75.44 MiB. Tetratbu/four-worker JSON parent peak falls 403.0 → 325.8 MiB (19%), while serial compact rises 416.8 → 447.0 MiB (7%). Serial in-memory also rises 387.2 → 409.5 MiB (6%). Thus this is not a universal RSS reduction; pooling/serialization have bookkeeping costs on smaller inputs. Parallel parent RSS excludes workers; sampled tree measurements are retained separately.

The previously saved pathological case-833 core probe showed 866.95 → 558.32 MiB at cap 2000, with exact graph equality. That was a separate single-repeat measurement, not a replacement for these repeated tests. See `../golden_memoryfix833_20260907`. This campaign does not establish a new 48-worker peak or eliminate the previously observed high memory during parallel verification.

## Reproduction and retained artifacts

Harness: `bench/aam_resource_regression.py`. Summaries and every repetition's timing/count/hash records are in `initial/`, `final/`, and `confirmation/`. The confirmation report intentionally contains only three executed configurations out of the common 83-task manifest.

Full frozen sources, inputs, logs, stage timings, sampled memory, raw cuts and complete result archives remain on the cluster:

```
/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_resource_regression_v2_20260907
/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_resource_regression_final_20260907
/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_resource_regression_confirm_20260907
```

The earlier first harness run used Python container equality after JSON round trips; tuple/list normalization caused assertions in both versions. Those pipeline runs are excluded. All reported comparisons use canonical full-graph equality.
