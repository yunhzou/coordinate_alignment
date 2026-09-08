# Fixed-policy Golden benchmark

Status: complete, with timeouts/unresolved checks explicitly retained. This is a fresh benchmark, not a union with earlier diagnostic runs.

## Results

All percentages use all 1,851 reactions as denominator.

| Mode | Mapping-family recovery | Verified absent | Unresolved | Representative top-1 | Representative top-5 |
|---|---:|---:|---:|---:|---:|
| Smaller to larger | 1821 (98.38%) | 17 | 13 | 1437 (77.63%) | 1786 (96.49%) |
| Bidirectional | 1839 (99.35%) | 9 | 3 | 1434 (77.47%) | 1789 (96.65%) |

Bidirectional search adds 18 certified recoveries. Verified misses: 7, 19, 871, 986, 1228, 1285, 1377, 1475, 1553. Unresolved: 590, 1358, 1799. These are not claimed to be impossible under further search.

All 3,702 directional attempts are accounted for: 3,689 completed searches and 13 five-minute watchdog timeouts. All 1,851 mode-comparison workers completed successfully. One additional whole-family verification subprocess timed out. No unresolved record is dropped from accuracy denominators.

Among completed directional searches, elapsed time including persistence has median 4.156 s, mean 15.956 s, and p95 75.861 s. Net search CPU totals 325,633.278 CPU-seconds over those 3,689 searches; timed-out work is additional and separately visible in accounting. CPU time is not elapsed latency. The frozen benchmark does not measure concurrent bidirectional latency under equal per-reaction resources.

`analytics.json`, `per_case.csv`, and `per_direction.csv` contain full metrics and timings. `summary.json` is the compact overview. Top-five compressed-family recovery is distinct from representative top-five: 1797/1851 single-direction and 1795/1851 bidirectional, with 8 and 14 unresolved respectively. The separately reported bond-event-window metrics are not top-k accuracy.

Full artifacts:
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_publication_20260908`

Frozen search/analysis engine commit: `bb0a7d7`.

## Protocol

- All 1,851 Golden reactions; 3,702 fresh directional searches.
- Single mode: smaller explicit-atom endpoint to larger; R to P on ties.
- Bidirectional mode: union of both directions, without splicing branches.
- Ten seeds per cut, branch cap 100, tolerance 1.0, mature cut sweep and symmetry finalization, explicit hydrogens.
- Fixed root seed 42 and deterministic cut-derived seed orders; hash seed 0. Actual seed-order digests and compressed search contexts are saved.
- Reference-blind ranking by product-heavy coverage, explicit-atom coverage, then bond events. No reference-guided selection or additional energy/chirality postfilter.
- Report representative top-k, event-window recovery, and compressed-family reference recovery separately. Reference accuracy compares heavy-atom relations modulo endpoint chemical symmetry; it does not score hydrogen identity.
- All records stay in accuracy denominators. Unknown verification, timeouts, cap hits, and incomplete rankings are explicit.

## Timing and resources

Search tasks use 16 CPUs and 32 GiB on `cpunodes_nia`, up to 160 simultaneous tasks (2,560 CPUs), subject to availability. Ranking and verification run in separately scheduled single-CPU tasks. Search subprocess watchdog: 300 seconds; job limit: 10 minutes.

CPU time excluding measured checkpoint encoding/compression/writing and loading/decoding is the primary additive compute measure. Raw elapsed time, persistence, ranking, and reference verification are logged separately. Summed worker saving time is **not** subtracted from parallel elapsed time. IPC, process setup, profiling, and metadata overhead remain in compute accounting. Summed direction elapsed times are not measured concurrent bidirectional latency. Timed-out tasks are censored, not counted as zero-time successes.

The artifact root contains pinned inputs, source/native binary hashes, versions, configuration, hardware, Slurm jobs, compressed intermediate/final archives, witnesses, and per-process timing logs. Final postprocessing writes `publication/analytics.json`, `publication/per_case.csv`, `publication/per_direction.csv`, and `publication/README.md` without rerunning AAM.

Search arrays: 443233, 443235, 443237, 443239.
Analysis arrays: 443234, 443236, 443238, 443240.
Mode-combination arrays: 443241, 443242.

These settings are not a claim that this fresh run reproduces the earlier 1840/1851 diagnostic-union result. Its accuracy will be measured independently.

## Scheduler recovery

Sixteen original search tasks remained in Slurm `CONFIGURING` for over an hour, without creating search directories. Only those unstarted tasks were cancelled and replaced (arrays 447646, 447647, 447648), excluding the affected nodes. The existing analysis dependencies were redirected to the replacements. Completed searches were not rerun and the frozen engine/configuration were not modified. The full action log is `reporting/recovery.json` in the artifact root; replacement allocations are included in final resource accounting.

Thirteen other searches reached the fixed five-minute subprocess watchdog. Their partial cut checkpoints remain saved; they must not be described as completed negative searches. Their original outcomes remain visible in the fixed-budget benchmark.

Further scheduler-only recovery is recorded in the same action log: unstarted analysis/comparison tasks were replaced, and concurrency adjusted to available capacity. Once all 3,702 analysis completion records existed, comparison dependencies were released despite stale Slurm cleanup states. Recurrent short-job startup stalls prompted dispatching the remaining 557 unstarted comparisons through eight unchanged workers in one persistent allocation (452121). No completed comparison was selected for this batch. The final report is gated on actual worker completion records rather than scheduler cleanup alone.

The persistent batch completed in 9m31s, and report job 446266 completed in 1m58s. These are scheduling/postprocessing durations, not AAM search latency. The final summary classifies missing top-five rankings from completed, incomplete searches as `unknown`, not `pending`; this reporting-label correction changes no search or mapping result.
