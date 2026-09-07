# Golden AAM benchmark — measured campaign

Run: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_full_20260906`. Search revision: `9dd6d8d2bb4ae20ec6e6482743d8727655bb724d`. Frozen source hashes and inputs: `manifest.json`.

Completion campaign: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_completion_20260906`; 152 records completed from saved evidence. These are not cold rerun timings.

## Recovery and coverage

| Measure | Count |
|---|---:|
| Records / evaluated | 1851 / 1851 |
| Complete heavy-atom reference annotations | 1760 |
| Complete-reference family recovery | 1595 / 1760 |
| Selected-representative correct / scored eligible | 1302 / 1760 |
| Partial-reference consistency (not full mapping accuracy) | 88 / 91 |
| One mapping covers all P heavy atoms | 1758 / 1851 |
| One mapping covers all P atoms, including explicit H | 1551 / 1851 |
| Incomplete searches | 0 |

Unresolved scoring and unscored records remain in `summary.json` and `cases.csv`; they are not counted as successes. A positive recovery from a completed cut is sound even when later search timed out. Its coverage is a lower bound. Cap-limited and seed-limited searches do not prove global exhaustiveness. Coverage is not a union of incompatible branches and does not establish chemical correctness.

## Measured performance

| Phase | Completed measurements | Median (s) | p95 (s) | Maximum (s) |
|---|---:|---:|---:|---:|
| Original completed searches including checkpoints | 1677 | 5.036 | 111.233 | 229.008 |
| Original additional gzip archive | 1677 | 2.062 | 64.620 | 136.328 |
| Checkpoint completion, not cold search | 152 | 129.523 | 255.433 | 284.549 |
| Evaluation | 1851 | 0.638 | 19.370 | 120.869 |

Original search watchdog: 300 s; outer task watchdog: 600 s. Original CPU allocation ceiling: 256; one AAM worker per reaction. 14 searches reused from the saved pilot.

Completion attempts use a 300-second finish watchdog and 270-second scoring watchdog within a 600-second outer watchdog. Frozen completion engine manifests record each attempt’s cut-worker and node allocation. Completion CPU statistics are parent-process only; scheduler accounting includes child processes. Completion latencies describe the successful resumption, not the cumulative cost of previous attempts.

Completed-search latency excludes timed-out searches and additional gzip serialization; report those separately, not as a cold full-run speed claim. Phase metrics are diagnostic, can overlap, and should not be blindly added. Slurm measured CPU and allocated CPU-seconds are different quantities. Peak RSS in the case table is the search process high-water mark; raw scheduler accounting also retains scoring/task memory.

## Evaluation scope

Full explicit-H graphs are mapped, with original reference map labels removed before canonical input ordering. Heavy-atom reference agreement is assessed modulo joint, chemically colored endpoint automorphisms. Family queries preserve correlated group actions and do not enumerate all bijections. Top-1 refers to the reference-blind selected terminal representative, not the best possible member of every compressed family. Reference annotations do not establish identities for automatically added H atoms; all-atom coverage is not H-mapping accuracy. This is graph AAM with discrete bond orders, not a 3D stereo-validation or mechanism-grouping benchmark.

No external baseline has been rerun here. Published Golden subsets/reference corrections differ; these results are not yet a like-for-like paper comparison.

## Saved evidence

- Every input, separate reference, status, search log and available complete compressed AAM archive.
- `.pkl.gz` files are trusted internal typed-graph snapshots: load with `read_aam_checkpoint` and the frozen producing engine, never from untrusted sources. They retain the same AAMResult; JSON export remains available.
- Atomic per-cut checkpoints, including for interrupted searches; partial manifests identify their exact files.
- Per-case evaluations and symbolic witness actions when required; source/configuration hashes and Slurm accounting.
- `cases.csv`: every dataset row, including failures. `summary.json`: stratified denominators and timings.

Recreate this report without remapping: `python bench/report_golden_campaign.py --run RUN --output OUTPUT`.
Rescore a completed archive without remapping: `PYTHONPATH=RUN/engine/src python RUN/engine/bench/golden_campaign.py score --run RUN --index INDEX`. For completion snapshots use the producing engine’s `complete_golden_campaign.py score` command instead. Preserve the previous evaluation before intentionally rescoring; use the frozen engine for exact reproducibility.
