# Opposite-direction search on all 21 unresolved Golden cases

**10 newly recovered; accumulated recovery 1,840 / 1,851 = 99.41%.**

New recoveries: **603, 780, 833, 845, 865, 867, 1380, 1568, 1574, 1793**.

Not recovered in the returned opposite-direction families:
7, 19, 590, 871, 986, 1228, 1285, 1358, 1377, 1475, 1553.

No pending cases, unknown verification results, watchdog terminations, or failed
processes remain. Twenty cases hit the branch cap (all except 986); negatives
therefore do not establish exhaustive unrepresentability.

## Scope and verification

- Twenty new searches; the already completed case-590 reverse run was reused.
- Opposite of each saved default direction, not uniformly P→R. The default was
  smaller explicit-atom endpoint first, R first on ties.
- Same 10 random seeds, cap 100, tolerance 1.0, explicit hydrogens, and ordinary
  single-edge sweep. Each direction cuts its own source edges.
- No core AAM changes, guided anchors, or reference-directed seed selection.
- These are ground-truth heavy-atom mapping recoveries, modulo endpoint chemical
  symmetry. They are not atom-coverage percentages or top-1 ranking scores.
- Nine recoveries are present among stored representatives; 865 was certified
  inside its compressed family with a full explicit-atom realization.
- All ten recovered mappings were independently rechecked after conversion to
  original R→P coordinates: injectivity, element correspondence, and reference
  chemical certificates. Certificates also match those derived from the original
  mapped reaction strings. See `validation.json`.
- Case 1793 now has a complete saved search and certified recovery: 71 mapped
  explicit atoms from 148 source atoms, preserving the benchmark's partial
  mapping semantics rather than forcing the extra atoms into the target.

The 99.41% figure is a **diagnostic union with previous saved experiments**;
it is not a new full-dataset score for one fixed bidirectional configuration.

## Resources and persistence

Slurm array **443210**, 16 CPUs / 48 GiB requested per case, at most 10 cases
concurrent (up to 160 CPUs). All 20 tasks completed successfully. Search watchdog
300 seconds; verification watchdog 260 seconds; each allocation limited to
10 minutes.

From first worker start to last worker finish: **279.63 seconds**. Across the
20 new cases, median search-through-archive time was **12.64 seconds**, maximum
**125.24 seconds** (1793), and sum **437.72 seconds**. These are per-case wall
times, not summed CPU time. Verification is additional and included in campaign
elapsed time. Do not compare these directly with runs on different allocations.

All raw cuts, finalized cuts, complete AAM archives, evaluation results, realized
mapping witnesses, job statuses, and the frozen engine are retained at:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_reverse_remaining21_20260908`

Per-case complete archive: `cases/<index>/cuts/aam.pkl.gz`.
Per-case reference result and witness: `cases/<index>/evaluation.json`.
Case 590 is a link to its previously completed archive; it was not recomputed.

Repo evidence: `summary.json`, `cases.json`, `validation.json`, `manifest.json`.
Runners: `bench/golden_reverse_campaign.py` and `bench/golden_forced_direction.py`.
Direction/evaluation regression tests: 15 passed before submission.
