# Full-dataset adaptive search confirmation

Status: full confirmation of the exact shared-policy implementation is running.
The earlier four-case pilot and faster but incomplete adaptive searches are
not confirmations of equivalence.

## Current confirmation

Frozen candidate `2d04f60`, original `98b01b1`. All 1,851 Golden and 140
holdout cases, both directions, identical original ten seed orders per cut,
full sweep, cap 100, explicit H, matching tolerance 1.0 and scoring 0.5.
Both implementations now receive eight CPUs per reaction. The current batch
campaign allocates at most 1,024 CPUs and retains independent 300-second
search and 240-second evaluation watchdogs.

The new backend shares identical conditional decisions across the original
policies. A necessary element-capacity check skips growth only when no unused
target atom of the seed element exists. Frontier admission and both cap stages
remain unchanged. Original native fragment growth is not replaced.

The diagnostic Golden 1636 R-to-P uncut profile fell from 21.483 to 0.857 CPU
seconds with 1,124 states and 1,123 transitions in both runs. This is a local
profile, **not** a claimed full-benchmark speedup. Focused tests compare complete
correlated fragment records, states and stop reasons against the original.
The complete repository suite passed: **518 tests in 203.61 seconds**, Slurm
job `466076`; the unabridged log is retained in the preceding campaign folder
as `capacity_full_tests_466076.out`. Saved-certificate verification tests include
both orientations, unequal composition, partial reference relations and
different explicit-H witnesses within one heavy-atom class.

Current full artifacts:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/shared_capacity_full_20260910`

The preceding shared per-cut campaign and its explicit checkpoint continuations
remain at `../shared_cut_workers_full_20260910`. No completed cut was rematched
to regenerate reporting. Saved full-atom ranking witnesses can now certify
representative recovery without loading and rescoring the whole archive; absent
representative certificates still require the unchanged compressed-family query.
Search never receives reference labels. Fresh paired timing excludes continued
archive construction, failures, loading and persistence; queue time is separate.

## Earlier adaptive campaign (retained experiment)

Every proposed search optimization must be checked against the original on
**all 1,851 Golden reactions and all 140 XYZ/WBO holdout steps**. Small examples
are debugging tests only. Do not claim no regression from those examples.

## Frozen protocol

- Original: commit `98b01b1`, rebuilt native extension, mature reused-native
  sequential pipeline, ten seeds per cut, full cut sweep, branch cap 100.
- Experimental: balanced adaptive seed agenda from `c5b74d1`, no sweep,
  saved budgets of 400, 1,600 and 6,400 growth operations; 30 CPU-second
  search budget. It retains the native growth cap of 100, but does **not**
  impose the original synchronized live-frontier cap. This is an explicitly
  different search policy, not an assumed output-equivalent optimization.
- Both: explicit hydrogens, matching tolerance 1.0, scoring tolerance 0.5,
  fixed root seed 42, hash seed 0. Search never reads reference labels.
- Both directions for every record: 3,982 directional comparisons,
  7,964 search calls. Smaller-to-larger and bidirectional unions reported
  separately. Branches are never spliced between directions.
- Original uses eight CPUs per call; adaptive one CPU. Same partition,
  hardware recorded per call. Compare CPU cost; elapsed timings are **not**
  an equal-resource speedup comparison. Maximum aggregate allocation is
  1,024 CPUs, subject to cluster availability.
- Five-minute external search watchdog, four-minute analysis watchdog,
  five-second hard-kill grace, ten-minute Slurm allocation. No automatic
  retries, budget escalation or substitution of a different algorithm.

## Confirmation metrics

Golden uses the unchanged original heavy-atom reference/symmetry evaluator,
including partial annotations. Report compressed-family recovery and
representative bond-event windows separately, with all 1,851 in denominators.
Unknown verification and timed-out search are not successful negatives.

The holdout has **no annotated ground truth**. Compare full mapping feasibility,
full-H event counts, and original score-response heavy-atom pattern classes.
Missing representative patterns require symbolic compressed-family queries;
they are not automatically missing from the returned compressed results.
The pilot's alternative element-pair equivalence is not substituted for the
original comparison metric.

Search CPU excludes measured checkpoint encoding/compression/writes and
loading. Evaluation is separate; parallel elapsed includes persistence but
not queue time. Failed/censored work stays in status/accounting records and
is never silently counted as zero-time success. All intermediate cuts,
adaptive snapshots, final AAM graphs, evaluations and hashes are retained.

Full artifacts:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/adaptive_full_20260910`

Driver: `bench/adaptive_full_benchmark.py`. `compare` reads saved artifacts
without rerunning AAM. The final report must identify recovery losses as well
as speedups before any decision to replace the mature pipeline.

Original search arrays: `457271`–`457274`; adaptive arrays: `457275`–`457278`.
Search engine/driver freeze: `0ce0f17`; subsequent reporting and monitor code
does not alter that frozen search or its evaluations.

A separate bounded campaign monitor refreshes the full comparison from saved
outputs. It replaces only allocations stuck in `CONFIGURING` for ten minutes
that have never started a worker, using the same frozen worker on `cpunodes`.
Started/failed calculations are not silently rerun. Scheduler actions are
journaled separately. The monitor's two-hour allocation is not an AAM run;
each AAM subprocess retains its independent five-minute watchdog.
