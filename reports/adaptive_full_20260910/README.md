# Full-dataset adaptive search confirmation

Status: running; the four-case pilot is not a confirmation of equivalence.

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
