# Full-AAM core regression: measured fix

The refactor's approximately 2.3× core slowdown was Python bookkeeping, not
additional chemical matching. The fix restores old-core throughput while keeping
the complete fragment-decision DAG and the reusable Python fragment API.

## Controlled comparison

59 saved, actual augmented inputs from precursor `INVENTORY-000400` against the
BIAN target were replayed unchanged. This is full augmented AAM, **not** the old
approximate retro scan. Each version received identical input bytes, explicit
hydrogens, branch limit 100, tolerance 0.5, anchors, seeds, and no cut sweep.

Each case was measured three times in a fresh process, rotating version order on
`bosque20`. Each call used one CPU; eight independent calls ran concurrently in
a 16-CPU allocation. Timers surround `find_islands` only. Input loading, output
validation, persistence, symmetry finalization, and assembly are outside these
core timings. Every Slurm job had a ten-minute limit and a 580-second watchdog;
individual benchmark subprocesses had a five-minute limit.

| Version | Sum of per-case median core time | Slowest case median | Maximum process peak during core |
| --- | ---: | ---: | ---: |
| Old full-AAM core, `4a7b57f` | 88.78 s | 4.35 s | 83.0 MiB |
| Pre-fix refactor, `4dbb493` | 204.70 s | 7.46 s | 198.8 MiB |
| Fixed core | 87.84 s | 4.35 s | 106.5 MiB |

This is **2.33× faster than the pre-fix refactor**, approximately the old core's
speed, while retaining substantially more history than the old terminal-branch
list. The sum is **not** a parallel wall time, a full-bank scan time, or a full
retro pipeline time. The fixed DAG still has a memory cost relative to the old
terminal-only representation; it is not an identical-output memory comparison.

One cap-only case (45; one state, no transitions) initially measured 6.6%
slower. An isolated nine-repeat paired check did not reproduce that difference:
old median 0.37537 s, fixed 0.37555 s. The original outlier and all follow-up
measurements are retained rather than excluded from the totals. Timing noise
does not justify claiming every individual observation must be faster.

The same existing C++ binaries were used for all versions. No native engine,
growth rule, seed schedule, tolerance, branch admission rule, symmetry domain,
or cut-sweep policy was changed by this fix.

## Causes and changes

Profiling case 2 found exactly 1,319 growth calls in both versions, producing
7,740 fragment commits. The refactor added repeated whole-graph edge scans,
deep copies, partition reconstruction, and redundant cumulative state copies.

- `FragmentPlacement.from_match` obtains induced bonds from fragment adjacency
  once. The DAG consumes these bonds instead of scanning every source edge again.
- Internal growth records transfer ownership into the DAG. They are not mutated
  afterward; the public caller-owned placement conversion still copies at its
  boundary. Compressed symmetry is retained, not expanded into bijections.
- Sibling target placements reuse an exactly keyed source-partition merge.
  Their distinct mappings and history transitions remain distinct.
- Cumulative snapshots share unchanged immutable atom pairs. One canonical
  snapshot serves both DAG storage and frontier deduplication, with a cached hash.
- Product-island labels are derived from source labels and the witness. The
  unused eagerly maintained duplicate table is now an inspection property.
- Optional event records retain their previous contents and ordering. Their
  construction no longer runs when tracing is disabled.

There are no new search caps, approximation rules, fallback paths, or discarded
capped histories. Public result fields and persistence schemas are unchanged.

## Correctness checks

- All 59 cases, across all three repetitions: exact equality of every saved
  context, root, state, transition, fragment/symmetry record, preserved bond,
  and stop/cap record against the pre-fix DAG.
- All 59 terminal mapping relations, source-fragment partitions, and deferred
  cuts agree with the old full-AAM core after normalizing arbitrary island IDs.
- A fresh fixed-core result for case 2 was saved and passed through the actual
  fragment consumer. All **26,487 occupations** agree with the saved reference,
  including correlated hierarchy/actions, derivations, exact groups, caps and
  history structure. This check does not substitute counts for content equality.
- Full multi-seed/cut-sweep workflows on TEMPO, tetraphenylmethane, and the
  symmetry-heavy tetra-tert-butyl example retain identical graphs and grouped
  mechanisms. Serial AAM-wrapper-plus-grouping medians were respectively
  17.91 → 17.42 s, 3.32 → 3.23 s, and 73.40 → 71.08 s. These broader timings
  include optional grouping and must not be confused with core-only timings.
- **250 tests pass**, including native-versus-Python growth comparisons,
  symmetry, partial composition, sparse indices, caps, parallel/serial equality,
  exact assembly, and persistence. Added tests independently check 600 commits
  across 20 randomized sequences against literal pre-fix partition/event logic, fragment bonds, and
  sibling snapshot ownership. Pytest discovery is scoped to `tests` so saved
  baseline worktrees cannot shadow the installed package during collection.

These checks establish no detected correctness regression on the tested corpus,
not a proof covering every possible molecule. Per-case timings and all repeats
are included in [the machine-readable evidence](aam_core_regression_2026-09-05.json).

## Reproduction and saved evidence

Use `bench/core_aam_regression.py suite` with the captured input directory,
`--baseline` pointing to `4a7b57f`, and `--before` pointing to `4dbb493`.
Both baseline worktrees must use the same native extensions as the current
package. `replay --profile` isolates a case; `workflows` exercises the broader
multi-seed/cut-sweep examples.

Cluster-relative evidence under `data/retro_runs/`:

- Original inputs and immutable raw-graph references:
  `aam_memory_boundary_20260905/inputs/`.
- Final three-version measurements, individual logs, and input hashes:
  `aam_core_regression_20260905/verified/`.
- Profiles and intermediate optimization trials:
  `aam_core_regression_20260905/{profile,fix1,fix1_profile,fix2_profile,fix3,fix4,final_profile}/`.
- Full workflow graphs and mechanism snapshots:
  `aam_core_regression_20260905/workflows/`.
- Fresh core graph, complete consumer result, and content comparison:
  `aam_core_regression_20260905/consumer/`.

Large intermediate files remain on the cluster; the compact measurements,
replay tool, tests, and explanation are versioned with the fix.
