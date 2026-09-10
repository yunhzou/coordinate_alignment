# Adaptive sequential fragment choices: bounded pilot

**Outcome: useful prototype, not a demonstrated replacement for sweep-cut.**
Native growth can expose earlier fragment boundaries cheaply and continue from
saved conditional states. It found both historical vanadium event patterns
without external cuts. Across the four test cases, however, no tested fixed
priority recovered every saved full-sweep best score. The difficult Ni case
remained at 11 events versus the saved sweep result's 8.

Default AAM search, seed selection and sweep-cut policy are unchanged. The
experimental scheduler is isolated on `experiment_adaptive_fragment_choices`.
Baseline before this work: `98b01b1`; initial implementation: `e914b3a`.

## What changed

1. One normal native fragment-growth call retains compressed frontiers before
   extensions. With the observer disabled, it follows the original path.
2. A lazy `close` finalizes one earlier frontier using the existing fragment
   finalizer. It does not regrow the prefix or delete a source bond. Boundary
   bonds become deferred constraints, not automatically counted bond breaks.
3. A persistent sequential scheduler continues under the locked mapping,
   fragment partition and boundary constraints. Exact equal states can share
   continuation; incoming correlated symmetry paths remain in the search DAG.
4. Earlier closures are explored in a bounded agenda. Full completed results
   and pending-work descriptors are saved at each checkpoint.

No search-family pruning uses a representative's event score. No explicit
atom-bijection space is enumerated by this search. The saved representative
mappings used for scoring are not the entire compressed result.

## Initial experiment

One CPU per case/policy on bosque6; fixed order generated with RNG seed 42;
explicit H; matching tolerance 1.0; event threshold 0.5; native growth branch
cap 100; no external cuts. The table gives best **representative event score**
and cumulative **compute CPU seconds** after 1,600 grow/close operations.
CPU includes setup, search and checkpoint symmetry finalization; it excludes
persistence and separate diagnostic scoring.

| Case | Saved 10-seed full-sweep best | Largest-first: events / CPU s | Smallest-first: events / CPU s |
|---|---:|---:|---:|
| 25, Pd, R→P | 5 | 5 / 2.328 | 5 / 3.219 |
| 76, Ni TS11, P→R | 3 | 21 / 1.565 | 3 / 4.809 |
| 77, Ni TS14, P→R | 8 | 11 / 1.702 | 19 / 2.305 |
| 114, V, R→P | 4 | 4 / 0.688 | 4 / 0.585 |

Running **both** initial policies would reach the saved best score on 3/4
cases, but its cost is the sum of both runs. Selecting the best policy per
case after seeing results is not an independently validated algorithm.
Equal event scores do not establish equal mapping-family coverage.

For case 114, smallest-first's saved representatives include **both** historical
4-event signatures in `mechanism_check_114.json`. These are historical outputs,
not independently annotated chemical ground truth. Their witness ordinals and
the raw checkpoint paths are recorded in `results.json`.

## Bounded follow-up tests

- Raising the two original policies to 25,600 operations / 30 seconds of
  compute produced no improvement in best scores. The largest Pd run reached
  119,143 pending choices and about 2.1 GiB peak RSS. More choices are not free.
- The original depth priority starved deeper revisions: Ni TS14 had 19,535
  second-level closure choices waiting behind first-level work.
- Event-guided ordering alone did not solve this. It orders closures by bond
  discrepancies touching released atoms; this is a heuristic, not an event
  lower bound or a proof that a fragment is wrong.
- Fair-depth scheduling removes that strict level starvation, but the final
  6,400-operation test still gave scores **5, 21, 11, 4**, costing respectively
  **6.601, 6.558, 4.913, 2.623 CPU seconds**. It can also delay completing the
  initial greedy path. It did not improve recovery.
- Fifteen-second, correlation-preserving family diagnostics on initial and
  fair-depth final checkpoints did not improve these minima. They inspected
  selected paths, not all paths, so absence is not a proof of unrepresentability.

Tuning stopped after these tests. No larger-budget loop, automatic retry or
production-default change was made. The evidence does not justify an
equal-accuracy speedup claim against either full AAM or SLAP.

## Correctness and resource checks

- **483 tests passed**, 26 existing multiprocessing/fork deprecation warnings,
  97.18 seconds. No test failure. Full-suite Slurm job: `457020`.
- **80 saved graph checkpoints validated**: acyclicity, injective and
  element-compatible states, monotone locked mapping, and preserved-bond
  compatibility at tolerance 1.0. This is not an exhaustive generator proof;
  correlated-family diagnostics and symmetry regression tests are separate.
- Default single-seed/no-cut runs compared against the frozen pre-experiment
  engine on the same pinned CPU, three trials per case: **all 12 graph digests
  exactly equal**. Median current/baseline CPU ratios were 0.992, 1.004, 1.016,
  1.015. No material slowdown observed in this small regression check; it is
  not a rerun of the entire benchmark. Job: `457039`.
- Every pilot worker used external `timeout --kill-after=5s 300`, plus a
  ten-minute Slurm allocation. Search additionally had a finite operation
  budget and a 30-second compute soft limit. Snapshot finalization can extend
  beyond the soft limit; the external watchdog still applies.
- Watchdog smoke test: a subprocess ignored SIGTERM and attempted a ten-second
  sleep. With 0.2-second timeout plus 0.2-second kill grace, it was killed by
  SIGKILL in 0.404 seconds. Native code cannot bypass this external watchdog.
- New submissions disable Slurm requeue. Audits and the full regression suite
  also have five-minute external watchdogs and ten-minute Slurm limits.

The native **per-growth** branch cap remains 100 and capped results are marked.
This experiment does **not** impose the mature scheduler's synchronized global
live-frontier cap on its agenda. Time/work budgets stop agenda processing and
expose pending work. All tested final snapshots had work pending. One deeper
Pd/smallest-first run hit the native growth cap; see the raw rows.

Other explicit limitations: one fixed seed order; anchors and external growth
replay are unsupported and rejected; no symmetry-repair pass; snapshots save
results and pending descriptors but not native sessions for cross-process
resume. In-process continuation is supported. Early closures do not constitute
an exhaustive alternative to all cuts, seeds or possible growth decisions.

## Evidence and reproduction

`results.json` collects all policy/budget rows, graph validations, limited family
diagnostics, default-regression results, test counts, input/native hashes and
frozen implementation hashes. Full checkpoints and witness mappings remain at:

```
/project/yunhengzou/coordinate_alignment/aam_benchmarks/adaptive_fragment_20260910_lRLOXa
/project/yunhengzou/coordinate_alignment/aam_benchmarks/adaptive_fragment_deeper_20260910_mzYH3n
/project/yunhengzou/coordinate_alignment/aam_benchmarks/adaptive_fragment_feedback_20260910_I8Kjaa
/project/yunhengzou/coordinate_alignment/aam_benchmarks/adaptive_fragment_fair_20260910_fKvIsb
```

Mapping jobs: `456959`, `456973`, `456996`, `457015`. Audit jobs: `456988`,
`457022`, `457023`, `457021`. Each policy/case is an independent one-CPU task;
no case consumed the entire node. No watchdog kill occurred in these jobs.

Create a new empty run directory, then use `bench/adaptive_fragment_pilot.py`
`prepare` with explicit `--policies`, `--work-budgets` and `--search-seconds`,
followed by `submit`. It freezes implementation and inputs. Do not reuse an
existing output directory. Audit saved outputs with
`bench/audit_adaptive_fragments.py`; no AAM rerun is needed. The separate
`bench/adaptive_default_regression.py` checks default-path equality and timing.
`bench/publish_adaptive_fragments.py` regenerates this evidence bundle solely
from saved artifacts.

Saved full-sweep score source:
`/project/yunhengzou/coordinate_alignment/aam_benchmarks/fragment_native_final_20260909_JGOrLU/results/{case}/{direction}/independent_native_dependency/summary.json`.
The previous full-sweep results are historical comparisons, not freshly timed
adaptive-equivalent coverage measurements.
