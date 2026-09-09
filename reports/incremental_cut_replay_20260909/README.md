# Exact incremental cut replay: experimental results

**Outcome:** replay preserves the tested compressed results, but is not a major
speed breakthrough. The same-node comparison gives **1.20x** compute-CPU speedup
when all cuts share seed orders; preserving the existing cut-specific orders is
**4% slower**. Production search defaults are unchanged.

Branch: `experiment_incremental_cut_replay`, based on `3e9a70a`.
This is an implementation/performance experiment, not a new paper accuracy run.

## What is reused

For a fixed source/target, a native growth session records the uncut trajectory
in compressed form. For a new cut:

1. Require the same seed, locked mapping, island partition, previous boundaries,
   tolerance and branch cap.
2. Find the first source-topology read affected by the cut.
3. Reuse the entire growth result if no read changes; otherwise resume the exact
   saved state immediately before the first affected growth iteration.
4. Continue ordinary AAM growth with the cut graph.

Candidates retain correlated symmetry blocks throughout. No bijections are
enumerated, and no growth decision or branch-cap rule is replaced. An exact
per-target symmetry workspace also reuses finalization calculations across cuts.

This implementation reuses **conditional fragment-growth calls**, not a complete
outer AAM search subtree. Different locked mappings remain different queries.
It does not implement the previously discussed family-exclusion search.

Code boundaries:

- `src/rxn_core/cut_replay.py`: opt-in session and cut views; no seed policy.
- `native/src/growth_checkpoint.h`: compressed state and topology dependencies.
- `native/src/growth_replay.h`: exact query cache and checkpoint selection.
- `native/src/engine.cpp`: existing growth function with optional checkpoints.
- `src/rxn_core/search_symmetry.py`: independently reusable `SymmetryWorkspace`.

The native trace-cache residency budget is 64 MiB per worker. Eviction affects
performance only, never search coverage. This is **not** a process-memory limit:
live traces, graph results and the symmetry workspace use additional memory.
Sessions/workspaces belong to one fixed graph pair and one serial worker; they
are not shared concurrently between threads.

## Same-node paired result

Each row compares four directional cases: 25 R-to-P, 76 P-to-R, 77 P-to-R and
114 R-to-P. These include three long tails and the vanadium example with two
previously saved mechanism patterns. Each policy uses ten parallel seed workers,
all sweep cuts, cap 100, matching tolerance 1.0, event tolerance 0.5, bond floor
0.2 and explicit hydrogen. Seed generation is fixed and reproducible.

Each fresh/reference/rolling triplet ran sequentially in the same ten-CPU Slurm
allocation. Method order was counterbalanced. Hosts were bosque7 and bosque11,
both Intel Xeon Gold 6248R. Eight allocations completed successfully; each had a
ten-minute limit and each method a 300-second watchdog. No failed attempts are
included in these timings.

| Seed policy | Reuse policy | Fresh compute CPU seconds | Reused compute CPU seconds | Speedup |
|---|---|---:|---:|---:|
| Same orders across cuts | Uncut-reference replay + symmetry workspace | 1497.92 | 1252.58 | 1.196x |
| Same orders across cuts | Rolling replay + symmetry workspace | 1497.92 | 1373.06 | 1.091x |
| Existing cut-specific orders | Uncut-reference replay + symmetry workspace | 1344.96 | 1397.00 | 0.963x |
| Existing cut-specific orders | Rolling replay + symmetry workspace | 1344.96 | 1479.04 | 0.909x |

CPU seconds sum actual process CPU time across all seed workers. Included:
graph/orbit setup, search and its graph construction/combination, symmetry
finalization, and terminal-witness scoring. Excluded and separately recorded:
audit JSON encoding, archive persistence, queue, imports and input loading.
This is not the entire downstream publication candidate-collection pipeline.

For the best shared-order variant, search CPU changed **835.22 -> 810.17 s**;
symmetry finalization changed **654.87 -> 434.35 s**. Most net savings therefore
came from finalization, not growth replay.

Per-case compute CPU seconds for that comparison:

| Case / direction | Fresh | Uncut-reference replay + shared symmetry |
|---|---:|---:|
| 25 / R-to-P | 604.96 | 480.20 |
| 76 / P-to-R | 554.39 | 485.15 |
| 77 / P-to-R | 335.53 | 285.07 |
| 114 / R-to-P | 3.04 | 2.16 |

Maximum single-worker high-water RSS rose from **1927 to 3008 MiB** for shared
orders and **1576 to 2335 MiB** for cut-specific orders. These are not whole-job
RSS numbers. Slurm's sampled peak job-step RSS across the paired allocations was
16067 MiB. Caching trades memory for reuse and is not ready for default use.

Actual per-method ten-worker job elapsed time, including saved artifacts, is in
each task summary (`job_elapsed_including_io`). It is deliberately not presented
as clean compute latency. Summed per-phase worker wall times are also **not**
parallel job latency. These are single paired observations, not repeated timing
distributions or a whole-holdout speed claim.

## Correctness and alternative coverage

- The paired run compared **1656 directional cut digests**, each aggregating ten
  individually hashed seed graphs: **16,560 exact seed/cut graph checks**, zero
  mismatches. Digests include states, transitions, complete compressed symmetry,
  terminals and cap/stop information, not merely one witness or minimum score.
- The preceding twelve-case, both-direction, ten-seed experiment compared
  **32,320 seed/cut graphs**, also with zero mismatches against saved fresh
  controls. Timing from that experiment is secondary because hosts differed.
- Both historical vanadium event patterns are present among saved terminal
  witnesses in every paired policy: breaking (8,15)/(12,22), or breaking
  (11,12)/(15,27), with their respective two formations. This check used only
  saved outputs; reference patterns were not supplied during search.
- Full test suite: **429 passed** after the final invalid-cut session-state test.
  Differential tests cover tolerances 0.5/1.0, caps 1/3/100, locked mappings,
  island partitions, multiple cuts, changed source symmetry, whole-result hits,
  prefix replay, tiny cache residency and complete final symmetry generators.
- A separate 21-call native-kernel probe against the pre-change compiled engine
  produced identical outputs. Its CPU totals were 1.815 versus 1.753 seconds;
  this small diagnostic does not establish a universal default-path speedup.

**Exactness applies within the same seed policy.** Sharing seed orders across
cuts changes the finite search and is not guaranteed to preserve the old set of
alternatives. In the initial one-seed pilot, shared orders missed a full mapping
for case 123 that cut-specific orders found. At ten seeds, shared orders improved
that case's P-to-R best score (19 -> 7); the other 23 directional best scores
were unchanged. Neither observation proves equivalent alternative coverage
between those seed policies. The native replay experiment itself never chooses
seed orders.

The 140 elementary steps have no independently curated atom-mapping ground
truth. The vanadium patterns are previous saved examples, not proof of chemical
mechanism correctness. Exact replay also does not make capped AAM exhaustive.

## What limits the gain

For shared seed orders, uncut-reference replay reused 848,395 of 4,180,349 logical
extension attempts, and 652,540 of 2,518,502 canonicalization calls. Costs are not
uniform per operation, and dependency checks, checkpoint copies and conversion
still cost time. Once a cut changes the locked fragment history, subsequent
conditional queries no longer match the uncut trace.

With existing cut-specific orders, only 716 of 928,772 growth calls found a
reusable whole/prefix trace. The native cache cannot identify unrelated
conditional histories as equivalent without additional proof. Rolling caches
record more traces but suffered 0.91-1.16 million evictions in the paired run.

**Decision:** retain this isolated experiment and its evidence; do not promote it
as the default or claim it closes the SLAP performance gap. The next genuinely
larger change would need to share conditioned search states/suffixes or reusable
local constraint results beyond identical growth-call histories. Such sharing
must retain full lock/boundary/symmetry constraints. It is not implemented here.

## Saved artifacts and reproduction

The compact summaries, logs, per-task/per-seed statistics, timing provenance and
selected vanadium witnesses are in `evidence.tar.gz`. Aggregate paired metrics
are also directly available in `paired_summary.json`; `mechanism_check_114.json`
contains source-index event patterns and concrete witnesses. Complete compressed
graphs are retained on the cluster for every cut and seed, so figures and
alternative checks do not require rerunning AAM.

Full run folders:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/incremental_cut_replay_20260909_RTPXuu
/project/yunhengzou/coordinate_alignment/aam_benchmarks/cut_replay_seed10_20260909_urKm12
/project/yunhengzou/coordinate_alignment/aam_benchmarks/cut_replay_reference_20260909_vr5WWI
/project/yunhengzou/coordinate_alignment/aam_benchmarks/cut_replay_paired_20260909_aihvmc
```

The first folder holds the pre-change native binary and kernel comparison. The
second captures the earlier rolling-cache implementation. The third/fourth use
the uncut-reference option and reuse the active native source between growth
calls for the same cut. Frozen drivers, engines, manifests and submission commands
record the implementation used for each run. The final source additionally
validates invalid cuts before mutating session state; valid-cut search is
unchanged. Final build/test logs record that post-benchmark correction.

Graphs: `results/<index>/<direction>/<policy>/seed_XX/cut_NNNN.pkl.gz`.
Each seed directory also contains `summary.json`, `progress.json` and
`witnesses.json`. The one-seed initial pilot omits the `seed_XX` layer.

Rebuild the analysis from saved data (no search):

```bash
PYTHONPATH=bench .venv/bin/python bench/cut_replay_pilot.py report --run /project/yunhengzou/coordinate_alignment/aam_benchmarks/cut_replay_paired_20260909_aihvmc
PYTHONPATH=bench .venv/bin/python bench/cut_replay_pilot.py mechanisms --run /project/yunhengzou/coordinate_alignment/aam_benchmarks/cut_replay_paired_20260909_aihvmc --index 114 --mechanism-reference /h/399/yunhengzou/appendix_final/bgcp_rerank_latest/stages/pr7.V.dodh_ts910/rp_stage.json
```

For a fresh experiment, create a **new** run directory, use `prepare` to freeze
the implementation, then `submit` or `submit_paired`. Never resubmit completed
tasks into these saved run folders. `tasks.json` in the paired run contains 48
prepared specifications, but `pairs.json` schedules only the 24 reported tasks.

Full report path:
`/h/399/yunhengzou/coordinate_alignment/reports/incremental_cut_replay_20260909/README.md`
