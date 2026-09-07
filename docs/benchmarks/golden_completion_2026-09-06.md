# Golden completion: cap 100 and archive overhead

The initial Golden campaign evaluated all 1,851 records, but 152 searches had
timed out before finishing all configured work. That is not a completed AAM
benchmark. The original evidence is retained; a separate completion campaign
resumes its saved cuts with identical inputs and configuration.

## What was slow

- Branch cap 100 applies separately to each seed/cut search, not to the sum of
  three seeds and every one-edge sweep cut. It does not bound archive size,
  serialization, symmetry finalization, or reference scoring.
- JSON serialization recursively copied shared dataclass data. Loading rebuilt
  millions of repeated atom-pair tuples and repeatedly scanned a growing heap
  with cyclic garbage collection.
- Cached cuts were restored before missing cuts ran. A timeout could therefore
  repeat expensive restoration without making search progress.
- Cached graph restoration and exact symmetry finalization were serial even
  when cut search workers were available.

A profile of case 824 recorded 2.316 million pair-restoration calls and 97.26
million pair-generator iterations. Of the interrupted 267.6-second profile,
saved graph restoration took 82.65 seconds, JSON decoding 48.99 seconds, and
seven new cuts took 47.4 seconds. Native growth accounted for 13.3 seconds.
These are profiling measurements, not a fresh end-to-end speed comparison.

## Changes, without changing matching semantics

1. Borrow graph data for immediate JSON encoding; keep detached records as the
   public default. Pause cyclic GC only around owned archive operations.
2. Intern identical immutable atom-pair and mapping sequences when restoring
   JSON. Logical alternatives, branches, and correlated symmetry groups remain.
3. Run missing cuts before loading saved graphs. Validate input/configuration
   identity before resuming; preserve original checkpoints.
4. Restore and finalize independent cuts in worker processes, then combine in
   the original order. Save each finalized cut atomically for further resumption.
5. Use trusted internal compressed pickle checkpoints to retain shared Python
   graph objects. Keep the producing source snapshot; never load untrusted
   pickle. JSON remains the interchange format.
6. Deduplicate identical heavy-atom reference queries after projecting away
   H-only actions. Preserve correlated group actions and their transition IDs.

No C++ growth, fragment locking, seed selection, tolerance, or branch-limit
algorithm was changed. Tests compare complete serial/parallel graph records,
including symmetry, and checkpoint round-trips and resumption.

The last five cases completed archive construction in 140–184 seconds each
with eight cut/finalization workers per reaction. All their cuts were already
saved; these are completion timings, not fresh search timings. Their recorded
live-branch maxima remained at or below 100. The relevant regression suite
passes 84 tests.

## Evidence and interpretation

Original: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_full_20260906`

Completion: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_completion_20260906`

Each completion engine version is frozen and hashed. Prior logs/status files,
including failed attempts, are preserved under each case's `attempts` directory.
Some original cut paths refer to the saved pilot campaign; retain it too.

Early high-concurrency recovery attempts included process crashes and an
out-of-memory Slurm job. An isolated debugger run did not reproduce a native matcher
crash. Reduced concurrency and removal of periodic stack-dump timers stopped
the observed crashes; this does not establish their precise root cause.

Report completed configured searches separately from globally exhaustive
searches (cap/seed limits still apply), and distinguish scorer `unknown` from
search timeout. Successful resume latency is not cold search latency or total
cost. Parent-process CPU alone excludes cut workers. Archive and scoring time
must be included in workflow costs.

## Final audit

All 152 recovery cases finished, leaving **zero incomplete configured searches
out of 1,851**. Inputs, references, configuration, expected cut counts, archive
presence, and recorded live-branch limits were checked. No previously confirmed
reference was lost; seven additional references were recovered. The complete
reference stratum recovered 1,595/1,760 (90.625%); 23 remain scoring-unknown and
are not counted as successes. Selected-representative top-1 is 1,302/1,760
(73.98%). Partial annotations are reported separately.

See [report](../../reports/golden_aam_completed_20260906/README.md),
[integrity audit](../../reports/golden_aam_completed_20260906/completion_integrity.json),
and [raw completion accounting](../../reports/golden_aam_completed_20260906/completion_slurm_accounting.psv).
Accounting includes failed attempts and diagnostics; the interrupted profiler
job under-reports CPU, so it is not an exact total-cost measurement.
