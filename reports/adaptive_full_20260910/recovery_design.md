# Recovery investigation — before changing matching/search code

Baseline and first candidate are frozen at the full-benchmark artifact root.
Do not overwrite those outputs, redefine equivalence, feed reference mappings
to blind searches, or append original answers to a candidate to claim recovery.

## Confirmed structural gap

Golden 1740, P→R: 12→60 explicit atoms. Adaptive exhausted after 27 growth
calls, with two full terminals, no pending work and no cap. The original's
reference-equivalent witness uses cut (1, 2).

Golden 1321, P→R: 15→16 explicit atoms. Adaptive exhausted after 55 growth
calls, with 40 full terminals, no pending work and no cap. The original's
reference-equivalent witness uses cut (1, 8).

These demonstrate missing cut-conditioned search states, not just an
insufficient numerical seed/work budget. Seed-only exhaustive continuation
from an uncut graph is not equivalent to allowing changed fragment boundaries.
The successful contexts are diagnostic evidence only, not production inputs.

## Revision to test

Restore the general source-edge cut conditions to the adaptive search domain.
Share prepared target data, conditioned symmetry and dependency-validated
fragment results across those conditions. Interleave the conditional agendas
so one uncut subtree cannot consume the entire budget before another boundary
condition is considered. Keep every correlated placement from each growth.

This restores missing types of search choices; it does not by itself prove
that a finite agenda budget equals the original ten-seed coverage. Validate
against **all Golden and all 140 holdout records**. Report unknown family
verification, lower-budget snapshots, capped calls and unvisited agenda work.
Retain original comparison metrics and explicit-H event scoring.

No changes to the native fragment-growth kernel are needed for this revision.

## Complete cut-interleaved run and next correction

Candidate `386f688` completed all 3,982 directional calls and all 496 tests.
Holdout saved best-event scores have no regressions. Golden bidirectional
reference recovery is 1,825 versus 1,836 original: eight previously recovered
cases are now negative and four are unresolved. Equal single-direction totals
hide different recovered cases and are not equivalence.

Saved-path diagnosis (no new mapping calculations) confirms:

- Case 21, P-to-R: the successful cut is visited, but the original root seed
  16 has not been explored before work exhaustion.
- Case 1246, R-to-P: the original 29-atom prefix is present, but its subsequent
  singleton decision is not completed before the work budget.
- Case 1801, P-to-R: only one root seed is explored in the successful cut
  before the global CPU budget. The original successful seed is unvisited.

Additionally, the adaptive scheduler uses one guide after every fragment,
whereas the original keeps ten distinct shuffled continuation orders. The
finite adaptive agenda has no guarantee of covering those original paths.

Next experiment: preserve the original cut and ten-order search domain and
its synchronized frontier cap, but memoize exact conditional decisions and
intern their fragment DAG online across seed policies. A policy carries only
its live frontier/cursor; repeated `(state, seed)` decisions share computation,
and equal fragment transitions share representation. Full compressed matches
are retained. No ground-truth-derived guides and no original-result fallback.
This changes scheduling/storage, not growth, matching tolerance, or scoring.
Measure the cost honestly; restoring equivalence takes precedence over calling
a smaller, incomplete search an acceleration. Use a fresh full benchmark.

## Shared-policy full-run resource finding

The exact shared-policy differential tests pass (weighted graphs, unequal
composition, both cap stages and complete correlated fragment records).
Previously missing cases 21 and 1246 now recover after policy exhaustion.
The first full run nevertheless exposed an implementation/resource failure:
round-robin retained every cut's branch dictionaries and decision caches at
once. Slurm job 465803 reached 65,068,852 KiB, with individual workers above
6 GiB; workers were killed in its 64-GiB allocation. Further expensive calls
hit the 300-second one-CPU watchdog. These are not successful comparisons.

Move the same shared policy scheduler into the established per-cut worker
pipeline. Each worker finishes, checkpoints and releases a cut; multiple cuts
run concurrently with an explicit CPU budget. Reuse existing producer-side
checkpoints, independent finalization, timing and bounded worker orchestration.
Do not alter fragment growth or retry via another mapper. Preserve this failed
attempt and its counters, and run the complete benchmark again under the
corrected, frozen execution backend.
