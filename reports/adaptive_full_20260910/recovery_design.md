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
