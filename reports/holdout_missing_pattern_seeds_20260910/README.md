# Can additional seeds recover the three missing AAM patterns?

This follow-up tests **only cases 11, 64 and 101**, the three known alternative
omissions in the forward one-seed/cap-1000 benchmark. It does not rerun the 140-case
benchmark. Each selected case uses 1, 3, 10 and 30 seed orders per cut, with cap
1000, tolerance 1.0, explicit H, the original uncut/single-edge sweep and
**reactant → product only**. The stable frozen engine remains unchanged.

Open [the original-style 3D comparison](viewer.html) for all three cases. Yellow
mechanism buttons select the one-seed AAM alternatives, the missing SLAP witness,
and the actual recovered AAM witnesses for case 101. The R-order/alignment toggle
switches between native P and display-only rigid fitting. Expand the event details
to inspect exact bonds and mapping provenance. The nine saved witnesses were
independently rescored; producing this viewer does not start another search.

## Recovery of the specific missing patterns

| Case | Name | Target event score | 1 seed | 3 seeds | 10 seeds | 30 seeds |
| --- | --- | ---: | --- | --- | --- | --- |
| 11 | pr10.carborane_ts3 | 4 | Absent from saved families | Absent from saved families | Absent from saved families | Absent from saved families |
| 64 | pr16.carbocation_ts5 | 5 | Absent from saved families | Absent from saved families | Absent from saved families | Absent from saved families |
| 101 | pr7.V.dodh_ts19 | 4 | Absent from saved families | Absent from saved families | Recovered | Recovered |

“Recovered” requires a validated mapping with the same joint symmetry-normalized
broken/formed/strengthened/weakened bond pattern. “Absent from saved families”
requires scanning all full terminal representatives and completing the correlated
family membership checks, or a complete necessary-condition exclusion certificate.
It concerns the saved output for that budget, not every
possible seed or every possible mapping. An unresolved query never proves absence.

The target event score is held fixed even if a new search finds a lower score.
`per_case.csv` records each run's best saved representative score separately.
The first successful tested budget does not identify the exact minimum seed count
between tested budgets. This is a selected-case recovery diagnostic, not a new
whole-dataset coverage estimate or a validation of physical reaction pathways.

## Mapping CPU seconds

| Case | 1 seed | 3 seeds | 10 seeds | 30 seeds |
| --- | ---: | ---: | ---: | ---: |
| 11 | 5.434 | 10.411 | 18.871 | 42.933 |
| 64 | 0.768 | 0.834 | 1.689 | 2.679 |
| 101 | 0.855 | 1.005 | 1.623 | 3.809 |

All four fresh budgets for each case run sequentially on the same host with eight
workers. CPU includes the parent and all workers; **do not multiply it by eight**.
The frozen profiler excludes measured checkpoint persistence/loading. Mapping
CPU and observed wall time, which includes archive I/O, are recorded separately.
Analysis is also separate: exact target membership can cost more than mapping.
These are single measurements on three selected cases, not a full timing benchmark.

The 12 fresh mapping runs consumed
90.912 CPU s in total. Completed timed analysis
stages, including the three cached checks and the certificates, consumed
390.300 CPU s. The interrupted general-query retry is
separately recorded in Slurm accounting and is excluded from this completed-stage
subtotal. Its inclusive CPU total is unavailable after cancellation; the recorded
Slurm counters are retained verbatim. The completed-stage subtotal is not the
whole campaign cost.
Searches with cap flags: [].

The initial 30-seed case-11 family check reached its 240-second analysis limit.
Its unresolved result is preserved; `extended_submission.json` records the longer
check on the same archive, with no mapping rerun. That redundant query was stopped
after an independent certificate excluded the target from all 5,637 full-terminal
families. Per-case records retain the initial unresolved status and its analysis cost.

The certificate checks every recorded family action against the product's binary
heavy-atom graph. These actions preserve each terminal's broken/formed heavy-bond
set, and none of the 5,637 sets matches the target under the benchmark's reactant
symmetry. This is sufficient to exclude the full target regardless of hydrogen
or bond-order choices. The certificate was checked on all 12 fresh outputs: it
agrees with the completed queries and does not exclude either known positive
case-101 output. `coarse_validation.json` records 366 literal invariance checks.

## Existing ten-seed, cap-100 outputs

The earlier saved ten-seed run was inspected without rerunning its search:

- Case 11: Absent from saved families; capped=True.
- Case 64: Absent from saved families; capped=False.
- Case 101: Recovered; capped=False.

That historical run changes both the seed count and the cap relative to the
one-seed/cap-1000 benchmark. Its results are contextual evidence; the controlled
seed comparison uses the fresh cap-1000 rows above.

## Verification and reproduction

The three targets are frozen before the diagnostic. The mapping subprocess uses
only endpoint inputs, seed configuration and the original cut policy. Target
patterns are read solely by the post-search analysis. Input hashes, configuration
checks, manifests and commands establish that only these three forward cases ran.
The one-seed control reproduces the three original omissions. Every published
minimum witness and every positive target witness was independently rescored
(18 scalar/canonical checks).

`targets.json` and `per_case.json` include the bond edits, actual mappings and
family-query provenance. Compressed search records, source hashes, frozen analysis
scripts and Slurm accounting accompany this report. Large search archives remain
at `/project/yunhengzou/coordinate_alignment/aam_benchmarks/holdout_missing_pattern_seeds_20260910` and retain their hashes. Existing full benchmark reports are
unchanged. Run `bench/holdout_missing_pattern_seeds.py prepare` with a fresh
`--run`, then `submit`; its task list is restricted to these three cases. After
the case analyses finish, run `bench/certify_heavy_event_exclusion.py --run` on
that directory. This supplies the independent exclusion certificates used here.
Finally publish with `bench/publish_missing_pattern_seeds.py --run` and a fresh
`--destination`. The submission record supplies the frozen Python/source paths.
