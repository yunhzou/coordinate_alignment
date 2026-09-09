# SLAPMapper: expanded alternative-search budget

**Best combined reference recovery: 1,691 / 1,851 = 91.36%.**
All planned searches have a saved outcome. Returned candidate lists were never
top-K truncated. This experiment does not establish exhaustive search coverage.

## Recovery versus search budget

Each ordering is searched in **both directions**, using the native heavy-atom
symmetry-breaking setting. Both modes still add explicit hydrogens internally.

| Input orderings per direction | Binary mode | Weighted mode | Union of both modes |
|---|---:|---:|---:|
| 1 | 1,620 | 1,576 | 1,661 |
| 2 | 1,639 | 1,583 | 1,675 |
| 5 | 1,658 | 1,593 | 1,689 |
| 10 | 1,664 (89.90%) | 1,593 (86.06%) | 1,689 (91.25%) |

The previous default, single-direction binary run recovered 1,591 (85.95%).
The expanded binary run therefore adds 73 recovered reactions. The union of
binary and weighted results gains no further reference recoveries between five
and ten orderings, although binary alone gains six. This is an observed plateau,
not a proof that a still larger search could never help.

The main experiment comprises **74,040 searches**: 1,851 reactions × 10 input
orderings × 2 directions × 2 bond modes. Ordering zero is the original input;
the other nine use fixed, reference-blind RDKit atom permutations. Equivalent
SMILES can still arise for symmetric/small molecules. Seed: `20260908 + case`.

| Main mode | Returned candidates | Distinct heavy-mapping classes, summed across cases | Timed-out attempts |
|---|---:|---:|---:|
| Binary | 57,495 | 4,572 | 54 / 37,020 |
| Weighted | 48,201 | 3,096 | 46 / 37,020 |

The distinct-class counts use the same chemical symmetry equivalence as the
reference evaluator. They are not counts of raw bijections or mere atom labels.

## Full atom symmetry-breaking check

We additionally ran `break_sym='all'` in both directions and both bond modes on
every Golden case: **7,404 additional searches**, without random-order repeats.
These are separate experiments, not fallback replacements for the heavy mode.

| Mode | Recovered cases | Mapping errors | Timed-out attempts |
|---|---:|---:|---:|
| All-atom binary | 1,537 | 167 | 15 |
| All-atom weighted | 1,491 | 183 | 11 |

The native full-symmetry path can raise errors, including an empty-result
`min()` error on case 590. We did not edit SLAP to suppress those errors.
Despite lower aggregate recovery, these variants add **cases 1304 and 1504** to
the expanded heavy-mode union, producing the final **1,691 / 1,851** result.

## What was—and was not—expanded

- All returned alternatives were retained, with their native cost and compressed
  labeled graphs. No user-facing top-K output cap exists in this SLAP API.
- We expanded input-order and direction diversity, and tested both native bond
  modes and native symmetry-breaking settings. No reference-guided constraints.
- The upstream algorithm itself was unchanged. Its minimum-cost retention,
  symmetry heuristics, and internal LAP-loop guard remain in place. This is not
  an invented near-optimal-cost enumeration mode or a claim to exhaust every
  mapping permitted by the endpoint graphs.
- Scoring uses the same exact heavy-atom mapping relation modulo endpoint
  chemical symmetry as our previous competitor benchmark. Reversed output is
  restored to the original R/P orientation before scoring. No mapping is repaired.

## Resources, failures and reproducibility

Each search request has a hard **300-second process watchdog**. Initial worker
memory was 6 GB. All 82 observed worker-crash attempts were preserved and retried
with **32 GB**; those retries reached the watchdog rather than finishing. There
are no remaining worker-crash or missing statuses in the final evaluations.
Five cases unrecovered by the main two-mode union have at least one timed-out
attempt, so their unsearched possibilities must not be represented as excluded.

Measured worker CPU, including startup, serialization, initial failed attempts
and retries: **16.65 CPU-hours** for the expanded heavy-mode experiment and
**4.09 CPU-hours** for the all-atom check. Evaluation is excluded from these
worker totals. The lower `cpu_hours` fields in per-mode summaries cover only
calls that returned a CPU reading; killed workers have no final process-time
reading. Slurm accounting and request timings preserve their cost separately.

Native source pin: `ea248fd9494f52f4865193e87a98cc92c62b5f9e`.
The same pinned competitor environment and original 1,851-record Golden input
were reused. Main searches used 128 shards per bond mode; the all-atom check
used 32 per mode. Every mapping call is single-threaded. Scoring used 32 shards
per mode, checkpointed to avoid recomputing completed cases. Node-setup-stalled
jobs were relocated; all submissions and retries are recorded. **19 tests pass.**

## Saved results

Full raw mappings, graphs, timings, errors, scripts and logs:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_slap_expanded_20260908`

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_slap_all_atoms_20260908`

This report includes both summaries, the cross-mode analysis, and a compressed
archive of the case-level evaluations and provenance. The much larger raw
candidate banks remain in the cluster directories above. No rerun is needed to
inspect or rescore any returned mapping.

Driver: `bench/golden_slap_budget.py`; mapper adapter:
`bench/golden_competitors.py`. Exact commands, input plans and hashes are stored
with each run.
