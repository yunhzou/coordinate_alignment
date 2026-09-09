# Do the native methods return different mappings?

Comparison of saved AAM and native-XYZ SLAP outputs on all 140 elementary steps.
**No new mapping searches and no ground-truth correctness labels.**

After scoring both on the same original WBO matrices, including explicit H:

| Best saved-result event count | Cases |
|---|---:|
| Equal | 133 |
| Our AAM fewer | 6 |
| SLAP fewer | 1 |

SLAP includes a heavy-atom pattern equivalent to at least one of our minimum-event
saved witnesses on **128/140** cases. On 12 cases, no such shared minimum-witness
pattern was found. This compares against all returned SLAP families, not just its
first candidate. Equal event counts do not imply the same atom mapping.

## The seven differing event counts

| Index | Reaction | Our AAM | SLAP |
|---|---|---:|---:|
| 31 | Pd hydroamination, alkene insertion | 7 | 9 |
| 59 | Carbocation TS11 | 2 | 4 |
| 64 | Carbocation TS5 | 5 | 9 |
| 69 | Carbene insertion TS8 | 5 | 6 |
| 96 | V DODH TS13 | 1 | 3 |
| 104 | V DODH TS32 | 5 | 7 |
| 135 | TS_11 | 7 | 5 |

Fewer events is a graph-based preference, not proof of correct chemistry. In
particular, the TS_11 result does not establish that our compressed archive lacks
SLAP's better-scoring pattern: unexplored symmetry alternatives are not excluded.

## Exact comparison scope

- AAM: score all saved full terminal witnesses from both directions, retaining
  every minimum-score heavy-index pattern for symmetry comparison. Case 123 uses
  the explicitly separate cap-200 follow-up; other cases use cap 100. No new seeds,
  cuts, or search were run. These are minima among saved witnesses, not proven
  minima over the whole compressed AAM output space.
- SLAP: all 176 saved native candidate families. Every heavy-atom label is unique
  within each endpoint. Some hydrogen labels are grouped. Their native costs are
  not used as our WBO event scores.
- Common score: count broken edges, formed edges, and retained edges with WBO
  change greater than 0.5, with bond presence defined by WBO greater than 0.2.
  This reuses the existing scorer's convention. All H atoms participate. No
  separate metal threshold or energy/stereo filter is introduced in this analysis.
- Heavy-pattern equivalence uses exact colored-graph certificates of both original
  endpoints and the mapping relation. Explicit-H endpoint graphs remain present,
  but mapping links are projected to heavy atoms. Bond colors encode their exact
  response to the opposite endpoint's WBO values at the scoring threshold. Thus
  allowed endpoint automorphisms preserve event classification; no arbitrary
  WBO rounding or raw-index equality is substituted for symmetry equivalence.
  This is **score-preserving graph symmetry**, not a stereochemical equivalence
  certificate. No geometry-based chirality information was added.

## Why hydrogen-group handling matters

77 of SLAP's 176 families were not certified score-invariant under their native
hydrogen-label permutations. We optimized those assignments symbolically with
fixed heavy-atom mappings and within-label injectivity, without enumerating
bijections. All family minima were certified: **zero unresolved optimizations**.
Original representatives and optimized witnesses are both saved.

For example, case 90 (`pr5.Noyori_ts65`) initially scored 83 events for an arbitrary
within-label H assignment. Another assignment in the same saved SLAP family has
only 5 events, tying our result. Counting only the arbitrary representative would
have created a misleading apparent difference. Native SLAP search was not rerun
or modified; this is transparent postprocessing of its returned alternatives.

SLAP H-group optimization took approximately 3.21 summed worker seconds. The
initial all-witness comparison took approximately 122.03 summed case seconds,
parallelized by case. Neither figure is new AAM search cost or campaign wall time.
The original exact heavy-index overlap is 107 cases; symmetry-aware comparison
is required and cannot be replaced by that raw-index count.

## Saved evidence and reproduction

- `per_case.csv`: all 140 case names, event totals and broken/formed/order-change
  breakdowns, and shared-pattern indicators.
- `per_case_mappings.tar.gz`: original comparisons and refined comparisons,
  concrete mappings, source terminals/directions, SLAP hydrogen optimality bounds.
- `summary.json`: aggregate results, explicitly null accuracy.
- Frozen comparison/refinement scripts, submissions and Slurm accounting are
  included. Native search inputs/archives remain in the original holdout run.
- Three comparison tests pass: batched event-score agreement with the existing
  explicit-atom scorer, symmetry-equivalent vs distinct mappings, and exact
  hydrogen refinement with an unchanged heavy-atom assignment.

Full analysis directory:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_output_comparison_20260908
```

Use `bench/compare_elementary_outputs.py report` on that directory to regenerate
the summary without rereading the full AAM archives or rerunning mapping.
