# Reference-blind direct-anchor pilot

## Outcome

Direct anchors are useful as a complementary search action, but these policies do not replace cuts or reach 100%.

The change-oriented policy recovered **case 1739** with one atom anchor and no cut. The unanchored run, random-anchor policy, stable-environment policy, and the sixteen tested single cuts did not recover it. The anchor was selected before reading reference labels: proposal 7, planned-source N13 → planned-target N16. Search took 0.108 seconds, returned 140 terminals, and was uncapped. The reference is the top representative of this anchored task, certified by the existing symmetry-aware evaluator.

This source nitrogen has one triple-bonded nitrogen neighbor; the selected target nitrogen has single bonds to N and P. Favoring a different local environment supplies a chemically changed atom assignment instead of demanding preservation of its original local motif. The core matcher then determines the remaining assignments. This is a graph heuristic, not a validated chemical reaction rule.

The other unresolved cases in this pilot—986, 1285 and 1553—were not recovered by any tested policy. These outcomes concern the returned families under finite budgets; they do not prove no other anchor can work.

## Equal task-budget comparison

Every policy gets sixteen independent hypotheses per case. Each hypothesis uses ten seeds, branch cap 100, tolerance 1, explicit H and one CPU. An anchor hypothesis contains exactly one same-element atom pair and no cuts. A cut hypothesis contains exactly one source edge cut and no anchors. This is not a full sweep per anchor.

| Case | Unanchored, no cut | Random anchors | Similar-environment anchors | Different-environment anchors | Random single cuts |
|---|---:|---:|---:|---:|---:|
| 986 | No | 0/16 | 0/16 | 0/16 | 0/16 |
| 1285 | No | 0/16 | 0/16 | 0/16 | 0/16 |
| 1553 | No | 0/16 | 0/16 | 0/16 | 0/16 |
| 1739 | No | 0/16 | 0/16 | **1/16** | 0/16 |
| 300, successful control | Yes | 11/16 | 9/16 | 1/16 | 16/16 |
| 600, successful full-sweep control | No | 0/16 | 0/16 | 0/16 | **1/16** |

Counts are successful hypotheses, not probabilities estimated over repeated RNG streams. Case 600 is a control previously recovered with full sweep, not an uncut-success control. Wrong anchors can exclude valid mappings, as case 300 shows. Keep independent baseline and cut results; never replace the whole search with one chosen hard anchor.

## Selection rules

`bench/golden_anchor_policy.py` is diagnostic-only. Its `proposals(problem, budget)` function reads only endpoint elements and WBO matrices—not references, saved mappings or evaluation outcomes.

- Random: sample same-element pairs without replacement.
- Similar environment: prioritize rare element-pair classes, then descending multiset-Jaccard similarity of neighboring element/bond-order descriptors.
- Different environment: same rarity ordering, but ascending neighborhood similarity.
- Random cut control: sample sixteen source edges without replacement.

All explicit H pairs remain eligible. Rarity is the product of the source and target population of that element. Ties use a seeded shuffle. The finite budget and rarity ordering are explicitly experimental hypotheses, not necessary chemical filters. The experiment does not enumerate full atom bijections or expand automorphism groups. Shared proposals between policies reuse the same saved calculation.

## Execution and verification

305 distinct tasks cover 390 policy/baseline slots after sharing duplicate proposals. All 305 searches and evaluations completed with no unknowns or watchdog terminations. 120 tasks recorded branch-cap hits; absence of the reference in those returned families is not an exhaustive negative. Every returned representative was checked to preserve its requested anchor.

Summed measured core search time: **80.234 seconds**; median per task **0.166 seconds**, maximum **3.360 seconds**. These are single-run search-stage wall timings, excluding process startup, scheduling, verification and archive finalization. They are not complete campaign elapsed time and are not isolated-node speed comparisons.

Each task had a 90-second search watchdog, 45-second score watchdog and five-minute Slurm limit. Arrays: 442866, 442931, 443171. Eleven unstarted tasks stalled in node configuration and were canceled/resubmitted without rerunning completed tasks. Raw accounting is retained.

Two policy tests pass: deterministic/reference-metadata-independent selection, element compatibility, uniqueness, budget behavior and explicit-H eligibility. Runtime anchor-preservation assertions pass for all completed searches. No production files under `src/rxn_core` changed; production baseline remains `29b0722` for this experiment. The experiment harness was checkpointed separately as `44c5135`.

## Artifacts and next decision

`manifest.json` records selections, task sharing and frozen source hashes. `cases.json` preserves every outcome and witness, and `summary.json` groups them by policy. Full AAM graphs, logs, status and input-source paths are saved under:

```
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_anchor_policy_20260908
```

The positive case is task `174`; its `aam.pkl.gz` and `evaluation.json` contain the complete result and verification. Selection was reference-blind, but the pilot cases were chosen from known failures and controls, so this is not a new full-dataset accuracy estimate.

Recommendation: retain direct anchors as one optional candidate-generation action alongside unanchored/cut runs. Before promoting a selector, test broader cases and multiple RNG streams at fixed total cost. Coordinated anchors remain untested; their potential should not be confused with the demonstrated single-anchor result.
