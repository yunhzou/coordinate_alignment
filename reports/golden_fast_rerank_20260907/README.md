# Fast bond-cost and compactness reranking: exploratory ablation

**Neither simple heuristic improved top-1. The production ranker is unchanged.** All candidates and original archives remain available; nothing was filtered out of the saved AAM result.

## Results

All-record percentages use all 1,851 Golden cases, including the 19 incomplete original searches with unknown rankings. Compactness was evaluated on all 1,832 completed archives.

| Ranking rule | Correct / 1,851 | Accuracy | Gains / losses versus baseline |
|---|---:|---:|---:|
| Existing event-count ranking | 1,435 | **77.53%** | — |
| Compactness as tie-break | 1,432 | 77.36% | 8 / 11 |
| Compactness soft penalty, weight 0.5 | 1,432 | 77.36% | 8 / 11 |
| Compactness soft penalty, weight 1 | 1,432 | 77.36% | 8 / 11 |
| Compactness soft penalty, weight 2 | 1,429 | 77.20% | 8 / 14 |

Energy-based rules were evaluated on **1,586 identical supported cases**. For 246 completed cases, at least one candidate in the best coverage tier needed a bond type/order absent from the energy table. Those cases were explicitly marked unsupported for energy scoring, not assigned guessed values or silently switched to another ranker.

| Energy-supported subset only | Correct / 1,586 | Accuracy | Gains / losses |
|---|---:|---:|---:|
| Existing baseline on exactly these cases | 1,262 | **79.57%** | — |
| Energy-weighted edit cost as tie-break | 1,260 | 79.45% | 33 / 35 |
| Energy-weighted edit cost as primary score | 1,260 | 79.45% | 33 / 35 |
| Compactness + energy tie-break | 1,258 | 79.32% | 35 / 39 |
| Combined soft score | 1,258 | 79.32% | 35 / 39 |

Energy weighting changed the selected class in 89 cases as a tie-break (97 as a primary score). Compactness changed it in 24 cases. Thus these are not merely inactive heuristics: they fix some cases and damage slightly more. This experiment does not establish that every possible chemically informed score will fail.

## Definitions and scope

- Every method preserves the existing priority: maximize mapped heavy atoms, then total mapped atoms including H.
- Candidates are the **same 91,571 cached class representatives** from the previous top-k study. This isolates reranking; it does not optimize a new witness within each compressed class or enumerate additional alternatives.
- Event definitions exactly match the existing ranker: broken, formed and bond-order-changed bonds, including H. An assertion checks this for every candidate.
- Compactness counts connected components of the edited-atom induced subgraph in the union of R and mapped P bonds. Unmatched P atoms get distinct nodes. The penalty is the number of components beyond one; edited-atom count is a secondary tie-break. This is a local connectivity heuristic, not a reaction-mechanism rule.
- Energy cost is the sum of absolute changes in generic bond energies for the edited bonds, divided by the C–C single-bond value (346 kJ/mol) to give a convenient scale. Broken and formed bonds both contribute positively. Aromatic order 1.5 uses the midpoint of listed single/double values. No unknown bond type is imputed.
- Values come from the [CDK 2.8 BondEnergies documentation](https://cdk.github.io/cdk/2.8/docs/api/org/openscience/cdk/smsd/tools/BondEnergies.html). This is **not an activation energy, a computed molecular energy or a quantitative barrier estimate**. Signed net bond-energy changes alone can be mapping-independent for complete endpoint assignments, hence the unsigned edit-cost experiment.
- The combined soft score is event count + extra center components + 0.25 × normalized energy edit cost. The nine fixed formulas are listed in the script. They were tested as an exploratory Golden ablation, not held-out validation or an optimized new model.

No rule uses the reference to choose candidates. Reference equivalence is read only to evaluate the selected class. Original candidate-set coverage is unchanged because no candidates were removed. The claim is not that the top-five shortlist is unchanged under a different ranking.

## Speed and saved data

No AAM reruns, quantum calculations, conformer generation, chirality filtering or compressed-family queries were performed.

- 32 one-CPU workers; all exited zero.
- Longest allocated task: **90 seconds**, including Python startup, archive loading and file writing.
- Summed archive-loading time: **1,209.77 seconds** across workers.
- Summed feature extraction time: **508.67 seconds**, about **5.6 ms per candidate** or **0.28 seconds per completed reaction**. These are summed measured elapsed times, not parallel wall time.
- The features are now cached, so further score-only experiments do not require loading AAM archives again.
- **151 regression tests pass**, including explicit-H event accounting, partial mappings, unsupported energies and reference-independent ranking keys.

Full frozen source snapshot, per-case feature tables, logs and statuses:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_fast_rerank_20260907`

Slurm array: **442443**. One CPU and 16 GiB requested per worker, five-minute process watchdog and six-minute Slurm limit. Source: `bench/golden_fast_rerank.py`.

Local `features.json.gz` stores all candidate features; `cases.json` stores selected terminals, scores and correctness under each rule. `summary.json` includes every gained/lost case index. `manifest.json` stores executed source hashes, and `slurm_accounting.psv` stores allocation accounting.
