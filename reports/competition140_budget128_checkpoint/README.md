# B-priority competition: optimized 140-case experiment

All 140 cases completed with four workers and five-minute per-case watchdogs. The added competition pass improves saved SLAP-sweep minimum-class coverage from **163/169 to 167/169**, and cases with all such classes from **136/140 to 139/140**. Only case 64 still has missing reference alternatives. This measures coverage of SLAP output, not annotated mapping accuracy.

| Reference | Baseline classes | Baseline + competition | Cases with every reference minimum class |
|---|---:|---:|---:|
| native_slap | 151/160 | 155/160 | 136 → 139 / 140 |
| slap_sweep | 163/169 | 167/169 | 136 → 139 / 140 |

## Improvements and controlled pilot

The implementation batches terminal scoring, avoids scoring duplicate mapping vectors, caches dependent-component release calculations, uses a deque, and replaces repeated full progress serialization with compact throttled checkpoints. At 320 repairs, the four-case pilot preserves exactly the old parents, operation/rejection counts, output classes, representative mappings and remaining queue length. Recorded CPU fell from 47.25 to 15.36 seconds (3.08×); batch wall time fell from 33.91 to 9.88 seconds. These are selected-case timings, not a measured whole-dataset speedup.

| Fixed repair budget | Recovered missing sweep classes (four-case pilot) | Recorded CPU seconds | Batch wall seconds |
|---:|---:|---:|---:|
| 64 | 3/6 | 1.50 | 1.76 |
| 128 | 4/6 | 3.07 | 2.99 |
| 320 | 4/6 | 15.36 | 9.88 |

The 128-repair budget was fixed before the full 140-case run. Its choice used the four known gap cases, so the larger run is development-set validation, not an independent test set. No SLAP label or witness is read during proposal search.

## Larger-run cost

Competition: **79.89 s wall**, 252.77 recorded CPU seconds including initial archive loading; peak combined sampled worker RSS **1.84 GiB**.
Full accepted-family decoding: **43.56 s wall**, 70.48 recorded CPU seconds; peak combined sampled worker RSS 0.40 GiB.
Sequential pass wall total: **123.45 s**. Recorded worker CPU total: **5.39 minutes**. All processes passed; no watchdog or memory-limit stops. Longest competition case was case 35 at 30.79 s including startup, closely followed by case 52.

These are additional costs after the existing one-seed, cap-2,000 baseline search. The older successful baseline searches recorded approximately 2.75 CPU minutes, so this extra competition search (4.21 CPU minutes) is substantial relative to matching alone, despite modest wall time. The existing baseline decoder was not rerun. CPU counters exclude interpreter imports before measurement, final output serialization after measurement, parent orchestration, and killed-process costs (none here); development and this separate proof/report audit are additional. RSS is sampled at 0.5 s intervals and may miss short peaks.

## Coverage and complete decoding

| Case | Sweep baseline | Sweep with competition | Native-SLAP baseline | Native-SLAP with competition |
|---:|---:|---:|---:|---:|
| 6 | 3/4 | 4/4 | 2/2 | 2/2 |
| 11 | 1/2 | 2/2 | 0/1 | 1/1 |
| 59 | 1/1 | 1/1 | 0/1 | 1/1 |
| 64 | 2/5 | 3/5 | 0/6 | 1/6 |
| 101 | 1/2 | 2/2 | 1/2 | 2/2 |

Every one of **118,082 accepted full paths** was decoded completely in its baseline event window; 276 solver queries sufficed, with most paths handled by exact invariance or reachable-image lower-bound certificates. There was no class-count cap. All 140 windows completed. No permutation-group or full-bijection enumeration was introduced.

The pass adds **36 event classes across 9 cases** within those windows. Six classes in case 129 require decoding the compressed families and are absent from the repair representatives. These extra classes do not lower any case's previously found minimum event count. All 36 new window classes passed independent exact family-membership queries, anchor checks, and raw event-class recomputation; each newly recovered reference-minimum representative also passed a replay of its local B collision and preservation of unaffected A assignments.

Case 64 still misses two SLAP-sweep minimum classes at four events and five native-SLAP classes at seven events. Those are different reference minima, so their missing counts should not be combined. Full saved-family decoding does not prove that the pending priority search or omitted parent partitions cannot recover them.

## Search limits and kernel finding

The run made 7,401 native proposal calls and 17,737 repair attempts, with cap 2,000 per native call, up to eight baseline parent partitions, and two collision levels. 138 cases retained pending queued states; 13 native repair calls hit their branch caps. Baseline classes are retained by construction, so zero lost baseline classes is a property of union, not a test of replacing the original search.

Independent validation rejected **99 whole repair graphs** (out of 17,737 attempts), containing 1155 invalid paths across 37 cases. They were neither expanded nor credited; otherwise-valid paths in the same rejected graphs were also excluded. Thus accepted-witness claims are sound, but this is still an experimental integration, not a production-ready complete search.

A traced case-64 rejection from the 320-repair pilot reproduces in both Python growth and the native engine. Before extending atom 17, source atom 16 retains target pool {12,14} while anchored source atom 1 maps to target 3. The current representative uses 16→12. Extending 17 admits 16→14, making previously preserved source bond (1,16), WBO 1.10987949312885, map to target bond (3,14), WBO 0. The exact final family is UNSAT. This localizes the problem to retained support in a compressed pool during growth, rather than merely export or the later event thresholds. The safe repair needs to retain the earlier bond-support constraint when refining the pool; changing the event threshold alone would not repair it. No broad kernel change was made on the strength of this single trace.

## Artifacts

- `summary.json`: all per-case coverage, timing, proof checks, configuration, and source hashes.
- `viewer.html`: all 140 cases using the existing white R/P viewer style. New competition classes appear first, with baseline and both SLAP references available.
- `case*/full_decode.json`: complete accepted-path certificates in the same per-case windows as the baseline.
- `case*/takeover*.pkl.gz`: accepted compressed repair archives; `rejected*.json` records excluded graphs.
- `scripts/`: exact experiment, decoder, audit, and diagnostic sources. The production repository was not modified by this experiment.
- `../competition_optimization/invalid_growth_trace.json`: diagnostic Python growth trace; the companion probe script checks the corresponding family with the exact solver.
