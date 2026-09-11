# Binary metal connections for search; raw-WBO events downstream

Binary metal search recovers the previously missing PR7 four-event alternative with one seed. The lowest saved event count is unchanged in all 140 holdout cases, and Golden reference recovery is unchanged. The observed mapping CPU reduction is small (2.6%). This supports a targeted benefit for metal-sensitive matching, without establishing general chemical accuracy or a consistent speed advantage.

The frozen engine `98b01b1` searches original weights or weights in which every metal-containing pair is 1 when raw WBO >= 0.2 and 0 otherwise. Nonmetal weights remain unchanged. One seed, tolerance 1.0, explicit H, and the original uncut/single-edge sweep are retained. The holdout uses R→P and cap 1000; Golden uses both directions and cap 100.

**Both arms use the original raw WBOs for event scoring:** decrease >= 0.5 means breaking/weakening and increase >= 0.5 means forming/strengthening; pairs involving a metal use 0.3. This is the existing `classify_bonds` rule, with inclusive thresholds and no extra 0.2 event-presence gate. The older holdout benchmark helper used a different uniform-0.5/floor-0.2 rule, so its event totals cannot be directly compared to these.

| 140-case forward holdout | Original WBO search | Binary metal search |
|---|---:|---:|
| Full saved mappings | 140/140 | 140/140 |
| Mapping CPU, seconds | 941.484 | 917.013 |
| Median mapping CPU, seconds/case | 2.229 | 2.170 |
| Cases with cap flags | 13 | 13 |
| Unique full representatives, summed | 79601 | 88364 |

CPU ratio binary/original: 0.9740. Fresh pairs use the same CPU allocation and alternate order. CPU includes parent plus workers and excludes measured archive persistence/loading. Conversion and analysis are separate. Single paired measurements do not establish a statistical speed guarantee.

Best representative event comparison (binary versus original): {'equal': 140}. After targeted saved-family checks, best verified scores: {'equal': 140}. These are achievable scores, not proofs of the global minimum. The holdout has no annotated chemical ground truth.

## Previously missing alternatives

| Index | Target events at 0.5/0.3 | Original membership | Binary membership |
|---|---:|---|---|
| 101 | 4 | excluded_from_saved_families | represented |
| 11 | 4 | excluded_from_saved_families | excluded_from_saved_families |
| 64 | 4 | excluded_from_saved_families | excluded_from_saved_families |

PR7 gained witness: two breaking/weakening and two forming/strengthening events. The table uses zero-based input atom indices; all weights are the original WBOs.

| Event | R pair | P pair | Raw R WBO | Raw P WBO | ΔWBO |
|---|---|---|---:|---:|---:|
| broken | 0–2 | 8–9 | 2.191932 | 0.695110 | -1.496822 |
| broken | 19–20 | 10–15 | 0.993625 | 0.000000 | -0.993625 |
| formed | 0–20 | 8–15 | 0.000000 | 1.216286 | +1.216286 |
| formed | 2–19 | 9–10 | 0.000000 | 1.053606 | +1.053606 |

Targets are the historical missing witness mappings, reclassified under the requested delta thresholds. Membership includes compressed families, with search constraints checked on search weights and events evaluated on raw WBO. Unresolved checks are not counted as exclusions. Pattern identity uses signed event pairs modulo exact raw-WBO score-response symmetry. This can distinguish patterns merged by an unweighted connectivity convention.

Cases with differing minimum representative pattern sets: [{'index': 101, 'shared': 1, 'original_only': 0, 'binary_only': 1}]

Targeted family comparison: [{'index': 101, 'shared': 1, 'binary_gain': 1}, {'index': 11, 'shared': 1, 'neither': 1}, {'index': 64, 'shared': 2, 'neither': 1}]

The query set includes both arms’ minimum representative patterns wherever their sets differ, plus the three historical targets. Identical representative sets do not prove equality of the entire compressed-family event spectrum; unseen patterns are not enumerated.

## Golden

Normalization changes only 11 Golden inputs. Those are rerun in both directions; the other 1840 search inputs are exactly unchanged and retain their archived reference checks. This is not a fresh timing run of all 1851 cases.

Reference recovery including unchanged archives: {'original': {'recovered': 1833, 'not_recovered': 17, 'unknown': 1}, 'metal_binary': {'recovered': 1833, 'not_recovered': 17, 'unknown': 1}}

| Index | Archived original | Fresh original | Binary |
|---|---|---|---|
| 691 | recovered | recovered | recovered |
| 692 | recovered | recovered | recovered |
| 904 | recovered | recovered | recovered |
| 937 | recovered | recovered | recovered |
| 959 | recovered | recovered | recovered |
| 1093 | recovered | recovered | recovered |
| 1312 | recovered | recovered | recovered |
| 1418 | recovered | recovered | recovered |
| 1447 | recovered | recovered | recovered |
| 1474 | recovered | recovered | recovered |
| 1475 | not_recovered | not_recovered | not_recovered |

## Validation and reproduction

Small exhaustive tests cover pair-specific inclusive thresholds, input immutability, symmetry-equivalent events, and a compressed family whose raw-WBO event scores differ despite identical binary metal weights. Every emitted event witness is independently rescored with `classify_bonds`. All saved states are checked for element preservation and injectivity; preserved transition bonds and stored-generator target adjacency are checked. Frozen source and original input hashes are verified. Engine defaults and earlier reports remain unchanged.

All 55 holdout inputs unaffected by normalization reproduce the same checked search counts and minimum representative patterns. The initial auxiliary audit failed on 48 analyses because discarded branches legitimately lack finalized symmetry groups. Only those analyses were retried after correcting the ancestry check; every returned path still requires finalized groups. No mapping searches were rerun. Original failure statuses, repair hashes and scheduler accounting are retained.

Input conversion: 0.386 CPU seconds for 140 cases (separate post-run measurement). Successful holdout event analysis CPU: original 577.155 s; binary 614.502 s. These analysis totals exclude archive loading, auxiliary structural audits, failed attempts and optional family queries.

Raw outputs and logs: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/metal_binary_20260911`. Drivers: `bench/metal_binary_benchmark.py`, `bench/metal_binary_events.py`, `bench/publish_metal_binary.py`. The run manifest records source hashes, configurations and scoring semantics.
