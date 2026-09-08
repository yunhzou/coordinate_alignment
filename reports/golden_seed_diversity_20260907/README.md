# Distance-biased AAM seed experiment

## Findings

All 46 paired jobs have finished (Slurm array 442802); no results are pending.

| Unresolved baseline cases (33) | Random | Distance |
|---|---:|---:|
| Verified recovered | 9 | 10 |
| Verified not recovered in searched families | 17 | 17 |
| Unknown / interrupted verification | 7 | 6 |

Both policies recover all 13 controls. Case 1793 remains unknown with incomplete search under both policies; the other unknown cases have saved completed searches but inconclusive verification. Successful Slurm completion means the watchdog workflow finished, not that every AAM/evaluation finished exhaustively.

Adding the two new certified cases (576 and 1033) to the previously accumulated 1827 recoveries gives **1829/1851 = 98.811% archival union coverage**. This is not 100%, not top-one accuracy, and not a fresh full-benchmark ten-seed result. This experiment alone adds eleven distinct recoveries across the two policies to the original 1818-case archive.

Distance bias is complementary, not a strict improvement. Case 576 is recovered only by distance bias (random is a verified miss); case 964 is recovered only by random (distance is a verified miss). Case 1033 is recovered by distance bias while random verification is unknown. All 13 successful controls are recovered under both policies.

Archive inspection locates the actual successful search contexts: distance case 576 uses trial 6, anchor 19, cut (0, 1); distance case 1033 uses trial 2, anchor 23, cut (29, 31). Random case 964 uses trial 10, anchor 5, cut (0, 1). Atom indices are in each planned search's source space. These are ordinary unguided search results, not reference-guided seeds or injected mappings.

On 45 pairs with completed searches, summed per-case search wall time is **1832.93 seconds random versus 1826.69 seconds distance** (-0.34%). The median paired ratio is 1.0008: effectively speed-neutral. These are sums across parallel cases, not campaign elapsed time or CPU-hours. Case 1793 hits the 300-second search watchdog under both policies and is excluded from this completed-pair timing comparison; its partial archives are retained for verification.

Do not promote this policy to the default based on this diagnostic. It produces different useful alternatives at approximately the same cost, but loses a random-policy success. A union of both runs consumes two search budgets and must not be advertised as one ten-seed policy. A single random stream per case also cannot establish average recovery probability across RNG streams.

## Policy

Optional `AAMSearchConfig(seed_selection='distance')`; the default remains `random`.

Each trial is a complete source-atom ordering. Preserve the original random first trial. For later trial anchors, sample without replacement with weight `1 + nearest graph distance` to previously chosen anchors. Disconnected regions receive distance equal to the graph's atom count. This is a soft preference: nearby eligible atoms remain selectable. Retain the existing handling of ambiguous isolated atoms and the original random ordering of all remaining atoms.

Distances are recomputed on each cut graph. Every seed ordering is prepared before matching, without consulting search results, reference atom labels, or coordinates. This adapts to earlier selected anchors, not discovered fragments. Matching, compressed symmetry, branch caps, tolerance and sweep-cut definitions are unchanged. The existing scheduler parallelizes cuts; this change introduces no additional dependency between matching tasks.

Selecting ten orders is O(10 × (vertices + edges)). A local seven-repeat microbenchmark on these 46 uncut source graphs measured median 0.03369 s for random selection versus 0.04295 s for distance selection, about 0.20 ms additional per graph. This is preprocessing only, not an estimate of full-search speedup.

## Experiment

- All 33 unresolved cases in the saved 1,851-record baseline: 32 verified misses and one unknown.
- Thirteen fixed index-spaced previously successful controls; no reference mapping inspected to choose seeds.
- Ten seeds, cap 100, original tolerance 1.0, original smaller-first orientation and full single-edge sweep for both policies.
- Eight workers per policy. Both policies run concurrently on the same 16-CPU Slurm allocation, sharing the same node and filesystem. These are single paired search runs per case, not repeated timing confidence intervals.
- Up to 24 allocations, 96 GiB memory reservation per paired job, search watchdog 300 seconds and scoring watchdog 240 seconds, ten-minute Slurm limit. Incomplete verification is unknown, never a negative.
- Ground-truth recovery uses the existing exact compressed-family verifier, including one-sided reference semantics. This is coverage of the reference mapping, not top-one accuracy or atom coverage.

All initial inputs, frozen source hashes, complete search archives, cuts and verification records are retained under:

```
/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_seed_diversity_20260907
```

`summary.json` and `cases.json` contain the comparison; separate policy manifests fingerprint the executed source. This is an unresolved-case diagnostic plus a small control set, not a full 1,851-record new-policy evaluation. Combining independently recovered cases is an archival union, not one ten-seed run.

## Software checks

Native-enabled full suite: 378 tests passed in 90.58 seconds. Subsequently added checkpoint-policy identity coverage; all five focused distance-policy tests passed. Tests cover reproducibility, unchanged first trial, stable seed-count prefixes, complete orderings, anchor uniqueness, disconnected graphs, isolated atoms, more trials than atoms, invalid policies, and statistical separation across 200 RNG streams. Default random checkpoint identity remains compatible; distance checkpoints cannot be silently reused as random-policy results.

Final combined seed-policy/checkpoint regression check: 17 passed in 1.43 seconds. `git diff --check` passes.
