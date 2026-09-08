# Compressed-path event scoring experiment

Outcome: promising per-path scoring speed, not a full Golden reranking or a universal latency guarantee. No AAM searches were rerun. Production scoring does not enumerate complete bijections or group elements.

## Method

1. Keep the stored path, coupled group factors and native domain constraints.
2. Check whether all recorded generating actions preserve the labels needed by the event scorer. If yes, certify the representative's score for the whole path family.
3. Otherwise compile the existing full-explicit-atom feasibility constraints shared with ground-truth verification.
4. Express bond-event cost using local bond tables; remove constant terms. Minimize with bounded SMT queries, retaining an achievable mapping and a lower bound.
5. Report `optimal` only when the bounds coincide. Other results remain bounded and unproven.

Matching tolerance is 1.0; scoring counts broken, formed and bond-order-changed events at tolerance 0.5, including explicit H. Both original R-to-P orientation and reversed search on incomplete compositions are tested. The native AAM search algorithm is unchanged. `family_query.py` only extracts its existing constraint builder for reuse; mechanism grouping is not involved.

## Targeted checks

Five Golden reactions, six directional archives. The selected paths were already known to contain ground-truth-compatible mappings, so this first experiment is explicitly ground-truth-selected. Optimization itself receives no reference mapping. Three repeats per path, all with proven minima; all 18 resulting mappings independently match the heavy-atom ground-truth certificate modulo exact endpoint symmetry.

| Case | Search direction | Representative events | Minimum events | Scoring time range |
|---|---|---:|---:|---:|
| 0 | P to R | 11 | 11 | 8.85–9.97 ms |
| 44 | R to P | 6 | 4 | 187–216 ms |
| 44 | P to R | 6 | 4 | 178–201 ms |
| 865 | P to R | 13 | 11 | 90.8–116 ms |
| 1033 | P to R | 8 | 6 | 176–228 ms |
| 1786 | P to R | 5 | 3 | 234–305 ms |

These results verify that tolerance-based matching equivalence can contain distinct 0.5 event scores, and that finding a minimum does not require exhaustive mapping enumeration.

## Reference-blind stress sample

From each of the same six archives: first eight reference-blind ranked classes plus 24 random terminals (fixed random seed 42), one recorded path per terminal. These are 192 path families, **not 192 independent reactions**.

| Subset | Proven minima | Median scoring time | Maximum scoring time |
|---|---:|---:|---:|
| Top-ranked classes | 48/48 | 0.122 s | 0.397 s |
| Random terminals | 135/144 | 0.184 s | 6.388 s |
| Combined | 183/192 | 0.179 s | 6.388 s |

27 paths had an invariance certificate; 40 required no solver query after invariance/constant-term analysis. Nine remain bounded, not proven optimal. The three-second setting is a soft encoding/solver budget, not a hard API wall-clock guarantee. Completed solver calls and witness decoding/validation can overrun it. Each case also had an external 240-second watchdog; none hit it.

The 6.388-second outlier did not reproduce in a separate instrumented run: it took 3.001 s, including 0.840 s constraint encoding, 0.068 s objective encoding, 1.947 s solver checks and 0.022 s witness decoding/validation. Its bounds remained 40–41. Therefore the original outlier cannot be attributed specifically to decoding from these measurements.

## Negative optimization result

Arc-consistency/forced-image support pruning was tested on exactly the same 192 paths. It produced 182 proven minima and a 0.184 s median, versus 183 and 0.179 s without it. All paired bounds were consistent. No benefit was demonstrated, so this extra pruning was removed from the retained prototype. The archived experiment remains available.

## Timing and scope

Each batch used six concurrent single-process cases on six allocated CPUs of bosque79 (64 GiB allocation). Initial targeted batch: 42 s elapsed; blind batch: 76 s; pruning experiment: 80 s. These batch times include loading/startup. Archives took 1.5–14.5 s to load per case; peak process RSS of roughly 0.5–3.7 GiB includes the complete pre-existing AAM archive.

The per-path timers include validation, invariance testing, constraint/objective construction, solving and witness checks. They exclude archive loading, path indexing and initial library import. Summed per-path wall times are not parallel batch latency. This experiment establishes neither a full-reaction/global minimum over all paths nor a new 99% bond-window metric. A full ranking would still need to account for all competitive families and unresolved searches. No speedup multiplier against a like-for-like previous optimizer has been measured.

28 relevant verification, scoring and benchmark tests pass. Tiny exhaustive enumeration is used only as a unit-test oracle, including coupled explicit-H symmetry and incomplete compositions in both orientations.

## Artifacts

`summary.json` and `per_family.csv` hold measurements and proof bounds. Source snapshots, settings, concrete witnesses and logs:

- `/project/yunhengzou/coordinate_alignment/aam_benchmarks/family_event_scoring_20260908` — targeted, job 452124
- `/project/yunhengzou/coordinate_alignment/aam_benchmarks/family_event_scoring_blind_20260908` — blind baseline, job 452125
- `/project/yunhengzou/coordinate_alignment/aam_benchmarks/family_event_scoring_pruned_20260908` — rejected pruning, job 452126
- `/project/yunhengzou/coordinate_alignment/aam_benchmarks/family_event_scoring_profile_20260908` — outlier profile, job 452127
