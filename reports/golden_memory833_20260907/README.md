# Case833 memory diagnosis — no core changes

The evidence identifies substantial history-storage overhead amplified by parallel workers. It does not show that the compressed automorphism representation itself exploded.

## Controlled core comparison

Case833 has 82 explicit atoms on each side. Replay the same uncut input and first ten seed orders, tolerance1.0, in one process with the same C++ binaries. Old full AAM: `4a7b57f`; current: `0fa0e17`. These are core comparisons, not old retro comparisons. One repetition per setting; timing is descriptive, not a precise speed benchmark.

| Identical completed work | Old core RSS | Current core RSS | Old/current time |
| --- | ---: | ---: | ---: |
| cap100, all ten seeds | 101MiB | 173MiB | 4.96 / 4.77s |
| cap2000, first seed | 458MiB | 867MiB | 28.68 / 26.82s |
| cap2000, first five seeds | 1041MiB | 1963MiB | 86.13 / 75.76s |

First-seed growth-profile counts match (41,210), as do cumulative counts after five seeds (128,251). This diagnostic does not claim independent verification of equality of every output mapping. The data-structure versions deliberately retain different history detail.

The 90-second cap2000 probes stopped at different work positions: current finished six seeds, old five. Their final RSS values must not be treated as equal-work comparisons; the table uses matching completed seed prefixes instead. Probe instrumentation retains builder references for inspection, adding list/reference overhead; it does not duplicate their state or fragment payloads.

## Findings

1. **Dead/capped histories remain resident.** At cap100, all ten seeds retain 18,827 transitions, of which only1,380 reach returned terminals (92.7% do not). At cap2000, the first completed seed retains112,906 transitions but only6,771 reach returned terminals (94.0% do not). It returns2,000 terminals and retains84,401 cap stops. The live frontier cap is not a bound on accumulated history size.
2. **Deferred-bond payloads are repeatedly copied.** A separate traced cap2000 partial first-seed run records362.6MB live Python allocations. `alignment/branch.py:140`, which constructs a fresh list of lists from `iso.deferred_edges` for each transition, accounts for166.6MB (46.0%). Mapping tuple snapshots at line134 add35.2MB. This is storage overhead, not evidence of group-element enumeration. Tracemalloc substantially increases time and RSS; its process RSS is not compared with untraced runs.
3. **Parallelism multiplies those working sets.** The failed production allocation had48 cut workers, each accumulating ten seed graphs. Slurm reported an OOM at a184000MiB node allocation. The probes demonstrate roughly2GB per worker after five uncut seeds; individual cuts differ, so this is not an exact reconstruction of the production180GB peak. The old core also grows substantially at cap2000.
4. **A further orchestration risk exists.** `_search_cut` retains all seed graphs and then combines them. The parent uses ordered `Pool.imap`: a slow earlier cut can delay consuming/checkpointing completed later results, while their large returned graphs remain buffered. This buffering contribution was identified in code, not quantified in the killed run.
5. **Global semantic deduplication is not the first established culprit.** The first current cap2000 seed has zero join transitions. The dominant measured costs already exist inside one seed: accumulated discarded histories and repeated payloads. Cross-seed/cut reuse may help, but merging merely similar states would require proof that their future continuation domains agree.

## Principled next changes — proposed, not implemented

- Keep full histories needed by surviving solutions; separate discarded-search trace storage from resident solution state. Preserve cap-hit evidence in compact records or an optional streamed trace. Merely deleting all history would violate the intended result/animation model.
- Share immutable deferred-edge collections and other repeated payloads. Avoid copying whole edge lists into every transition. Validate exact fragment domains, mappings and cap semantics before/after.
- Persist completed seed/cut results promptly and consume cut results as they complete; transport archive references rather than large graphs where possible.
- Size concurrent workers against measured memory, not just CPU count. Full-node CPU saturation is unsafe when the per-worker footprint exceeds available memory divided by concurrency.

No production search logic was changed during this investigation.

## Reproducibility

- Jobs442571 (four bounded old/current probes) and442575 (allocation tracing).
- Each probe: one CPU,16GiB Slurm allocation,90-second diagnostic search limit (60seconds traced), three-minute overall limit. Signal handling can be delayed by native calls; Slurm provides the outer bound.
- Inputs, old source snapshot and logs: `/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_memory833_20260907`.
- Input SHA256: `324d8debe660b7d5d14036083769cfe8b9b3854a352064080b2f29b66b90b836`.
- Old/current use the same copied native binaries, as in the previous full-AAM core regression protocol.
- `probe_0/1`: current cap100/2000; `probe_2/3`: old cap100/2000; `trace`: separately instrumented current cap2000. Partial-graph terminal-ancestor counts must not be used to declare live unfinished paths dead.
