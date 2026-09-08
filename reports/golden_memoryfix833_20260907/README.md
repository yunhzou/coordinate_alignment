# AAM memory fix: shared records and producer-side checkpoints

No search semantics were changed: same seeds, branch limits, tolerance, fragment domains, mappings, histories and cap diagnostics. Discarded histories remain available; this change reduces their storage cost rather than silently removing them.

## Changes

- DAG-owned deferred-bond records share identical collections and individual edge rows. JSON restoration applies the same sharing. Public record conversion remains detached; no branches or alternative assignments are merged by this storage optimization.
- Cut workers reduce profiling events to the requested maximum counter online instead of retaining every profiling dictionary.
- Workers persist raw cuts themselves and return archive references. Results are consumed as they complete; final graph ordering remains canonical. A slow cut no longer prevents other completed cuts from being saved or causes the parent to buffer their raw graph objects.
- `archive_format='checkpoint'` now selects compact, sharing-preserving raw cuts as well as the final checkpoint. JSON is still the default interchange format. Raw binary and finalized binary records have distinct schemas. Common readers support either explicit raw format and reject duplicate cut IDs.
- Resume, finalization and benchmark consumers use the same raw-cut API. Existing JSON cuts remain reusable; binary deserialization is restricted to trusted internally produced files.

## Controlled measurements

Case833, same uncut input, same C++ binaries, one CPU per probe. Before: `c2e7efe`. Timing is a single repetition; graph equality is exhaustive for each tested returned graph.

| Work | Before RSS | After RSS | Before time | After time |
| --- | ---: | ---: | ---: | ---: |
| cap100, ten seeds | 173.43MiB | 125.41MiB | 4.74s | 4.02s |
| cap2000, first seed | 866.95MiB | 558.32MiB | 26.75s | 21.89s |

The larger probe reduces RSS by35.6%. Its complete canonical graph SHA256 is identical before/after, including112,906 transitions,112,907 states and86,401 stops. All ten cap100 graph hashes also match. Fingerprinting happens after the timed search and RSS measurement. These core probes deliberately retain identical diagnostic profile lists, so they do not include the worker-side profiling-counter saving.

`probe_0/1.json` are before/after cap100; `probe_2/3.json` are before/after cap2000. Probe jobs442578. Full source baseline and outputs reside at:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_memoryfix833_20260907`

## Full-node verification

The first retest (job442582) used shared records and producer-side JSON. It avoided OOM, with Slurm MaxRSS107428452KiB (about102.45GiB), but JSON serialization limited completed work to nine of86 cuts within the five-minute search watchdog. Its verified outcome was unknown, not a miss.

The compact-cut continuation (job442587) reuses those nine cuts. It is a continuation, not a fresh-runtime comparison. The exact frozen source, raw/finalized checkpoints and checks are retained at:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_memoryfix833_checkpoint_20260907`

Final compact-run result: no OOM;80/86 raw cuts saved (nine reused,71 newly saved),46 cuts verified negatively. Search and verification watchdogs both fired, so reference recovery remains **unknown**. Slurm recorded551seconds overall and peak163172012KiB (155.61GiB), reached during parallel verification. During search, observed memory was approximately82GiB near four minutes; that observation is not a separately measured search-stage peak. This demonstrates more retained work without OOM, not that high-cap full-history verification has become lightweight. The `fullnode_json` and `fullnode_compact` folders contain the immutable manifests, outcomes and Slurm accounting.

Both runs retain48 workers, full-node memory, a300-second search watchdog, a240-second verification watchdog, and a ten-minute overall limit. Compared with the original OOM, source history and elapsed work differ; the controlled table above is the equal-work memory comparison.

## Validation and limits

All374 tests pass. New tests cover exact shared-bond records, detached public copies, raw checkpoint transport, both explicit formats, duplicate rejection and mixed-format resume without recomputing saved cuts. Existing serial/parallel, symmetry, native engine and full-graph tests pass.

This is not a claim of constant-memory search: retaining all histories still scales with generated work. The patches do not prune additional paths, loosen chemistry conditions, or claim new benchmark recovery without verification.
