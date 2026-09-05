# AAM RAM regression: stage-isolated diagnosis

## Finding

The main blow-up begins in **shared AAM symmetry finalization**, then grows
further when the fragment consumer materializes occupations and collects
parallel results. It is not explained by core growth exceeding its branch cap.
The in-memory symmetry representation repeats huge numbers of identical dense
generator arrays; the sharing applied later by archive serialization is too late.

No production algorithm was changed during this investigation.

## Baseline and identical inputs

- Old **core AAM**: `4a7b57f`, before the search-graph refactor. No old retro code
  was used as a baseline.
- Current core AAM: `9e2221a`.
- Same native growth binary; `native/src` has no changes between those commits.
- Source `INVENTORY-000400`, one of the actual 100-GB OOM failures; 76 explicit
  atoms, bulky BIAN target. All 59 actual augmented inputs were captured directly
  at the current consumer's `find_islands` call boundary.
- Both versions received identical serialized graphs, WBO matrices, anchors,
  seed order and settings. Cap 100, tolerance 0.5, explicit H, no sweep.
- Native graph caches were omitted from the input archive and rebuilt equally
  by both cores; no molecular input or search setting was changed.

All **59** comparisons produced identical terminal mappings, matched-fragment
partitions and deferred cuts, after normalizing arbitrary island labels.
This checks terminal search relations, not a claim that the two versions have
identical history schemas or identical generator encodings.

| Measurement across 59 inputs | Old core | Current core |
| --- | ---: | ---: |
| Process RSS at core return | 55.3–82.9 MB | 55.3–198.1 MB |
| Slowest individual core call | 4.30 s | 7.40 s |
| Maximum terminal branches | 100 | 100 |

There is a smaller core representation/time regression: the current graph keeps
the explored history and capped subtrees, not just surviving branch records.
That alone did not create gigabytes per worker in this test.

## Where RAM grows

For captured input **family 2**, the current saved AAM graph was replayed without
running AAM again. Measurements are process RSS, not recursive object-size
estimates:

| Boundary | RSS |
| --- | ---: |
| Loaded current core result | 154 MB |
| After `finalize_graph_symmetry` | 842 MB |
| After occupation generation and projection | 1,509 MB |
| Peak while pickling the result | 2,308 MB |

The finalizer took 26.49 seconds; the complete post-core consumer took 58.97
seconds. It generated 26,487 pre-merge candidates from 100 returned paths.
This experiment used one worker, so the increase is not solely an aggregate-RSS
artifact from shared pages across processes.

### AAM result-construction defect

`src/rxn_core/search_symmetry.py:17` processes every recorded transition,
including the history of discarded/capped subtrees. Line 41 creates fresh
`list(g)` arrays for every transition's generators.

The saved graph audit found:

- 7,740 transitions were finalized, but only **148** reach returned terminals.
- **218,432** stored generator occurrences contain only **49 distinct values**.
- Those occurrences store **29,706,752 atom-index slots**; the unique generators
  need only **6,664 slots** plus references/group membership.
- All 7,740 finalization requests missed its conditioned-state/coloring cache.

These are generator arrays, not enumeration of the whole permutation group.
Nevertheless, copying the same dense arrays this many times defeats the memory
benefit of the compressed representation. This function is also called by
`search_aam`, so it is not merely a retro-only serializer problem.

### Additional consumer/execution amplification

`fragment_matching/augmentation.py` walks returned paths and eagerly retains
every materialized occupation. `fragment_matching/parallel.py:140` collects all
family results into a tuple before the final merge. Parallel transport and
parent-side collection add more resident data. Increasing per-source concurrency
to up to 47 workers made this representation problem more dangerous; four real
jobs exceeded their 100-GB allocation.

The measured single-family sizes explain why high concurrency can cross that
limit. They do not prove every worker reached the same peak simultaneously,
nor identify the exact allocation that triggered each historical OOM.

## Fix direction — not implemented by this diagnosis

1. Keep exact generators interned/shared **in memory**, not only in saved JSON.
2. Avoid eagerly finalizing symmetry for histories with no returned outcome;
   retain their cap/provenance records and compute groups only when requested.
3. Preserve compressed correlated occupation families through detection; avoid
   eagerly expanding and collecting all family payloads before consumption.

No new chemical filter, branch pruning, witness sampling, or weakening of
hydrogen/symmetry semantics is needed to address these representation defects.

## Reproduction and evidence

`bench/aam_memory_boundary.py` supplies prepare/core/post/audit/compare modes.
The post mode reuses a saved core result; it does not repeat AAM. The Slurm
launcher `hpc/diagnose_aam_memory.sbatch` has a ten-minute hard allocation limit;
individual core diagnostic subprocesses are bounded to five minutes.

Evidence root:
`/h/399/yunhengzou/coordinate_alignment/data/retro_runs/aam_memory_boundary_20260905/`

- `inputs/inputs.json`: exact input hashes and family metadata.
- `inputs/*.old.metrics.json`, `*.current.metrics.json`: isolated core results.
- `inputs/core_comparison.json`: 59 matching terminal-relation comparisons.
- `inputs/2.graph.pkl`: saved current core result.
- `post_2.out`: memory at symmetry/occupation/result boundaries.
- `inputs/2.post.pkl`: saved post-core result for further analysis.
- `inputs/2.storage_audit.json`: generator duplication and reachable-history counts.
- `old_core/`: detached old-core worktree using the identical native binary.

Successful diagnostic jobs: `432051` (core pairs), `432052` (post-core replay),
`432055` (old terminal-state recording), `432056` (saved storage audit).
