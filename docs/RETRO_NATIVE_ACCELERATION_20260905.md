# Exact bookkeeping acceleration — 5 September 2026

The baseline here is commit `608af96`, not the historical 6.77-second workflow.
Explicit H, tolerance 0.5, branch cap 100, orbit-representative seed policy,
every initial family, full augmented AAM, no candidate cap and no sweep are
unchanged. The native AAM growth/search engine itself is unchanged.

## Implemented modules

1. **Symmetry finalization:** one prepared target topology, sequential recoloring,
   and exact ordered-coloring caching within each graph. Atom-generator projection
   moves to C++, preserving generator order, atom labels and full frames.
2. **Hierarchy representation:** `AAMHierarchyChain` composes independently
   transformed immutable segments. It does not transform the locked prefix with
   the residual action. Archives preserve segments and shared generators without
   eager conjugation. Explicit `.fragments` / `.to_record()` still materialize the
   same typed/legacy evidence. Conjugation runs in C++ when actually requested.
3. **Occupation processing:** native arrays replace Python nested sorting,
   permutation application and relation-key lookup during stage closure. The key
   still includes target atoms, attachments, labeled fragment sets and preserved
   bonds. Chronological stages, discovery order and first witnesses are unchanged.
   This is not a coverage-only quotient or enumeration of all atom bijections.
4. **Scheduling:** large jobs start early; CPU-budgeted admission immediately
   reuses freed slots. There is no ordinary/outlier phase-wide barrier.
5. **Persistence:** workers encode/compress complete records. The coordinator
   copies gzip members instead of serially unpickling and recompressing archives.

No seed removal, placement sampling, chemical rules, new search caps or fallback
paths were added. C++ entry points validate array frames/indices for memory safety.
Python remains the public abstraction. Install with `pip install -e .`; builds
require C++17 and use pybind11 through normal build isolation.

## Same-host molecular checks

Bulky BIAN target, one CPU per detector; fresh runs, with full evidence saved.
These are individual measurements, not statistically stable throughput estimates.

| Source | Before | After native + sharing + coloring | Speedup | Candidates |
| --- | ---: | ---: | ---: | ---: |
| Iodobenzene, INVENTORY-001161 | 1.100 s | 0.852 s | 1.29× | 404 |
| Aniline, INVENTORY-001301 | 4.691 s | 3.498 s | 1.34× | 896 |
| Long tail, INVENTORY-000435 | 56.892 s | 30.962 s | 1.84× | 1,303 |

All three expanded-evidence SHA-256 hashes match their baseline: candidate
mappings, hierarchies, symmetry, complete search graphs, derivations and cap
diagnostics. Aniline and the long-tail result remain explicitly capped. The
iodobenzene result remains marked rough under the existing seed policy.

The long-tail native-only intermediate took 38.611 s; adding shared hierarchy
chains took 30.323 s. Coloring reuse gave no additional measurable improvement
on that individual pair. Do not multiply these measurements into a larger claim.

Reproduction and saved files:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python \
  bench/native_retro_acceleration.py --output data/retro_runs/NEW_RUN
```

`data/retro_runs/native_exact_20260905/` holds `before`, `native_first`,
`shared_chain`, `reused_coloring`, full gzip records, metrics, and the final
profile. The benchmark watchdog bounds each fresh source run to 180 seconds.
`methanol/methanol.html` is a regenerated end-to-end assembly/viewer smoke test.

## Remaining limits

The final long-tail profile spends 14.62 of 41.88 instrumented detection seconds
in augmented `find_islands`, 11.62 in exact symmetry finalization, and 1.67 in
the native occupation closure. These are a completed-pair profile, not the old
first-20-second sample. Profile overhead is not a production runtime estimate.

Cross-family continuation memoization is **not implemented**: matching atom
coverage is not a sufficient state key. Anchors, cuts, seed position, branch-cap
accounting and the future conditioned search must agree before states can share
work. No unsafe merge was substituted to claim a speedup.

Full-bank scan job `432000`, with replacement `432029` for two nodes stuck in
configuration, uses the same 1,919-entry inventory and 28 × 48 CPUs. Results are
under `full_bank/`. Its completion and timing must be measured separately; pair
speedups do not establish full-bank or assembly throughput.
