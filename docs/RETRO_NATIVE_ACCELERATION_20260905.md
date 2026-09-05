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
under `full_bank/`. I stopped the remaining jobs at 9 min 25 s (replacement jobs
at about 8 minutes), before the ten-minute watchdog. **1,901 of 1,919 molecules
were saved; 18 are unfinished.** Fifteen complete shards took 106–485 seconds
of scanner wall time. The other shards retain their completed records but are
explicitly incomplete. There is no completed full-bank runtime or speedup claim.

`full_bank/progress_audit.json` lists all unfinished source IDs/SMILES and the
slowest completed sources; regenerate it with `bench/catalog_scan_progress.py`.
At the tail, one inspected node had three busy augmented-search workers while
most catalog workers were idle. Budgeted job admission removes the old phase
barrier, but does not redistribute an already-running molecule's work across
newly idle CPUs. A shared fragment-task scheduler is a remaining opportunity;
it has not been substituted or claimed as implemented here.

Verification: 206 full-suite tests passed, followed by four scheduler/archive
audit tests including one newly added test (207 tests total). The methanol
saved-record assembly and standalone viewer also completed successfully.

## Completing the unfinished experiment

At the user's request, job `432033` resumes only the 18 unfinished sources,
one 48-CPU node per source. The detector and configuration are unchanged; the
original row indices are preserved. The 1,901 saved records are not recomputed.
`hpc/resume_catalog_tail.sbatch` runs `bench/resume_catalog_source.py`, which
uses the existing catalog detector and records source-level start/completion
checkpoints, full compressed AAM evidence, and timings under
`data/retro_runs/native_exact_20260905/tail_completion/`.

The allocation watchdog is one hour, with inspection rather than automatic
cancellation at ten minutes. This redistributes CPUs compared with the original
scan. The resumed completion time must therefore be reported separately from a
fresh full-bank benchmark. A regression test verifies configuration/row-index
preservation and refusal to rerun an already saved source.
