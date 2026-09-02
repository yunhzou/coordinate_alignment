# AAM Acceleration Analysis (100-1000x target)

Static analysis of `rxn_core` on branch `organic_single_step_retro_syntehsis`
(commit 23b1e8b), 2026-09-01. No code was executed; every claim below cites
the source it was read from. Numbers marked *measured* come from
ALGORITHM.md section 18 and committed run artifacts; everything else is an
estimate from operation counts and should be confirmed by the instrumentation
in section 4 before large engineering commitments.


> **Status (branch `Fable_AAM_Opt`):** the Layer 0 engine now runs as a
> compiled extension and the Layer 1 tails were trimmed; measured,
> output-identical results are in `LAYER0_CORE_AAM_ACCELERATION.md`, section 6
> (12 workers: TEMPO 72x, tetraphenylmethane 124x, tetra-tert-butylmethane
> 123x versus the original serial time; serial 20x, 35x, 27x).  The tiers below are the original
> analysis and ranking.

## 1. Where the time actually goes

### 1.1 Correct the baseline first

`benchmarks/aam_regression_contract.json` gives 380/650 s (Fe TS8) and
30/900 s (Fe TS2) for AAM/post-AAM. Those are **ceilings** set from the
pre-July code (commit ecca6a0), not measurements. Measured values
(ALGORITHM.md:927-931):

| Case | Atoms | Workers | AAM search | post-AAM | Total |
|---|---:|---:|---:|---:|---:|
| Fe TS8 reductive elimination | 82 | 16 | ~299 s | 13.7 s | 313 s (was 823 s) |
| Fe TS2' redPath | 83 | 16 | ~22.9 s | 0.98 s (was 792 s) | 23.9 s |
| Pd TS12 (post-AAM only, old pipeline path) | 133 | 8 | n/a | 11.5 s | n/a |
| PR9 carbene ts41a | 40 | ? | | | 3.47 s |
| TEMPO ts3 R/P stage | 57 | ? | | | 15.6 s |
| 140-case manifest (post-AAM reprocessing) | mixed | ? | | 9,165 s summed, ~66 s/case | |

Two conclusions follow. AAM search is ~95% of wall time on the medium cases
today. And per-unit cost is governed by symmetry-driven branching, not atom
count: TS8 and TS2 have the same size and the same caps but differ 13x.

Catalog fragment detection (this branch's main workload):

| Run | Target atoms | Rows | Workers | Wall | Notes |
|---|---:|---:|---:|---:|---|
| Bulky BIAN inventory, rough mode | 78 | 1919 | 28 shards x 64 | 14.8 / 42.6 / 204.7 s (min/median/max per shard) | ~1 row per worker per shard, so shard wall = slowest row |
| Vancomycin gap, exact mode, region-directed | 132-176 | 412 | 48 | 958 s | ~100 CPU-s per row; median 76 live branches per row |
| Large-star stress, single process | 183 | 2 | 1 | 14.6 s (12-atom source), 133 s (79-atom source) | no cap hit: pure constant-factor cost |

### 1.2 Anatomy of one AAM search work unit

Work list = (no-cut + one item per R edge with WBO >= 0.2) x 3 seed orders,
dispatched with ordered `imap`, chunksize 1 (sweep.py:1699-1726). For an
82-atom complex that is ~270 units; TS8 spends ~18 CPU-s per unit, TS2 ~1.4.

Per unit, in cost order (all pure Python except the nauty C call):

1. **Extension scan** (extend.py:426-459). For every live candidate C and every
   popped R edge, the code loops over *all* product atoms
   (`for v in ctx.g_P.nodes()`), rebuilding `_sym_block_indexes` per v
   (extend.py:413, support.py:133), running the support check with three
   Python calls per WBO read (primitives.py:7-29), and computing a dense
   relation signature of length |fragment| per v (dedupe.py:72-87).
   Estimated 40-75% of per-step time. The active-edge rule
   (support.py:157: `w_P >= graph_floor`) means only P-neighbours of the
   already-mapped images can ever pass, so the admissible set is ~1-4 atoms,
   not N.
2. **Per-child pynauty certificate** (dedupe.py:185, canonical.py:284-293).
   Every child of every extension gets a fresh `pynauty.Graph` over the whole
   product graph (N+E ~ 170 vertices at N=82, ~380 at N=183) plus two passes
   of `_colored_vertices` and a repr-sorted cell list, then a saturation
   re-dedupe (island.py:373). It runs even when there is exactly one child.
   Python marshalling is most of the ~0.5-2 ms per call; nauty itself is
   ~20-100 us on these near-discrete partitions. Estimated 20-40% of
   per-step time.
3. **Per-branch scoring** (sweep.py:1013-1066 -> branch.py:354-571).
   `classify_bonds` is an O(N^2) Python double loop run 3-4 times per
   completed branch; `symmetry_repair_mapping` rebuilds the full colored P
   graph (O(N^2) `_wbo_tolerance_bucket_lookup`, orbits.py:82-124) and runs
   `autgrp` once per touched orbit group per branch, then up to 20,000
   vectorised evaluations. Identical complete mappings recur across cuts and
   seeds (TEMPO: 166 hits of one mechanism) and are repaired independently.
   Estimated 5-15% of unit CPU, up to ~1 s per branch when the cap binds.
4. **Per-unit setup** (sweep.py:1176-1189): event canonicaliser, `build_graph`,
   `_nauty_orbits` on the cut graph, seed orders, all rebuilt per (cut, seed)
   although two are cut-invariant and one is fully constant. ~25 ms per unit,
   <1% on TS8, ~2% on TS2. `r_orbits_cut` is a provable no-op: it only shifts
   boundary-signature key values inside one dedupe call.
5. **Parent-side merge** (sweep.py:1527-1598): every result is pickled to
   bucket files, a second pool reduces them (parallelism = number of distinct
   mechanisms, typically 1-3), then reloaded; hierarchy dicts are deep-copied
   2-4 times on the way. Instrumented (`parent_merge_seconds`) but never
   persisted; estimated 5-20 s serial on TS8.
6. **Fragment-group finalisation** (aam.py:21-42 -> sweep.py:934-983): serial
   in the parent, over *all* mechanisms including non-minimum ones, one
   `autgrp` per cache miss. Its output is not read by the R/P or TS path
   (rp.py:64 uses the compiled family's generators).

### 1.3 Anatomy of post-AAM

- `compile_mapping_families` runs `_masked_relation_data` twice per unique
  payload and 12 nauty runs per family (index_chirality.py:897-956,
  1201-1225). With the default empty anchors, graph A is byte-identical to
  graph B, so 9 of the 12 runs are redundant. The pool is capped at 8 workers
  unless there are >= 128 payloads and is skipped below 16 (analytical.py:115-122).
- `select_rp_mappings` is **fully serial** (rp.py:34-81). The 48-worker
  (branch, coset) scheduler described in ALGORITHM.md:818-823 lived in the
  deleted `pipeline.py`; the typed API has no `workers` parameter here.
  `--post-workers` reaches only family compilation.
- Direct chirality route (114/166 mechanisms): after per-factor frame
  filtering it enumerates the Cartesian product of surviving constrained
  factors (index_chirality.py:1733-1751) and runs a full covariance
  branch-and-bound per valid base (1815-1819), rebuilding identical factor
  matrices and ball trees each time. Chirality-neutral rotors (CH3, CF3, tBu)
  multiply the product by 3 each: TS04 has 2187 bases.
- Compiled route (52/166): 5 + N_c individualised `autgrp` runs plus 4 nauty
  runs per constraint unit on graphs that grow by 120 vertices per accepted
  simplex; N_c ~ N/2. Consistent with ~60 s on 133-atom TS11.
- Per (family, coset), `analyze_group_chirality_branch` rebuilds both endpoint
  graphs with O(N^2) `build_graph` (1316-1317) although `static_context`
  already holds them; sympy Schreier-Sims runs even when |G| == |K|.

### 1.4 Anatomy of one catalog row

`detect_fragments` (fragment_matching/detection.py) grows an island from
*every* explicit atom including H (mode `all`) against the full target, pays a
joint source+target certificate per placement (deep-copying a ~430-600 vertex
base each time, canonical.py:94-101), then runs competitive augmentation which
recomputes target orbits for a fresh available-target subgraph per residual
component (augmentation.py:85-100 -> branch.py:616-624). Two full
`_PartialMappingCanonicalizer` bases are rebuilt per row (O(N_target^2)
Python) although the target never changes. Catalog bond orders are formal
integers (rdkit_adapter.py:23), so all target-side buckets are trivially
precomputable. The driver has no necessary-condition prefilter (every row is
searched) and no per-row timing.

### 1.5 Documentation versus code

- docs/ARCHITECTURE.md:32-36 says cut-sweep work is grouped per cut with one
  orbit computation per cut; sweep.py:1699-1703 dispatches (cut, seed) pairs
  and recomputes everything per pair.
- ALGORITHM.md:818-823 says (branch, coset) units run on up to 48 workers;
  rp.py is serial.
- The contract's post-AAM ceilings are 47-900x above measured values.

## 2. Ranked proposals

Exactness classes: **E0** bit-identical outputs (pure implementation change);
**E1** identical selected mappings and mechanism certificates, but internal
metrics or provenance labels may change; **E2** changes the search space or
the computed object and needs owner sign-off plus 140-case re-baselining
against certificate digests.

### Tier A: exact, small, high yield (1-5 days each)

| ID | Change | Where | Est. gain | Class |
|---|---|---|---|---|
| A1 | Skip the certificate when there is one child; bucket children by the color profile and by the sorted multiset of (target orbit, role) before calling nauty; compute `_colored_vertices` once; one canonicaliser per `grow_island` | dedupe.py:177-209, canonical.py:253-293 | 3-10x fewer nauty calls; 1.4-2x per step | E0 |
| A2 | Replace the full-target scan with the exact admissible set (intersection of P-neighbourhoods of the fixed bonded images, union over open-block pools, same-element index); hoist per-candidate structures out of the per-v loop; use the compact relation signature; vectorise the fixed-image predicate with numpy | extend.py:426-459, support.py:97-206 | 10-100x fewer inner evaluations; 1.6-3x per step on AAM, 3-5x per row on 183-atom targets | E0 |
| A3 | Vectorised `classify_bonds`; reuse base events in repair; memoise `_nauty_atom_generators` by target set and repair by mapping; cache the bucket lookup per graph; build the event canonicaliser once per worker; drop `r_orbits_cut` | frag.py:199-262, branch.py:354-571, sweep.py:1176-1189 | 10-50x on scoring kernels; 1.2-3x on branch-heavy units | E0 |
| A4 | Stream-merge results in the parent (no bucket files, no second pool, no deep copies); carry a precomputed branch key; restrict fragment-group finalisation to minimum-event mechanisms and run it in the worker pool | sweep.py:1527-1598, 700-737, aam.py:21-42 | removes est. 5-20 s serial tail on TS8; needed before anything else matters at the 3-10 s scale | E0 |
| A5 | Parallelise `select_rp_mappings` over (mechanism, family, coset) with the existing total-order reduce; compute `static_context` once; pass it to `analyze_group_chirality_branch`; skip sympy when structural and conservative orders are equal; 12 -> 3 nauty runs per compile and one relation build; lift the compile-pool caps; membership shortcut for payloads with equal fragment sets | rp.py:19-93, analytical.py:99-195, index_chirality.py:897-1142 | post-AAM 13.7 s -> ~1-3 s on TS8; 60-90 s -> ~10-20 s on 133-atom compiled-route cases | E0 |
| A6 | Fold chirality-filtered constrained factors into the single covariance branch-and-bound (kill product enumeration and per-base searches); pass generators, not closures; batch ordinary units into one trial and use certificate-only existence tests on the compiled route; pair-orbit BFS instead of N_c individualised `autgrp` | index_chirality.py:1696-1820, 2024-2308 | up to V_b-fold (2187x on TS04-like) on direct-route selection; 3-10x on compiled-route chirality | E0 |
| A7 | Catalog: move the target half of the joint canonicaliser into `FragmentTargetContext`; thread target orbits into augmentation; pre-bucket joint certificates by (source orbit, target orbit) multisets; per-row dispatch with largest-first ordering; exact prefilters when `minimum_fragment_size > 1`; canonical-structure memo | detection.py, augmentation.py, tools/search_mcule_retro.py | per-row 2-5x (with A1+A2 10-30x); 3-10x fewer node-hours from scheduling alone | E0 (prefilter/memo: E1 diagnostics) |

Expected after Tier A alone: TS8 313 s -> ~40-60 s at 16 workers (5-8x),
TS2 23.9 s -> ~5-7 s, catalog exact rows ~100 CPU-s -> ~5-10 CPU-s.

### Tier B: exact, larger engineering (2-8 weeks)

| ID | Change | Est. gain | Class |
|---|---|---|---|
| B1 | Native growth kernel (Rust/PyO3 or C++/pybind11): compact integer `_SymCand`, bitset adjacency, identical predicates and cap semantics, direct nauty C API with one persistent sparsegraph per target and recolouring per call | 40-100x on unit growth; combined with A3/A4 unit CPU ~18 s -> 0.25-0.5 s on TS8 | E0 (requires a record/replay differential corpus) |
| B2 | Fork-join inside a unit: dispatch the live branches of each seed step to a flat pool, then replay `_admit_subtree` in branch order; per-branch scoring as tasks | breaks the slowest-unit wall bound; makes 150-500 cores useful per reaction (today ~270 units, tail-bound at 48+) | E0 |
| B3 | Exact branch-and-bound on a robust lower bound of broken+formed events, seeded by the no-cut units (prune live branches/candidates whose bound exceeds the no-cut minimum) | 3-6x on branch-dominated cases like TS8; ~1x on TS2 | E1 (superset under the cap; min-event output unchanged) |
| B4 | Deferred-atom short-circuit; skip redundant saturation dedupe; memoise `grow_island` per (branch state, seed) to remove the guaranteed no-progress final pass; incremental frontier bookkeeping | 1.1-1.5x | E0 |
| B5 | Quotient completed branches into families inside the sweep workers by relation transport, so the parent, pickles, finalisation and compile see U' families instead of B branches | 5-50x on everything downstream of growth on symmetric cases | E1 (pool schema change) |

### Tier C: change what is computed (owner decisions)

| ID | Change | Est. gain | Risk |
|---|---|---|---|
| C1 | Assignment (1-hop environment + Hungarian) lower bound on event count; when it equals the no-cut minimum, run one core-relaxed multi-cut search instead of the full sweep, else fall back | 10-45x fewer work units; the only route to ~1000x on search | family antichain may differ within the same certificate; TS2/TS13/TS41a contract entries pin atom labels, not digests |
| C2 | Compile one coarsest-partition family per mechanism instead of B per-branch hierarchies (the compiled family depends only on mapping, fragment partition, anchors) | post-AAM compile B -> 1 | RMSD can only decrease; degeneracy-group labels may change |
| C3 | Proposer -> exact completion fast path (VF2/RI/FMCS proposes, LAP bound proves optimality, existing compile/chirality/RMSD verifies) | search -> ms when one proven-optimal mechanism suffices | loses the complete mechanism set unless guarded |
| C4 | `heavy_cuts_only`, adaptive cut ordering with early stop, fewer seed orders | 2-30x fewer units | can miss mechanisms only reachable through a skipped unit (pr14 bH-elimination is H-transfer chemistry) |
| C5 | Catalog: rough MCS-proposer screen with exact finalist refinement; region-anchored seeding for gap searches; one seed per source orbit | 100-1000x on the screening pass | changes the discovered family set; RETROSYNTHESIS_DESIGN.md sections 2 and 10 already intend a two-tier design |

## 3. Composite arithmetic (Amdahl, honest)

**Single medium reaction (Fe TS8, 313 s at 16 workers).**
Tier A: search 299 -> ~40-50 s, post-AAM 13.7 -> ~2 s, merge ~1 s: **~45-55 s (6x)**.
Plus B1 (native kernel) and B3/B4: unit CPU ~0.3-0.5 s x 270 / 16 = ~6-9 s
search, ~2 s post-AAM: **~8-12 s (25-40x)**. Plus B2 on 200-500 cores:
search ~2-3 s: **~4-5 s (60-80x)**. Reaching **100x (3.1 s)** additionally
requires C1 or C2-scale reductions in what is computed, or accepting the
E1 superset semantics of B3/B5. **1000x (0.3 s) is not reachable** while
computing the same object: a 48-process pool spawn plus per-worker graph and
orbit initialisation, one family compile (3 nauty runs on ~350 vertices), the
Schreier-Sims quotient and one covariance search already sum to ~1-3 s.

**TS2-like cases (23.9 s)** are fixed-cost bound after Tier A: ~2-3 s (8-12x).

**Large reactions (133 atoms)** gain more from A2/A3/A5/A6 because the O(N^2)
Python loops scale worse; no measured AAM time exists for 133-atom cases in
the current code, so instrument first.

**Catalog scan.** Exact mode against a 130-180-atom target: ~100 CPU-s/row
-> ~5-10 (Tier A) -> ~1-3 (B1) per row, i.e. 30-100x CPU. Node-hours drop a
further 3-10x from scheduling fixes (today 1792 cores served 1919 rows with
~1 row per worker; shard wall equals the slowest row). Small-target scans
with `minimum_fragment_size >= 6` gain 2-4x more from exact prefilters.
1000x in exact mode is not reachable because the answer itself is large
(median 76 live partial embeddings per row on vancomycin); the design
document's two-tier plan (rough screen, exact finalists) is the legitimate
route.

**Resource scaling.** Today: the sweep saturates at ~3N work units and is
straggler-bound earlier; the post-AAM stage is serial, so adding cores to the
83-atom case buys nothing past 16. After A5 and B2: 150-500 cores useful per
reaction. Catalog rows are embarrassingly parallel already; the constraint is
CPU-hours per row, not parallelism. GPU is unsuitable (backtracking over
pointer-heavy state, serial nauty, sub-launch-latency kernels); threads are
GIL-bound until B1.

## 4. First week: measure, then quick wins

1. Persist the timing the code already computes but drops: `seed_orders`,
   `growth_calls`, `parent_route/reduce/load/stream_seconds`
   (sweep.py:1484-1501 vs domain.py:147-166); return the `cut_end` per-unit
   split (search/expand/score) even without `trace_path`; add counters
   `children_before_dedupe`, `certificates_computed`, `certificate_collisions`,
   `scan_visits` to the growth profile dict (island.py:65-86); add per-row
   `elapsed_seconds` and `getrusage` deltas to catalog records.
2. Run TS8, TS2 and one 133-atom case at 16 and 48 workers with tracing; run
   one catalog shard. This decides between branch-level (B3) and per-step
   (A1/A2) priorities and sizes the merge tail (A4).
3. Add `mechanism_event_certificate_digests` for TS2, TS13 and TS41a to the
   regression contract so later family-level changes can be validated by
   certificate rather than by atom labels.
4. Land A1 (singleton early-out is a two-line change), A3's memoisation of
   repair generators, A5's `select_rp_mappings` pool and `static_context`
   hoisting, and A7's per-row dispatch. All are E0 and each is under a week.
5. Build the record/replay differential harness (inputs and outputs of
   `_extend_sym_cands`, `_dedup_sym_cands`, `grow_island`,
   `symmetry_repair_mapping` on the contract cases and a few catalog rows).
   It is the gate for A2, B1 and B5.

## 5. What not to do

- Do not merge candidates, placements or live branches by orbit labels,
  hashes or event summaries; pre-grouping is allowed only as a
  necessary-condition split before the exact certificate (ALGORITHM.md
  sections 4.3, 7.2, 16). The earlier automorphic live-state quotient lost 2
  of 4 TS04 mechanisms.
- Do not change child ordering, heap tie-breaks, which edges are recorded as
  deferred, or the 4096-state support cap when porting to native code; they
  determine retained witnesses, branch keys and pool representatives.
- Do not switch the sweep to `imap_unordered` or larger chunks: pool insertion
  order selects the retained representative (tests pin this).
- Do not touch `_MechanismEventCanonicalizer`'s graph or cell order: its
  certificate hex is persisted in the regression contract.
- Do not treat `heavy_cuts_only`, seed skipping or cut-orbit dedupe as
  optimisations; they are search-space changes.
- Do not grow heavy atoms first and attach H afterwards in exact mode: the
  greedy interleaving decides which partial embeddings survive.
- Do not use the contract ceilings as measurements when prioritising.
