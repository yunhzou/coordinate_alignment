# Layer 0: Core AAM Engine Acceleration

Scope note (2026-09-02). Work proceeds in modules, innermost first. Only exact
changes are considered: for identical inputs the engine must return an
identical list of completed branches (mappings, island partitions, deferred
edges, fragment hierarchies, in the same order). Search-space changes (fewer
cuts, fewer seed orders, pruning by event bounds) are out of scope.

## 1. Module map

| Layer | Package / files | Role | Status |
|---|---|---|---|
| 0 | `rxn_core.matcher` (`state`, `extend`, `support`, `dedupe`, `canonical`, `orbits`, `primitives`), `rxn_core.growth.island`, `alignment/branch.py::find_islands` | the matching engine: symmetry-compressed candidate growth of one fragment, multi-island branch construction | **now** |
| 1 | `alignment/sweep.py` | cut enumeration, per-unit setup, in-sweep branch scoring and symmetry repair, worker pool, merge, fragment-group finalisation | after Layer 0 |
| 2 | `analytical.py`, `rp.py`, `alignment/index_chirality.py` | families, chirality, RMSD | separate track; RMSD may be dropped |
| 3 | `fragment_matching`, `retrosynthesis`, `tools/search_mcule_retro.py` | catalog detection and assembly | later; reuses Layer 0 directly |

Layer 0 contract: `find_islands(g_R, g_P, seed_order, ...)` and
`grow_island(g_R, g_P, seed, mapping, ...)` are deterministic functions of
their arguments. Acceptance for every change below is a bit-identical replay
on recorded inputs, plus the existing test suite and unchanged mechanism
certificate digests in `benchmarks/aam_regression_contract.json`.

## 2. What one extension step costs today

### 2.1 Measured Python primitive costs (this machine, CPython 3.9, shapes for an 82-atom molecule with a 60-atom fragment)

| Primitive mirrored | Where in code | Cost |
|---|---|---:|
| WBO read via dict of dicts, 3 nested calls | `primitives._edge_wbo` shape | 0.09 us (the real path indexes a numpy scalar and calls `float()`; expect roughly 0.3 to 0.6 us) |
| support predicate per bonded pair | `support._pair_ok` | 0.21 us (dict form) |
| rebuild block index maps, 3 blocks | `state._sym_block_indexes` | 0.6 us |
| dense relation signature, F = 60 | `dedupe._p_relation_signature_from_parts` dense branch | 15 us |
| `_SymCand.__init__` shape, F = 60 | `state.py:49-85` | 4.5 us |
| `_colored_vertices`, N_P = 82, repr-sorted cells | `canonical.py:253-269` | 52 us (done twice per certificate today) |
| pynauty `Graph` Python-side construction and validation | `canonical.py:276-282` | 27 us |
| branch state key, three sorted 60-item dicts | `branch._progress_key` | 3.7 us |

Not measurable here: nauty's C time (literature and experience put it at 20
to 100 us for a 170-vertex near-discrete partition) and pynauty's dict-to-C
conversion (tens of microseconds). A candidate certificate is therefore about
0.2 to 0.4 ms at N = 82 and roughly 2x that at N = 183, of which more than half
is Python.

### 2.2 Structural counts on the committed TEMPO example (57 atoms, 58 active edges, 30 of them X-H)

| Quantity | Value |
|---|---|
| same-element atoms scanned per extension today | C: 19, H: 30, O: 6 |
| admissible atoms under the neighbourhood restriction (product neighbours of the mapped image with the right element) | C: mean 1.5 max 3; H: mean 2.4 max 3; O: mean 2.7 max 4 |
| scan reduction over all extensions | **11x** (2260 evaluations down to 206) |
| extension steps for one full-molecule growth | 58 (one pop per active edge) |
| cut work units, seed-order tasks | 59, 177 |

### 2.3 Cost model per extension step

Let C be the live candidate count, S the scanned same-element atoms per
candidate (today ~20 to 40 on organics), F the fragment size, Ch the number
of children after extension.

| Term | Cost today | Notes |
|---|---|---|
| scan: C x S x (block index rebuild + support check + dense signature) | C x S x 20 to 50 us | the dense signature is the largest piece and grows with F |
| certificates: Ch x (2 colourings + Graph build + conversion + nauty) | Ch x 0.2 to 0.4 ms | runs even when Ch = 1; repeated at saturation |
| child construction and bookkeeping | Ch x 5 to 15 us + ~0.1 ms | set copies, frontier rescan, state key |

Calibration against the two measured cases: a TS2-like unit with C of 1 to 3
costs 1 to 3 ms per step, times ~60 to 90 steps, times a few growths per
unit, which reproduces the measured ~1.4 CPU-s per unit. A TS8-like unit with
C of 10 to 40 costs 10 to 65 ms per step and forks many more growths, which
reproduces the measured ~18 CPU-s per unit.

## 3. Optimisation strategy for Layer 0

Ordered by gain per unit of change. Each entry states why the output cannot
differ.

### L0-1. Singleton dedupe early-out
`_dedup_sym_cands` returns its input unchanged when it receives one candidate,
both per step and at saturation.
Exactness: with one candidate the function computes a certificate, sees a
count of one, builds no boundary signature, and returns that same object.
Effect: removes one full certificate per step wherever growth is
unambiguous, which is most steps on asymmetric regions. Per-step factor
1.2 to 1.6x when C is small; 1.05 to 1.2x when C is large.

### L0-2. Neighbourhood-restricted candidate generation
In `_collect_free_target_entries`, iterate only over
`same_element(n) AND (intersection over fixed bonded u of N_P(m(u))) AND
(for each open-block bonded u: union over pool atoms p of N_P(p))`, in
ascending atom order; hoist the block index maps, the fixed-image set and the
edge-weight table out of the per-atom loop; group with the compact relation
signature (O(degree) instead of O(F)) and order groups by the dense key
computed once per group. Fall back to the full scan for callable node policies.
Exactness: every atom outside that set is rejected today by
`_growth_edge_supported`, which requires a product bond at or above the graph
floor to every mapped bonded neighbour, and graph edges are exactly the pairs
above that floor. Ascending order keeps the retained witness and the child
order identical. The compact signature is documented lossless when the orbit
map carries a structural zero bucket; equality is what grouping uses, and
ordering keeps the dense key.
Effect: scan term shrinks by the measured 11x from the restriction and a
further 5 to 10x from the compact signature. The scan stops being the
dominant term. Per-step factor 2 to 4x.

### L0-3. Exact pre-bucketing before nauty
Group children by the colour profile that the certificate key already
contains, then by the sorted multiset of (target orbit id, role) over
role-carrying atoms; run nauty only inside buckets of size two or more; build
the colouring once per certificate; construct one canonicaliser per
`grow_island` because the locked mapping is constant within it.
Exactness: the merge key today is (certificate, colour profile), so children
in different profile buckets never merge. Any colour-preserving isomorphism
between two role-coloured graphs is an automorphism of the underlying
WBO-coloured target graph at the same tolerance, hence preserves nauty orbit
ids; equal certificates therefore imply equal orbit multisets, and splitting
by the multiset never separates candidates that would have merged. Boundary
signatures still run inside collision classes as today. Fall back to full
certificates when the orbit map is a plain dict or its tolerance differs.
Effect: 3 to 10x fewer nauty calls on asymmetric targets, 2 to 5x on
complexes with equivalent ligands. After L0-1 and L0-2 the certificate term
is what remains, so this is worth another 1.5 to 3x per step.

### L0-4. Skip repeated growth attempts on an unchanged branch state
In `find_islands`, remember (branch state key, seed) pairs already attempted;
skip the call when the same pair recurs.
Exactness: `grow_island` is deterministic in inputs that are all captured by
the existing branch state key; a repeated key means the earlier attempt
produced no state change, so the repeat cannot either. Only the
`growth_calls` metric changes.
Effect: removes the guaranteed no-progress final pass and repeated failing
seeds. 1.0 to 1.3x per unit.

### L0-5. Deferred-atom short-circuit
In `grow_island`, when an atom already failed extension for the current
candidate list, record later edges to it as deferred without re-running the
extension, unless the earlier failure involved the support-state cap.
Exactness argument: within one growth the fragment only grows, locked images
are constant, candidates only refine, and block closure only adds rejections,
so an image rejected earlier stays rejected; the cap is the one non-monotone
path and is excluded by a flag. This needs a written proof and a shadow run
(compute both ways, assert equality) before landing.
Effect: 5 to 30% of growth on multi-island units, none on single-island ones.

### L0-6. Bookkeeping
Compute the branch state key once per child and cache it on the child; cache
frozen path keys for history merges; maintain the frontier edge set
incrementally; avoid the three fragment set copies per pop.
Exactness: same values, computed fewer times.
Effect: under 5%, but it makes profile attribution clean.

Not proposed: skipping the saturation dedupe after a final commit. The
saturation call uses a coarser boundary (deferred edges only, no frontier
edges) than the per-step call and can merge more candidates, so it is not
redundant.

## 4. Estimated acceleration, Layer 0 alone

| Regime | Per step today | After L0-1 to L0-6 | Factor |
|---|---:|---:|---:|
| low symmetry, C of 1 to 3 (TS2-like, most catalog rows) | 1 to 3 ms | 0.2 to 0.5 ms | 4 to 6x |
| high symmetry, C of 10 to 40 (TS8-like) | 10 to 65 ms | 1.5 to 8 ms | 6 to 10x |
| large asymmetric target, N_P = 183 | 20 to 25 ms per (candidate, step) | 0.5 to 1 ms | 20 to 40x |

Central estimate for the engine in pure Python: **4 to 8x on kernel time**,
with the residual dominated by nauty calls on genuine symmetry collisions and
by candidate object construction.

Native kernel (Rust with PyO3 or C++ with pybind11, vendored nauty, one
persistent target graph recoloured per call, integer-array candidate state,
bitset adjacency, identical predicates and identical support-cap semantics):
per step floor is roughly 0.05 to 0.5 ms, set by nauty on real collisions.
That is a further **10 to 30x** over the optimised Python, for a Layer 0
total of **50 to 200x on kernel time** relative to today.

What Layer 0 cannot change: the number of work units (cuts times seed orders),
the number of live branches the cap admits, the number of genuine symmetry
collisions that need nauty, and everything in Layers 1 to 3. On the whole
search stage the Python-only Layer 0 work translates to roughly 3 to 6x today
because scoring and the serial tail in Layer 1 stay as they are.

## 5. Measurement plan and acceptance

Counters to add to the existing growth profile dict (`island.py:65-86`),
all off by default:
`children_before_dedupe`, `certificates_computed`, `certificate_collisions`,
`boundary_signature_calls`, `scan_visits`, `admissible_visits`,
`single_child_steps`, `dedupe_elapsed_sec`, `scan_elapsed_sec`.
The ratio `certificates_computed / children_before_dedupe` after L0-3 and
`admissible_visits / scan_visits` after L0-2 are the two numbers that confirm
or refute the model above.

Replay corpus: record inputs and outputs of `_extend_sym_cands`,
`_dedup_sym_cands` and `grow_island` on the TEMPO example, the four contract
cases, and a few catalog rows (large asymmetric target, small source).
Every Layer 0 change is accepted only when the recorded outputs are
reproduced exactly, including order, and the 154 tests pass.

Cases to time: TEMPO locally (inputs are committed under
`docs/example_runs/pr1.tempo_ts3/work/endpoints`), Fe TS8 and Fe TS2 on the
cluster with `return_trace` enabled so the per-unit search/score split is
recorded.

## 6. Results so far (branch `Fable_AAM_Opt`)

Measured with `bench/replay_harness.py` (serial `search_aam`, one process,
Python 3.9, this Mac).  Every step was accepted only with an identical pool
(mechanisms, branches, hierarchies, encounter counts, cuts) and identical
running digests of every extension, dedupe and growth result against the
recording made with the original code.

| Step | TEMPO 57 atoms | tetraphenylmethane H-shift 45 atoms, 4 equivalent rings |
|---|---:|---:|
| original code | 18.8 s, 31,749 certificates | 16.0 s, 53,598 certificates |
| L0-1 singleton dedupe early-out | 17.9 s, 28,430 | |
| L0-3 orbit-role bucketing before nauty, one colouring per certificate | 13.3 s, 8,679 | |
| L0-2 admissible-neighbourhood candidate generation, hoisted block indexes | 10.7 s | |
| boundary-signature context hoisting and memo, hashed role key | 8.3 s | |
| reusable pynauty graph, shared role dictionaries, colour-order memo | 7.8 s | 4.9 s, 12,385 |
| one canonicalizer per island, WBO row cache, hoisted witness view in the scan | **7.55 s (2.5x)** | **4.71 s (3.4x)** |

Inside the engine (extension plus dedupe) the factor is 3.6x on TEMPO and
4.0x on the symmetric case; the remainder of the wall time is Layer 1 work
(in-sweep branch scoring and symmetry repair ~2.5 s, fragment-group
finalisation ~0.5 s, pool merging ~0.5 s on TEMPO), which is untouched by
design.  The measured engine residual is now: nauty itself ~0.7 s (the
floor), Python around the remaining certificates ~0.5 s, role keys ~0.6 s,
boundary signatures ~0.7 s, the restricted scan ~1.4 s, candidate object
construction ~0.8 s.  Further Python-level work on the engine is worth at
most another ~1.3x; the next step for Layer 0 is the native kernel, whose
acceptance test is the same replay corpus.

Environment note: the package declares Python >= 3.10, and the only
interpreter here is 3.9.  Two module-level type aliases were rewritten with
`typing.Union`/`Optional` (annotation-only) so the engine imports; the Layer 3
retrosynthesis tests need `int.bit_count` and could not be run here.  One
pre-existing Layer 2 test (`test_index_chirality_scores_every_valid_atom_action_before_rmsd_choice`)
fails on this numpy/LAPACK build with an RMSD of 2e-8 instead of 0 within 1e-12;
it does not touch the engine.
