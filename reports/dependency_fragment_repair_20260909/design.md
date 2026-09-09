# Experimental dependency repair and native scheduling

These are two separable components. Neither changes seed selection, the cut
sweep, matching tolerance, element compatibility, branch admission or mechanism
selection. Neither is dispatched by the production public AAM workflow.

## Conditional fragment repair

The old cut replay required an identical complete locked mapping and previous
island history. `FragmentRepair` records a smaller sufficient dependency set for
the **whole compressed growth result**, allowing an unaffected later fragment
to be reused after earlier mapping decisions change.

Session invariants: source elements/full WBO matrix, target graph and orbit
model stay fixed. A query fixes seed, graph floor, matching tolerance, minimum
size, branch cap, mapped-seed mode and the complete occupied-target mask.

Recorded dependencies:

1. Every source-topology row/pair consulted, including failed extension tests.
2. Every source mapping value consulted, including currently unmapped atoms.
3. Ordered membership of every consulted locked island, including all its
   mapped values and the absence of an island when that is observed.
4. The union of source atoms used in boundary comparisons during all attempted
   extensions and final saturation.
5. Boundary edges generated locally during the call, separately from inherited
   boundary edges.

For a later query with the same invariant key:

```text
if a consulted topology entry, mapping value, or island membership changed:
    recompute the conditional fragment result
else if an inherited boundary change touches the observed boundary region:
    recompute the conditional fragment result
else:
    reuse every candidate, correlated symmetry block, multiplicity and cap stop
    combine the new inherited boundary list with the recorded local boundaries
```

This is dependency validation, not acceptance of a single still-valid witness.
Newly unlocked source atoms are detected by the recorded unmapped-value reads.
Newly available target atoms change the occupied-target mask and force a new
query. Thus the implementation does not silently exclude newly enabled growth
or additional placements just because an old representative still works.

### Why target occupancy is a global dependency

Native candidate deduplication conditions target automorphisms on individually
locked target atoms. With the same occupied target set, this is the same
pointwise stabilizer. Renaming the distinct source labels on untouched locked
target vertices changes color names, not which candidate colorings are
equivalent. Any mapping actually consulted by growth, or belonging to a consulted
island, must still agree exactly. Final AAM transition groups are finalized
under the actual new parent state, not copied from an unrelated parent.

The invariant mask is deliberately global: a local geometric neighborhood alone
does not prove that the conditioned symmetry and allowed target domains agree.

### What this does not implement

This is a **dependency-aware conditional-result reuse primitive**, invoked by the
existing chronological search. It does not yet replace seed-times-cut scheduling
with dynamic backtracking or a new alternative-discovery policy. It does not
freeze a collection of independently matched fragments and assemble it blindly.
All existing conditional calls and branch admission decisions remain in place;
calls with changed relevant inputs recompute normally.

Native cache residency is bounded at 64 MiB per worker in these experiments.
Eviction changes time, not answers. The budget does not include live computation,
native target canonicalization tables, returned graphs or Python objects. Sessions
are serial, per fixed graph pair; workers do not share mutable native sessions.

## Native branch scheduler

`native/src/fragment_scheduler.h` ports the scheduling portion of
`alignment.branch.find_islands`, calling the existing `grow_island` function.
It can use the preceding repair primitive without passing each fragment through
Python. It retains:

- The same ordered seed traversal and repeated passes.
- The same source-island merging and exact live-state deduplication.
- Atomic subtree admission at the same branch cap.
- All state/transition/stop records, including rejected capped subtrees.
- Anchors, requested-core stopping and incomplete endpoint composition.
- Compressed symmetry records and the same later symmetry finalization.

The optimized frontier stores cached key hashes and references immutable keys;
it does not copy long serialized keys into every temporary map. A seed with no
eligible live branch can advance without rebuilding the already unique, sorted
frontier. Requested-core stops and trace positions still advance exactly. A
specialized one-branch carry operation applies the same join/cap rules without
allocating a temporary subtree container.

The Python boundary remains concrete and unchanged:

Build the native extension after checking out this experimental branch:

```bash
.venv/bin/python native/build_engine.py
```

```python
from rxn_core.native_search import find_islands_native

graph = find_islands_native(source_graph, target_graph, seed_order,
                           iso_tol=1.0, max_branches=100)
# graph is the existing AAMSearchGraph, not a new result representation.
```

The native entry point supports prepared WBO graphs with native element/orbit
matching. Python event-tracing and custom node policies are not silently
substituted; this separate experimental API does not accept those options.
The original Python scheduler and public production dispatch remain intact.

## Validation boundary

Acceptance compares the entire serialized graph, not only best scores or sampled
assignments. This covers representative mappings, fragment boundaries, compressed
symmetry, all transitions, terminal states and cap evidence. Final group
generation is also included in chemical-case comparisons.

Finite seed/cap search is still finite seed/cap search. Exact reuse of its results
does not establish exhaustive AAM discovery or chemical mechanism correctness.
The elementary-step set used here has no independent curated atom-map references.
