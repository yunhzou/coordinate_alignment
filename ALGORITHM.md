# R↔P atom alignment — current algorithm

WBO-graph subgraph-isomorphism atom mapper with from-scratch re-enumeration,
heap-driven propagation, and chemistry-signature branching at saturation.

## Inputs

- Two xyz geometries: reactant `R` and product (or TS / IG) `T`.
- xtb GFN2 single-point on each → element list, coords, and Wiberg bond-order
  matrix (`wboR`, `wboT`).
- Identical composition: `Counter(elR) == Counter(elT)` (required).

## Output

A bijection `mapping: R-atom-index → T-atom-index` plus a list of *broken*
and *formed* bonds (the chemistry).

---

## Core principles

1. **Graphs only.** All growth decisions use the WBO graph; Cartesian coords
   are only used for the final spectator-chirality tiebreak.

2. **From-scratch re-enumeration at every grow step.** When the heap pops an
   edge that extends the fragment by one atom `n` (or absorbs a whole island
   into the fragment), we re-enumerate **all** valid subgraph isomorphisms of
   the new fragment against `g_P` via NetworkX `GraphMatcher`. The new `cands`
   list is whatever that enumeration returns — it does **not** depend on
   which `cands` we held one step earlier. The previous step's canonical
   choices cannot lock us into a wrong branch.

3. **Lock only at set-unique or saturation.** During the heap loop:
   - `set_unique(cands)` (every cand is the same bijection) **and** `len(cands)==1`
     → lock that single iso as an island, exit the heap loop.
   - Heap empty with `>1` non-unique cands → saturation. Lock by forking one
     branch per *chemistry-distinct* iso (see point 6).
   The algorithm never picks a canonical "winner" mid-propagation. Every
   surviving cand stays live until one of these two terminating conditions.

4. **Edge consumption rule.** Each `(frag_atom, ext_atom)` edge pops the heap
   exactly once. The re-enumeration either succeeds (cands replaced, fragment
   grown by `n` and possibly its whole island) or returns empty (`cands` would
   collapse to 0). On empty, the edge is **consumed** — marked used,
   fragment is unchanged, heap continues popping the next edge.

5. **Reseed across islands.** `find_islands_pq` iterates seeds in a given
   ordering. After one seed locks an island, the next seed grows from a still-
   unmapped atom; if that seed's edges later reach an already-mapped atom,
   the from-scratch re-enumeration must respect the existing mapping
   (`forced` images) — that's the island-merge path.

6. **Chemistry-signature dedup at iso forking.** When `grow_island_pq` returns
   K saturating isos, `find_islands_pq` does **not** fork K branches. It
   computes the `(broken_edges, formed_edges)` tuple for each iso (treating
   any R-edge whose endpoints' images aren't connected in P as broken,
   symmetrically for formed). Isos with identical `(broken, formed)` are
   spectator-permutations of the same mechanism — they collapse to one
   representative. Branch count grows only on chemistry-distinct alternatives
   (e.g. Mech A vs Mech B), not on spectator symmetries.

7. **Chirality tiebreak.** After all islands lock, compute chirality flips
   at spectator sp³ stereocenters (atoms whose neighbors are all mapped and
   not in `broken ∪ formed`). Score each branch:
   `(n_broken + n_formed,  chirality_violations,  -n_mapped)` — lex-min wins.

---

## State transitions during growth

```
cand-count transition           action
─────────────────────────       ─────────────────────────────────────────
>1 → >1                         commit, keep popping
>1 →  1                         commit, keep popping
 1 →  1                         commit, keep popping
 1 →  1   (set_unique met)      LOCK as island, exit heap loop, reseed
>1 →  0   (extension cuts all)  edge CONSUMED, fragment unchanged,
                                pop next heap entry
heap empty, set_unique           lock the unique iso
heap empty, >1 non-unique cands  saturation: dedup by chemistry signature,
                                FORK one branch per distinct (broken, formed)
```

---

## Grow algorithm — `grow_island_pq(g_R, g_P, seed, mapping, inv, …)`

```
if seed already in mapping: return []
cands = [{seed: v} for v in g_P
         if v not in inv and g_P.nodes[v].element == g_R.nodes[seed].element]
if not cands: return []

fragment = {seed}
heap     = pq of (-wbo, seed, neighbor) for every R-edge out of seed with WBO ≥ graph_floor
used_edges = ∅

while heap not empty:
    if set_unique(cands) and len(cands) == 1 and len(fragment) >= min_lock_size:
        return [cands[0]]

    (-wbo, u, n) = heap.pop()
    if frozenset({u, n}) in used_edges: continue
    used_edges.add(frozenset({u, n}))
    if n in fragment: continue

    candidate_fragment = fragment | {n}
    if n in mapping and islands_R has n:
        candidate_fragment ∪= whole island of n   # absorb the whole locked island

    new_cands = _compute_all_isos_FROM_SCRATCH(
        candidate_fragment, g_R, g_P, mapping, inv, iso_tol)

    if new_cands:
        cands     = new_cands             # full replacement, no inheritance from previous cands
        fragment  = candidate_fragment
        for r in fragment - old_fragment:
            push every outgoing R-edge of r into heap (if WBO ≥ graph_floor and not in used_edges)
    else:
        # edge consumed; fragment unchanged

# heap empty
if not cands or len(fragment) < min_lock_size:
    return []
if set_unique(cands):
    return [cands[0]]

# >1 non-unique saturation: return all isos, dedup happens in find_islands_pq
return list_of_distinct_isos_by_full_bijection
```

`_compute_all_isos_FROM_SCRATCH(fragment, g_R, g_P, mapping, inv, iso_tol)`:

```
sub_R   = g_R.subgraph(fragment)
forced  = {r: mapping[r]  for r in fragment if r in mapping}
node_match: element equality
edge_match: |wboR - wboP| ≤ iso_tol
for every subgraph isomorphism iso : g_P[?] ↔ sub_R via NetworkX GraphMatcher:
    rev = inverted iso (R-atom → P-atom)
    reject if any forced[r] != rev[r]                     # respect already-locked islands
    reject if any non-forced rev[r] is in global inv      # respect P-atoms already used by other islands
    keep
return the kept list
```

---

## Multi-island driver — `find_islands_pq(g_R, g_P, seed_order, …)`

```
branches = [empty _Branch()]
while progressed:
    for seed in seed_order:
        new_branches = []
        for b in branches:
            if seed in b.mapping: continue (carry b forward unchanged)
            isos = grow_island_pq(g_R, g_P, seed, b.mapping, b.inv, islands_R=b.islands_R)

            if not isos: continue (carry b forward unchanged)

            # chemistry-signature dedup: collapse spectator-permutation siblings
            seen_chem = {}
            for iso in isos:
                full_m = b.mapping | iso
                br_edges = {(u,v) for (u,v) in g_R.edges(full_m.keys())
                            if not g_P.has_edge(full_m[u], full_m[v])}
                fm_edges = {(inv[u], inv[v]) for (u,v) in g_P.edges(inv.keys())
                            if not g_R.has_edge(inv[u], inv[v])}
                key = (sorted(br_edges), sorted(fm_edges))
                seen_chem.setdefault(key, iso)

            for iso in seen_chem.values():
                b2 = b.fork(); b2.commit(iso); new_branches.append(b2)

        branches = dedup_by_full_bijection(new_branches)[:max_branches]

return branches
```

`_Branch.commit(iso, g_R)`:
- For each `(r, p)` in iso, set `mapping[r] = p`, `inv[p] = r`.
- All atoms in iso share an `island_id`. If iso atoms were already in islands
  (via prior island-merge paths), merge those islands into one.

---

## Outer multi-trial driver — `align_from_arrays(…)`

```
generate n_seeds random seed orderings (one per heavy atom + padding)
for each ordering:
    branches = find_islands_pq(g_R, g_P, ordering)
    for each branch b:
        full = expand_mapping(b.mapping, g_R, g_P)        # greedy fill of unmapped
        broken, formed = classify_bonds(full, wboR, wboP)
        chir = chirality_violations(full, xyzR, xyzP, broken, formed)
        score = (len(broken)+len(formed), chir, -len(full))
        all_results.append((score, full, broken, formed, chir))
return min(all_results, key=score)
```

---

## Optional perturbation layer — `cut_sweep` (in `build_bgcp_views_v2.py`)

For tracking *multiple* chemistry-distinct mechanisms (e.g. Mech A and Mech B
on `pr7.V.dodh_ts910`) the alignment is run repeatedly:

- once with `g_R` unchanged (baseline)
- once for each strong R-bond `(i, j)` with `WBO ≥ 0.5`, removing that edge
  from `g_R` before alignment

Each run can produce a different `(broken, formed)` chemistry. The pool is
the union of all distinct chemistries seen. With from-scratch re-enumeration
in `grow_island_pq`, the **baseline alone** already finds all min-bond
mechanisms in practice (`ts910`: 2 mechs in 0.5 s). `cut_sweep` is mostly
redundant under from-scratch and can be disabled when only the min-bond set
matters.

`expand_chemistry_relevant_atoms(mapping, core_atoms, g_R, g_P)`: a post-hoc
swap-based enumeration that, given a canonical bijection, generates
alternative bijections by swapping each core atom's image with a
locally-equivalent P-target (same element, same WBO to fragment neighbors,
same multiset of unmapped-neighbor (element, WBO)). Also largely redundant
under from-scratch.

---

## Bond classification — `classify_bonds(mapping, wboR, wboP, …)`

For each `(r1, r2)` pair with `r1 < r2`:

- Let `wR = wboR[r1, r2]`, `wP = wboP[mapping[r1], mapping[r2]]`.
- **Bond in R**: `wR ≥ bond_high` (default 0.5).
- **Bond in P**: `wP ≥ bond_high`.
- **Broken**: bond in R, not in P, and `|wR - wP| ≥ dwbo_threshold` (default 0.5).
- **Formed**: bond in P, not in R, and `|wR - wP| ≥ dwbo_threshold`.
- Same in P-frame symmetrically for atoms with `inv[p1], inv[p2]` defined.

`core_atoms_in_R_frame(mapping, broken, formed)`: union of R-atoms appearing
in any broken or formed bond. These are the chemistry-relevant atoms.

---

## Knobs (defaults in `src/rxn_core/pq.py`)

| name | default | meaning |
|---|---|---|
| `graph_floor` | 0.2 | edge inclusion in `g_R` / `g_P` |
| `iso_tol` | 1.0 | edge-WBO match tolerance during subgraph iso |
| `min_lock_size` | 1 | fragment must be ≥ this to lock |
| `dwbo_threshold` | 0.5 | classify_bonds: |ΔWBO| threshold |
| `max_branches` | 10⁶ | cap on concurrent live branches per seed ordering |
| `n_seeds` | 10 | random seed orderings in `align_from_arrays` |
| `chirality` | True | spectator-stereocenter tiebreak |
| `max_isos` (from-scratch) | 10⁹ | cap on isos returned per grow step (effectively no cap) |

---

## Verification

- `ts910` (pr7.V.dodh, 28 atoms): baseline finds Mech A and Mech B at 2/2
  each in ~0.5 s, 3 seed orderings.
- `pr14` (Pd hydroamination, 105 atoms): expected to find 2/3 mechanisms;
  NetworkX subgraph-iso cost grows with fragment size, so runtime depends
  on g_P's symmetry. Spectator-permutation isos are absorbed by the chemistry-
  signature dedup at iso forking, not by per-step grouping.
