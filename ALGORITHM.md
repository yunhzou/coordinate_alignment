# Priority-queue alignment algorithm

## High-level principles

1. **Two thresholds, three uses.**
   - `graph_floor = 0.2` — "what edges enter the search at all" (lowered from 0.5; lets partial bonds in TS-like geometries propagate).
   - `iso_tol = 0.5` — "what counts as a matched edge during ISO" (a fragment edge is preserved iff `|WBO_R − WBO_P| ≤ iso_tol`).
   - `bond_high = 0.5`, `dwbo_thr = 0.5` — "what counts as a real bond and a real change" at classification time.

2. **Fragment.** A connected, ISO-matched subgraph being grown from a seed atom. Internal edges all satisfy `iso_tol`. Cartesian coords play no role during growth — pure graph + WBO.

3. **Island.** A locked fragment that has been committed to `mapping`. Every island has a unique R-set ↔ P-set correspondence at the moment of locking (set-unique, sequence may be symmetric).

4. **Priority-queue propagation.** Each fragment maintains a global heap of `(−wbo, frag_atom, ext_atom)` candidate edges. The frontier is dynamic: when an atom is added to the fragment, its outgoing edges enter the heap. Each heap entry is *consumed exactly once* — popping pops it permanently, regardless of whether the extension succeeded. Termination = heap empty.

5. **Hits-free vs hits-island unified.** Same pop step, two outcomes:
   - `ext_atom` is unmapped (free) → try to extend cands. Success: add to fragment + push neighbors. Fail (`>1 → 0`): consumed.
   - `ext_atom` is mapped (island) → merge test: do all current cands have a consistent assignment of `ext_atom` to its island image? Success: fold in. Fail: consumed (this edge is a chemistry boundary).

6. **Branching at non-set-unique saturation.** When the heap empties with `>1` cands and the cands cover *different* P-atom sets, we don't pick arbitrarily — we fork into one branch per unique P-atom-set. Branches run the rest of `find_islands` independently; final scoring picks the best.

7. **Chirality as a 3D tiebreak.** For each branch's final mapping, count chirality violations at *spectator* sp³ stereocenters (atoms whose neighbors are all mapped and whose bonds are not in `broken ∪ formed`). Score = `(n_broken + n_formed, chirality_violations, −n_mapped)`.

## State transitions during growth

```
cand-count transition          action
─────────────────────────      ─────────────────────────────────
>1 → >1                        commit, keep growing
>1 →  1                        commit, keep growing (1→1 regime)
 1 →  1                        commit, keep growing
 1 →  0  (heap empty)          LOCK as island, reseed
>1 →  0  (extension cuts all)  pop is consumed; pick next heap entry
>1 →  0  (heap empty)          force lock; if set-unique pick one,
                               else BRANCH on different P-atom sets
```

## Why this generalizes

- **TS partial bonds:** floor 0.2 catches bonds in the WBO 0.3–0.6 partial-bond range that 0.5 misses.
- **Pseudo-symmetric atoms 2+ bonds out (the README known bug):** branching keeps both atom-set permutations alive; chirality scoring picks the chemistry-correct one.
- **Multi-mechanism handling:** unified hits-free / hits-island rule means there's only one place where chemistry boundaries are detected.
- **Deterministic given seed:** priority-queue order is fixed by WBO; only seed-choice randomness drives the multi-trial outer loop.

## Pipeline

```
xtb on R, P                    →  wboR, wboP
build_graph (floor 0.2)        →  g_R, g_P
                                                                
for seed_order in N random orderings:
    branches = [empty mapping]
    for seed in seed_order:
        new_branches = []
        for b in branches:
            isos = grow_island_pq(g_R, g_P, seed, b.mapping, …)
            # isos is [] (skip), [single_iso], or [iso_a, iso_b, …] (branch)
            if not isos: new_branches.append(b); continue
            for iso in isos:
                b' = b.fork(); b'.commit(iso)
                new_branches.append(b')
        branches = prune_to_top_k(new_branches, K=8)   # by partial score
    for b in branches:
        b.mapping = expand_mapping(b.mapping, g_R, g_P)
        b.broken, b.formed = classify_bonds(b.mapping, …)
        b.chir = chirality_violations(b.mapping, coords_R, coords_P,
                                      b.broken, b.formed)
        b.score = (len(b.broken)+len(b.formed), b.chir, -len(b.mapping))
    keep best b across all (seed_order, branches)
```

## Knobs (defaults in `rxn_core_pq.py`)

| name | default | meaning |
|---|---|---|
| `graph_floor` | 0.2 | edge inclusion in g_R / g_P |
| `iso_tol` | 0.5 | bond-WBO match within fragment ISO |
| `min_lock_size` | 2 | fragment must be ≥ this to lock |
| `bond_high` | 0.5 | classify_bonds: WBO is a bond |
| `dwbo_thr` | 0.5 | classify_bonds: change is real |
| `max_branches` | 8 | concurrent live branches per seed_order |
| `n_seeds` | 10 | random seed orderings |
| `chirality` | True | spectator-stereocenter tiebreak |
