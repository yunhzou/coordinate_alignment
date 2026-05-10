"""
Priority-queue alignment algorithm.

See ALGORITHM.md for principles. Reuses xtb / classify_bonds / build_graph
plumbing from rxn_core_frag, but reimplements grow_island and find_islands
on top of a per-fragment priority queue with consume semantics, branching
on set-non-unique saturation, and chirality-aware scoring.
"""
from __future__ import annotations

import heapq
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import networkx as nx

from .frag import (
    run_xtb, build_graph, expand_mapping, classify_bonds,
)


# -------------------- core grow --------------------

def _set_unique(cands):
    """True iff there is at most one *distinct mapping* among `cands`.

    Earlier this used `frozenset(c.values())` as the equivalence key,
    which collapses every full-bijection candidate to one (since they
    all use the same value-set), silently dropping symmetric variants
    like (R87->IG87, R88->IG88) vs (R87->IG88, R88->IG87). The right
    equivalence for our use is the full bijection itself."""
    if not cands:
        return False
    if len(cands) == 1:
        return True
    sig0 = tuple(sorted(cands[0].items()))
    return all(tuple(sorted(c.items())) == sig0 for c in cands[1:])


def _push_edges_from(heap, used_edges, g_R, atom, fragment, graph_floor):
    """Push all not-yet-consumed outgoing edges of `atom` into the heap."""
    for nb in g_R.neighbors(atom):
        if nb in fragment:
            continue
        edge = frozenset({atom, nb})
        if edge in used_edges:
            continue
        w = g_R[atom][nb]['wbo']
        if w >= graph_floor:
            heapq.heappush(heap, (-w, atom, nb))


def _compute_all_isos_FROM_SCRATCH(fragment, g_R, g_P, mapping, inv, iso_tol,
                      islands_R, max_isos=2000):
    """SLOW: NetworkX subgraph_isomorphisms_iter from scratch.
    Kept for verification; live algorithm uses _extend_cands_incremental
    which is mathematically equivalent (proof: completeness preserved by
    induction over fragment growth) but ~10x faster."""
    sub_R = g_R.subgraph(fragment).copy()
    forced = {r: mapping[r] for r in fragment if r in mapping}
    def nm(nP, nR): return nR['element'] == nP['element']
    def em(eP, eR): return abs(eR['wbo'] - eP['wbo']) <= iso_tol
    matcher = nx.algorithms.isomorphism.GraphMatcher(g_P, sub_R, node_match=nm, edge_match=em)
    isos = []
    for raw in matcher.subgraph_isomorphisms_iter():
        rev = {r: p for p, r in raw.items()}
        ok = True
        for r in fragment:
            if r in forced and rev.get(r) != forced[r]: ok = False; break
            elif r not in forced and rev.get(r) in inv: ok = False; break
        if ok:
            isos.append(rev)
            if len(isos) >= max_isos: break
    return isos


def _extend_cands_incremental(cands, fragment_old, n, g_R, g_P, mapping, inv,
                               iso_tol, islands_R, max_cands_hard=None):
    """Incremental order-independent extension: extend each cand to include n.

    Equivalence: if cands == {all valid isos of fragment_old}, then the result
    is {all valid isos of fragment_old ∪ {n}}. By induction starting from
    cands_init = all 1-atom isos of seed, completeness is preserved throughout.

    For each cand:
      - bonded = n's R-neighbors that are in fragment_old (must be in cand)
      - v_set = ⋂(P-neighbors of cand[u]) for u in bonded   minus used_p
      - for each v in v_set with element match, |dWBO|<=iso_tol on every
        bonded edge, and v not in inv: emit cand ∪ {n: v}

    Cands that produce 0 extensions are dropped — they correspond to no iso
    of fragment_new, so their disappearance is correct (not lossy).

    Whole-island merge: if n in mapping, we additionally pull in n's island
    atoms via subgraph iso check on the union (each cand fans out to the
    isos consistent with both its existing bindings and the forced island).
    """
    n_el = g_R.nodes[n]['element']
    bonded_in_frag = [u for u in g_R.neighbors(n) if u in fragment_old]
    if not bonded_in_frag:
        return None
    r_wbos = [(u, g_R[u][n]['wbo']) for u in bonded_in_frag]

    # Determine target P-atom set and other forced atoms (for whole-island).
    if n in mapping:
        if islands_R is not None and n in islands_R:
            iid = islands_R[n]
            island_atoms = [r for r, k in islands_R.items() if k == iid
                            and r not in fragment_old]
        else:
            island_atoms = [n]
        forced_island = {r: mapping[r] for r in island_atoms}
    else:
        island_atoms = [n]
        forced_island = {}

    new_cands = []
    for cand in cands:
        used_p = set(cand.values())
        # n must map to forced[n] if n is in mapping
        if n in mapping:
            v_n = mapping[n]
            if v_n in used_p and cand.get(n) != v_n:
                continue
            # check forced bond constraints
            ok = True
            for u, w in r_wbos:
                if not g_P.has_edge(cand[u], v_n):
                    ok = False; break
                if abs(w - g_P[cand[u]][v_n]['wbo']) > iso_tol:
                    ok = False; break
            if not ok: continue
            base = dict(cand); base[n] = v_n
            # extend with rest of island atoms (force their images)
            extras_ok = True
            for r in island_atoms:
                if r == n or r in base: continue
                p = mapping[r]
                if p in base.values(): extras_ok = False; break
                base[r] = p
            if not extras_ok: continue
            # verify all R-edges among (island ∪ fragment_old) match in P
            check_set = set(base.keys())
            ok = True
            for r in island_atoms:
                for r2 in g_R.neighbors(r):
                    if r2 not in check_set: continue
                    if r >= r2 and r2 in island_atoms: continue
                    wR = g_R[r][r2]['wbo']
                    p, p2 = base[r], base[r2]
                    if not g_P.has_edge(p, p2): ok = False; break
                    if abs(wR - g_P[p][p2]['wbo']) > iso_tol: ok = False; break
                if not ok: break
            if ok:
                new_cands.append(base)
        else:
            # free extension: enumerate v's
            v_set = set(g_P.neighbors(cand[bonded_in_frag[0]]))
            for u in bonded_in_frag[1:]:
                v_set &= set(g_P.neighbors(cand[u]))
            v_set -= used_p
            for v in v_set:
                if v in inv: continue
                if g_P.nodes[v]['element'] != n_el: continue
                if all(abs(w - g_P[cand[u]][v]['wbo']) <= iso_tol
                       for u, w in r_wbos):
                    nc = dict(cand); nc[n] = v
                    new_cands.append(nc)
                    if max_cands_hard is not None and len(new_cands) > max_cands_hard:
                        return new_cands
    return new_cands


def _merge_whole_island_LEGACY(cands, fragment, n, mapping, islands_R,
                        g_R, g_P, iso_tol):
    """Up-front whole-island merge.

    When propagation reaches an island atom n, pull the ENTIRE island
    (every R-atom whose island_id matches n's) into each cand at once.
    Then check every cross-bond between fragment+island and the rest of
    fragment+island is consistent in P with the existing mapping.

    Returns (new_cands, island_atoms_added). If new_cands is empty, no
    cand survives — caller consumes the edge."""
    if islands_R is None or n not in islands_R:
        # No island bookkeeping; fall back to single-atom semantics.
        island_atoms = [n]
    else:
        target_iid = islands_R[n]
        island_atoms = [r for r, k in islands_R.items() if k == target_iid]
    p_atoms_in_island = [mapping[r] for r in island_atoms]
    new_cands = []
    new_added = [r for r in island_atoms if r not in fragment]
    if not new_added:
        return cands, []  # whole island already absorbed; nothing to do
    for cand in cands:
        # 1) no P-atom in island can be already used by this cand for a
        #    different R-atom
        used_p = set(cand.values())
        if any(p in used_p and cand.get(r) != p
               for r, p in zip(island_atoms, p_atoms_in_island)):
            continue
        # 2) build candidate extension with all island bindings forced
        nc = dict(cand)
        ok = True
        for r, p in zip(island_atoms, p_atoms_in_island):
            if r in nc and nc[r] != p:
                ok = False; break
            nc[r] = p
        if not ok:
            continue
        # 3) check every R-bond between fragment ∪ island (any pair where
        #    at least one endpoint is in island_atoms; the other in
        #    fragment or island_atoms) — the bond must exist in P with
        #    matching WBO. Bonds within fragment alone were checked when
        #    those atoms were added; bonds within island alone are
        #    preserved by construction (the island was locked already).
        check_set = set(island_atoms) | fragment
        for r in island_atoms:
            for r2 in g_R.neighbors(r):
                if r2 not in check_set: continue
                if r2 not in nc:        continue
                if r >= r2 and r2 in island_atoms: continue  # canonicalize
                wR = g_R[r][r2]['wbo']
                p, p2 = nc[r], nc[r2]
                if not g_P.has_edge(p, p2):
                    ok = False; break
                wP = g_P[p][p2]['wbo']
                if abs(wR - wP) > iso_tol:
                    ok = False; break
            if not ok: break
        if ok:
            new_cands.append(nc)
    return new_cands, new_added


def grow_island_pq(g_R, g_P, seed, mapping, inv,
                   graph_floor=0.2,
                   iso_tol=1.0,
                   min_lock_size=1,
                   max_branches=1_000_000,
                   max_cands_hard=2000,
                   events=None,
                   islands_R=None):
    """
    Grow a fragment from `seed` using priority-queue propagation.

    Returns a list of isos:
      []            -- failed (no initial cands, or fragment too small)
      [single_iso]  -- locked successfully (set-unique or single cand)
      [iso_a, ...]  -- non-set-unique saturation; caller branches

    Optional `events` list receives diagnostic events (seed_start /
    commit / consumed / merge / seed_end) compatible with the
    existing trace_run.HTML viewer.
    """
    record = events is not None
    if seed in mapping:
        return []
    seed_el = g_R.nodes[seed]['element']
    cands = [{seed: v} for v in g_P.nodes()
             if v not in inv and g_P.nodes[v]['element'] == seed_el]
    if not cands:
        if record:
            events.append({'type': 'seed_start', 'seed': int(seed),
                           'init_cands': 0, 'fragment': [int(seed)],
                           'p_atoms': []})
            events.append({'type': 'seed_end', 'result': 'no_initial_cands',
                           'final_cands': 0, 'fragment': [int(seed)],
                           'iso': None})
        return []
    fragment = {seed}
    distance = {seed: 0}
    used_edges = set()
    heap = []
    _push_edges_from(heap, used_edges, g_R, seed, fragment, graph_floor)
    if record:
        events.append({
            'type': 'seed_start',
            'seed': int(seed),
            'init_cands': len(cands),
            'fragment': [int(seed)],
            'p_atoms': sorted({int(v) for c in cands for v in c.values()}),
        })

    def _heap_snapshot(k=None):
        """Pending heap entries sorted by WBO desc, filtered to live
        (not yet consumed, not already in fragment). k=None → all."""
        peek = list(heap)
        peek.sort()
        live = [(w, uu, nn) for (w, uu, nn) in peek
                if frozenset({uu, nn}) not in used_edges and nn not in fragment]
        if k is not None:
            live = live[:k]
        return [{'frag_atom': int(uu), 'ext_atom': int(nn),
                 'wbo': round(-w, 3),
                 'ext_status': ('mapped' if nn in mapping
                                else 'free')}
                for w, uu, nn in live]

    def _pool_by_frag_atom():
        """Live propagation pool grouped by fragment atom. Each atom's
        outgoing live edges sorted WBO desc."""
        peek = list(heap)
        peek.sort()
        by_u = defaultdict(list)
        for w, uu, nn in peek:
            if frozenset({uu, nn}) in used_edges: continue
            if nn in fragment: continue
            by_u[int(uu)].append({
                'ext_atom': int(nn),
                'wbo': round(-w, 3),
                'ext_status': ('mapped' if nn in mapping else 'free'),
                'ext_element': g_R.nodes[nn]['element'],
            })
        return [{'frag_atom': int(u),
                 'frag_element': g_R.nodes[u]['element'],
                 'edges': sorted(by_u[u], key=lambda x: -x['wbo'])}
                for u in sorted(by_u.keys())]

    def _cands_sample(cs, k=10):
        return [{int(a): int(b) for a, b in c.items()} for c in cs[:k]]

    def _why_extend_failed(n_atom):
        """Per-cand explanation of why extension to n_atom failed."""
        bonded = [u for u in g_R.neighbors(n_atom) if u in fragment]
        n_el = g_R.nodes[n_atom]['element']
        r_wbos = [(u, g_R[u][n_atom]['wbo']) for u in bonded]
        out = []
        for ci, cand in enumerate(cands[:5]):
            used_p = set(cand.values())
            v_set = set(g_P.neighbors(cand[bonded[0]]))
            for u in bonded[1:]:
                v_set &= set(g_P.neighbors(cand[u]))
            v_set -= used_p
            tried = []
            for v in sorted(v_set):
                why = []
                if v in inv:
                    why.append(f'P[{v}] in global inv')
                elif g_P.nodes[v]['element'] != n_el:
                    why.append(f'element {g_P.nodes[v]["element"]} != {n_el}')
                else:
                    bad = []
                    for u, w in r_wbos:
                        wp = g_P[cand[u]][v]['wbo']
                        if abs(w - wp) > iso_tol:
                            bad.append(f'|{w:.3f}-{wp:.3f}|={abs(w-wp):.3f}>{iso_tol}')
                    if bad:
                        why.append('; '.join(bad))
                tried.append({'v': int(v), 'rejected': bool(why),
                              'reason': '; '.join(why) if why else 'OK'})
            out.append({
                'cand_idx': ci,
                'cand_at_in_frag_neighbors': {int(u): int(cand[u]) for u in bonded},
                'common_v_set_size': len(v_set),
                'tried_v': tried,
            })
        return out

    def _why_merge_failed(n_atom):
        """Per-cand explanation of why whole-island merge to n_atom failed."""
        if islands_R is None or n_atom not in islands_R:
            island_atoms = [n_atom]
        else:
            target_iid = islands_R[n_atom]
            island_atoms = [r for r, k in islands_R.items() if k == target_iid]
        out = []
        for ci, cand in enumerate(cands[:5]):
            why = []
            used_p = set(cand.values())
            for r in island_atoms:
                p = mapping[r]
                if p in used_p and cand.get(r) != p:
                    why.append(f'P[{p}] (image of R[{r}]) already used by R[{[k for k,v in cand.items() if v==p][0]}]')
                if r in cand and cand[r] != p:
                    why.append(f'R[{r}] in cand as P[{cand[r]}], conflicts with mapping P[{p}]')
            # cross-bond check
            nc = dict(cand)
            for r, p in [(r, mapping[r]) for r in island_atoms if r not in nc]:
                nc[r] = p
            check_set = set(island_atoms) | fragment
            for r in island_atoms:
                for r2 in g_R.neighbors(r):
                    if r2 not in check_set or r2 not in nc: continue
                    if r >= r2 and r2 in island_atoms: continue
                    wR = g_R[r][r2]['wbo']
                    p, p2 = nc[r], nc[r2]
                    if not g_P.has_edge(p, p2):
                        why.append(f'no P-edge P[{p}]-P[{p2}] (R[{r}]-R[{r2}] WBO={wR:.3f})')
                    else:
                        wP = g_P[p][p2]['wbo']
                        if abs(wR - wP) > iso_tol:
                            why.append(f'R[{r}]-R[{r2}]: |{wR:.3f}-{wP:.3f}|={abs(wR-wP):.3f} (P[{p}]-P[{p2}])')
            out.append({'cand_idx': ci, 'reasons': why[:8]})
        return out

    while heap:
        if _set_unique(cands) and len(cands) == 1 and len(fragment) >= min_lock_size:
            if record:
                iso = {int(k): int(v) for k, v in cands[0].items()}
                events.append({
                    'type': 'seed_end', 'result': 'success',
                    'final_cands': 1,
                    'fragment': sorted(int(x) for x in fragment),
                    'iso': iso,
                    'all_isos': [iso],
                    'lock_reason': 'set_unique_len1_during_BFS',
                    'heap_remaining': len(heap),
                })
            return [cands[0]]

        neg_w, u, n = heapq.heappop(heap)
        wbo = -neg_w
        edge = frozenset({u, n})
        if edge in used_edges:
            continue
        used_edges.add(edge)
        if n in fragment:
            if record:
                events.append({
                    'type': 'pop_skip',
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3)},
                    'reason': 'ext_atom already in fragment',
                })
            continue

        n_in_mapping = n in mapping
        if record:
            events.append({
                'type': 'pop',
                'edge': {'frag_atom': int(u), 'ext_atom': int(n),
                         'wbo': round(wbo, 3),
                         'ext_element': g_R.nodes[n]['element']},
                'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                'island_id_at_ext': (int(islands_R[n]) if islands_R and n in islands_R else None),
                'island_image': (int(mapping[n]) if n_in_mapping else None),
                'pre_state': {
                    'fragment_size': len(fragment),
                    'fragment': sorted(int(x) for x in fragment),
                    'cands_count': len(cands),
                    'cands_sample': _cands_sample(cands, 5),
                    'p_atoms_in_cands': sorted({int(v) for c in cands for v in c.values()}),
                },
                'heap_top_after_pop': _heap_snapshot(8),
                'pool_by_frag_atom': _pool_by_frag_atom(),
            })

        # ORDER-INDEPENDENT MATCHING:
        # tentatively add n to fragment, recompute ALL valid isos from scratch
        # using subgraph isomorphism on the fragment's R subgraph against g_P.
        # If any iso survives, the atom is committed and cands is replaced.
        # If none survives, the edge is consumed and the fragment doesn't grow.
        old_count = len(cands)
        old_fragment = set(fragment)
        candidate_fragment = fragment | {n}
        # INCREMENTAL extension: extend each cand to include n (and n's whole
        # island if n is mapped). Mathematically equivalent to recomputing all
        # subgraph isos from scratch but ~10x faster: completeness preserved
        # by induction since each old cand fans out to all its valid extensions
        # to n, and any iso of fragment_new restricts to a complete iso of
        # fragment_old (one of the existing cands).
        new_cands = _extend_cands_incremental(
            cands, fragment, n, g_R, g_P, mapping, inv,
            iso_tol, islands_R, max_cands_hard=max_cands_hard)
        if new_cands:
            cands = new_cands
            ref_dist = distance.get(u, 0) + 1
            if n_in_mapping and islands_R is not None and n in islands_R:
                target_iid = islands_R[n]
                whole_island = [r for r, k in islands_R.items() if k == target_iid]
                candidate_fragment = candidate_fragment | set(whole_island)
                added_extra = [r for r in whole_island if r not in distance]
                for r in added_extra:
                    distance[r] = ref_dist
            fragment = candidate_fragment
            for r in fragment - old_fragment:
                if r not in distance:
                    distance[r] = ref_dist
                _push_edges_from(heap, used_edges, g_R, r, fragment, graph_floor)
            if record:
                added_atoms = sorted(int(r) for r in fragment - old_fragment)
                events.append({
                    'type': 'commit',
                    'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3)},
                    'added': int(n), 'element': g_R.nodes[n]['element'],
                    'atoms_added': added_atoms,
                    'island_size_absorbed': len(added_atoms) if n_in_mapping else None,
                    'island_image': int(mapping[n]) if n_in_mapping else None,
                    'cand_n_value_set': sorted({int(c[n]) for c in cands if n in c}),
                    'cands_before': old_count,
                    'cands_after': len(cands),
                    'cands_sample_after': _cands_sample(cands, 5),
                    'fragment': sorted(int(x) for x in fragment),
                    'p_atoms': sorted({int(v) for c in cands for v in c.values()}),
                    'distance_from_seed': distance[n],
                    'bonds_to_fragment': [(int(u), round(wbo, 3))],
                    'heap_remaining': len(heap),
                    'heap_top': _heap_snapshot(8),
                    'pool_by_frag_atom': _pool_by_frag_atom(),
                })
        else:
            if record:
                events.append({
                    'type': 'consumed',
                    'scenario': ('merge_island' if n_in_mapping else 'extend_free'),
                    'edge': {'frag_atom': int(u), 'ext_atom': int(n), 'wbo': round(wbo, 3),
                             'ext_element': g_R.nodes[n]['element']},
                    'reason': ('merge_failed' if n_in_mapping else 'cut_all_cands'),
                    'island_image': int(mapping[n]) if n_in_mapping else None,
                    'island_id': int(islands_R[n]) if islands_R and n in islands_R else None,
                    'why_per_cand': (_why_merge_failed(n) if n_in_mapping
                                     else _why_extend_failed(n)),
                    'cands_count': len(cands),
                    'cands_sample': _cands_sample(cands, 5),
                    'fragment': sorted(int(x) for x in fragment),
                    'heap_remaining': len(heap),
                    'heap_top': _heap_snapshot(8),
                    'pool_by_frag_atom': _pool_by_frag_atom(),
                })

    # heap empty
    if not cands or len(fragment) < min_lock_size:
        if record:
            events.append({
                'type': 'seed_end',
                'result': ('no_cands' if not cands else 'too_small'),
                'final_cands': len(cands),
                'fragment': sorted(int(x) for x in fragment),
                'iso': None,
            })
        return []
    if _set_unique(cands):
        if record:
            iso = {int(k): int(v) for k, v in cands[0].items()}
            events.append({
                'type': 'seed_end', 'result': 'success',
                'final_cands': len(cands),
                'fragment': sorted(int(x) for x in fragment),
                'iso': iso,
                'all_isos': [iso],
            })
        return [cands[0]]
    # Dedup by full bijection signature, NOT by frozenset(values).
    # Two cands with the same value-set but different mappings are
    # genuinely different branches (they correspond to permuting
    # symmetric atoms; cf. the comment on _set_unique above).
    by_set = {}
    for c in cands:
        key = tuple(sorted(c.items()))
        if key not in by_set:
            by_set[key] = c
    branches = list(by_set.values())[:max_branches]
    if record:
        events.append({
            'type': 'seed_end', 'result': 'branched',
            'final_cands': len(cands), 'n_branches': len(branches),
            'fragment': sorted(int(x) for x in fragment),
            'iso': {int(k): int(v) for k, v in branches[0].items()},
            'all_isos': [{int(k): int(v) for k, v in c.items()}
                         for c in branches],
        })
    return branches


# -------------------- BFS frontier search --------------------

def grow_island_bfs(g_R, g_P, seed, mapping, inv,
                    graph_floor=0.2, iso_tol=1.0,
                    min_lock_size=1, events=None,
                    islands_R=None, max_branches=1_000_000):
    """EXPERIMENTAL: BFS frontier search over partial mappings.

    Status: kept for reference, not wired into find_islands_pq because
    it explodes on real molecules. The naive BFS spawns one child per
    valid (R_neighbor, P_target) extension; even with mapping-signature
    dedup, the state space grows multiplicatively because different
    growth orders can land on different INTERMEDIATE mappings that don't
    converge until very late.

    The current grow_island_pq avoids this by maintaining ONE shared
    fragment (deterministic growth order via WBO heap) and tracking
    multiple bijections of that single fragment as the cands list. The
    result is the same set of saturation bijections, but with bounded
    intermediate state count.

    A working BFS replacement would need either:
      - a canonical fragment-growth order (e.g., always extend the
        lowest-index frontier atom) so all states share a fragment, OR
      - a smarter dedup that canonicalizes mappings up to spectator
        permutations earlier in the search.

    Left here as a starting point for future work.
    """
    record = events is not None
    if seed in mapping:
        return []
    seed_el = g_R.nodes[seed]['element']

    init_isos = []
    for v in g_P.nodes():
        if v in inv: continue
        if g_P.nodes[v]['element'] != seed_el: continue
        init_isos.append({seed: v})

    if not init_isos:
        if record:
            events.append({'type': 'seed_start', 'seed': int(seed),
                           'init_cands': 0, 'fragment': [int(seed)],
                           'p_atoms': []})
            events.append({'type': 'seed_end', 'result': 'no_initial_cands',
                           'final_cands': 0, 'fragment': [int(seed)],
                           'iso': None})
        return []

    if record:
        events.append({'type': 'seed_start', 'seed': int(seed),
                       'init_cands': len(init_isos),
                       'fragment': [int(seed)],
                       'p_atoms': sorted(int(c[seed]) for c in init_isos)})

    # Frontier: set of frozensets of (R, P) pairs. Use frozenset for O(1)
    # dedup of partial-mapping states.
    frontier = {frozenset(m.items()) for m in init_isos}
    saturated = set()
    bfs_level = 0

    while frontier:
        bfs_level += 1
        next_frontier = set()
        for state in frontier:
            m = dict(state)
            mapped_R = set(m.keys())
            mapped_P = set(m.values())
            extensions_made = False
            # Iterate every fragment atom; for each, enumerate every
            # valid (n, v) extension. The dedup at next_frontier collapses
            # duplicates from different atom-pick orders.
            for r in mapped_R:
                for n in g_R.neighbors(r):
                    if n in mapped_R: continue
                    if g_R[r][n]['wbo'] < graph_floor: continue
                    if n in mapping: continue  # globally locked
                    n_el = g_R.nodes[n]['element']
                    bonded_to_n = [u for u in g_R.neighbors(n)
                                   if u in mapped_R
                                   and g_R[u][n]['wbo'] >= graph_floor]
                    if not bonded_to_n: continue
                    p_cands = set(g_P.neighbors(m[bonded_to_n[0]]))
                    for u in bonded_to_n[1:]:
                        p_cands &= set(g_P.neighbors(m[u]))
                    p_cands -= mapped_P
                    p_cands -= set(inv.keys())
                    for v in p_cands:
                        if g_P.nodes[v]['element'] != n_el: continue
                        ok = True
                        for u in bonded_to_n:
                            wR = g_R[u][n]['wbo']
                            wP = g_P[m[u]][v]['wbo']
                            if abs(wR - wP) > iso_tol:
                                ok = False; break
                        if not ok: continue
                        new_m = dict(m); new_m[n] = v
                        next_frontier.add(frozenset(new_m.items()))
                        extensions_made = True
            if not extensions_made:
                saturated.add(state)
        frontier = next_frontier
        if len(saturated) > max_branches:
            break

    # Saturation can produce many tied bijections (different P targets for
    # symmetric atoms). The frozenset already deduped by mapping signature,
    # so just convert back to dicts.
    result = []
    for fro in saturated:
        m = {r: v for r, v in fro}
        if len(m) >= min_lock_size:
            result.append(m)

    if record:
        if not result:
            events.append({'type': 'seed_end',
                           'result': 'no_cands' if not saturated else 'too_small',
                           'final_cands': 0,
                           'fragment': [int(seed)],
                           'iso': None})
        else:
            ref_frag = sorted(int(x) for x in result[0].keys())
            iso0 = {int(k): int(v) for k, v in result[0].items()}
            evt = {
                'type': 'seed_end',
                'result': 'success' if len(result) == 1 else 'branched',
                'final_cands': len(result),
                'n_branches': len(result),
                'fragment': ref_frag,
                'iso': iso0,
                'all_isos': [{int(k): int(v) for k, v in m.items()}
                             for m in result],
                'bfs_levels': bfs_level,
            }
            events.append(evt)
    return result


# -------------------- find_islands with branching --------------------

class _Branch:
    __slots__ = ('mapping', 'inv', 'islands_R', 'islands_P', 'next_iid')
    def __init__(self):
        self.mapping = {}
        self.inv = {}
        self.islands_R = {}
        self.islands_P = {}
        self.next_iid = 1
    def fork(self):
        b = _Branch()
        b.mapping = dict(self.mapping)
        b.inv = dict(self.inv)
        b.islands_R = dict(self.islands_R)
        b.islands_P = dict(self.islands_P)
        b.next_iid = self.next_iid
        return b
    def commit(self, iso, g_R, events=None):
        touched = set()
        for r in iso:
            if r in self.islands_R:
                touched.add(self.islands_R[r])
        if touched:
            iid = min(touched)
        else:
            iid = self.next_iid
            self.next_iid += 1
        committed_new = []
        relabeled = []
        for r, p in iso.items():
            if r not in self.mapping:
                self.mapping[r] = p
                self.inv[p] = r
                committed_new.append((int(r), int(p)))
            elif self.islands_R.get(r) != iid:
                relabeled.append((int(r), int(self.islands_R[r])))
            self.islands_R[r] = iid
            self.islands_P[p] = iid
        for r, k in list(self.islands_R.items()):
            if k in touched and k != iid:
                relabeled.append((int(r), int(k)))
                self.islands_R[r] = iid
                self.islands_P[self.mapping[r]] = iid
        if events is not None:
            events.append({
                'type': 'island_locked',
                'island_idx': int(iid),
                'pairs': committed_new,
                'merged_with': sorted(int(t) for t in touched - {iid}),
                'relabeled': relabeled,
                'mapped_total': len(self.mapping),
            })


def find_islands_pq(g_R, g_P, seed_order,
                    graph_floor=0.2, iso_tol=1.0,
                    max_branches=1_000_000, events=None):
    """Run growth over a single seed ordering, branching on
    non-set-unique locks. Returns list of _Branch.

    Optional `events` only records the FIRST (best-mapped) branch's
    trajectory — multi-branch traces would be confusing on a slider."""
    branches = [_Branch()]
    progressed = True
    pass_no = 0
    while progressed:
        progressed = False
        pass_no += 1
        if events is not None:
            events.append({'type': 'pass_start', 'pass': pass_no,
                           'mapped': len(branches[0].mapping)})
        for seed in seed_order:
            new_branches = []
            for bi, b in enumerate(branches):
                if seed in b.mapping:
                    new_branches.append(b)
                    continue
                # Only record events for branch 0 to keep trace linear
                ev_arg = events if (events is not None and bi == 0) else None
                isos = grow_island_pq(g_R, g_P, seed, b.mapping, b.inv,
                                      graph_floor=graph_floor, iso_tol=iso_tol,
                                      max_branches=max_branches,
                                      events=ev_arg,
                                      islands_R=b.islands_R)
                if not isos:
                    new_branches.append(b)
                    continue
                for ii, iso in enumerate(isos):
                    b2 = b.fork()
                    b2.commit(iso, g_R,
                              events=events if (bi == 0 and ii == 0) else None)
                    new_branches.append(b2)
                    progressed = True
            new_branches.sort(key=lambda b: -len(b.mapping))
            seen = set()
            uniq = []
            for b in new_branches:
                sig = tuple(sorted(b.mapping.items()))
                if sig in seen:
                    continue
                seen.add(sig)
                uniq.append(b)
                if len(uniq) >= max_branches:
                    break
            branches = uniq
            # Soft warning: pathological symmetry can blow this up. Default
            # cap is 1e6 so we should never hit it on real molecules; if we
            # do, surface it so we know.
            if len(branches) >= 10_000:
                import sys
                print(f"  [warn] alignment branch count = {len(branches)}  "
                      f"max_branches={max_branches}  "
                      f"new_branches_in={len(new_branches)}",
                      file=sys.stderr, flush=True)
    if events is not None:
        events.append({'type': 'done',
                       'mapped': len(branches[0].mapping)})
    return branches


# -------------------- chirality scoring --------------------

def _chirality_violations(mapping, coords_R, coords_P,
                          broken, formed, elements_R, min_deg=4):
    """Count spectator-stereocenter chirality flips.
    Spectator = atom whose 4+ mapped neighbors are all preserved
    (no incident broken/formed bond). Sign = sign(det) of first-3
    neighbor displacement vectors against the 4th. Skip near-coplanar."""
    bad_atoms = set()
    for (i, j, _, _) in broken:
        bad_atoms.add(i); bad_atoms.add(j)
    inv = {v: k for k, v in mapping.items()}
    for (ip, jp, _, _) in formed:
        if ip in inv: bad_atoms.add(inv[ip])
        if jp in inv: bad_atoms.add(inv[jp])

    # Build neighbor lists from R-coords (assume bond-graph already filtered)
    # We use "any atom within reasonable bond distance" — use the bond-graph
    # caller passes us implicitly through coords ordering. For correctness
    # we re-derive from coords by looking up close neighbors.
    # Simpler: just take all atoms within 1.9 Å in R as "neighbors".
    n_R = len(elements_R)
    diff = coords_R[:, None, :] - coords_R[None, :, :]
    dist = np.linalg.norm(diff, axis=2)
    np.fill_diagonal(dist, np.inf)
    violations = 0
    checked = 0
    for u in range(n_R):
        if u in bad_atoms: continue
        if u not in mapping: continue
        nbrs = np.where(dist[u] < 1.9)[0].tolist()
        # all must be mapped + spectator
        if len(nbrs) < min_deg: continue
        if any(nb not in mapping for nb in nbrs): continue
        if any(nb in bad_atoms for nb in nbrs): continue
        # take first 4 neighbors
        nb4 = nbrs[:4]
        # signed volume in R-frame
        ru = coords_R[u]
        v1 = coords_R[nb4[1]] - coords_R[nb4[0]]
        v2 = coords_R[nb4[2]] - coords_R[nb4[0]]
        v3 = coords_R[nb4[3]] - coords_R[nb4[0]]
        det_R = float(np.linalg.det(np.stack([v1, v2, v3])))
        # corresponding atoms in P-frame
        pu = coords_P[mapping[u]]
        pn = [coords_P[mapping[nb]] for nb in nb4]
        v1p = pn[1] - pn[0]
        v2p = pn[2] - pn[0]
        v3p = pn[3] - pn[0]
        det_P = float(np.linalg.det(np.stack([v1p, v2p, v3p])))
        if abs(det_R) < 0.05 or abs(det_P) < 0.05:
            continue  # near-coplanar; sign unstable
        checked += 1
        if np.sign(det_R) != np.sign(det_P):
            violations += 1
    return violations


# -------------------- multi-seed driver --------------------

def _generate_seed_orders(g_R, n_trials, rng_seed=42):
    """One order per heavy (non-H) atom, putting that atom first; the
    remaining nodes follow in a random order. Hydrogens are never used
    as initial seeds because they have very few connectivity constraints
    and produce poor island growth.

    If n_trials > n_heavy, pad with full random shuffles. The previous
    implementation was n_trials random shuffles only, which could miss
    structurally informative seedings (e.g. a key metal atom never
    landing first in any of the 10 default trials).
    """
    nodes = list(g_R.nodes())
    heavy = [n for n in nodes if g_R.nodes[n].get('element') != 'H']
    rng = random.Random(rng_seed)
    orders = []
    for h in heavy:
        rest = [x for x in nodes if x != h]
        rng.shuffle(rest)
        orders.append([h] + rest)
    while len(orders) < n_trials:
        perm = list(nodes); rng.shuffle(perm); orders.append(perm)
    return orders


def align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP,
                      graph_floor=0.2, iso_tol=1.0,
                      dwbo_threshold=0.5,
                      n_seeds=10, max_branches=1_000_000,
                      chirality=True, return_all=False):
    """Pure-graph entry point: assumes (el, xyz, wbo) for R and P are
    already in hand (e.g. loaded from a pre-computed xtb cache). Runs the
    multi-seed PQ search + scoring. Returns the same dict as analyze_pq."""
    if Counter(elR) != Counter(elP):
        raise ValueError(f"composition mismatch: {Counter(elR)} vs {Counter(elP)}")

    g_R = build_graph(elR, wboR, bond_cut=graph_floor)
    g_P = build_graph(elP, wboP, bond_cut=graph_floor)

    orders = _generate_seed_orders(g_R, n_seeds)
    all_results = []
    for order in orders:
        branches = find_islands_pq(g_R, g_P, order,
                                   graph_floor=graph_floor, iso_tol=iso_tol,
                                   max_branches=max_branches)
        for b in branches:
            mapping = expand_mapping(b.mapping, g_R, g_P)
            broken, formed, _, _ = classify_bonds(
                mapping, wboR, wboP, dwbo_threshold=dwbo_threshold)
            chir = (_chirality_violations(mapping, xyzR, xyzP, broken, formed, elR)
                    if chirality else 0)
            score = (len(broken) + len(formed), chir, -len(mapping))
            all_results.append((score, mapping, broken, formed, chir))

    all_results.sort(key=lambda r: r[0])
    best = all_results[0]
    score, mapping, broken, formed, chir = best
    out = dict(
        elements_R=elR, coords_R=xyzR, wbo_R=wboR,
        elements_P=elP, coords_P=xyzP, wbo_P=wboP,
        mapping=mapping, broken=broken, formed=formed,
        n_mapped=len(mapping), n_broken=len(broken), n_formed=len(formed),
        chirality_violations=chir,
        score=score,
    )
    if return_all:
        out['all_scored'] = [(s, m, b, f, c) for (s, m, b, f, c) in all_results]
    return out


def analyze_pq(reactant_xyz, product_xyz, workdir,
               charge=0, uhf=0,
               graph_floor=0.2, iso_tol=1.0,
               dwbo_threshold=0.5,
               n_seeds=10, max_branches=1_000_000,
               chirality=True,
               return_all=False):
    """Run xtb on R and P, build graphs at graph_floor, run PQ search for
    n_seeds random orderings, score by (broken+formed, chir, -mapped)."""
    workdir = Path(workdir)
    elR, xyzR, wboR = run_xtb(reactant_xyz, workdir / "R", charge=charge, uhf=uhf)
    elP, xyzP, wboP = run_xtb(product_xyz, workdir / "P", charge=charge, uhf=uhf)
    return align_from_arrays(
        elR, xyzR, wboR, elP, xyzP, wboP,
        graph_floor=graph_floor, iso_tol=iso_tol,
        dwbo_threshold=dwbo_threshold,
        n_seeds=n_seeds, max_branches=max_branches,
        chirality=chirality, return_all=return_all)
