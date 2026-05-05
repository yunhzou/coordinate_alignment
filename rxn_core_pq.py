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

from rxn_core_frag import (
    run_xtb, build_graph, expand_mapping, classify_bonds,
)


# -------------------- core grow --------------------

def _set_unique(cands):
    if not cands:
        return False
    if len(cands) == 1:
        return True
    s0 = frozenset(cands[0].values())
    return all(frozenset(c.values()) == s0 for c in cands[1:])


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


def _extend_cands_free(cands, fragment, n, g_R, g_P, iso_tol, max_cands_hard,
                       inv):
    """ext_atom n is unmapped. Extend each cand to include n.
    Require: element match, every fragment-bond |dWBO| <= iso_tol,
    and v not already mapped by another island (`inv`)."""
    bonded_in_frag = [u for u in g_R.neighbors(n) if u in fragment]
    if not bonded_in_frag:
        return None
    n_el = g_R.nodes[n]['element']
    r_wbos = [(u, g_R[u][n]['wbo']) for u in bonded_in_frag]
    new_cands = []
    for cand in cands:
        used_p = set(cand.values())
        v_set = set(g_P.neighbors(cand[bonded_in_frag[0]]))
        for u in bonded_in_frag[1:]:
            v_set &= set(g_P.neighbors(cand[u]))
        v_set -= used_p
        for v in v_set:
            if v in inv:
                continue
            if g_P.nodes[v]['element'] != n_el:
                continue
            if all(abs(w - g_P[cand[u]][v]['wbo']) <= iso_tol
                   for u, w in r_wbos):
                nc = dict(cand); nc[n] = v
                new_cands.append(nc)
                if len(new_cands) > max_cands_hard:
                    return None
    return new_cands


def _merge_island(cands, fragment, n, mapping, g_R, g_P, iso_tol):
    """ext_atom n is mapped (island). Each cand must accept n -> mapping[n]
    AND the bond from n's in-fragment R-neighbors must match in P."""
    p_n = mapping[n]
    new_cands = []
    for cand in cands:
        if p_n in cand.values():
            continue  # already used by this cand
        ok = True
        for u in g_R.neighbors(n):
            if u in fragment and u in cand:
                if not g_P.has_edge(cand[u], p_n):
                    ok = False; break
                if abs(g_R[u][n]['wbo'] - g_P[cand[u]][p_n]['wbo']) > iso_tol:
                    ok = False; break
        if ok:
            nc = dict(cand); nc[n] = p_n
            new_cands.append(nc)
    return new_cands


def grow_island_pq(g_R, g_P, seed, mapping, inv,
                   graph_floor=0.2,
                   iso_tol=0.5,
                   min_lock_size=2,
                   max_branches=8,
                   max_cands_hard=2000):
    """
    Grow a fragment from `seed` using priority-queue propagation.

    Returns a list of isos:
      []            -- failed (no initial cands, or fragment too small)
      [single_iso]  -- locked successfully (set-unique or single cand)
      [iso_a, ...]  -- non-set-unique saturation; caller branches
    """
    if seed in mapping:
        return []
    seed_el = g_R.nodes[seed]['element']
    cands = [{seed: v} for v in g_P.nodes()
             if v not in inv and g_P.nodes[v]['element'] == seed_el]
    if not cands:
        return []
    fragment = {seed}
    used_edges = set()
    heap = []
    _push_edges_from(heap, used_edges, g_R, seed, fragment, graph_floor)

    while heap:
        # Early-lock check before each pop (saves work on big symmetric scaffolds)
        if _set_unique(cands) and len(cands) == 1 and len(fragment) >= min_lock_size:
            return [cands[0]]

        neg_w, u, n = heapq.heappop(heap)
        edge = frozenset({u, n})
        if edge in used_edges:
            continue
        used_edges.add(edge)
        if n in fragment:
            continue  # already in (added through another edge)

        if n in mapping:
            new_cands = _merge_island(cands, fragment, n, mapping, g_R, g_P, iso_tol)
            if new_cands:
                cands = new_cands
                fragment.add(n)
                _push_edges_from(heap, used_edges, g_R, n, fragment, graph_floor)
            # else: edge consumed, continue
        else:
            new_cands = _extend_cands_free(
                cands, fragment, n, g_R, g_P, iso_tol, max_cands_hard, inv)
            if new_cands:
                cands = new_cands
                fragment.add(n)
                _push_edges_from(heap, used_edges, g_R, n, fragment, graph_floor)
            # else: edge consumed, continue

    # heap empty
    if not cands or len(fragment) < min_lock_size:
        return []
    if _set_unique(cands):
        return [cands[0]]
    # branch on distinct P-atom sets
    by_set = {}
    for c in cands:
        key = frozenset(c.values())
        if key not in by_set:
            by_set[key] = c
    branches = list(by_set.values())
    return branches[:max_branches]


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
    def commit(self, iso, g_R):
        # touched islands
        touched = set()
        for r in iso:
            if r in self.islands_R:
                touched.add(self.islands_R[r])
        if touched:
            iid = min(touched)
        else:
            iid = self.next_iid
            self.next_iid += 1
        for r, p in iso.items():
            if r not in self.mapping:
                self.mapping[r] = p
                self.inv[p] = r
            self.islands_R[r] = iid
            self.islands_P[p] = iid
        # transitive merge
        for r, k in list(self.islands_R.items()):
            if k in touched and k != iid:
                self.islands_R[r] = iid
                self.islands_P[self.mapping[r]] = iid


def find_islands_pq(g_R, g_P, seed_order,
                    graph_floor=0.2, iso_tol=0.5,
                    max_branches=8):
    """Run growth over a single seed ordering, branching on
    non-set-unique locks. Returns list of _Branch."""
    branches = [_Branch()]
    progressed = True
    while progressed:
        progressed = False
        for seed in seed_order:
            new_branches = []
            for b in branches:
                if seed in b.mapping:
                    new_branches.append(b)
                    continue
                isos = grow_island_pq(g_R, g_P, seed, b.mapping, b.inv,
                                      graph_floor=graph_floor, iso_tol=iso_tol,
                                      max_branches=max_branches)
                if not isos:
                    new_branches.append(b)
                    continue
                for iso in isos:
                    b2 = b.fork()
                    b2.commit(iso, g_R)
                    new_branches.append(b2)
                    progressed = True
            # cap by simple heuristic: prefer branches with more mapped atoms
            new_branches.sort(key=lambda b: -len(b.mapping))
            # dedupe by mapping signature
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
    nodes = list(g_R.nodes())
    rng = random.Random(rng_seed)
    orders = []
    while len(orders) < n_trials:
        perm = list(nodes); rng.shuffle(perm); orders.append(perm)
    return orders


def analyze_pq(reactant_xyz, product_xyz, workdir,
               charge=0, uhf=0,
               graph_floor=0.2, iso_tol=0.5,
               bond_high=0.5, dwbo_threshold=0.5,
               n_seeds=10, max_branches=8,
               chirality=True,
               return_all=False):
    """Run xtb on R and P, build graphs at graph_floor, run priority-queue
    grow with branching for n_seeds random orderings, score by
    (broken+formed, chirality_violations, -mapped). Returns best."""
    workdir = Path(workdir)
    elR, xyzR, wboR = run_xtb(reactant_xyz, workdir / "R", charge=charge, uhf=uhf)
    elP, xyzP, wboP = run_xtb(product_xyz, workdir / "P", charge=charge, uhf=uhf)
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
                mapping, wboR, wboP,
                bond_high=bond_high, dwbo_threshold=dwbo_threshold)
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
