"""Boundary-aware dedupe signatures for compressed candidates."""
from __future__ import annotations

from collections import Counter, defaultdict

from .orbits import _cand_canon_signature, _orbit_wbo_bucket
from .primitives import _edge_wbo, _orbit_id, _wbo_bucket
from .state import _SymCand, _cand_map, _cand_possible_p_atoms


def _p_relation_signature(cand, v, g_P, p_orbits):
    cm = _cand_map(cand)
    rel = []
    for r, p in sorted(cm.items()):
        if p == v:
            continue
        w = _edge_wbo(g_P, p, v)
        rel.append((r, _orbit_wbo_bucket(p_orbits, p, v, w)))
    block_rel = []
    if isinstance(cand, _SymCand):
        for i, b in enumerate(cand.blocks):
            edge_wbos = []
            for p in b.p_atoms:
                if p == v:
                    continue
                w = _edge_wbo(g_P, p, v)
                edge_wbos.append(_orbit_wbo_bucket(p_orbits, p, v, w))
            block_rel.append((i, v in b.p_atoms, tuple(sorted(edge_wbos))))
    return (g_P.nodes[v].get('element'), _orbit_id(p_orbits, v),
            tuple(rel), tuple(block_rel))


def _boundary_signature(cand, g_R, g_P, fragment=None, deferred_edges=(),
                        r_orbits=None, p_orbits=None, locked_mapping=None):
    if not fragment or not deferred_edges:
        return ()
    fragment = set(fragment)
    cm = _cand_map(cand)
    used_possible = _cand_possible_p_atoms(cand)
    locked_p_atoms = set((locked_mapping or {}).values())

    boundary = set()
    deferred_by_outside = defaultdict(list)
    for raw in deferred_edges or ():
        edge = tuple(raw)
        if len(edge) != 2:
            continue
        a, b = edge
        a_in = a in fragment
        b_in = b in fragment
        if a_in == b_in:
            continue
        inside, outside = (a, b) if a_in else (b, a)
        boundary.add(outside)
        deferred_by_outside[outside].append(inside)

    mapped_rs = sorted(r for r in fragment if r in cm)
    out = []
    for x in sorted(boundary):
        r_vec = tuple(
            (r, _orbit_id(r_orbits, r), _wbo_bucket(_edge_wbo(g_R, x, r)))
            for r in mapped_rs
        )
        target_sigs = []
        x_el = g_R.nodes[x].get('element')
        for v in g_P.nodes():
            if v in locked_p_atoms or v in used_possible:
                continue
            if g_P.nodes[v].get('element') != x_el:
                continue
            target_sigs.append(_p_relation_signature(cand, v, g_P, p_orbits))
        counts = Counter(target_sigs)
        p_vec = tuple(sorted(counts.items(), key=lambda item: str(item[0])))
        deferred = tuple(sorted(
            (r, _orbit_id(r_orbits, r), _wbo_bucket(_edge_wbo(g_R, x, r)))
            for r in deferred_by_outside.get(x, [])
        ))
        out.append((_orbit_id(r_orbits, x), x_el, r_vec, deferred, p_vec))
    return tuple(out)


def _dedup_sym_cands(cands, g_R, g_P, r_orbits=None, p_orbits=None,
                     fragment=None, deferred_edges=(), locked_mapping=None):
    if not cands:
        return cands
    seen = {}
    for cand in cands:
        if isinstance(cand, _SymCand):
            internal = cand.structural_signature(g_R, g_P, r_orbits, p_orbits)
        elif p_orbits is not None:
            internal = _cand_canon_signature(cand, p_orbits)
        else:
            internal = tuple(sorted(cand.items()))
        boundary = _boundary_signature(
            cand, g_R, g_P, fragment=fragment,
            deferred_edges=deferred_edges, r_orbits=r_orbits,
            p_orbits=p_orbits, locked_mapping=locked_mapping)
        sig = (internal, boundary)
        if sig not in seen:
            seen[sig] = cand
        elif isinstance(seen[sig], _SymCand) and isinstance(cand, _SymCand):
            seen[sig] = seen[sig].with_added_alternate(cand)
    return list(seen.values())
