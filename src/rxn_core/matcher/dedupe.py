"""Boundary-aware dedupe signatures for compressed candidates."""
from __future__ import annotations

from collections import Counter, defaultdict

from .orbits import _cand_canon_signature, _orbit_wbo_bucket
from .policy import as_node_match_policy
from .primitives import _edge_wbo, _orbit_id, _wbo_bucket
from .state import _SymCand, _cand_map, _cand_possible_p_atoms


def _p_relation_signature_from_parts(cand, v, g_P, p_orbits,
                                     cm_items=None, blocks=None,
                                     node_policy=None):
    """Signature of one possible target atom against a candidate state.

    Boundary dedupe calls this many times for the same candidate.  Accepting
    precomputed mapping/block views keeps the signature exact while avoiding
    repeated candidate materialization and sorting in the hot loop.
    """
    if cm_items is None:
        cm_items = tuple(sorted(_cand_map(cand).items()))
    if blocks is None:
        blocks = cand.blocks if isinstance(cand, _SymCand) else ()
    node_policy = as_node_match_policy(node_policy)

    rel = []
    for r, p in cm_items:
        if p == v:
            continue
        w = _edge_wbo(g_P, p, v)
        rel.append((r, _orbit_wbo_bucket(p_orbits, p, v, w)))
    block_rel = []
    for i, b in enumerate(blocks):
        edge_wbos = []
        for p in b.p_atoms:
            if p == v:
                continue
            w = _edge_wbo(g_P, p, v)
            edge_wbos.append(_orbit_wbo_bucket(p_orbits, p, v, w))
        block_rel.append((i, v in b.p_atoms, tuple(sorted(edge_wbos))))
    return (node_policy.key(g_P, v), _orbit_id(p_orbits, v),
            tuple(rel), tuple(block_rel))


def _p_relation_signature(cand, v, g_P, p_orbits, node_policy=None):
    return _p_relation_signature_from_parts(
        cand, v, g_P, p_orbits, node_policy=node_policy)


def _boundary_signature(cand, g_R, g_P, fragment=None, deferred_edges=(),
                        r_orbits=None, p_orbits=None, locked_mapping=None,
                        node_policy=None):
    if not fragment or not deferred_edges:
        return ()
    node_policy = as_node_match_policy(node_policy)
    fragment = set(fragment)
    cm = _cand_map(cand)
    cm_items = tuple(sorted(cm.items()))
    used_possible = _cand_possible_p_atoms(cand)
    locked_p_atoms = set((locked_mapping or {}).values())
    blocks = cand.blocks if isinstance(cand, _SymCand) else ()

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

    # For custom compatibility rules, two nodes with the same display element
    # can still have different target pools.  Cache by the concrete boundary
    # node to keep this signature exact for arbitrary policies.
    p_vec_by_node = {}

    def p_vec_for_node(x):
        cached = p_vec_by_node.get(x)
        if cached is not None:
            return cached
        target_sigs = []
        for v in g_P.nodes():
            if v in locked_p_atoms or v in used_possible:
                continue
            if not node_policy.compatible(g_R, x, g_P, v):
                continue
            target_sigs.append(_p_relation_signature_from_parts(
                cand, v, g_P, p_orbits, cm_items=cm_items, blocks=blocks,
                node_policy=node_policy))
        counts = Counter(target_sigs)
        p_vec = tuple(sorted(counts.items(), key=lambda item: str(item[0])))
        p_vec_by_node[x] = p_vec
        return p_vec

    for x in sorted(boundary):
        r_vec = tuple(
            (r, _orbit_id(r_orbits, r), _wbo_bucket(_edge_wbo(g_R, x, r)))
            for r in mapped_rs
        )
        x_key = node_policy.key(g_R, x)
        p_vec = p_vec_for_node(x)
        deferred = tuple(sorted(
            (r, _orbit_id(r_orbits, r), _wbo_bucket(_edge_wbo(g_R, x, r)))
            for r in deferred_by_outside.get(x, [])
        ))
        out.append((_orbit_id(r_orbits, x), x_key, r_vec, deferred, p_vec))
    return tuple(out)


def _dedup_sym_cands(cands, g_R, g_P, r_orbits=None, p_orbits=None,
                     fragment=None, deferred_edges=(), locked_mapping=None,
                     node_policy=None):
    if not cands:
        return cands
    internal_keys = []
    internal_counts = Counter()
    for cand in cands:
        if isinstance(cand, _SymCand):
            internal = cand.structural_signature(g_R, g_P, r_orbits, p_orbits)
        elif p_orbits is not None:
            internal = _cand_canon_signature(cand, p_orbits)
        else:
            internal = tuple(sorted(cand.items()))
        internal_keys.append(internal)
        internal_counts[internal] += 1

    seen = {}
    for cand, internal in zip(cands, internal_keys):
        if internal_counts[internal] == 1:
            # Boundary can only distinguish candidates that already share the
            # same internal symmetry signature.  Singleton internal classes are
            # necessarily unique under the full (internal, boundary) key.
            boundary = ()
        else:
            boundary = _boundary_signature(
                cand, g_R, g_P, fragment=fragment,
                deferred_edges=deferred_edges, r_orbits=r_orbits,
                p_orbits=p_orbits, locked_mapping=locked_mapping,
                node_policy=node_policy)
        sig = (internal, boundary)
        if sig not in seen:
            seen[sig] = cand
        elif isinstance(seen[sig], _SymCand) and isinstance(cand, _SymCand):
            seen[sig] = seen[sig].with_added_alternate(cand)
    return list(seen.values())
