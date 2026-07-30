"""Boundary-aware dedupe signatures for compressed candidates."""
from __future__ import annotations

from collections import Counter, defaultdict

from .canonical import _CandidateAutomorphismCanonicalizer
from .orbits import _cand_canon_signature, _orbit_wbo_bucket
from .policy import (
    AttributeNodeMatchPolicy,
    ElementNodeMatchPolicy,
    as_node_match_policy,
)
from .primitives import _edge_wbo, _orbit_id, _wbo_bucket
from .state import _SymCand, _cand_map, _cand_possible_p_atoms


def _p_relation_signature_from_parts(cand, v, g_P, p_orbits,
                                     cm_items=None, blocks=None,
                                     node_policy=None, compact=False):
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

    # Exact nauty orbit maps give structural zero its own bucket.  In that
    # representation all absent graph edges contribute the same zero value,
    # so recording them one by one is redundant.  Store the domain plus only
    # nonzero neighbor relations; this is losslessly equivalent to the former
    # dense vector and changes O(fragment_size) work per target atom into
    # O(target_degree) work.  Plain user-supplied orbit dicts lack this
    # structural-zero guarantee and retain the dense path below.
    structural_zero = getattr(p_orbits, 'zero_bucket', None)
    if compact and structural_zero is not None:
        active_cm = tuple((r, p) for r, p in cm_items if p != v)
        mapped_r = tuple(r for r, _ in active_cm)
        r_by_p = {p: r for r, p in active_cm}
        rel = []
        for p in g_P.neighbors(v):
            r = r_by_p.get(p)
            if r is None:
                continue
            w = _edge_wbo(g_P, p, v)
            bucket = _orbit_wbo_bucket(p_orbits, p, v, w)
            if bucket != structural_zero:
                rel.append((r, bucket))
        rel = ('sparse', mapped_r, tuple(sorted(rel)))

        neighbor_buckets = {}
        for p in g_P.neighbors(v):
            bucket = _orbit_wbo_bucket(
                p_orbits, p, v, _edge_wbo(g_P, p, v))
            if bucket != structural_zero:
                neighbor_buckets[p] = bucket
        block_rel = []
        for i, b in enumerate(blocks):
            member = v in b.p_atoms
            domain_size = len(b.p_atoms) - int(member)
            edge_wbos = tuple(sorted(
                neighbor_buckets[p] for p in b.p_atoms
                if p != v and p in neighbor_buckets
            ))
            block_rel.append(
                (i, member, domain_size, edge_wbos))
    else:
        rel = []
        for r, p in cm_items:
            if p == v:
                continue
            w = _edge_wbo(g_P, p, v)
            rel.append((r, _orbit_wbo_bucket(p_orbits, p, v, w)))
        rel = tuple(rel)
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
            rel, tuple(block_rel))


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

    # Element/attribute policies define compatibility entirely by ``key``.
    # Boundary nodes with the same key therefore have exactly the same target
    # pool, so compute that pool's candidate-dependent vector only once.  A
    # callable/custom policy may distinguish two same-key query nodes; retain
    # concrete-node caching for those policies to preserve exact semantics.
    compatibility_is_key_equality = isinstance(
        node_policy, (ElementNodeMatchPolicy, AttributeNodeMatchPolicy))
    p_vec_by_node = {}

    def p_vec_for_node(x):
        cache_key = (
            ('key', node_policy.key(g_R, x))
            if compatibility_is_key_equality else ('node', x)
        )
        cached = p_vec_by_node.get(cache_key)
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
                node_policy=node_policy, compact=True))
        counts = Counter(target_sigs)
        # This is a multiset used only for hashing/equality.  Sorting its very
        # large nested signatures by their string representation was the main
        # cost in large fragments and carries no semantic information.
        p_vec = frozenset(counts.items())
        p_vec_by_node[cache_key] = p_vec
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
    canonicalizer = _CandidateAutomorphismCanonicalizer(
        g_P, p_orbits=p_orbits, locked_mapping=locked_mapping,
        node_policy=node_policy)
    certificates = [canonicalizer.certificate(cand) for cand in cands]
    certificate_counts = Counter(certificates)
    kept = []
    representatives = []
    indices_by_signature = defaultdict(list)
    for cand, certificate in zip(cands, certificates):
        # Exact active-graph automorphism is the primary hierarchy.  Preserve
        # legacy deferred/full-WBO evidence only inside an automorphic class;
        # this catches sub-floor distinctions without comparing every pair of
        # unrelated candidates.
        boundary = ()
        if certificate_counts[certificate] > 1 and deferred_edges:
            boundary = _boundary_signature(
                cand, g_R, g_P, fragment=fragment,
                deferred_edges=deferred_edges, r_orbits=r_orbits,
                p_orbits=p_orbits, locked_mapping=locked_mapping,
                node_policy=node_policy)
        sig = (certificate, boundary)
        equivalent_index = None
        transporter = None
        for index in indices_by_signature[sig]:
            try:
                transporter = canonicalizer.transporter(
                    cand, representatives[index])
            except ValueError:
                # A pynauty certificate is a coarse partition certificate:
                # entire same-profile role cells can be canonically renamed.
                # Exact semantic-color transport is authoritative.
                continue
            equivalent_index = index
            break
        if equivalent_index is None:
            equivalent_index = len(kept)
            kept.append(cand)
            representatives.append(cand)
            indices_by_signature[sig].append(equivalent_index)
        elif (isinstance(kept[equivalent_index], _SymCand)
              and isinstance(cand, _SymCand)):
            kept_cand = kept[equivalent_index]
            # Exact automorphism transport guarantees that every continuation
            # of ``cand`` has an equivalent continuation of ``kept``.  Count
            # the represented states, but never enumerate alternate witnesses.
            # ``kept`` accumulates quotient metadata and consequently need
            # not retain the original certificate.  Transport every member
            # to the immutable first representative of this equivalence
            # class.
            kept[equivalent_index] = kept_cand.with_automorph_equivalent(
                cand, transporter=transporter)
    return kept
