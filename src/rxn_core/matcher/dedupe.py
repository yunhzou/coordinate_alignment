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


class _BoundaryContext:
    """Candidate-independent part of :func:`_boundary_signature`.

    Everything here depends only on the fragment, the deferred edges, the
    graphs, the R orbit map, the locked mapping and the node policy, all of
    which are shared by every candidate of one ``_dedup_sym_cands`` call.
    Building it once per call and reusing it changes nothing in the signature
    values; it only stops recomputing them per candidate.
    """

    def __init__(self, g_R, g_P, fragment, deferred_edges, r_orbits,
                 p_orbits, locked_mapping, node_policy):
        node_policy = as_node_match_policy(node_policy)
        self.g_R = g_R
        self.g_P = g_P
        self.r_orbits = r_orbits
        self.p_orbits = p_orbits
        self.node_policy = node_policy
        self.fragment = set(fragment)
        locked_p_atoms = set((locked_mapping or {}).values())
        boundary = set()
        deferred_by_outside = defaultdict(list)
        for raw in deferred_edges or ():
            edge = tuple(raw)
            if len(edge) != 2:
                continue
            a, b = edge
            a_in = a in self.fragment
            b_in = b in self.fragment
            if a_in == b_in:
                continue
            inside, outside = (a, b) if a_in else (b, a)
            boundary.add(outside)
            deferred_by_outside[outside].append(inside)
        # Element/attribute policies define compatibility entirely by
        # ``key``: boundary atoms with the same key share one target pool.  A
        # callable/custom policy may distinguish two same-key query atoms, so
        # it keeps a pool per concrete atom.
        key_equality = isinstance(
            node_policy, (ElementNodeMatchPolicy, AttributeNodeMatchPolicy))
        pools = {}
        self.entries = []
        for x in sorted(boundary):
            x_key = node_policy.key(g_R, x)
            pool_key = ('key', x_key) if key_equality else ('node', x)
            if pool_key not in pools:
                pools[pool_key] = tuple(
                    v for v in g_P.nodes()
                    if v not in locked_p_atoms
                    and node_policy.compatible(g_R, x, g_P, v))
            deferred = tuple(sorted(
                (r, _orbit_id(r_orbits, r), _wbo_bucket(_edge_wbo(g_R, x, r)))
                for r in deferred_by_outside.get(x, [])
            ))
            self.entries.append(
                (x, _orbit_id(r_orbits, x), x_key, pool_key, deferred))
        self.pools = pools
        self.structural_zero = getattr(p_orbits, 'zero_bucket', None)
        # Per target atom: policy key, orbit id and the nonzero-bucket
        # neighbour list, all candidate independent.
        self.target_static = {}
        if self.structural_zero is not None:
            for v in g_P.nodes():
                neighbours = []
                for p in g_P.neighbors(v):
                    bucket = _orbit_wbo_bucket(
                        p_orbits, p, v, _edge_wbo(g_P, p, v))
                    if bucket != self.structural_zero:
                        neighbours.append((p, bucket))
                self.target_static[v] = (
                    node_policy.key(g_P, v), _orbit_id(p_orbits, v),
                    tuple(neighbours))
        self._r_vec_cache = {}

    def r_vec(self, x, mapped_rs):
        key = (x, mapped_rs)
        cached = self._r_vec_cache.get(key)
        if cached is None:
            g_R, r_orbits = self.g_R, self.r_orbits
            cached = tuple(
                (r, _orbit_id(r_orbits, r), _wbo_bucket(_edge_wbo(g_R, x, r)))
                for r in mapped_rs)
            self._r_vec_cache[key] = cached
        return cached


def _boundary_signature(cand, g_R, g_P, fragment=None, deferred_edges=(),
                        r_orbits=None, p_orbits=None, locked_mapping=None,
                        node_policy=None, context=None, memo=None):
    """Deferred/frontier boundary evidence of one candidate.

    ``context`` may carry the shared :class:`_BoundaryContext` of the current
    dedupe call and ``memo`` a dict keyed by the candidate's literal state;
    two candidates with the same witness mapping and blocks have the same
    signature, so it is computed once.
    """
    if not fragment or not deferred_edges:
        return ()
    if context is None:
        context = _BoundaryContext(
            g_R, g_P, fragment, deferred_edges, r_orbits, p_orbits,
            locked_mapping, node_policy)
    node_policy = context.node_policy
    cm = _cand_map(cand)
    cm_items = tuple(sorted(cm.items()))
    blocks = cand.blocks if isinstance(cand, _SymCand) else ()
    memo_key = None
    if memo is not None:
        memo_key = (
            cm_items,
            tuple((b.r_atoms, b.p_atoms, bool(b.extendable)) for b in blocks))
        cached = memo.get(memo_key)
        if cached is not None:
            return cached
    used_possible = _cand_possible_p_atoms(cand)
    mapped_rs = tuple(sorted(r for r in context.fragment if r in cm))
    # Every pool atom is outside ``used_possible`` (mapped images and block
    # pools), so the compact relation signature of a pool atom sees the whole
    # witness mapping: hoist its inverse and domain out of the per-atom loop.
    mapped_r = tuple(r for r, _p in cm_items)
    inverse = {p: r for r, p in cm_items}
    target_static = context.target_static
    p_vec_by_pool = {}

    def p_vec_for_pool(pool_key):
        cached = p_vec_by_pool.get(pool_key)
        if cached is not None:
            return cached
        target_sigs = []
        for v in context.pools[pool_key]:
            if v in used_possible:
                continue
            static = target_static.get(v)
            if static is None or v in inverse:
                # Dense orbit map or unexpected image: exact general path.
                target_sigs.append(_p_relation_signature_from_parts(
                    cand, v, g_P, p_orbits, cm_items=cm_items, blocks=blocks,
                    node_policy=node_policy, compact=True))
                continue
            v_key, v_orbit, neighbours = static
            rel = []
            neighbor_buckets = {}
            for p, bucket in neighbours:
                neighbor_buckets[p] = bucket
                r = inverse.get(p)
                if r is not None:
                    rel.append((r, bucket))
            block_rel = []
            for i, b in enumerate(blocks):
                member = v in b.p_atoms
                edge_wbos = tuple(sorted(
                    neighbor_buckets[p] for p in b.p_atoms
                    if p != v and p in neighbor_buckets))
                block_rel.append(
                    (i, member, len(b.p_atoms) - int(member), edge_wbos))
            target_sigs.append((
                v_key, v_orbit, ('sparse', mapped_r, tuple(sorted(rel))),
                tuple(block_rel)))
        # This is a multiset used only for hashing/equality.  Sorting its very
        # large nested signatures by their string representation was the main
        # cost in large fragments and carries no semantic information.
        p_vec = frozenset(Counter(target_sigs).items())
        p_vec_by_pool[pool_key] = p_vec
        return p_vec

    out = []
    for x, x_orbit, x_key, pool_key, deferred in context.entries:
        out.append((x_orbit, x_key, context.r_vec(x, mapped_rs), deferred,
                    p_vec_for_pool(pool_key)))
    result = tuple(out)
    if memo is not None:
        memo[memo_key] = result
    return result


def _dedupe_certificates(canonicalizer, cands, p_orbits):
    """Certificate keys for ``cands``; nauty runs only where a merge is possible.

    Two candidates merge only when their exact coloured certificates and colour
    profiles agree.  Candidates are first grouped by the automorphism-invariant
    orbit-role key (``role_key``): equal certificates imply equal keys, so a
    class of size one cannot merge with anything, and a class whose role atoms
    all lie in singleton orbits has identical colourings throughout.  Neither
    needs a nauty call; the class key stands in for the certificate.  All other
    classes are certified exactly as before.  The stand-in is a tuple headed by
    a string and can never equal a real ``(bytes, profile)`` certificate, so
    the per-certificate counts used for boundary signatures are unchanged.
    """
    if not canonicalizer.role_keys_applicable(p_orbits):
        return [canonicalizer.certificate(cand) for cand in cands]
    keys = [canonicalizer.role_key(cand, p_orbits) for cand in cands]
    classes = defaultdict(list)
    for index, (key, _singleton) in enumerate(keys):
        classes[key].append(index)
    certificates = [None] * len(cands)
    for key, members in classes.items():
        if len(members) == 1 or keys[members[0]][1]:
            stand_in = ('orbit_role_key', key)
            for index in members:
                certificates[index] = stand_in
        else:
            for index in members:
                certificates[index] = canonicalizer.certificate(cands[index])
    return certificates


def _dedup_sym_cands(cands, g_R, g_P, r_orbits=None, p_orbits=None,
                     fragment=None, deferred_edges=(), locked_mapping=None,
                     node_policy=None):
    if not cands:
        return cands
    if len(cands) == 1:
        # One candidate cannot merge with anything.  The general path below
        # would compute a certificate, see a class of size one, skip the
        # boundary signature, and return this same object; skip straight to
        # that result.
        return list(cands)
    canonicalizer = _CandidateAutomorphismCanonicalizer(
        g_P, p_orbits=p_orbits, locked_mapping=locked_mapping,
        node_policy=node_policy)
    certificates = _dedupe_certificates(canonicalizer, cands, p_orbits)
    certificate_counts = Counter(certificates)
    seen = {}
    boundary_context = None
    boundary_memo = {}
    for cand, certificate in zip(cands, certificates):
        # Exact active-graph automorphism is the primary hierarchy.  Preserve
        # deferred/full-WBO evidence only inside an automorphic class;
        # this catches sub-floor distinctions without comparing every pair of
        # unrelated candidates.
        boundary = ()
        if certificate_counts[certificate] > 1 and deferred_edges:
            if boundary_context is None and fragment:
                boundary_context = _BoundaryContext(
                    g_R, g_P, fragment, deferred_edges, r_orbits, p_orbits,
                    locked_mapping, node_policy)
            boundary = _boundary_signature(
                cand, g_R, g_P, fragment=fragment,
                deferred_edges=deferred_edges, r_orbits=r_orbits,
                p_orbits=p_orbits, locked_mapping=locked_mapping,
                node_policy=node_policy, context=boundary_context,
                memo=boundary_memo)
        sig = (certificate, boundary)
        if sig not in seen:
            seen[sig] = cand
        elif isinstance(seen[sig], _SymCand) and isinstance(cand, _SymCand):
            # The colored pynauty certificate proves local automorphic
            # equivalence.  Keep one witness and its symbolic variation
            # domain; exact group generators are reconstructed only from the
            # completed AAM fragment relation.
            seen[sig] = seen[sig].with_automorph_equivalent(cand)
    return list(seen.values())
