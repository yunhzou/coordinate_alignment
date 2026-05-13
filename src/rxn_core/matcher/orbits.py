"""Orbit grouping for WBO graphs."""
from __future__ import annotations

from collections import defaultdict

import networkx as nx

from .primitives import _edge_wbo, _wbo_bucket


class _OrbitMap(dict):
    """dict-like orbit map with optional graph-specific WBO buckets."""

    def __init__(self, *args, wbo_buckets=None, zero_bucket=None,
                 wbo_tol=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.wbo_buckets = dict(wbo_buckets or {})
        self.zero_bucket = zero_bucket
        self.wbo_tol = wbo_tol


def _orbit_wbo_bucket(orbits, a, b, w):
    """WBO bucket consistent with an orbit backend when available."""
    if hasattr(orbits, "wbo_buckets"):
        key = (a, b) if a <= b else (b, a)
        if key in orbits.wbo_buckets:
            return orbits.wbo_buckets[key]
        if orbits.zero_bucket is not None:
            return orbits.zero_bucket
    return _wbo_bucket(w)


def _color_refine_orbits(g, iters=20):
    """Iterated color refinement (1-WL / Morgan extended connectivity) on
    an UNLABELED graph. Returns dict node -> int orbit_id such that two
    nodes have the same orbit_id iff they're indistinguishable under
    iterated neighbor-multiset refinement.

    Used to detect symmetric P-atom orbits: bijections that differ only
    by swapping P-atoms in the same orbit produce identical chemistry
    signatures and can be safely collapsed in the cand list.

    For molecules with little symmetry, every atom ends up in its own
    orbit (no collapse). For benzene rings, methyl Hs, B12 carborane,
    etc., symmetric atoms share an orbit and the K-factorial cand
    explosion is avoided.

    Algorithmic correctness: this is the standard 1-WL test used in
    InChI / RDKit canonical SMILES. It's exact for almost all molecular
    graphs and known to fail only on a small class of regular graphs
    (where 2-WL or full automorphism check is needed). For chemistry
    we're safe.
    """
    elements = nx.get_node_attributes(g, 'element')
    # Initial color: just element
    colors = {v: elements.get(v, '') for v in g.nodes()}
    # WBO bucket width: 0.2 so atoms whose WBO differs by ≤0.1 from xtb
    # noise (or geometry distortion at TS) end up in the same orbit.
    # round(x / 0.2) * 0.2 → bucket centers at 0.0, 0.2, 0.4, ...
    def _bucket(w):
        return int(round(w * 5))   # 0.2-wide buckets, integer-keyed
    for _ in range(iters):
        new = {}
        for v in g.nodes():
            nbr = tuple(sorted(
                (colors[w], _bucket(_edge_wbo(g, v, w)))
                for w in g.nodes()
                if w != v
            ))
            new[v] = (colors[v], nbr)
        # Compact to small ints (cheap downstream hashing)
        unique = {c: i for i, c in
                  enumerate(sorted(set(new.values()), key=str))}
        new_int = {v: unique[new[v]] for v in g.nodes()}
        if new_int == colors:
            break
        colors = new_int
    return colors


def _wbo_tolerance_bucket_lookup(g, tolerance):
    """Return a pair->bucket lookup for the complete WBO graph.

    Buckets are formed by greedy tolerance clustering of the WBO values present
    in this graph. This avoids a hard grid boundary, so values such as 1.0 and
    1.1 land in the same bucket when ``tolerance=0.2``.
    """
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    nodes = sorted(g.nodes())
    values = {0.0}
    pair_values = {}
    for i, a in enumerate(nodes):
        for b in nodes[i + 1:]:
            w = round(float(_edge_wbo(g, a, b)), 12)
            pair_values[(a, b)] = w
            values.add(w)

    reps = []
    value_to_bucket = {}
    eps = 1e-12
    for value in sorted(values):
        for idx, rep in enumerate(reps):
            if abs(value - rep) <= tolerance + eps:
                value_to_bucket[value] = idx
                break
        else:
            value_to_bucket[value] = len(reps)
            reps.append(value)
    return {pair: value_to_bucket[value]
            for pair, value in pair_values.items()}, value_to_bucket[0.0]


def _nauty_orbits(g, wbo_tol=0.2):
    """Exact automorphism orbits for a tolerance-bucketed WBO graph.

    Returns the same shape as :func:`_color_refine_orbits`: ``dict[node] ->
    orbit_id``.  Atom elements are vertex colors.  WBO edge colors are encoded
    by subdivision vertices colored by WBO bucket; the baseline 0-WBO bucket is
    represented by absence of a subdivision vertex.

    ``pynauty`` does the exact automorphism calculation on that colored graph.
    The WBO tolerance only controls bucket construction; the matcher still
    performs its normal complete-WBO ``iso_tol`` validity checks.
    """
    try:
        import pynauty
    except ImportError as exc:
        raise RuntimeError(
            "pynauty is required for _nauty_orbits; install pynauty or use "
            "_color_refine_orbits"
        ) from exc

    nodes = sorted(g.nodes())
    atom_index = {node: idx for idx, node in enumerate(nodes)}
    pair_buckets, zero_bucket = _wbo_tolerance_bucket_lookup(g, wbo_tol)

    adjacency = defaultdict(set)
    vertex_colors = defaultdict(set)
    elements = nx.get_node_attributes(g, 'element')
    for node, idx in atom_index.items():
        vertex_colors[('atom', elements.get(node, ''))].add(idx)

    next_idx = len(nodes)
    for (a, b), bucket in sorted(pair_buckets.items()):
        if bucket == zero_bucket:
            continue
        edge_idx = next_idx
        next_idx += 1
        ai = atom_index[a]
        bi = atom_index[b]
        adjacency[ai].add(edge_idx)
        adjacency[bi].add(edge_idx)
        adjacency[edge_idx].update((ai, bi))
        vertex_colors[('wbo', bucket)].add(edge_idx)

    adjacency_dict = {
        idx: sorted(adjacency.get(idx, ()))
        for idx in range(next_idx)
    }
    coloring = [
        set(vertices)
        for _, vertices in sorted(vertex_colors.items(), key=lambda item: str(item[0]))
    ]
    nauty_graph = pynauty.Graph(
        next_idx, directed=False,
        adjacency_dict=adjacency_dict,
        vertex_coloring=coloring)
    _, _, _, raw_orbits, _ = pynauty.autgrp(nauty_graph)

    atom_orbit_groups = defaultdict(list)
    for node, idx in atom_index.items():
        atom_orbit_groups[raw_orbits[idx]].append(node)
    compact = {}
    for orbit_id, (_, group) in enumerate(
            sorted(atom_orbit_groups.items(), key=lambda item: min(item[1]))):
        for node in group:
            compact[node] = orbit_id
    return _OrbitMap(
        compact, wbo_buckets=pair_buckets,
        zero_bucket=zero_bucket, wbo_tol=wbo_tol)


def _cand_canon_signature(cand, p_orbits):
    """Canonical signature for a cand under unlabeled-g_P orbits.

    Two cands with the same signature differ only by permuting P-atoms
    within their unlabeled orbits — they're spectator-equivalent and
    produce identical chemistry signatures.
    """
    return tuple(sorted((r, p_orbits[p]) for r, p in cand.items()))


def _dedup_cands_by_orbit(cands, p_orbits):
    """Collapse cands that differ only by g_P orbit permutations."""
    if not cands or p_orbits is None or len(cands) == 1:
        return cands
    seen = {}
    for c in cands:
        sig = _cand_canon_signature(c, p_orbits)
        if sig not in seen:
            seen[sig] = c
    return list(seen.values())


def _group_nodes_by_signature(nodes, sig_fn):
    groups = defaultdict(list)
    for node in sorted(nodes):
        groups[sig_fn(node)].append(node)
    return [tuple(vs) for _, vs in sorted(groups.items(), key=lambda kv: str(kv[0]))]
