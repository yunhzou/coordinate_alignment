"""Low-level WBO and orbit-id primitives for graph matching."""
from __future__ import annotations

SYM_SUPPORT_MAX_STATES = 4096


def _edge_wbo(g, a, b):
    """WBO for any atom pair in the complete weighted graph."""
    if a == b:
        return 0.0
    mat = g.graph.get("wbo_matrix")
    if mat is not None:
        return float(mat[a, b])
    if g.has_edge(a, b):
        return float(g[a][b].get("wbo", 0.0))
    return 0.0


def _growth_edge_supported(w_R, w_P, iso_tol):
    """Compatibility for the popped growth edge.

    The popped edge is not special chemically; it is one entry in the same
    complete weighted WBO vector checked for every fragment atom.  Keep this
    helper only so trace/debug paths share the exact same tolerance policy.
    """
    return abs(float(w_R) - float(w_P)) <= iso_tol


def _wbo_bucket(w):
    return int(round(float(w) * 5))


def _orbit_id(orbits, node):
    return orbits[node] if orbits is not None else node
