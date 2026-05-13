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


def _growth_edge_supported(w_R, w_P, iso_tol, graph_floor=0.2):
    """Compatibility for one active R-side growth/vector edge.

    The popped edge is not special chemically; it is one active R-side pair in
    the weighted WBO vector checked during extension.  If an R pair is active,
    the target pair must also be an active graph edge before the WBO tolerance
    is considered.  This prevents a loose ``iso_tol`` from treating a real
    bond as matching a target nonbond.
    """
    w_P = float(w_P)
    return w_P >= graph_floor and abs(float(w_R) - w_P) <= iso_tol


def _wbo_bucket(w):
    return int(round(float(w) * 5))


def _orbit_id(orbits, node):
    return orbits[node] if orbits is not None else node
