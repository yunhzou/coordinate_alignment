"""Low-level WBO and orbit-id primitives for graph matching."""
from __future__ import annotations

import os

SYM_SUPPORT_MAX_STATES = 4096


def _load_fast_kernels():
    """Compiled ``_fast`` kernels when ``RXN_CORE_FAST=1`` and they import.

    The pure-Python functions remain the default.  A matcher module calls
    this once at import time and rebinds its leaf functions only when the
    optional extension (built by ``bench/build_fast.py``) is present and the
    variable is set; otherwise it returns None and nothing changes.
    """
    if os.environ.get("RXN_CORE_FAST") != "1":
        return None
    try:
        from . import _fast
    except ImportError:
        return None
    return _fast


def _edge_wbo(g, a, b):
    """WBO for any atom pair in the complete weighted graph.

    The matrix is read through a per-graph list-of-lists view built once
    with ``ndarray.tolist()``, whose Python floats equal ``float(mat[a, b])``
    exactly; this avoids a numpy scalar allocation and a ``float()`` call per
    lookup.  The view is keyed on the matrix object so a replaced matrix is
    re-read.
    """
    if a == b:
        return 0.0
    graph = g.graph
    mat = graph.get("wbo_matrix")
    if mat is not None:
        cached = graph.get("_wbo_rows")
        if cached is None or cached[0] is not mat:
            rows = None
            if hasattr(mat, "tolist") and getattr(
                    getattr(mat, "dtype", None), "kind", None) == "f":
                rows = mat.tolist()
            if rows is None:
                rows = [[float(value) for value in row] for row in mat]
            cached = (mat, rows)
            graph["_wbo_rows"] = cached
        return cached[1][a][b]
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
