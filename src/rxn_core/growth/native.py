"""Bridge to the optional native growth engine (``rxn_core._engine``).

The native engine is a C++ port of :func:`rxn_core.growth.island.grow_island`
and the matcher step it drives.  It is used only when the extension is
importable, ``RXN_CORE_NATIVE`` is not ``"0"``, and the call is one the port
covers (default element policy, exact nauty orbit map with a structural zero
bucket, no trace events).  Its outputs are the same ``_IsoResult`` objects the
Python engine returns; bench/compare_grow_calls.py replays recorded Python
calls through it to prove that.
"""
from __future__ import annotations

import os
import time

from ..matcher.orbits import _OrbitMap
from ..matcher.policy import ElementNodeMatchPolicy
from .result import IslandBranchLimitExceeded, _IsoResult

try:  # pragma: no cover - depends on the build
    from .. import _engine
except ImportError:  # pragma: no cover
    _engine = None


def available():
    return _engine is not None and os.environ.get("RXN_CORE_NATIVE", "1") != "0"


def _graph_edges(g):
    return sorted((min(a, b), max(a, b)) for a, b in g.edges())


def _wbo_rows(g):
    matrix = g.graph.get("wbo_matrix")
    if matrix is None:
        return None
    rows = matrix.tolist() if hasattr(matrix, "tolist") else [list(r) for r in matrix]
    return [[float(x) for x in row] for row in rows]


def source_graph(g_R):
    """Native view of the (possibly cut) reactant graph, cached on the graph."""
    cached = g_R.graph.get("_native_source")
    if cached is not None and cached[0] == g_R.number_of_edges():
        return cached[1]
    nodes = sorted(g_R.nodes())
    if nodes != list(range(len(nodes))):
        return None
    rows = _wbo_rows(g_R)
    if rows is None:
        return None
    elements = [str(g_R.nodes[v].get("element")) for v in nodes]
    native = _engine.SourceGraph(elements, rows, float(g_R.graph.get("bond_cut", 0.2)),
                                 _graph_edges(g_R))
    g_R.graph["_native_source"] = (g_R.number_of_edges(), native)
    return native


def target_graph(g_P, p_orbits):
    """Native view of the product graph plus its orbit map, cached on the map."""
    cache = p_orbits.__dict__.setdefault("_native_target", {})
    key = id(g_P)
    entry = cache.get(key)
    if entry is not None and entry[0] is g_P:
        return entry[1]
    nodes = sorted(g_P.nodes())
    if nodes != list(range(len(nodes))):
        return None
    rows = _wbo_rows(g_P)
    if rows is None or set(p_orbits) != set(nodes):
        return None
    elements = [str(g_P.nodes[v].get("element")) for v in nodes]
    pair_buckets = [(int(a), int(b), int(bucket))
                    for (a, b), bucket in p_orbits.wbo_buckets.items()]
    # the bucket table must describe exactly this graph's edges
    zero = p_orbits.zero_bucket
    for a, b, bucket in pair_buckets:
        if (bucket != zero) != g_P.has_edge(a, b):
            return None
    native = _engine.TargetGraph(
        elements, rows, float(g_P.graph.get("bond_cut", 0.2)), _graph_edges(g_P),
        [int(p_orbits[v]) for v in nodes], pair_buckets, int(zero))
    cache[key] = (g_P, native)
    return native


def applicable(g_R, g_P, p_orbits, node_policy, events, defer=False):
    if not available() or events is not None or defer:
        return False
    if not isinstance(node_policy, ElementNodeMatchPolicy):
        return False
    if not isinstance(p_orbits, _OrbitMap) or p_orbits.zero_bucket is None:
        return False
    if p_orbits.wbo_tol is None:
        return False
    return True


def grow_island(g_R, g_P, seed, mapping, *, graph_floor, iso_tol, min_lock_size,
                max_branches, islands_R, p_orbits, prior_deferred_edges,
                allow_mapped_seed, profile, profile_context):
    """Run the native engine; returns None when the inputs are not covered."""
    source = source_graph(g_R)
    target = target_graph(g_P, p_orbits)
    if source is None or target is None:
        return None
    n_r = g_R.number_of_nodes()
    image = [-1] * n_r
    for r, p in mapping.items():
        image[int(r)] = int(p)
    islands = None
    if islands_R is not None:
        islands = [(int(r), int(k)) for r, k in islands_R.items()]
    deferred = [(int(a), int(b)) for a, b in (prior_deferred_edges or ())]
    started = time.perf_counter()
    out = _engine.grow_island(
        source, target, int(seed), image, float(graph_floor), float(iso_tol),
        int(min_lock_size), int(max_branches), islands, deferred,
        bool(allow_mapped_seed))
    elapsed = time.perf_counter() - started
    if profile is not None:
        prof = {'seed': int(seed)}
        prof.update(out["profile"])
        prof['extend_elapsed_sec'] = elapsed
        prof['max_extend_elapsed_sec'] = 0.0
        prof['slowest_extend'] = None
        prof['elapsed_sec'] = elapsed
        if profile_context:
            prof.update(profile_context)
        profile.append(prof)
    if out["capped"]:
        raise IslandBranchLimitExceeded(
            out["cap_count"], out["cap_limit"], seed=int(seed))
    results = []
    for iso in out["isos"]:
        results.append(_IsoResult(
            iso["mapping"], deferred_edges=iso["deferred_edges"],
            fragment=iso["fragment"], symmetry=iso["symmetry"]))
    return results
