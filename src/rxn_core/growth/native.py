"""Bridge to the optional native growth engine (``rxn_core._engine``).

The native engine is a C++ port of :func:`rxn_core.growth.island.grow_island`
and the matcher step it drives.  It is used only when the extension is
importable, ``RXN_CORE_NATIVE`` is not ``0`` (enabled by default), and the call is one the port
covers (default element policy, exact nauty orbit map with a structural zero
bucket, no trace events).  Its outputs are the same ``_IsoResult`` objects the
Python engine returns. ``tests/test_native_engine.py`` compares growth calls,
symmetry state, cap behavior, and typed AAM outputs against the Python engine.
"""
from __future__ import annotations

import os
import time
from dataclasses import dataclass

from ..matcher.orbits import _OrbitMap
from ..matcher.policy import ElementNodeMatchPolicy
from .result import IslandBranchLimitExceeded, _IsoResult

try:  # pragma: no cover - depends on the build
    from .. import _engine
except ImportError:  # pragma: no cover
    _engine = None


def built():
    return _engine is not None


def available():
    return built() and os.environ.get("RXN_CORE_NATIVE", "1") != "0"


@dataclass(frozen=True)
class _NativeGraphView:
    graph: object
    nodes: tuple[int, ...]
    index: dict[int, int]


def _graph_edges(g, index):
    return sorted(
        (min(index[a], index[b]), max(index[a], index[b]))
        for a, b in g.edges()
    )


def _wbo_rows(g, nodes):
    matrix = g.graph.get("wbo_matrix")
    if matrix is None:
        return None
    return [[float(matrix[a][b]) for b in nodes] for a in nodes]


def source_graph(g_R):
    """Native view of the (possibly cut) reactant graph, cached on the graph."""
    cached = g_R.graph.get("_native_source")
    if cached is not None and cached[0] is g_R:
        return cached[1]
    nodes = tuple(sorted(map(int, g_R.nodes())))
    index = {atom: position for position, atom in enumerate(nodes)}
    edges = tuple(_graph_edges(g_R, index))
    rows = _wbo_rows(g_R, nodes)
    if rows is None:
        return None
    elements = [str(g_R.nodes[v].get("element")) for v in nodes]
    native = _engine.SourceGraph(elements, rows, float(g_R.graph.get("bond_cut", 0.2)),
                                 edges)
    view = _NativeGraphView(native, nodes, index)
    g_R.graph["_native_source"] = (g_R, view)
    return view


def target_graph(g_P, p_orbits):
    """Native view of the product graph plus its orbit map, cached on the map."""
    cache = p_orbits.__dict__.setdefault("_native_target", {})
    key = id(g_P)
    entry = cache.get(key)
    if entry is not None and entry[0] is g_P:
        return entry[1]
    nodes = tuple(sorted(map(int, g_P.nodes())))
    index = {atom: position for position, atom in enumerate(nodes)}
    rows = _wbo_rows(g_P, nodes)
    if rows is None or set(p_orbits) != set(nodes):
        return None
    elements = [str(g_P.nodes[v].get("element")) for v in nodes]
    pair_buckets = [
        (index[int(a)], index[int(b)], int(bucket))
        for (a, b), bucket in p_orbits.wbo_buckets.items()
    ]
    # the bucket table must describe exactly this graph's edges
    zero = p_orbits.zero_bucket
    for (a, b), bucket in p_orbits.wbo_buckets.items():
        if (bucket != zero) != g_P.has_edge(a, b):
            return None
    native = _engine.TargetGraph(
        elements, rows, float(g_P.graph.get("bond_cut", 0.2)),
        _graph_edges(g_P, index),
        [int(p_orbits[v]) for v in nodes], pair_buckets, int(zero))
    view = _NativeGraphView(native, nodes, index)
    cache[key] = (g_P, view)
    return view


def _translate_block(block, source_nodes, target_nodes):
    translated = dict(block)
    translated["r_atoms"] = [source_nodes[int(atom)]
                             for atom in block.get("r_atoms", ())]
    translated["p_atoms"] = [target_nodes[int(atom)]
                             for atom in block.get("p_atoms", ())]
    return translated


def _translate_symmetry(symmetry, source_nodes, target_nodes):
    symmetry = dict(symmetry)
    return {
        **symmetry,
        "witness": {
            source_nodes[int(source)]: target_nodes[int(target)]
            for source, target in dict(symmetry.get("witness", {})).items()
        },
        "blocks": [
            _translate_block(block, source_nodes, target_nodes)
            for block in symmetry.get("blocks", ())
        ],
        "exact_fixed": [
            source_nodes[int(atom)]
            for atom in symmetry.get("exact_fixed", ())
        ],
        "automorph_blocks": [
            _translate_block(block, source_nodes, target_nodes)
            for block in symmetry.get("automorph_blocks", ())
        ],
    }


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
    source_view = source_graph(g_R)
    target_view = target_graph(g_P, p_orbits)
    if source_view is None or target_view is None:
        return None
    source_index = source_view.index
    target_index = target_view.index
    if seed not in source_index:
        return None
    n_r = len(source_view.nodes)
    image = [-1] * n_r
    for r, p in mapping.items():
        if r not in source_index or p not in target_index:
            return None
        image[source_index[int(r)]] = target_index[int(p)]
    islands = None
    if islands_R is not None:
        if any(r not in source_index for r in islands_R):
            return None
        islands = [(source_index[int(r)], int(k))
                   for r, k in islands_R.items()]
    deferred = []
    for a, b in prior_deferred_edges or ():
        if a not in source_index or b not in source_index:
            return None
        deferred.append((source_index[int(a)], source_index[int(b)]))
    started = time.perf_counter()
    out = _engine.grow_island(
        source_view.graph, target_view.graph, source_index[int(seed)], image,
        float(graph_floor), float(iso_tol),
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
        translated_mapping = {
            source_view.nodes[int(source)]: target_view.nodes[int(target)]
            for source, target in dict(iso["mapping"]).items()
        }
        translated_deferred = [
            (source_view.nodes[int(a)], source_view.nodes[int(b)])
            for a, b in iso["deferred_edges"]
        ]
        translated_fragment = [
            source_view.nodes[int(atom)] for atom in iso["fragment"]
        ]
        results.append(_IsoResult(
            translated_mapping,
            deferred_edges=translated_deferred,
            fragment=translated_fragment,
            symmetry=_translate_symmetry(
                iso["symmetry"], source_view.nodes, target_view.nodes),
        ))
    return results
