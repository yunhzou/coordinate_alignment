"""Mechanism-independent AAM search with persistent fragment-decision graphs."""
from __future__ import annotations

import json
import multiprocessing as mp
import time
from pathlib import Path

from .alignment.branch import _generate_seed_orders, find_islands
from .alignment.sweep import cut_sweep_items
from .domain import AAMProblem, AAMResult, AAMSearchConfig, AAMSearchMetrics
from .frag import build_graph
from .matcher import _nauty_orbits
from .search_graph import AAMSearchGraph
from .search_symmetry import finalize_graph_symmetry


_SEARCH_CONTEXT = None


def _initialize_search(problem, config):
    global _SEARCH_CONTEXT
    target = build_graph(problem.product.elements, problem.product.wbo,
                         bond_cut=config.graph_floor)
    _SEARCH_CONTEXT = (problem, config, target,
                       _nauty_orbits(target, wbo_tol=config.iso_tolerance))


def _search_cut(cut):
    problem, config, target, target_orbits = _SEARCH_CONTEXT
    source = build_graph(problem.reactant.elements, problem.reactant.wbo,
                         bond_cut=config.graph_floor)
    source.remove_edges_from(cut)
    source_orbits = _nauty_orbits(source, wbo_tol=config.iso_tolerance)
    graphs, profile = [], []
    for order in _generate_seed_orders(source, n_trials=config.seed_count):
        graphs.append(find_islands(source, target, order,
            graph_floor=config.graph_floor, iso_tol=config.iso_tolerance,
            max_branches=config.branch_limit, p_orbits=target_orbits,
            r_orbits=source_orbits, anchor_map=dict(config.anchors),
            profile=profile, cuts=cut))
    graph = AAMSearchGraph.combine(graphs)
    return graph, {
        'max_live_branches': max((len(g.terminals) for g in graphs), default=0),
        'max_growth_candidates': max((row.get('max_cands_before', 0)
                                      for row in profile), default=0),
    }


def search_aam(problem: AAMProblem, config: AAMSearchConfig | None = None,
               *, workers: int = 1, intermediate_dir=None) -> AAMResult:
    """Return raw matching histories; repair/grouping/ranking are separate calls.

    When supplied, intermediate_dir receives each completed cut graph before
    group finalization and a final reusable graph record. No path or bijection
    enumeration is needed to collect worker results.
    """
    if not isinstance(problem, AAMProblem):
        raise TypeError('search_aam requires an AAMProblem')
    config = config or AAMSearchConfig()
    started = time.perf_counter()
    cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor)
    directory = None if intermediate_dir is None else Path(intermediate_dir)
    if directory is not None:
        directory.mkdir(parents=True, exist_ok=True)
    graphs, metrics = [], {'max_live_branches': 0, 'max_growth_candidates': 0}

    def collect(payloads):
        for index, (graph, counts) in enumerate(payloads):
            graphs.append(graph)
            for key, value in counts.items():
                metrics[key] = max(metrics[key], value)
            if directory is not None:
                (directory / f'cut_{index:05d}.json').write_text(
                    json.dumps(graph.to_record()) + '\n')

    workers = min(max(1, int(workers)), len(cuts))
    if workers == 1:
        _initialize_search(problem, config)
        collect(map(_search_cut, cuts))
    else:
        with mp.get_context('fork').Pool(workers, initializer=_initialize_search,
                initargs=(problem, config)) as pool:
            collect(pool.imap(_search_cut, cuts, chunksize=config.task_chunksize))
    merge_started = time.perf_counter()
    graph = AAMSearchGraph.combine(graphs)
    metrics['parent_merge_seconds'] = time.perf_counter() - merge_started
    target = build_graph(problem.product.elements, problem.product.wbo,
                         bond_cut=config.graph_floor)
    graph, groups = finalize_graph_symmetry(graph, target, iso_tolerance=config.iso_tolerance)
    metrics.update(groups)
    metrics.update(cuts=len(cuts), raw_result_count=len(graph.terminals),
        retained_branch_count=len(graph.terminals),
        subtree_branch_cap_count=sum(stop.reason == 'capped' for stop in graph.stops))
    result = AAMResult(problem, config, graph,
                      AAMSearchMetrics.from_record(metrics, time.perf_counter()-started))
    if directory is not None:
        from .artifacts import aam_record
        (directory / 'aam.json').write_text(json.dumps(aam_record(result)) + '\n')
    return result
