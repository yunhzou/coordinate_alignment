"""Mechanism-independent AAM search with persistent fragment-decision graphs."""
from __future__ import annotations

import json
import gc
import multiprocessing as mp
import time
from dataclasses import asdict
from pathlib import Path

from .alignment.branch import _generate_seed_orders, find_islands
from .alignment.sweep import cut_sweep_items
from .domain import AAMProblem, AAMResult, AAMSearchConfig, AAMSearchMetrics
from .frag import build_graph
from .matcher import _nauty_orbits
from .search_graph import AAMSearchGraph
from .search_symmetry import finalize_graph_symmetry


_SEARCH_CONTEXT = None


def checkpoint_manifest(problem, config):
    """Identity required before reusing cut checkpoints (not reference labels)."""
    def endpoint(value):
        return dict(elements=list(value.elements), wbo=value.wbo.tolist(),
                    coordinates=value.coordinates.tolist())
    return json.loads(json.dumps(dict(schema='rxn_core.aam_checkpoints/v1',
        reactant=endpoint(problem.reactant),product=endpoint(problem.product),
        config=asdict(config))))


def _initialize_search(problem, config):
    global _SEARCH_CONTEXT
    target = build_graph(problem.product.elements, problem.product.wbo,
                         bond_cut=config.graph_floor)
    _SEARCH_CONTEXT = (problem, config, target,
                       _nauty_orbits(target, wbo_tol=config.iso_tolerance))


def _search_cut(cut):
    started = time.perf_counter()
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
        'search_seconds': time.perf_counter() - started,
        'max_live_branches': max((len(g.terminals) for g in graphs), default=0),
        'max_growth_candidates': max((row.get('max_cands_before', 0)
                                      for row in profile), default=0),
    }


def _initialize_finalization(problem,config):
    # Workers own acyclic archive data. Do not scan/copy a fork-inherited heap.
    gc.disable()
    _initialize_search(problem,config)


def _restore_finalized_cut(payload):
    from .artifacts import read_graph_checkpoint,write_graph_checkpoint
    index,data,checkpoint=payload
    if checkpoint is not None and Path(checkpoint).exists():
        return index,read_graph_checkpoint(checkpoint),{}
    graph=(data if isinstance(data,AAMSearchGraph) else
           AAMSearchGraph.from_record(json.loads(Path(data).read_bytes()),copy=False))
    _problem,config,target,_orbits=_SEARCH_CONTEXT
    graph,counts=finalize_graph_symmetry(graph,target,iso_tolerance=config.iso_tolerance)
    if checkpoint is not None:write_graph_checkpoint(graph,checkpoint)
    return index,graph,counts


def search_aam(problem: AAMProblem, config: AAMSearchConfig | None = None,
               *, workers: int = 1, intermediate_dir=None, resume=False,
               archive_format='json') -> AAMResult:
    """Return raw matching histories; repair/grouping/ranking are separate calls.

    When supplied, intermediate_dir receives each completed cut graph before
    group finalization and a final reusable graph record. No path or bijection
    enumeration is needed to collect worker results.

    resume=True reuses cuts only after verifying their input/configuration
    manifest. Invocation timings exclude the earlier checkpointed work;
    resulting graph and cut counts still describe the whole search.
    """
    if not isinstance(problem, AAMProblem):
        raise TypeError('search_aam requires an AAMProblem')
    config = config or AAMSearchConfig()
    requested_workers=max(1,int(workers))
    if archive_format not in ('json','checkpoint'):
        raise ValueError('archive_format must be json or checkpoint')
    started = time.perf_counter()
    cuts = cut_sweep_items(problem.reactant.wbo, config.cut_floor)
    directory = None if intermediate_dir is None else Path(intermediate_dir)
    if resume and directory is None:
        raise ValueError('Resuming requires an intermediate directory')
    if directory is not None:
        directory.mkdir(parents=True, exist_ok=True)
        manifest_path = directory / 'manifest.json'
        identity = checkpoint_manifest(problem,config)
        if resume:
            if json.loads(manifest_path.read_text()) != identity:
                raise ValueError('Checkpoint input or configuration differs from this search')
        else:
            if manifest_path.exists() or any(directory.glob('cut_*.json')):
                raise ValueError('Existing checkpoints require explicit resume=True')
            manifest_path.write_text(json.dumps(identity)+'\n')
    metrics = {'max_live_branches': 0, 'max_growth_candidates': 0}

    checkpoint_seconds = worker_search_seconds = 0.0
    def collect(payloads, indices):
        nonlocal checkpoint_seconds, worker_search_seconds
        for index, (graph, counts) in zip(indices,payloads):
            worker_search_seconds += counts['search_seconds']
            graphs[index] = graph
            for key in ('max_live_branches', 'max_growth_candidates'):
                metrics[key] = max(metrics[key], counts[key])
            if directory is not None:
                checkpoint_started = time.perf_counter()
                checkpoint = directory / f'cut_{index:05d}.json'
                temporary = checkpoint.with_suffix('.json.tmp')
                temporary.write_text(graph.to_json())
                temporary.replace(checkpoint)
                checkpoint_seconds += time.perf_counter() - checkpoint_started

    graphs = [None] * len(cuts)
    saved = {int(path.stem.split('_')[1]):path for path in directory.glob('cut_*.json')} if resume else {}
    missing = [index for index in range(len(cuts)) if index not in saved]
    # Finish independent missing cuts before rebuilding the large cached DAG.
    # This also prevents forked cut workers from inheriting that archive heap.
    workers = min(max(1, int(workers)), max(1,len(missing)))
    if workers == 1:
        _initialize_search(problem, config)
        collect(map(_search_cut, (cuts[index] for index in missing)),missing)
    else:
        with mp.get_context('fork').Pool(workers, initializer=_initialize_search,
                initargs=(problem, config)) as pool:
            collect(pool.imap(_search_cut, (cuts[index] for index in missing),
                             chunksize=config.task_chunksize),missing)
    # These JSON records are acyclic. Repeated cyclic-GC scans while restoring
    # millions of retained tuples needlessly revisit the entire growing graph.
    # Restore the caller's GC setting before any new matching begins.
    restore_started=time.perf_counter()
    if requested_workers>1:
        payloads=[(index,str(saved[index]) if index in saved else graphs[index],
                   str(directory/f'cut_{index:05d}.finalized.pkl.gz') if directory else None)
                  for index in range(len(cuts))]
        gc_enabled=gc.isenabled()
        try:
            gc.disable()
            with mp.get_context('fork').Pool(min(requested_workers,len(cuts)),
                    initializer=_initialize_finalization,initargs=(problem,config)) as pool:
                for index,graph,counts in pool.imap(_restore_finalized_cut,payloads,chunksize=1):
                    graphs[index]=graph
                    for key,value in counts.items():metrics[key]=metrics.get(key,0)+value
        finally:
            if gc_enabled:gc.enable()
        del payloads
    else:
        gc_enabled = gc.isenabled()
        tuple_pool = {}
        try:
            gc.disable()
            for index in sorted(saved):
                graphs[index] = AAMSearchGraph.from_record(json.loads(saved[index].read_bytes()),
                                                          copy=False,tuple_pool=tuple_pool)
        finally:
            if gc_enabled:gc.enable()
        del tuple_pool
    metrics['checkpoint_restore_and_finalize_seconds']=time.perf_counter()-restore_started
    for graph in graphs:
        metrics['max_live_branches']=max(metrics['max_live_branches'],
            max((sum(graph.states[t].context==context for t in graph.terminals)
                 for context in range(len(graph.contexts))),default=0))
    merge_started = time.perf_counter()
    gc_enabled = gc.isenabled()
    try:
        gc.disable()
        graph = AAMSearchGraph.combine(graphs)
    finally:
        if gc_enabled:gc.enable()
    metrics['parent_merge_seconds'] = time.perf_counter() - merge_started
    target = build_graph(problem.product.elements, problem.product.wbo,
                         bond_cut=config.graph_floor)
    symmetry_started = time.perf_counter()
    graph, groups = finalize_graph_symmetry(graph, target, iso_tolerance=config.iso_tolerance)
    metrics.update(symmetry_finalization_seconds=time.perf_counter()-symmetry_started,
                   worker_search_seconds=worker_search_seconds, checkpoint_seconds=checkpoint_seconds)
    for key,value in groups.items():metrics[key]=metrics.get(key,0)+value
    metrics.update(cuts=len(cuts), raw_result_count=len(graph.terminals),
        retained_branch_count=len(graph.terminals),
        subtree_branch_cap_count=sum(stop.reason == 'capped' for stop in graph.stops))
    result = AAMResult(problem, config, graph,
                      AAMSearchMetrics.from_record(metrics, time.perf_counter()-started))
    if directory is not None:
        from .artifacts import aam_json,write_aam_checkpoint
        if archive_format=='checkpoint':
            write_aam_checkpoint(result,directory/'aam.pkl.gz')
        else:
            temporary = directory / 'aam.json.tmp'
            temporary.write_text(aam_json(result))
            temporary.replace(directory / 'aam.json')
    return result
