"""Mechanism-independent AAM search with persistent fragment-decision graphs."""
from __future__ import annotations

import json
import hashlib
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


class _GrowthCounts:
    """Reduce worker profile events online instead of retaining every row."""
    def __init__(self):
        self.maximum = 0

    def append(self, row):
        self.maximum = max(self.maximum, row.get('max_cands_before', 0))


def cut_seed(cut, root_seed=42):
    """Distinct reproducible streams, independent of scheduling and cut order."""
    edges = tuple(sorted(tuple(sorted(edge)) for edge in cut))
    if not edges:
        return root_seed
    payload = json.dumps([root_seed, edges], separators=(',', ':')).encode('ascii')
    return int.from_bytes(hashlib.blake2b(payload, digest_size=16,
                                         person=b'AAM-cut-seeds-v1').digest(), 'big')


def checkpoint_manifest(problem, config):
    """Identity required before reusing cut checkpoints (not reference labels)."""
    def endpoint(value):
        return dict(elements=list(value.elements), wbo=value.wbo.tolist(),
                    coordinates=value.coordinates.tolist())
    return json.loads(json.dumps(dict(schema='rxn_core.aam_checkpoints/v2',
        seed_policy='independent_per_cut_blake2b_v1',
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
    graphs, profile = [], _GrowthCounts()
    for order in _generate_seed_orders(source, n_trials=config.seed_count, rng_seed=cut_seed(cut)):
        graphs.append(find_islands(source, target, order,
            graph_floor=config.graph_floor, iso_tol=config.iso_tolerance,
            max_branches=config.branch_limit, p_orbits=target_orbits,
            r_orbits=source_orbits, anchor_map=dict(config.anchors),
            profile=profile, cuts=cut))
    graph = AAMSearchGraph.combine(graphs)
    return graph, {
        'search_seconds': time.perf_counter() - started,
        'max_live_branches': max((len(g.terminals) for g in graphs), default=0),
        'max_growth_candidates': profile.maximum,
    }


def _search_cut_task(payload, *, in_process=False):
    """Checkpoint at the producer; a slow sibling cannot block persistence."""
    index, cut, checkpoint = payload
    graph, counts = _search_cut(cut)
    counts['checkpoint_seconds'] = 0.0
    if checkpoint is not None:
        from .artifacts import write_raw_cut
        started = time.perf_counter()
        path = Path(checkpoint)
        write_raw_cut(graph,path)
        counts['checkpoint_seconds'] = time.perf_counter()-started
        # No graph is pickled through IPC or buffered in the parent.
        return index, graph if in_process else str(path), counts
    return index, graph, counts


def _initialize_finalization(problem,config):
    # Workers own acyclic archive data. Do not scan/copy a fork-inherited heap.
    gc.disable()
    _initialize_search(problem,config)


def _restore_finalized_cut(payload):
    from .artifacts import read_graph_checkpoint,write_graph_checkpoint,read_raw_cut
    index,data,checkpoint=payload
    if checkpoint is not None and Path(checkpoint).exists():
        return index,read_graph_checkpoint(checkpoint),{}
    graph=(data if isinstance(data,AAMSearchGraph) else
           read_raw_cut(data))
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

    archive_format selects both raw-cut and final persistence: JSON interchange
    records, or trusted compact checkpoints that preserve shared storage.
    Completed cuts are persisted by workers and collected out of order, but
    the returned graph always retains canonical cut/context ordering.

    resume=True reuses cuts only after verifying their input/configuration
    manifest. Invocation timings exclude the earlier checkpointed work;
    resulting graph and cut counts still describe the whole search.
    """
    if not isinstance(problem, AAMProblem):
        raise TypeError('search_aam requires an AAMProblem')
    from .artifacts import raw_cut_paths,read_raw_cut
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
            if manifest_path.exists() or raw_cut_paths(directory):
                raise ValueError('Existing checkpoints require explicit resume=True')
            manifest_path.write_text(json.dumps(identity)+'\n')
    metrics = {'max_live_branches': 0, 'max_growth_candidates': 0}

    checkpoint_seconds = worker_search_seconds = 0.0
    def collect(payloads):
        nonlocal checkpoint_seconds, worker_search_seconds
        for index, graph, counts in payloads:
            worker_search_seconds += counts['search_seconds']
            checkpoint_seconds += counts['checkpoint_seconds']
            if isinstance(graph, str):
                saved[index] = Path(graph)
            else:
                graphs[index] = graph
            for key in ('max_live_branches', 'max_growth_candidates'):
                metrics[key] = max(metrics[key], counts[key])

    graphs = [None] * len(cuts)
    saved = {int(path.name.split('_')[1].split('.')[0]):path for path in raw_cut_paths(directory)} if resume else {}
    missing = [index for index in range(len(cuts)) if index not in saved]
    raw_suffix = '.raw.pkl.gz' if archive_format=='checkpoint' else '.json'
    tasks = [(index, cuts[index], str(directory/f'cut_{index:05d}{raw_suffix}') if directory else None)
             for index in missing]
    # Finish independent missing cuts before rebuilding the large cached DAG.
    # This also prevents forked cut workers from inheriting that archive heap.
    workers = min(max(1, int(workers)), max(1,len(missing)))
    if workers == 1:
        _initialize_search(problem, config)
        collect(_search_cut_task(task,in_process=requested_workers==1) for task in tasks)
    else:
        with mp.get_context('fork').Pool(workers, initializer=_initialize_search,
                initargs=(problem, config)) as pool:
            collect(pool.imap_unordered(_search_cut_task, tasks,
                                        chunksize=config.task_chunksize))
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
                for index,graph,counts in pool.imap_unordered(_restore_finalized_cut,payloads,chunksize=1):
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
                graphs[index] = read_raw_cut(saved[index],tuple_pool=tuple_pool)
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
