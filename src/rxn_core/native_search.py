"""Experimental native scheduling with the unchanged Python search-graph API.

No seed selection, cut policy, mechanism collection or production dispatch lives
here. Unsupported native graph inputs are errors, not another search method.
"""
from dataclasses import asdict
import time

from .alignment.branch import _normalize_anchor_map
from .growth import native
from .matcher import _nauty_orbits
from .search_graph import AAMSearchGraph, SearchContext


def find_islands_native(source, target, seed_order, *, graph_floor=.2, iso_tol=1.,
                        max_branches=100, p_orbits=None, r_orbits=None, cuts=(),
                        anchor_map=None, core_R=None, stop_when_core_mapped=False,
                        growth_replay=None, profile=None):
    order=list(dict.fromkeys(seed_order))
    anchors=_normalize_anchor_map(anchor_map,source,target)
    p_orbits=p_orbits if p_orbits is not None else _nauty_orbits(target,wbo_tol=iso_tol)
    if not native.available():
        raise ValueError('native scheduling requires the built native engine')
    repair=None;repair_cuts=()
    if growth_replay is not None:
        if growth_replay.source is not source or growth_replay.target is not target:
            raise ValueError('fragment repair view belongs to another graph pair')
        repair=growth_replay.engine
        if not isinstance(repair,native._engine.FragmentRepair):
            raise ValueError('native scheduling accepts a FragmentRepair session, not prefix replay')
        repair_cuts=growth_replay.cuts
    r=native.source_graph(source)
    p=growth_replay.target_view if growth_replay is not None else native.target_graph(target,p_orbits)
    if r is None or p is None:
        raise ValueError('native scheduling requires WBO graphs and exact target orbits')
    core=tuple(sorted(set(core_R or ())))
    context=SearchContext(tuple(sorted(source)),tuple(sorted(target)),tuple(order),tuple(cuts),
        core,tuple(sorted(anchors.items())),float(graph_floor),float(iso_tol),int(max_branches),
        objective='requested_core' if stop_when_core_mapped and core else 'full_source')
    cpu=time.process_time() if profile is not None else 0.
    wall=time.perf_counter() if profile is not None else 0.
    record=native._engine.search_fragments(r.graph,p.graph,[r.index[a] for a in order],
        float(graph_floor),float(iso_tol),int(max_branches),
        [(r.index[a],p.index[b]) for a,b in sorted(anchors.items())],
        [r.index[a] for a in core],bool(stop_when_core_mapped),r.nodes,p.nodes,repair,repair_cuts)
    native_cpu=time.process_time()-cpu if profile is not None else 0.
    native_wall=time.perf_counter()-wall if profile is not None else 0.
    growth_calls=record.pop('growth_calls')
    record.update(schema='rxn_core.aam_search_graph/v1',contexts=(asdict(context),))
    cpu=time.process_time() if profile is not None else 0.
    graph=AAMSearchGraph.from_record(record,copy=False)
    if profile is not None:
        profile.append(dict(growth_calls=growth_calls,native_schedule_and_export_cpu=native_cpu,
            python_graph_cpu=time.process_time()-cpu,extend_elapsed_sec=native_wall))
    return graph
