"""Exact differential tests, including correlated symmetry and cap stops."""
import numpy as np
import pytest

from rxn_core.cut_replay import CutReplay
from rxn_core.frag import build_graph
from rxn_core.growth import native
from rxn_core.matcher import _nauty_orbits
from rxn_core.alignment.branch import find_islands

pytestmark = pytest.mark.skipif(not native.built(), reason="native engine not built")


def graph(n, edges, elements=None):
    w = np.zeros((n,n))
    for edge in edges:
        a,b = edge[:2]
        w[a,b] = w[b,a] = edge[2] if len(edge)==3 else 1.
    return build_graph(elements or ('C',)*n,w,bond_cut=.2)


def raw(source,target,orbits,seed,mapping,cut,*,replay=None,cap=100,tol=1.,islands=None,deferred=()):
    if replay is None:
        cut_graph = source.copy()
        cut_graph.remove_edges_from(cut)
        view = native.source_graph(cut_graph)
        return native._engine.grow_island(view.graph,native.target_graph(target,orbits).graph,
            seed,mapping,.2,tol,1,cap,islands,deferred,False)
    result = replay.engine.grow(cut,seed,mapping,.2,tol,1,cap,islands,deferred,False)
    result.pop('reused_extensions')
    return result


@pytest.mark.parametrize('checkpoints',[False,True])
@pytest.mark.parametrize('tol',[.5,1.])
@pytest.mark.parametrize('reference_only',[False,True])
def test_every_cut_matches_fresh_native(checkpoints,tol,reference_only):
    rng = np.random.default_rng(819)
    for trial in range(4):
        edges = [(i,i+1,float(rng.choice([.8,1.2,1.8]))) for i in range(7)]
        edges += [(0,4,1.),(2,6,1.)]
        source=graph(8,edges)
        target=graph(8,edges if trial==0 else [(a,b,w) for a,b,w in edges[:-1]])
        orbits=_nauty_orbits(target,wbo_tol=tol)
        replay=CutReplay(source,target,orbits,checkpoints=checkpoints,reference_only=reference_only)
        cuts=[()]+[(e[:2],) for e in edges]+[((1,2),(4,5))]
        for seed in (0,4):
            for cap in (1,3,100):
                for cut in cuts:
                    expected=raw(source,target,orbits,seed,[-1]*8,cut,cap=cap,tol=tol)
                    actual=raw(source,target,orbits,seed,[-1]*8,cut,replay=replay,cap=cap,tol=tol)
                    assert actual==expected,(trial,seed,cap,cut)


def test_prefix_is_reused_before_distant_cut():
    source=graph(12,[(i,i+1) for i in range(11)])
    target=graph(12,[(i,i+1) for i in range(11)])
    orbits=_nauty_orbits(target,wbo_tol=1.)
    replay=CutReplay(source,target,orbits)
    for cut in ((),((8,9),),((4,5),),((8,9),)):
        assert raw(source,target,orbits,0,[-1]*12,cut,replay=replay)==raw(source,target,orbits,0,[-1]*12,cut)
    stats=replay.stats()
    assert stats['prefix_hits']>=1 and stats['full_hits']>=1
    assert stats['reused_extensions']>0


def test_unread_component_can_change_source_symmetry_and_still_reuse():
    source=graph(6,[(0,1),(1,2),(3,4),(4,5)])
    target=graph(6,[(0,1),(1,2),(3,4),(4,5)])
    orbits=_nauty_orbits(target,wbo_tol=1.)
    replay=CutReplay(source,target,orbits)
    for cut in ((),((4,5),)):
        assert raw(source,target,orbits,0,[-1]*6,cut,replay=replay)==raw(source,target,orbits,0,[-1]*6,cut)
    assert replay.stats()['full_hits']==1


def test_locked_mapping_islands_and_prior_boundaries_are_part_of_key():
    source=graph(8,[(i,i+1) for i in range(7)]+[(0,3)])
    target=graph(8,[(i,i+1) for i in range(7)]+[(0,3)])
    orbits=_nauty_orbits(target,wbo_tol=1.)
    replay=CutReplay(source,target,orbits)
    for mapped in ((0,1),(1,0)):
        mapping=[*mapped]+[-1]*6
        for islands in ([(0,1),(1,1)],[(0,1),(1,2)]):
            for deferred in ((),((1,2),)):
                for cut in ((),((5,6),)):
                    kw=dict(islands=islands,deferred=deferred)
                    assert raw(source,target,orbits,7,mapping,cut,replay=replay,**kw)==raw(source,target,orbits,7,mapping,cut,**kw)


def test_cache_budget_does_not_limit_search():
    source=graph(4,[(0,1),(1,2),(2,3)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    replay=CutReplay(source,source,orbits,cache_bytes=1)
    for cut in ((),((1,2),)):
        assert raw(source,source,orbits,0,[-1]*4,cut,replay=replay)==raw(source,source,orbits,0,[-1]*4,cut)
    assert replay.stats()['entries']==0


def test_invalid_cut_does_not_change_session_topology():
    source=graph(4,[(0,1),(1,2),(2,3)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    replay=CutReplay(source,source,orbits)
    for _ in range(2):
        with pytest.raises(ValueError,match='not an active edge'):
            raw(source,source,orbits,0,[-1]*4,((0,1),(0,3)),replay=replay)
    for cut in ((),((1,2),)):
        assert raw(source,source,orbits,0,[-1]*4,cut,replay=replay)==raw(source,source,orbits,0,[-1]*4,cut)


def test_full_search_graph_is_identical():
    source=graph(8,[(i,i+1) for i in range(7)]+[(0,3)])
    target=graph(8,[(i,i+1) for i in range(7)]+[(2,5)])
    orbits=_nauty_orbits(target,wbo_tol=1.)
    replay=CutReplay(source,target,orbits)
    for cut in ((),((3,4),),((0,1),),((0,3),)):
        view=replay.for_cut(cut)
        options=dict(cuts=cut,max_branches=100,p_orbits=orbits)
        expected=find_islands(view.source,target,list(range(8)),**options)
        actual=find_islands(view.source,target,list(range(8)),growth_replay=view,**options)
        assert actual==expected


def test_shared_finalization_preserves_all_generators_and_reports_reuse():
    from rxn_core.search_symmetry import finalize_graph_symmetry, SymmetryWorkspace
    source=graph(8,[(i,i+1) for i in range(7)]+[(0,3)])
    target=graph(8,[(i,i+1) for i in range(7)]+[(2,5)])
    workspace=SymmetryWorkspace(target,1.)
    orbits=_nauty_orbits(target,wbo_tol=1.)
    for cut in ((),((3,4),),((0,1),),((0,3),),()):
        r=source.copy();r.remove_edges_from(cut)
        raw=find_islands(r,target,list(range(8)),cuts=cut,max_branches=100,p_orbits=orbits)
        expected,_=finalize_graph_symmetry(raw,target,iso_tolerance=1.)
        actual,metrics=finalize_graph_symmetry(raw,target,iso_tolerance=1.,workspace=workspace)
        assert actual==expected
    assert metrics['completed_candidate_group_calculations']==0
