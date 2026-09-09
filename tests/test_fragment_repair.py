"""Differential tests for dependency-directed conditional fragment repair."""
import numpy as np
import pytest

from rxn_core.alignment.branch import find_islands
from rxn_core.cut_replay import FragmentRepair
from rxn_core.growth import native
from rxn_core.matcher import _nauty_orbits
from rxn_core.search_symmetry import finalize_graph_symmetry
from test_cut_replay import graph, raw

pytestmark=pytest.mark.skipif(not native.built(),reason='native engine not built')


def test_reuse_later_fragment_after_independent_locked_mapping_changes():
    source=graph(10,[(0,1),(1,2),(2,3),(4,5),(6,7),(7,8),(8,9)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    repair=FragmentRepair(source,source,orbits)
    for mapping,islands in (([-1]*4+[4,5]+[-1]*4,[(4,1),(5,1)]),
                            ([-1]*4+[5,4]+[-1]*4,[(4,1),(5,2)])):
        assert raw(source,source,orbits,0,mapping,(),islands=islands,replay=repair)==raw(
            source,source,orbits,0,mapping,(),islands=islands)
    assert repair.stats()['history_hits']==1


def test_rebase_unrelated_previous_boundaries():
    source=graph(10,[(0,1),(1,2),(2,3),(4,5),(6,7),(7,8),(8,9)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    repair=FragmentRepair(source,source,orbits)
    for deferred in (((4,5),),((6,7),),()):
        assert raw(source,source,orbits,0,[-1]*10,(),deferred=deferred,replay=repair)==raw(
            source,source,orbits,0,[-1]*10,(),deferred=deferred)
    assert repair.stats()['boundary_rebases']==2


def test_changed_occupancy_cannot_hide_new_seed_alternatives():
    source=graph(6,[(0,1),(1,2),(3,4),(4,5)])
    orbits=_nauty_orbits(source,wbo_tol=1.)
    repair=FragmentRepair(source,source,orbits)
    expected=[]
    for mapping in ([-1]*3+[3,4,5],[-1]*6):
        fresh=raw(source,source,orbits,0,mapping,())
        assert raw(source,source,orbits,0,mapping,(),replay=repair)==fresh
        expected.append(fresh)
    assert expected[0]!=expected[1]
    assert repair.stats()['full_hits']==0


@pytest.mark.parametrize('tol',[.5,1.])
@pytest.mark.parametrize('budget',[1,64*1024*1024])
def test_changed_contexts_and_cuts_match_fresh_results(tol,budget):
    rng=np.random.default_rng(93852)
    for trial in range(5):
        elements=['C']*6+['H']*3
        edges=[(i,i+1,float(rng.choice([.8,1.2,1.8]))) for i in range(5)]
        edges += [(0,6,1.),(2,7,1.),(4,8,1.)]
        if trial%2:edges.append((0,4,1.))
        source=graph(9,edges,elements)
        target=graph(9,edges[:-1]+[(1,5,.8)],elements)
        orbits=_nauty_orbits(target,wbo_tol=tol)
        repair=FragmentRepair(source,target,orbits,cache_bytes=budget)
        for _ in range(40):
            mapping=[-1]*9
            locked=list(map(int,rng.choice(6,3,replace=False)))
            for a,b in zip(locked,rng.permutation([1,3,5])):mapping[a]=int(b)
            islands=[(r,1+int(rng.integers(2))) for r in locked]
            cut=() if rng.integers(3)==0 else (edges[int(rng.integers(len(edges)))][:2],)
            deferred=() if rng.integers(2)==0 else (edges[int(rng.integers(len(edges)))][:2],)
            seed=int(rng.integers(9));cap=int(rng.choice([1,3,100]))
            kw=dict(cap=cap,tol=tol,islands=islands,deferred=deferred)
            assert raw(source,target,orbits,seed,mapping,cut,replay=repair,**kw)==raw(
                source,target,orbits,seed,mapping,cut,**kw),(trial,mapping,islands,cut)


def test_complete_compressed_search_graphs_and_generators():
    source=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(0,6),(2,7),(4,8)],['C']*6+['H']*3)
    target=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(1,6),(3,7),(5,8)],['C']*6+['H']*3)
    orbits=_nauty_orbits(target,wbo_tol=1.)
    repair=FragmentRepair(source,target,orbits)
    for order in (list(range(9)),list(reversed(range(9))),[4,2,6,0,1,3,5,7,8]):
        for cut in ((),((2,3),),((0,6),),((0,1),)):
            view=repair.for_cut(cut)
            options=dict(cuts=cut,max_branches=100,p_orbits=orbits,iso_tol=1.)
            expected=find_islands(view.source,target,order,**options)
            actual=find_islands(view.source,target,order,growth_replay=view,**options)
            expected,_=finalize_graph_symmetry(expected,target,iso_tolerance=1.)
            actual,_=finalize_graph_symmetry(actual,target,iso_tolerance=1.)
            assert actual==expected
