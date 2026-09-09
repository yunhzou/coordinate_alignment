"""Native scheduling must preserve the complete Python decision DAG."""
import networkx as nx
import numpy as np
import pytest

from rxn_core.alignment.branch import find_islands
from rxn_core.growth import native
from rxn_core.cut_replay import FragmentRepair
from rxn_core.native_search import find_islands_native
from rxn_core.matcher import _nauty_orbits
from rxn_core.search_symmetry import finalize_graph_symmetry
from test_cut_replay import graph

pytestmark=pytest.mark.skipif(not native.built(),reason='native engine not built')


def compare(source,target,order,**options):
    orbits=_nauty_orbits(target,wbo_tol=options.get('iso_tol',1.))
    expected=find_islands(source,target,order,p_orbits=orbits,**options)
    actual=find_islands_native(source,target,order,p_orbits=orbits,**options)
    assert actual==expected
    return expected


@pytest.mark.parametrize('tol',[.5,1.])
@pytest.mark.parametrize('cap',[0,1,3,100])
def test_graphs_cuts_hydrogen_and_caps(tol,cap):
    rng=np.random.default_rng(7244)
    for trial in range(12):
        elements=['C']*6+['H']*3
        edges=[(i,i+1,float(rng.choice([.8,1.2,1.8]))) for i in range(5)]
        edges += [(0,6,1.),(2,7,1.),(4,8,1.)]
        source=graph(9,edges,elements)
        target=graph(9,edges[:-1]+[(1,5,.8)],elements)
        order=list(map(int,rng.permutation(9)))
        cut=() if trial%3==0 else (edges[int(rng.integers(len(edges)))][:2],)
        source.remove_edges_from(cut)
        compare(source,target,order,cuts=cut,iso_tol=tol,max_branches=cap)


def test_empty_order_duplicate_seeds_and_partial_composition():
    source=graph(5,[(0,1),(1,2),(2,3)],['C','C','O','H','Cl'])
    target=graph(4,[(0,1),(1,2)],['C','C','O','H'])
    for a,b in ((source,target),(target,source)):
        for order in ([],list(a),[0,0,1,0,1]):
            compare(a,b,order,max_branches=100)


def test_anchors_and_requested_core():
    source=graph(8,[(i,i+1) for i in range(7)],['C']*6+['H']*2)
    target=graph(8,[(i,i+1) for i in range(7)],['C']*6+['H']*2)
    for anchors in ({0:0},{0:5,5:0},{1:6}):
        for core in (None,[0,1],[0,1,2,3]):
            for cap in (1,100):
                compare(source,target,list(range(8)),anchor_map=anchors,core_R=core,
                        stop_when_core_mapped=True,max_branches=cap)


def test_nonidentity_atom_labels_and_final_groups():
    source=graph(6,[(0,1),(1,2),(2,3),(3,4),(4,5)])
    target=graph(6,[(0,1),(1,2),(2,3),(3,4),(4,5)])
    # The existing weighted-graph matrix uses physical integer atom indices.
    # Select/relabel a graph while retaining a correctly indexed full matrix.
    r_ids={i:2*i+1 for i in range(6)};p_ids={i:3*i+2 for i in range(6)}
    source=nx.relabel_nodes(source,r_ids);target=nx.relabel_nodes(target,p_ids)
    for g in (source,target):
        matrix=np.zeros((max(g)+1,max(g)+1))
        for a,b,data in g.edges(data=True):matrix[a,b]=matrix[b,a]=data['wbo']
        g.graph['wbo_matrix']=matrix
    expected=compare(source,target,list(source),max_branches=100)
    actual=find_islands_native(source,target,list(source),max_branches=100)
    expected,_=finalize_graph_symmetry(expected,target,iso_tolerance=1.)
    actual,_=finalize_graph_symmetry(actual,target,iso_tolerance=1.)
    assert actual==expected


def test_native_scheduling_with_dependency_repair():
    source=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(0,6),(2,7),(4,8)],['C']*6+['H']*3)
    target=graph(9,[(0,1),(1,2),(2,3),(3,4),(4,5),(1,6),(3,7),(5,8)],['C']*6+['H']*3)
    orbits=_nauty_orbits(target,wbo_tol=1.)
    repair=FragmentRepair(source,target,orbits)
    for order in (list(source),list(reversed(list(source)))):
        for cut in ((),((0,1),),((2,7),)):
            view=repair.for_cut(cut)
            for cap in (1,3,100):
                options=dict(p_orbits=orbits,cuts=cut,max_branches=cap)
                expected=find_islands(view.source,target,order,**options)
                actual=find_islands_native(view.source,target,order,growth_replay=view,**options)
                assert actual==expected


def test_random_graphs_exercise_both_cap_stages():
    rng=np.random.default_rng(1719)
    cap_stages=set()
    for _ in range(120):
        matrices=[]
        for side in range(2):
            w=np.triu(rng.choice([0.,0.,0.,1.,1.5,2.],size=(8,8)),1)
            matrices.append(w+w.T)
        from rxn_core.frag import build_graph
        source,target=[build_graph(['C']*4+['O']*2+['H']*2,w,bond_cut=.2) for w in matrices]
        result=compare(source,target,list(map(int,rng.permutation(8))),max_branches=int(rng.choice([3,4,8,100])))
        cap_stages.update(stop.stage for stop in result.stops if stop.reason=='capped')
    assert cap_stages=={'fragment_growth','combined_live_leaves'}
