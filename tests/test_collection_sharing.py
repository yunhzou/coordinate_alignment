from dataclasses import replace

from rxn_core import AAMSearchConfig, search_aam
from rxn_core.frag import build_graph
from rxn_core.search_graph import AAMSearchGraph
from rxn_core.search_symmetry import finalize_graph_symmetry
from test_fragment_choices import problem


def old_combine(graphs):
    contexts, roots, states, edges, stops = [], [], [], [], []
    for graph in graphs:
        c,s,t = len(contexts),len(states),len(edges)
        contexts.extend(graph.contexts)
        roots.extend(r+s for r in graph.roots)
        states.extend(replace(v,id=v.id+s,context=v.context+c) for v in graph.states)
        edges.extend(replace(e,id=e.id+t,source=e.source+s,target=e.target+s) for e in graph.transitions)
        stops.extend(replace(v,state=v.state+s) for v in graph.stops)
    return AAMSearchGraph(tuple(contexts),tuple(roots),tuple(states),tuple(edges),tuple(stops))


def test_direct_offset_construction_preserves_every_field():
    graphs = [search_aam(problem(n), AAMSearchConfig(seed_count=2,branch_limit=cap)).graph
              for n,cap in [(4,1),(6,100),(3,3)]]
    for values in ([],graphs[:1],graphs,graphs+graphs):
        assert AAMSearchGraph.combine(values) == old_combine(values)


def test_pre_finalized_immutable_group_is_only_interned_once():
    class ObservedTuple(tuple):
        iterations = 0
        def __iter__(self):
            self.iterations += 1
            return super().__iter__()
    p = problem(4)
    graph = search_aam(p, AAMSearchConfig(seed_count=1,cut_floor=10,iso_tolerance=1.)).graph
    edge = next(e for e in graph.transitions if e.match is not None)
    raw = ObservedTuple(edge.match['symmetry']['automorph_generators'])
    match = {**edge.match,'symmetry':{**edge.match['symmetry'],'automorph_generators':raw}}
    repeated = replace(graph, transitions=tuple(replace(edge,id=i,match=match) for i in range(100)))
    target = build_graph(p.product.elements,p.product.wbo,bond_cut=.2)
    finalized,_ = finalize_graph_symmetry(repeated,target,iso_tolerance=1.)
    assert raw.iterations == 1
    assert repeated == finalized
    assert repeated.transitions[0].match['symmetry']['automorph_generators'] is raw
