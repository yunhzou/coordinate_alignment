import itertools
import random
import networkx as nx
import numpy as np
import pytest

from rxn_core import AAMProblem, AAMSearchConfig
from rxn_core.adaptive_search import AdaptiveFragmentSearch
from rxn_core.domain import MolecularEndpoint
from rxn_core.fragment import match_fragment, FragmentMatchConfig, FragmentMatchContext
from rxn_core.fragment_choices import FragmentChoices
from test_cut_replay import graph


@pytest.mark.parametrize('tol',[.5,1.])
@pytest.mark.parametrize('cap',[1,3,100])
def test_observer_preserves_normal_result_exactly(tol,cap):
    rng=np.random.default_rng(801019)
    for _ in range(8):
        source=graph(8,[(i,i+1) for i in range(5)]+[(0,6),(3,7)],['C']*6+['H']*2)
        target=graph(8,[(i,i+1) for i in range(5)]+[(2,6),(5,7)],['C']*6+['H']*2)
        seed=int(rng.integers(8))
        cfg=FragmentMatchConfig(iso_tolerance=tol,branch_limit=cap)
        expected=match_fragment(source,target,seed=seed,config=cfg)
        actual=FragmentChoices(source,target,seed=seed,config=cfg)
        assert actual.normal.matches==expected.matches
        assert actual.normal.capped==expected.capped
        assert actual.normal.branch_count==expected.branch_count
        for closure in actual.closures:
            for iso in actual.close(closure).matches:
                assert len(iso.fragment)==closure.size
                assert len(set(iso.values()))==len(iso)
                for a,b in iso.preserved_bonds:
                    assert abs(source[a][b]['wbo']-target.graph['wbo_matrix'][iso[a],iso[b]])<=tol+1e-12


def test_close_does_not_remove_bonds_or_lose_locked_mapping():
    source=graph(7,[(i,i+1) for i in range(6)])
    target=graph(7,[(i,i+1) for i in range(6)])
    edges=set(source.edges)
    choices=FragmentChoices(source,target,seed=1,
        context=FragmentMatchContext({0:0},{0:1}),config=FragmentMatchConfig(iso_tolerance=1.))
    for closure in choices.closures:
        for iso in choices.close(closure).matches:
            if 0 in iso:assert iso[0]==0
            assert set(source.edges)==edges
    assert choices.closures


def problem(size=7):
    w=np.zeros((size,size))
    for i in range(size-1):w[i,i+1]=w[i+1,i]=1.
    r=MolecularEndpoint(['C']*size,np.zeros((size,3)),w)
    p=MolecularEndpoint(['C']*size,np.zeros((size,3)),w.copy())
    return AAMProblem(r,p)


def test_sequential_agenda_resumes_keeps_full_graph_and_reports_pending():
    session=AdaptiveFragmentSearch(problem())
    session.advance()
    first=session.snapshot()
    assert first.pending and not first.exhausted
    assert first.aam.graph.terminals
    for _ in range(40):session.advance()
    result=session.snapshot()
    assert result.work>first.work
    assert set(first.aam.graph.terminals)<=set(result.aam.graph.terminals)
    assert len(result.aam.graph.roots)==1 and not result.aam.graph.contexts[0].cuts
    assert nx.is_directed_acyclic_graph(nx.DiGraph((e.source,e.target) for e in result.aam.graph.transitions))
    for path in itertools.islice(result.aam.graph.paths(),150):
        assert len(path.mapping)==7 and len(set(path.mapping.values()))==7
        realized=path.sample(random.Random(41))
        assert len({p for r,p in realized.mapping})==7


def test_partial_composition_and_visible_growth_caps():
    original=problem(7);other=problem(5)
    partial=AAMProblem(original.reactant,other.product)
    session=AdaptiveFragmentSearch(partial,AAMSearchConfig(seed_count=1,branch_limit=1))
    for _ in range(20):session.advance()
    result=session.snapshot()
    for path in itertools.islice(result.aam.graph.paths(),100):
        assert len(path.mapping)<=5 and len(path.mapping)==len(set(path.mapping.values()))
    assert result.aam.graph.capped or result.pending or result.aam.graph.terminals


def test_event_feedback_only_reorders_pending_choices():
    session=AdaptiveFragmentSearch(problem(),policy='event_guided')
    session.advance()
    assert session.best_score==0
    before={serial for _,serial,_ in session.agenda}
    session.best_score=None
    session._feedback({i:i for i in range(7)})
    assert {serial for _,serial,_ in session.agenda}==before
    for _ in range(30):session.advance()
    assert session.snapshot().aam.graph.terminals


def test_fair_agenda_does_not_starve_deeper_choices():
    from rxn_core.adaptive_search import _ChoiceAgenda
    agenda=_ChoiceAgenda(fair=True)
    for i in range(100):agenda.push((1,0,i),i)
    assert agenda.pop()[2]==0
    agenda.push((2,1,0),'deeper')
    assert agenda.pop()[2]=='deeper'
    assert len(agenda)==99


def test_fair_search_can_revise_more_than_one_fragment():
    session=AdaptiveFragmentSearch(problem(),policy='fair_depth')
    for _ in range(150):session.advance()
    assert session.agenda.served[2]>0
    g=session.snapshot().aam.graph
    assert nx.is_directed_acyclic_graph(nx.DiGraph((e.source,e.target) for e in g.transitions))


def test_experimental_unsupported_context_is_explicit():
    with pytest.raises(ValueError, match='anchored growth'):
        AdaptiveFragmentSearch(problem(), AAMSearchConfig(seed_count=1, anchors=((0,0),)))
    with pytest.raises(ValueError, match='external growth replay'):
        FragmentChoices(graph(3,[(0,1),(1,2)]), graph(3,[(0,1),(1,2)]), seed=0,
                        context=FragmentMatchContext(growth_replay=object()))
