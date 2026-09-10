import itertools
import networkx as nx

from rxn_core import AAMProblem, AAMSearchConfig
from rxn_core.adaptive_seed_search import AdaptiveSeedSearch
from test_fragment_choices import problem


def test_seed_alternatives_are_retained_and_no_seed_repeats_at_a_state():
    session=AdaptiveSeedSearch(problem())
    session.advance()
    first=session.snapshot()
    assert first.aam.graph.terminals and first.pending
    for _ in range(500):
        if not session.advance():break
    result=session.snapshot()
    assert result.closure_calls==0
    assert session.reused_states>0
    assert sum(len(v['tried']) for v in session.seed_trials.values())==session.growth_calls
    assert any(len(v['tried'])==7 for v in session.seed_trials.values())
    assert nx.is_directed_acyclic_graph(nx.DiGraph((e.source,e.target) for e in result.aam.graph.transitions))
    for path in itertools.islice(result.aam.graph.paths(),100):
        assert len(path.mapping)==7 and len(set(path.mapping.values()))==7


def test_partial_composition_eventually_stops_without_repeating_empty_growth():
    session=AdaptiveSeedSearch(AAMProblem(problem(5).reactant,problem(3).product),
                              AAMSearchConfig(seed_count=1,branch_limit=100))
    for _ in range(1000):
        if not session.advance():break
    result=session.snapshot()
    assert result.exhausted
    assert all(len(s.mapping)<=3 for s in result.aam.graph.states)
    assert all(len(t['tried'])<=5 for t in session.seed_trials.values())


def test_complete_normal_path_before_fairly_serving_deeper_seed_alternatives():
    from rxn_core.adaptive_seed_search import _SeedAgenda
    agenda=_SeedAgenda()
    agenda.push((1,1,0),'first alternative')
    agenda.push((0,0,0),'normal')
    assert agenda.pop()[2]=='normal'
    assert agenda.pop()[2]=='first alternative'
    agenda.push((1,1,0),'another shallow alternative')
    agenda.push((2,1,0),'deeper alternative')
    agenda.push((1,0,0),'normal continuation')
    assert agenda.pop()[2]=='normal continuation'
    assert agenda.pop()[2]=='deeper alternative'


def test_unexplored_regions_precede_redundant_region_seeds_without_dropping_them():
    from rxn_core.adaptive_seed_search import _SeedAgenda
    agenda=_SeedAgenda()
    agenda.push((1,1,1,0),'covered region')
    agenda.push((2,1,0,10),'new region')
    assert agenda.pop()[2]=='new region'
    assert agenda.pop()[2]=='covered region'


def test_completed_route_yields_before_all_normal_siblings_are_expanded():
    from rxn_core.adaptive_seed_search import _SeedAgenda
    agenda=_SeedAgenda()
    agenda.push((1,1,0,0),'alternative')
    for i in range(100):agenda.push((0,0,0,0),i)
    assert agenda.pop()[2]==99
    agenda.yield_to_alternatives=True
    assert agenda.pop()[2]=='alternative'
    assert len(agenda)==99


def test_known_region_seeds_are_not_starved_by_new_regions():
    from rxn_core.adaptive_seed_search import _SeedAgenda
    agenda=_SeedAgenda()
    for i in range(100):agenda.push((1,1,0,i),i)
    agenda.push((1,1,1,0),'known region')
    assert agenda.pop()[2]==0
    assert agenda.pop()[2]=='known region'
    assert len(agenda)==99
