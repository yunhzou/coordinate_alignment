from dataclasses import replace

import networkx as nx

from rxn_core import AAMSearchConfig
from rxn_core.adaptive_cut_search import AdaptiveCutSearch
from rxn_core.adaptive_seed_search import AdaptiveSeedSearch
from test_fragment_choices import problem


def test_cut_conditions_receive_work_before_one_condition_is_exhausted():
    session = AdaptiveCutSearch(problem(5))
    session.advance()
    first = session.snapshot()
    assert len(first.pending) > 0
    assert sum(p['kind'] == 'unvisited_cut' for p in first.pending) == len(session.cuts)-1
    for _ in range(len(session.cuts)-1):
        session.advance()
    assert set(session.sessions) == set(range(len(session.cuts)))
    assert all(child.work == 1 for child in session.sessions.values())
    result = session.snapshot()
    assert {context.cuts for context in result.aam.graph.contexts} == set(session.cuts)
    assert nx.is_directed_acyclic_graph(nx.DiGraph((e.source,e.target) for e in result.aam.graph.transitions))


def test_shared_repair_changes_neither_cut_condition_results_nor_correlations():
    session = AdaptiveCutSearch(problem(5), AAMSearchConfig(seed_count=1, branch_limit=100, iso_tolerance=1.))
    for _ in range(80):
        if not session.advance():
            break
    for child in session.sessions.values():
        # Same condition and decisions, but independently computed growth.
        plain = AdaptiveSeedSearch(session.problem, session.config,
            condition=replace(child.condition, growth_replay=None))
        for _ in range(child.work):
            assert plain.advance()
        assert plain.snapshot().aam.graph.to_record() == child.snapshot().aam.graph.to_record()
    assert session.repair.stats()['calls'] == session.work


def test_saved_prefixes_do_not_change_when_more_cut_conditions_are_explored():
    session = AdaptiveCutSearch(problem(4))
    for _ in range(8):
        session.advance()
    first = session.snapshot()
    record = first.aam.graph.to_record()
    for _ in range(50):
        if not session.advance():
            break
    session.snapshot()
    assert record == first.aam.graph.to_record()
