from dataclasses import replace

import networkx as nx
import numpy as np
import pytest

from rxn_core import AAMSearchConfig
from rxn_core.aam import cut_seed
from rxn_core.adaptive_cut_search import AdaptiveCutSearch
from rxn_core.alignment.branch import _generate_seed_orders
from rxn_core.domain import AAMProblem, MolecularEndpoint
from rxn_core.native_search import find_islands_native
from rxn_core.search_graph import frozen_value
from rxn_core.shared_seed_search import SharedSeedSearch


def state_key(state):
    return state.mapping, state.islands, state.deferred_edges


def contents(graph):
    # Ignore duplicated histories, node IDs and seed-provenance multiplicity.
    # Retain exact correlated fragment records, constraints and stop reasons.
    states = {state_key(s) for s in graph.states}
    edges = {(state_key(graph.states[e.source]), state_key(graph.states[e.target]),
              frozen_value(e.match), e.preserved_bonds) for e in graph.transitions if e.match is not None}
    stops = {(state_key(graph.states[s.state]), s.reason, s.stage) for s in graph.stops}
    return states, edges, stops


@pytest.mark.parametrize('tol', [.5, 1.])
@pytest.mark.parametrize('cap', [1, 3, 100])
def test_exact_union_of_original_policies_including_frontier_caps(tol, cap):
    rng = np.random.default_rng(82197)
    for trial in range(6):
        endpoints = []
        for side in range(2):
            n = 7 + (side if trial % 2 else 0)
            elements = tuple(['C']*4 + ['O']*2 + ['H']*(n-6))
            w = np.triu(rng.choice([0., 0., 1., 1.5, 2.], size=(n,n)), 1)
            endpoints.append(MolecularEndpoint(elements, np.zeros((n,3)), w+w.T, str(side)))
        problem = AAMProblem(*endpoints)
        config = AAMSearchConfig(seed_count=3, branch_limit=cap, iso_tolerance=tol)
        outer = AdaptiveCutSearch(problem, config, policy='shared')
        child = outer._session(trial % len(outer.cuts))
        while child.advance():
            pass
        actual = contents(child.builder.finish())
        reference = [set(), set(), set()]
        for order in child.orders:
            graph = find_islands_native(child.source, child.target, order,
                graph_floor=config.graph_floor, iso_tol=tol, max_branches=cap,
                p_orbits=child.condition.orbits, r_orbits=child.source_orbits,
                cuts=child.condition.cuts)
            for target, values in zip(reference, contents(graph)):
                target.update(values)
        assert actual == tuple(reference)
        assert nx.is_directed_acyclic_graph(nx.DiGraph((e.source,e.target)
            for e in child.builder.transitions))
        assert child.growth_calls == len(child.decisions)
        assert child.work == child.growth_calls + child.reused_states


def test_shared_growth_cache_and_saved_snapshot_preserve_groups():
    from test_fragment_choices import problem
    outer = AdaptiveCutSearch(problem(6), AAMSearchConfig(seed_count=3, iso_tolerance=1.), policy='shared')
    for _ in range(25):
        outer.advance()
    prefix = outer.snapshot()
    record = prefix.aam.graph.to_record()
    while outer.advance():
        pass
    outer.snapshot()
    assert prefix.aam.graph.to_record() == record
    for child in outer.sessions.values():
        plain = SharedSeedSearch(child.problem, child.config,
            condition=replace(child.condition, growth_replay=None))
        while plain.advance():
            pass
        assert plain.snapshot().aam.graph == child.snapshot().aam.graph
    assert sum(c.reused_states for c in outer.sessions.values()) > 0


@pytest.mark.parametrize('workers', [1, 2])
def test_cut_worker_backend_preserves_full_fragment_content_and_checkpoint_identity(tmp_path, workers):
    from rxn_core import search_aam
    from test_fragment_choices import problem
    config = AAMSearchConfig(seed_count=3, iso_tolerance=1., branch_limit=100)
    original = search_aam(problem(6), config, workers=workers, execution='reused_native')
    directory = tmp_path/'cuts'
    shared = search_aam(problem(6), config, workers=workers, execution='shared_policies',
                        intermediate_dir=directory, archive_format='checkpoint')
    # Conditioned symmetry is finalized too: compare complete records, not only scores.
    assert contents(shared.graph) == contents(original.graph)
    assert len(list(directory.glob('*.raw.pkl.gz'))) == shared.metrics.cut_count
    resumed = search_aam(problem(6), config, workers=workers, execution='shared_policies',
                         intermediate_dir=directory, archive_format='checkpoint', resume=True)
    assert resumed.graph == shared.graph
    with pytest.raises(ValueError, match='configuration differs'):
        search_aam(problem(6), config, workers=workers, execution='reused_native',
                   intermediate_dir=directory, archive_format='checkpoint', resume=True)
