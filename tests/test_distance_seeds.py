import networkx as nx
import pytest
from rxn_core import AAMSearchConfig
from rxn_core.alignment.branch import _generate_seed_orders


def orders(graph, count=10, seed=42, policy='distance'):
    nx.set_node_attributes(graph, 'C', 'element')
    return _generate_seed_orders(graph, count, rng_seed=seed, seed_selection=policy)


def test_distance_orders_are_reproducible_complete_and_prefix_stable():
    graph = nx.path_graph(20)
    result = orders(graph)
    assert result == orders(graph)
    assert result[:3] == orders(graph, 3)
    assert len({order[0] for order in result}) == 10
    assert all(sorted(order) == list(graph) for order in result)
    assert result[0] == orders(graph, policy='random')[0]
    assert result != orders(graph, seed=43)


def test_distance_policy_increases_average_initial_separation():
    graph = nx.path_graph(40)
    totals = {}
    for policy in ('random', 'distance'):
        totals[policy] = sum(abs(a[0]-b[0]) for seed in range(200)
            for a,b in [orders(graph, 2, seed, policy)])
    assert totals['distance'] > totals['random'] * 1.15


def test_disconnected_isolated_and_exhausted_anchors():
    graph = nx.disjoint_union(nx.path_graph(5), nx.path_graph(5))
    graph.add_nodes_from([10,11])
    result = orders(graph, 25)
    assert len(result) == 25
    assert {order[0] for order in result[:10]} == set(range(10))
    assert all(set(order) == set(graph) for order in result)
    assert len(orders(nx.empty_graph(2), 5)) == 5
    assert orders(nx.Graph()) == []


def test_invalid_policy_is_rejected():
    with pytest.raises(ValueError):
        AAMSearchConfig(seed_selection='wrong')
    with pytest.raises(ValueError):
        orders(nx.path_graph(4), policy='wrong')


def test_checkpoint_identity_distinguishes_policy_without_breaking_random():
    from dataclasses import replace
    from rxn_core.aam import checkpoint_manifest
    from rxn_core import AAMProblem
    from test_aam_search_policy import endpoint
    problem = AAMProblem(endpoint(2), endpoint(2))
    config = AAMSearchConfig()
    random = checkpoint_manifest(problem, config)
    distance = checkpoint_manifest(problem, replace(config, seed_selection='distance'))
    assert 'seed_selection' not in random['config']
    assert distance['config']['seed_selection'] == 'distance'
    assert random != distance
