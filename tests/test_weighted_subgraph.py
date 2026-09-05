import numpy as np

from rxn_core import WeightedGraph, WeightedNode, match_weighted_subgraph


def _edge_matrix(n, edges):
    wbo = np.zeros((n, n))
    for i, j, w in edges:
        wbo[i, j] = wbo[j, i] = w
    return wbo


def test_weighted_subgraph_default_policy_is_same_element():
    query = WeightedGraph(
        [
            WeightedNode(element="O", features={"outer_shell": 6}),
            WeightedNode(element="W", features={"outer_shell": 6}),
        ],
        _edge_matrix(2, [(0, 1, 1.2)]),
    )
    target = WeightedGraph(
        [
            WeightedNode(element="S", features={"outer_shell": 6}),
            WeightedNode(element="Mo", features={"outer_shell": 6}),
        ],
        _edge_matrix(2, [(0, 1, 1.1)]),
    )

    assert not match_weighted_subgraph(
        query, target, graph_floor=0.2, iso_tol=0.3)


def test_weighted_subgraph_can_match_by_custom_node_feature():
    query = WeightedGraph(
        [
            WeightedNode(element="O", features={"outer_shell": 6}),
            WeightedNode(element="W", features={"outer_shell": 6}),
        ],
        _edge_matrix(2, [(0, 1, 1.2)]),
    )
    target = WeightedGraph(
        [
            WeightedNode(element="S", features={"outer_shell": 6}),
            WeightedNode(element="Mo", features={"outer_shell": 6}),
            WeightedNode(element="H", features={"outer_shell": 1}),
        ],
        _edge_matrix(3, [(0, 1, 1.1), (1, 2, 0.9)]),
    )

    matches = match_weighted_subgraph(
        query,
        target,
        node_policy="outer_shell",
        graph_floor=0.2,
        iso_tol=0.3,
        orbit_dedup=False,
    )

    assert [m.mapping for m in matches] == [{0: 0, 1: 1}, {0: 1, 1: 0}]


def test_weighted_subgraph_respects_wbo_iso_after_custom_node_match():
    query = WeightedGraph(
        [
            {"element": "O", "features": {"outer_shell": 6}},
            {"element": "W", "features": {"outer_shell": 6}},
        ],
        _edge_matrix(2, [(0, 1, 1.2)]),
    )
    target = WeightedGraph(
        [
            {"element": "S", "features": {"outer_shell": 6}},
            {"element": "Mo", "features": {"outer_shell": 6}},
        ],
        _edge_matrix(2, [(0, 1, 0.4)]),
    )

    assert not match_weighted_subgraph(
        query,
        target,
        node_policy="outer_shell",
        graph_floor=0.2,
        iso_tol=0.3,
    )


def test_weighted_subgraph_can_match_by_multiple_node_features():
    query = WeightedGraph(
        [
            {"features": {"outer_shell": 6, "nuclear_charge": 8}},
            {"features": {"outer_shell": 6, "nuclear_charge": 74}},
        ],
        _edge_matrix(2, [(0, 1, 1.0)]),
    )
    target = WeightedGraph(
        [
            {"features": {"outer_shell": 6, "nuclear_charge": 74}},
            {"features": {"outer_shell": 6, "nuclear_charge": 8}},
        ],
        _edge_matrix(2, [(0, 1, 1.0)]),
    )

    matches = match_weighted_subgraph(
        query,
        target,
        node_policy=("outer_shell", "nuclear_charge"),
        graph_floor=0.2,
        iso_tol=0.1,
        orbit_dedup=False,
    )

    assert [m.mapping for m in matches] == [{0: 1, 1: 0}]


def test_weighted_subgraph_anchor_seed_can_grow_from_locked_pair():
    query = WeightedGraph(
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        _edge_matrix(2, [(0, 1, 1.0)]),
    )
    target = WeightedGraph(
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        _edge_matrix(4, [(0, 1, 1.0), (2, 3, 1.0)]),
    )

    matches = match_weighted_subgraph(
        query,
        target,
        node_policy="outer_shell",
        anchor_map={0: 2},
        seed_order=[0],
        graph_floor=0.2,
        iso_tol=0.1,
        orbit_dedup=False,
    )

    assert [m.mapping for m in matches] == [{0: 2, 1: 3}]


def test_weighted_subgraph_anchor_map_filters_to_exact_target():
    query = WeightedGraph(
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        _edge_matrix(2, [(0, 1, 1.0)]),
    )
    target = WeightedGraph(
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        _edge_matrix(4, [(0, 1, 1.0), (2, 3, 1.0)]),
    )

    matches = match_weighted_subgraph(
        query,
        target,
        node_policy="outer_shell",
        anchor_map={0: 2, 1: 3},
        graph_floor=0.2,
        iso_tol=0.1,
        orbit_dedup=False,
    )

    assert [m.mapping for m in matches] == [{0: 2, 1: 3}]


def test_weighted_subgraph_complete_anchor_map_still_checks_wbo_edges():
    query = WeightedGraph(
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        _edge_matrix(2, [(0, 1, 1.0)]),
    )
    target = WeightedGraph(
        [
            {"features": {"outer_shell": 6}},
            {"features": {"outer_shell": 6}},
        ],
        _edge_matrix(2, [(0, 1, 0.3)]),
    )

    assert not match_weighted_subgraph(
        query,
        target,
        node_policy="outer_shell",
        anchor_map={0: 0, 1: 1},
        graph_floor=0.2,
        iso_tol=0.1,
    )


def test_weighted_subgraph_incompatible_anchor_has_no_match():
    query = WeightedGraph(["O"], _edge_matrix(1, []))
    target = WeightedGraph(["S"], _edge_matrix(1, []))

    assert not match_weighted_subgraph(
        query, target, anchor_map={0: 0})
