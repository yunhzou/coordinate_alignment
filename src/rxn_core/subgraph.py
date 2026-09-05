"""General weighted-subgraph matching API."""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Sequence

import networkx as nx

from .alignment.branch import find_islands
from .frag import WeightedGraph
from .matcher import _edge_wbo, _growth_edge_supported
from .matcher.policy import as_node_match_policy
from .search_graph import AAMSearchGraph, SearchPath, frozen_value


@dataclass(frozen=True)
class SubgraphMatch:
    """One query->target weighted-subgraph placement."""

    mapping: dict[int, int]
    symmetry_fragments: tuple = ()
    deferred_edges: tuple = ()
    search_path: SearchPath | None = None

    @property
    def query_nodes(self):
        return tuple(sorted(self.mapping))

    @property
    def target_nodes(self):
        return tuple(self.mapping[i] for i in self.query_nodes)


@dataclass(frozen=True)
class SubgraphSearchResult(Sequence):
    """Validated query placements and all search evidence, including caps."""
    matches: tuple[SubgraphMatch, ...]
    graph: AAMSearchGraph

    def __len__(self):
        return len(self.matches)

    def __getitem__(self, index):
        return self.matches[index]

    @property
    def capped(self):
        return self.graph.capped


def _coerce_graph(graph, graph_floor):
    if isinstance(graph, WeightedGraph):
        return graph.to_networkx(bond_cut=graph_floor)
    if isinstance(graph, nx.Graph):
        if "wbo_matrix" not in graph.graph:
            raise ValueError("NetworkX graph must carry graph['wbo_matrix']")
        if "bond_cut" not in graph.graph:
            graph = graph.copy()
            graph.graph["bond_cut"] = float(graph_floor)
        return graph
    if hasattr(graph, "to_networkx"):
        return graph.to_networkx(bond_cut=graph_floor)
    raise TypeError("graph must be WeightedGraph or a NetworkX graph")


def _query_edges_supported(g_q, g_t, mapping, iso_tol):
    graph_floor = float(g_t.graph.get("bond_cut", 0.2))
    for i, j in g_q.edges():
        if i not in mapping or j not in mapping:
            return False
        if not _growth_edge_supported(
                _edge_wbo(g_q, i, j),
                _edge_wbo(g_t, mapping[i], mapping[j]),
                iso_tol,
                graph_floor):
            return False
    return True


def match_weighted_subgraph(query, target, *,
                            node_policy=None,
                            node_match=None,
                            node_key=None,
                            anchor_map=None,
                            graph_floor=0.2,
                            iso_tol=1.0,
                            symmetry_wbo_tol=0.2,
                            seed_order=None,
                            orbit_dedup=True,
                            target_orbits=None,
                            max_branches=1_000_000):
    """Find query placements inside target with a replaceable node rule.

    The node rule decides which query/target node pairs are allowed.  The edge
    verifier then applies the existing WBO ``iso_tol`` logic to active query
    edges.  The default node rule is same element, preserving the old sub-iso
    behavior.  ``anchor_map`` gives exact query->target constraints that are
    preloaded as locked single-node islands and can still seed growth.
    """
    policy = as_node_match_policy(
        node_policy, node_match=node_match, node_key=node_key)
    g_q = _coerce_graph(query, graph_floor)
    g_t = _coerce_graph(target, graph_floor)
    query_nodes = set(g_q.nodes())
    anchor_map = {int(r): int(p) for r, p in dict(anchor_map or {}).items()}
    order = list(seed_order) if seed_order is not None else sorted(g_q.nodes())
    search_graph = find_islands(
        g_q,
        g_t,
        order,
        graph_floor=graph_floor,
        iso_tol=iso_tol,
        symmetry_wbo_tol=symmetry_wbo_tol,
        max_branches=max_branches,
        orbit_dedup=orbit_dedup,
        core_R=tuple(sorted(query_nodes)),
        stop_when_core_mapped=True,
        p_orbits=target_orbits,
        node_policy=policy,
        anchor_map=anchor_map,
    )
    matches = []
    seen = set()
    for branch in search_graph.paths():
        mapping = {int(k): int(v) for k, v in branch.mapping.items()
                   if k in query_nodes}
        if set(mapping) != query_nodes:
            continue
        if not _query_edges_supported(g_q, g_t, mapping, iso_tol):
            continue
        key = (tuple(sorted(mapping.items())), frozen_value(branch.fragments))
        if key in seen:
            continue
        seen.add(key)
        matches.append(SubgraphMatch(
            mapping=mapping,
            symmetry_fragments=branch.fragments,
            deferred_edges=tuple(sorted(tuple(e) for e in branch.deferred_edges)),
            search_path=branch,
        ))
    return SubgraphSearchResult(tuple(matches), search_graph)
