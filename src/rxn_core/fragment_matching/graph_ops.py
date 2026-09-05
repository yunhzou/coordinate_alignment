"""Graph operations used by augmented fragment detection."""
from __future__ import annotations

import networkx as nx
import numpy as np

from ..frag import WeightedGraph


def weight_matrix(graph):
    return np.asarray(graph.graph["wbo_matrix"], dtype=float)


def weighted_graph_from_nx(graph, matrix):
    nodes = [dict(graph.nodes[index]) for index in sorted(graph.nodes())]
    return WeightedGraph(nodes, np.asarray(matrix, dtype=float))


def partition_at_retained_fragment(graph, retained_atoms):
    retained = set(retained_atoms)
    outside = set(graph) - retained
    boundary = tuple(sorted(
        tuple(sorted((int(left), int(right))))
        for left, right in graph.edges()
        if (left in retained) != (right in retained)
    ))
    fragments = tuple(sorted(
        (tuple(sorted(map(int, component)))
         for component in nx.connected_components(graph.subgraph(outside))),
        key=lambda component: (component[0], len(component)),
    )) if outside else ()
    return outside, boundary, fragments


def fragment_equivalence_classes(source, cuts, fragments, tolerance):
    """Prove interchangeable connected fragment units by exact source orbits."""
    from ..matcher import _nauty_orbits
    cut_graph = source.copy()
    cut_graph.remove_edges_from(cuts)
    matrix = weight_matrix(source).copy()
    for a, b in cuts:
        matrix[a, b] = matrix[b, a] = 0
    cut_graph.graph["wbo_matrix"] = matrix
    orbits = _nauty_orbits(cut_graph, wbo_tol=tolerance)
    signatures = [tuple(sorted(orbits[a] for a in fragment)) for fragment in fragments]
    classes = {signature: i for i, signature in enumerate(sorted(set(signatures)))}
    return tuple(classes[s] for s in signatures)
