"""Competitive augmented matching and placement projection."""
from __future__ import annotations

import networkx as nx
import numpy as np

from ..frag import WeightedGraph
from ..growth import IslandBranchLimitExceeded
from ..subgraph import match_weighted_subgraph
from .graph_ops import weight_matrix, weighted_graph_from_nx


def match_augmented_residuals(
        source, target, retained_mapping, outside, boundary, *,
        graph_floor, iso_tolerance, branch_limit):
    source_matrix = np.array(weight_matrix(source), copy=True)
    for left, right in boundary:
        source_matrix[left, right] = source_matrix[right, left] = 0.0
    query = weighted_graph_from_nx(source, source_matrix)

    target_atom_count = len(target)
    outside_order = tuple(sorted(outside))
    copied_index = {
        atom: target_atom_count + offset
        for offset, atom in enumerate(outside_order)
    }
    target_nodes = [dict(target.nodes[index]) for index in sorted(target)]
    target_nodes.extend(dict(source.nodes[atom]) for atom in outside_order)
    augmented_size = target_atom_count + len(outside_order)
    augmented_matrix = np.zeros((augmented_size, augmented_size), dtype=float)
    augmented_matrix[:target_atom_count, :target_atom_count] = weight_matrix(
        target)
    full_source_matrix = weight_matrix(source)
    for left in outside_order:
        for right in outside_order:
            augmented_matrix[copied_index[left], copied_index[right]] = (
                full_source_matrix[left, right])
    augmented_target = WeightedGraph(target_nodes, augmented_matrix)

    fixed_mapping = dict(retained_mapping)
    try:
        matches = match_weighted_subgraph(
            query,
            augmented_target,
            anchor_map=fixed_mapping,
            graph_floor=graph_floor,
            iso_tol=iso_tolerance,
            max_branches=branch_limit,
        )
    except IslandBranchLimitExceeded as exc:
        return (), True, int(exc.count), augmented_size

    mappings = tuple(
        tuple(sorted(
            (int(source_atom), int(target_atom))
            for source_atom, target_atom in match.mapping.items()
        ))
        for match in matches
        if all(
            match.mapping.get(source_atom) == target_atom
            for source_atom, target_atom in fixed_mapping.items()
        )
    )
    target_matrix = weight_matrix(target)

    def preserves_cut_attachments(mapping):
        placement = dict(mapping)
        return all(
            not (placement[left] < target_atom_count
                 and placement[right] < target_atom_count)
            or target_matrix[placement[left], placement[right]] >= graph_floor
            for left, right in boundary
        )

    mappings = tuple(filter(preserves_cut_attachments, mappings))
    if mappings:
        maximum_target_ownership = max(
            sum(image < target_atom_count for _source, image in mapping)
            for mapping in mappings
        )
        mappings = tuple(
            mapping for mapping in mappings
            if sum(image < target_atom_count for _source, image in mapping)
            == maximum_target_ownership
        )
    return mappings, False, len(matches), augmented_size


def project_augmented_placement(
        source, target_mapping, cut_boundary, augmented_mapping,
        target_atom_count):
    retained = set(target_mapping)
    residual = set(source) - retained
    cut_graph = source.copy()
    for left, right in cut_boundary:
        if cut_graph.has_edge(left, right):
            cut_graph.remove_edge(left, right)

    retained_fragments = tuple(sorted(
        (tuple(sorted(map(int, component)))
         for component in nx.connected_components(
             cut_graph.subgraph(retained))),
        key=lambda component: (component[0], len(component)),
    )) if retained else ()
    leftover_fragments = tuple(sorted(
        (tuple(sorted(map(int, component)))
         for component in nx.connected_components(
             cut_graph.subgraph(residual))),
        key=lambda component: (component[0], len(component)),
    )) if residual else ()
    copied_residual_placements = tuple(sorted(
        (int(source_atom), int(target_atom))
        for source_atom, target_atom in augmented_mapping.items()
        if int(target_atom) >= target_atom_count
    ))
    attachment_atoms = tuple(sorted({
        int(atom)
        for edge in cut_boundary for atom in edge if atom in retained
    }))
    return (
        retained_fragments,
        leftover_fragments,
        copied_residual_placements,
        attachment_atoms,
    )
