"""Competitive augmented matching and placement projection."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import networkx as nx

from ..alignment.post_aam import AAMHierarchy
from ..matcher import _PartialMappingCanonicalizer, _nauty_orbits
from ..subgraph import match_weighted_subgraph
from .graph_ops import weight_matrix


@dataclass(frozen=True)
class AugmentedFragmentPlacement:
    mapping: tuple[tuple[int, int], ...]
    hierarchy: AAMHierarchy


def _hierarchy_with(base, extension):
    return AAMHierarchy(base.fragments + extension.fragments)


def _component_graph(source, atoms, graph_floor):
    graph = source.subgraph(atoms).copy()
    graph.graph["wbo_matrix"] = weight_matrix(source)
    graph.graph["bond_cut"] = float(graph_floor)
    return graph


def _composition_fits(component, target):
    source_counts = Counter(
        component.nodes[atom].get("element") for atom in component)
    target_counts = Counter(
        target.nodes[atom].get("element") for atom in target)
    return all(count <= target_counts[element]
               for element, count in source_counts.items())


def _component_placements(
        component, target, available_atoms, available_target,
        available_target_orbits, retained_mapping, boundary, *, graph_floor,
        iso_tolerance, branch_limit):
    if not _composition_fits(component, target):
        return ()
    if len(component) > len(available_atoms):
        return ()
    if len(component) == 1:
        source_atom = int(next(iter(component)))
        element = component.nodes[source_atom].get("element")
        attached_retained = []
        for left, right in boundary:
            if left == source_atom and right in retained_mapping:
                attached_retained.append(right)
            elif right == source_atom and left in retained_mapping:
                attached_retained.append(left)
        placements = []
        for target_atom in available_atoms:
            if target.nodes[target_atom].get("element") != element:
                continue
            if not all(target.has_edge(
                    target_atom, retained_mapping[retained_atom])
                       for retained_atom in attached_retained):
                continue
            hierarchy = AAMHierarchy.from_record({
                "fragments": ({
                    "fragment_index": 0,
                    "island_idx": 0,
                    "fragment": [source_atom],
                    "deferred_edges": [],
                    "symmetry": {
                        "witness": {source_atom: int(target_atom)},
                        "blocks": [],
                    },
                },),
            })
            placements.append(AugmentedFragmentPlacement(
                mapping=((source_atom, int(target_atom)),),
                hierarchy=hierarchy,
            ))
        return tuple(placements)
    seed_order = sorted(component, key=lambda atom: (
        -component.degree(atom),
        str(component.nodes[atom].get("element")),
        int(atom),
    ))
    matches = match_weighted_subgraph(
        component,
        available_target,
        graph_floor=graph_floor,
        iso_tol=iso_tolerance,
        max_branches=branch_limit,
        seed_order=seed_order,
        target_orbits=available_target_orbits,
    )

    component_atoms = set(component)
    relevant_boundary = tuple(
        (left, right) for left, right in boundary
        if left in component_atoms or right in component_atoms)

    def attachment_ok(mapping):
        for left, right in relevant_boundary:
            if left in component_atoms:
                component_atom, retained_atom = left, right
            else:
                component_atom, retained_atom = right, left
            if retained_atom not in retained_mapping:
                continue
            if not target.has_edge(
                    mapping[component_atom],
                    retained_mapping[retained_atom]):
                return False
        return True

    placements = []
    for match in matches:
        mapping = {int(source): int(image)
                   for source, image in match.mapping.items()}
        if not attachment_ok(mapping):
            continue
        placements.append(AugmentedFragmentPlacement(
            mapping=tuple(sorted(mapping.items())),
            hierarchy=AAMHierarchy.from_record({
                "fragments": match.symmetry_fragments,
            }),
        ))
    return tuple(placements)


def match_augmented_residuals(
        source, target, retained_mapping, outside, boundary, *,
        graph_floor, iso_tolerance, branch_limit, target_region_atoms=None,
        retained_symmetry=None):
    target_atom_count = len(target)
    fixed_mapping = dict(retained_mapping)
    baseline_hierarchy = AAMHierarchy.from_record({
        "fragments": ({
            "fragment_index": 0,
            "island_idx": 0,
            "fragment": sorted(map(int, retained_mapping)),
            "deferred_edges": [list(map(int, edge)) for edge in boundary],
            "symmetry": dict(retained_symmetry or {
                "witness": dict(retained_mapping),
                "blocks": [],
            }),
        },),
    })
    states = (AugmentedFragmentPlacement(
        mapping=tuple(sorted(fixed_mapping.items())),
        hierarchy=baseline_hierarchy,
    ),)
    canonicalizer = _PartialMappingCanonicalizer(
        source, target, wbo_tol=iso_tolerance)
    cut_graph = source.subgraph(outside).copy()
    retained_images = set(map(int, retained_mapping.values()))
    available_atoms = tuple(
        atom for atom in target if atom not in retained_images)
    available_target = target.subgraph(available_atoms).copy()
    available_target.graph["wbo_matrix"] = weight_matrix(target)
    available_target.graph["bond_cut"] = float(graph_floor)
    available_target_orbits = _nauty_orbits(
        available_target, wbo_tol=iso_tolerance)
    components = sorted(
        (tuple(sorted(map(int, component)))
         for component in nx.connected_components(cut_graph)),
        key=lambda atoms: (-len(atoms), atoms),
    )
    maximum_branch_count = 1
    competitive_atom_count = 0
    for atoms in components:
        component = _component_graph(source, atoms, graph_floor)
        options = _component_placements(
            component,
            target,
            available_atoms,
            available_target,
            available_target_orbits,
            fixed_mapping,
            boundary,
            graph_floor=graph_floor,
            iso_tolerance=iso_tolerance,
            branch_limit=branch_limit,
        )
        maximum_branch_count = max(maximum_branch_count, len(options))
        if not options:
            continue
        competitive_atom_count += len(atoms)
        next_states = {}
        for state in states:
            next_states.setdefault(
                canonicalizer.certificate(dict(state.mapping)), state)
            state_mapping = dict(state.mapping)
            used_target = set(state_mapping.values())
            for option in options:
                option_mapping = dict(option.mapping)
                if used_target.intersection(option_mapping.values()):
                    continue
                combined_mapping = dict(state_mapping)
                combined_mapping.update(option_mapping)
                combined = AugmentedFragmentPlacement(
                    mapping=tuple(sorted(combined_mapping.items())),
                    hierarchy=_hierarchy_with(
                        state.hierarchy, option.hierarchy),
                )
                next_states.setdefault(
                    canonicalizer.certificate(combined_mapping), combined)
        maximum_branch_count = max(
            maximum_branch_count, len(next_states))
        if len(next_states) > branch_limit:
            return (
                (), True, len(next_states),
                target_atom_count + competitive_atom_count,
            )
        states = tuple(next_states.values())

    placements = states
    augmented_size = target_atom_count + competitive_atom_count
    if placements:
        region = (frozenset(map(int, target_region_atoms))
                  if target_region_atoms is not None else None)

        def ownership(candidate):
            target_images = {
                image for _source, image in candidate.mapping
                if image < target_atom_count
            }
            if region is None:
                return (len(target_images),)
            return (
                len(target_images & region),
                -len(target_images - region),
            )

        best_ownership = max(map(ownership, placements))
        placements = tuple(
            placement for placement in placements
            if ownership(placement) == best_ownership
        )
    return placements, False, maximum_branch_count, augmented_size


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
