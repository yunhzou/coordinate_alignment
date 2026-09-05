"""Competitive residual matching through the shared AAM search engine."""
from dataclasses import dataclass
import networkx as nx
import numpy as np

from ..alignment.branch import find_islands
from ..alignment.post_aam import AAMHierarchy, AAMHierarchyChain
from ..search_symmetry import finalize_graph_symmetry
from .graph_ops import weight_matrix, fragment_equivalence_classes


@dataclass(frozen=True)
class AugmentedFragmentPlacement:
    mapping: tuple[tuple[int, int], ...]
    hierarchy: AAMHierarchyChain
    search_paths: tuple = ()
    target_action: tuple = ()


@dataclass(frozen=True)
class AugmentedMatchResult:
    placements: tuple[AugmentedFragmentPlacement, ...]
    capped: bool
    maximum_branch_count: int
    augmented_target_atom_count: int
    search_graphs: tuple


def match_augmented_residuals(source, target, retained_mapping, outside, boundary, *,
        graph_floor, iso_tolerance, branch_limit, target_region_atoms=None,
        retained_symmetry=None):
    """Lock the initial fragment and search R against P plus residual copies.

    Copies are ordinary competing graph components. No attachment-neighbour
    rule, singleton special case, ownership-maximization, or manual product of
    component witnesses participates in matching.
    """
    from .models import FragmentCandidate, FragmentDerivation
    from .symmetry import materialize_target_coverage_orbit
    outside = tuple(sorted(outside))
    target_count = len(target)
    augmented_count = target_count + len(outside)
    augmented = target.copy()
    matrix = np.zeros((augmented_count, augmented_count))
    matrix[:target_count, :target_count] = weight_matrix(target)
    copied = {atom: target_count + i for i, atom in enumerate(outside)}
    for atom, image in copied.items():
        augmented.add_node(image, **dict(source.nodes[atom]))
    for a, b, attributes in source.edges(data=True):
        if a in copied and b in copied:
            augmented.add_edge(copied[a], copied[b], **dict(attributes))
            matrix[copied[a], copied[b]] = matrix[copied[b], copied[a]] = weight_matrix(source)[a, b]
    augmented.graph["wbo_matrix"] = matrix
    augmented.graph["bond_cut"] = graph_floor
    cut_source = source.copy()
    source_matrix = weight_matrix(source).copy()
    for a, b in boundary:
        if cut_source.has_edge(a, b):
            cut_source.remove_edge(a, b)
        source_matrix[a, b] = source_matrix[b, a] = 0
    cut_source.graph["wbo_matrix"] = source_matrix
    # No sweep: one conditional continuation of the selected initial fragment.
    seeds = sorted(outside, key=lambda a: (-cut_source.degree(a), a))
    graph = find_islands(cut_source, augmented, seeds, graph_floor=graph_floor,
        iso_tol=iso_tolerance, symmetry_wbo_tol=iso_tolerance,
        max_branches=branch_limit, anchor_map=retained_mapping)
    graph, _metrics = finalize_graph_symmetry(graph, augmented, iso_tolerance=iso_tolerance)
    baseline = AAMHierarchy.from_record({"fragments": ({
        "fragment_index": 0, "island_idx": 0, "fragment": sorted(retained_mapping),
        "deferred_edges": boundary,
        "symmetry": dict(retained_symmetry or {
            "witness": dict(retained_mapping), "blocks": []}),
    },)})
    placements = []
    for path in graph.paths():
        if set(path.mapping) != set(source):
            continue
        fragments = (tuple(sorted(retained_mapping)),) + tuple(
            fragment.r_atoms for fragment in path.hierarchy.fragments)
        # The quotient tracks complete fragment image sets, including competitor
        # positions, before projecting ownership to P. Dropping competitor atoms
        # earlier would not commute with the recorded generator action.
        candidate = FragmentCandidate(
            source_id="", mapping=tuple(sorted(path.mapping.items())),
            retained_atoms=tuple(sorted(source)),
            covered_target_atoms=tuple(sorted(path.mapping.values())),
            leftover_fragments=(), boundary_bonds=boundary,
            attachment_atoms_source=(), attachment_atoms_target=(),
            copied_residual_placements=(), augmented_target_atom_count=augmented_count,
            retained_fragments=fragments, aam_hierarchy=path.hierarchy,
            derivations=(FragmentDerivation((), (path,)),),
            preserved_source_bonds=tuple(sorted(tuple(sorted((a, b))) for a, b in cut_source.edges()
                if tuple(sorted((a, b))) not in path.deferred_edges)),
            fragment_classes=fragment_equivalence_classes(source,
                tuple(boundary) + tuple(path.deferred_edges), fragments, iso_tolerance))
        for variant in materialize_target_coverage_orbit(candidate, augmented,
                iso_tolerance=iso_tolerance, generators=()):
            # Record the actual action without editing the saved search graph.
            placements.append(AugmentedFragmentPlacement(
                variant.mapping, AAMHierarchyChain((baseline, variant.aam_hierarchy)),
                (path,), variant.derivations[0].target_action))
    caps = [s for s in graph.stops if s.reason == "capped"]
    return AugmentedMatchResult(tuple(placements), bool(caps),
        max([len(graph.terminals)] + [s.count for s in caps]), augmented_count, (graph,))


def project_augmented_placement(source, target_mapping, cut_boundary, augmented_mapping,
                                 target_atom_count):
    retained = set(target_mapping)
    residual = set(source) - retained
    cut_graph = source.copy()
    cut_graph.remove_edges_from(cut_boundary)
    def parts(atoms):
        return tuple(sorted(tuple(sorted(component))
                            for component in nx.connected_components(cut_graph.subgraph(atoms))))
    return (
        parts(retained), parts(residual),
        tuple(sorted((int(a), int(b)) for a, b in augmented_mapping.items()
                     if b >= target_atom_count)),
        tuple(sorted({a for edge in cut_boundary for a in edge if a in retained})),
    )
