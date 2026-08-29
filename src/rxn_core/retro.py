"""Fragment-unit search primitives for single-step retrosynthesis."""
from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np

from .frag import WeightedGraph
from .growth import IslandBranchLimitExceeded, grow_island
from .matcher import _nauty_orbits
from .subgraph import _coerce_graph, match_weighted_subgraph


@dataclass(frozen=True)
class RetroFragmentSearchConfig:
    """Bounds for one precursor-to-target retained-fragment search."""

    graph_floor: float = 0.2
    iso_tolerance: float = 0.5
    minimum_fragment_size: int = 1
    branch_limit: int = 100
    candidate_limit: int = 512
    maximum_boundary_bonds: int | None = None
    maximum_leftover_fragments: int | None = None

    def __post_init__(self):
        if self.graph_floor <= 0 or self.iso_tolerance <= 0:
            raise ValueError("graph floor and isomorphism tolerance must be positive")
        if self.minimum_fragment_size < 1:
            raise ValueError("minimum fragment size must be positive")
        if self.branch_limit < 1 or self.candidate_limit < 1:
            raise ValueError("branch and candidate limits must be positive")


@dataclass(frozen=True)
class RetroFragmentCandidate:
    """One retained precursor fragment placed inside the target."""

    precursor_id: str
    mapping: tuple[tuple[int, int], ...]
    retained_atoms: tuple[int, ...]
    covered_target_atoms: tuple[int, ...]
    leftover_fragments: tuple[tuple[int, ...], ...]
    boundary_bonds: tuple[tuple[int, int], ...]
    attachment_atoms_R: tuple[int, ...]
    attachment_atoms_P: tuple[int, ...]
    augmented_anchors: tuple[tuple[int, int], ...]
    augmented_target_atom_count: int
    retained_fragments: tuple[tuple[int, ...], ...] = ()

    @property
    def atom_mapping(self):
        return dict(self.mapping)

    @property
    def retained_size(self):
        return len(self.retained_atoms)


@dataclass(frozen=True)
class RetroFragmentSearchResult:
    """Candidates and completeness information for one precursor."""

    precursor_id: str
    candidates: tuple[RetroFragmentCandidate, ...]
    status: str
    complete: bool
    branch_limit: int
    maximum_branch_count: int
    capped_seed_count: int
    best_fragment_size: int


@dataclass(frozen=True)
class RetroAssembly:
    """A non-overlapping fragment-unit cover of the target."""

    candidates: tuple[RetroFragmentCandidate, ...]
    formed_bonds: tuple[tuple[int, int], ...]
    broken_bonds: tuple[tuple[str, int, int], ...]

    @property
    def precursor_ids(self):
        return tuple(candidate.precursor_id for candidate in self.candidates)


@dataclass(frozen=True)
class RetroAssemblySearchResult:
    assemblies: tuple[RetroAssembly, ...]
    status: str
    complete: bool
    assembly_limit: int


def _matrix(graph):
    return np.asarray(graph.graph["wbo_matrix"], dtype=float)


def _weighted_graph_from_nx(graph, matrix):
    nodes = [dict(graph.nodes[index]) for index in sorted(graph.nodes())]
    return WeightedGraph(nodes, np.asarray(matrix, dtype=float))


def _candidate_parts(graph_R, retained):
    retained = set(retained)
    outside = set(graph_R) - retained
    boundary = tuple(sorted(
        tuple(sorted((int(left), int(right))))
        for left, right in graph_R.edges()
        if (left in retained) != (right in retained)
    ))
    fragments = tuple(sorted(
        (tuple(sorted(map(int, component)))
         for component in nx.connected_components(graph_R.subgraph(outside))),
        key=lambda component: (component[0], len(component)),
    )) if outside else ()
    return outside, boundary, fragments


def _augmented_validation(graph_R, graph_P, retained_mapping,
                          outside, boundary, *, graph_floor, iso_tolerance,
                          max_branches):
    """Place residual fragments competitively in P or copied fallbacks.

    Only the initially retained fragment is anchored.  Residual components
    are copied into the augmented target as fallback locations, but remain
    free to occupy unused atoms of the original target when they fit there.
    """
    r_matrix = np.array(_matrix(graph_R), copy=True)
    for left, right in boundary:
        r_matrix[left, right] = r_matrix[right, left] = 0.0
    query = _weighted_graph_from_nx(graph_R, r_matrix)

    p_count = len(graph_P)
    outside_order = tuple(sorted(outside))
    copied_index = {
        atom: p_count + offset for offset, atom in enumerate(outside_order)
    }
    target_nodes = [dict(graph_P.nodes[index]) for index in sorted(graph_P)]
    target_nodes.extend(dict(graph_R.nodes[atom]) for atom in outside_order)
    target_matrix = np.zeros((p_count + len(outside_order),) * 2, dtype=float)
    target_matrix[:p_count, :p_count] = _matrix(graph_P)
    r_full = _matrix(graph_R)
    for left in outside_order:
        for right in outside_order:
            target_matrix[copied_index[left], copied_index[right]] = r_full[left, right]
    target = WeightedGraph(target_nodes, target_matrix)

    fixed_mapping = dict(retained_mapping)
    try:
        matches = match_weighted_subgraph(
            query,
            target,
            anchor_map=fixed_mapping,
            graph_floor=graph_floor,
            iso_tol=iso_tolerance,
            max_branches=max_branches,
        )
    except IslandBranchLimitExceeded as exc:
        return (), True, int(exc.count), len(target_nodes)
    mappings = tuple(
        tuple(sorted((int(source), int(image))
                     for source, image in match.mapping.items()))
        for match in matches
        if all(match.mapping.get(source) == image
               for source, image in fixed_mapping.items())
    )
    p_matrix = _matrix(graph_P)

    def preserves_cut_attachments(mapping):
        placement = dict(mapping)
        return all(
            not (placement[left] < p_count and placement[right] < p_count)
            or p_matrix[placement[left], placement[right]] >= graph_floor
            for left, right in boundary
        )

    mappings = tuple(filter(preserves_cut_attachments, mappings))
    if mappings:
        maximum_target_ownership = max(
            sum(image < p_count for _source, image in mapping)
            for mapping in mappings
        )
        mappings = tuple(
            mapping for mapping in mappings
            if sum(image < p_count for _source, image in mapping)
            == maximum_target_ownership
        )
    return mappings, False, len(matches), len(target_nodes)


def _competitive_parts(graph_R, target_mapping, cut_boundary,
                       augmented_mapping, target_atom_count):
    """Classify target-owned fragments and fallback-owned spectators."""
    retained = set(target_mapping)
    spectator = set(graph_R) - retained
    cut_graph = graph_R.copy()
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
             cut_graph.subgraph(spectator))),
        key=lambda component: (component[0], len(component)),
    )) if spectator else ()
    fallback_anchors = tuple(sorted(
        (int(source), int(image))
        for source, image in augmented_mapping.items()
        if int(image) >= target_atom_count
    ))
    attachments = tuple(sorted({
        int(atom)
        for edge in cut_boundary for atom in edge if atom in retained
    }))
    return (retained_fragments, leftover_fragments, fallback_anchors,
            attachments)


def discover_retained_fragments(
        precursor, target, *, precursor_id="",
        config: RetroFragmentSearchConfig | None = None):
    """Find maximum connected precursor fragments that can exist in target.

    Every precursor seed first grows one connected island.  Maximum-size
    placements are locked in the target, their boundary bonds are cut, and
    residual components compete for unused target atoms or copied fallback
    components in an augmented product.
    """
    config = config or RetroFragmentSearchConfig()
    graph_R = _coerce_graph(precursor, config.graph_floor)
    graph_P = _coerce_graph(target, config.graph_floor)
    r_orbits = _nauty_orbits(graph_R, wbo_tol=config.iso_tolerance)
    p_orbits = _nauty_orbits(graph_P, wbo_tol=config.iso_tolerance)
    raw = {}
    capped_seeds = 0
    maximum_branch_count = 0
    candidate_capped = False

    seed_order = sorted(
        graph_R,
        key=lambda atom: (-graph_R.degree(atom),
                          str(graph_R.nodes[atom].get("element")), int(atom)),
    )
    for seed in seed_order:
        try:
            placements = grow_island(
                graph_R,
                graph_P,
                seed,
                {},
                graph_floor=config.graph_floor,
                iso_tol=config.iso_tolerance,
                min_lock_size=config.minimum_fragment_size,
                max_branches=config.branch_limit,
                p_orbits=p_orbits,
                r_orbits=r_orbits,
            )
        except IslandBranchLimitExceeded as exc:
            capped_seeds += 1
            maximum_branch_count = max(maximum_branch_count, exc.count)
            continue
        maximum_branch_count = max(maximum_branch_count, len(placements))
        for placement in placements:
            retained = tuple(sorted(map(int, placement.fragment)))
            mapping = tuple(sorted(
                (int(source), int(image))
                for source, image in placement.items()
                if source in placement.fragment
            ))
            raw.setdefault((retained, mapping), (retained, mapping))
            if len(raw) >= config.candidate_limit:
                candidate_capped = True
                break
        if candidate_capped:
            break

    best_size = max((len(retained) for retained, _mapping in raw.values()),
                    default=0)
    candidates = []
    seen_units = set()
    for retained, mapping_pairs in raw.values():
        if len(retained) != best_size:
            continue
        outside, boundary, _fragments = _candidate_parts(graph_R, retained)
        if (config.maximum_boundary_bonds is not None and
                len(boundary) > config.maximum_boundary_bonds):
            continue
        initial_mapping = dict(mapping_pairs)
        (augmented_mappings, augmented_capped, augmented_branch_count,
         augmented_count) = _augmented_validation(
            graph_R, graph_P, initial_mapping, outside, boundary,
            graph_floor=config.graph_floor,
            iso_tolerance=config.iso_tolerance,
            max_branches=config.branch_limit,
        )
        maximum_branch_count = max(
            maximum_branch_count, augmented_branch_count)
        if augmented_capped:
            capped_seeds += 1
        if not augmented_mappings:
            continue
        for augmented_pairs in augmented_mappings:
            augmented_mapping = dict(augmented_pairs)
            target_mapping = {
                source: image for source, image in augmented_pairs
                if image < len(graph_P)
            }
            retained_all = tuple(sorted(target_mapping))
            (retained_fragments, leftover_fragments, anchors,
             attachment_R) = _competitive_parts(
                graph_R, target_mapping, boundary, augmented_mapping,
                len(graph_P))
            if (config.maximum_leftover_fragments is not None and
                    len(leftover_fragments)
                    > config.maximum_leftover_fragments):
                continue
            attachment_R = tuple(sorted(set(attachment_R) | {
                atom for atom in retained_all
                if int((graph_R.nodes[atom].get("features") or {}).get(
                    "formal_charge", 0) or 0) != 0
            }))
            attachment_P = tuple(sorted({
                target_mapping[atom] for atom in attachment_R
                if atom in target_mapping
            }))
            target_pairs = tuple(sorted(target_mapping.items()))
            unit_key = (
                retained_all,
                tuple(sorted(target_mapping.values())),
                retained_fragments,
                leftover_fragments,
                boundary,
                attachment_P,
            )
            if unit_key in seen_units:
                continue
            seen_units.add(unit_key)
            candidates.append(RetroFragmentCandidate(
                precursor_id=str(precursor_id),
                mapping=target_pairs,
                retained_atoms=retained_all,
                covered_target_atoms=tuple(sorted(target_mapping.values())),
                leftover_fragments=leftover_fragments,
                boundary_bonds=boundary,
                attachment_atoms_R=attachment_R,
                attachment_atoms_P=attachment_P,
                augmented_anchors=anchors,
                augmented_target_atom_count=augmented_count,
                retained_fragments=retained_fragments,
            ))
            if len(candidates) >= config.candidate_limit:
                candidate_capped = True
                break
        if len(candidates) >= config.candidate_limit:
            break

    candidates.sort(key=lambda candidate: (
        candidate.covered_target_atoms,
        candidate.attachment_atoms_P,
        candidate.mapping,
    ))
    best_size = max((candidate.retained_size for candidate in candidates),
                    default=best_size)
    incomplete = bool(capped_seeds or candidate_capped)
    if candidates:
        status = "capped" if incomplete else "matched"
    elif incomplete:
        status = "capped"
    else:
        status = "no_match"
    return RetroFragmentSearchResult(
        precursor_id=str(precursor_id),
        candidates=tuple(candidates),
        status=status,
        complete=not incomplete,
        branch_limit=config.branch_limit,
        maximum_branch_count=maximum_branch_count,
        capped_seed_count=capped_seeds + int(candidate_capped),
        best_fragment_size=best_size,
    )


def assemble_fragment_cover(
        target, candidates, *, maximum_precursors=2,
        assembly_limit=1_000, require_attachment_bonds=False,
        allow_repeated_precursors=True):
    """Combine fragment candidates into non-overlapping complete target covers.

    Attachment compatibility is opt-in because a coarse or tautomeric AAM can
    place the correct retained fragment while misidentifying its attachment
    atom.  Strict attachment validation belongs in finalist refinement.
    """
    if maximum_precursors < 1 or assembly_limit < 1:
        raise ValueError("precursor and assembly limits must be positive")
    graph_P = _coerce_graph(target, 0.2)
    target_atoms = frozenset(map(int, graph_P.nodes()))
    ordered = sorted(candidates, key=lambda candidate: (
        -candidate.retained_size,
        candidate.precursor_id,
        candidate.covered_target_atoms,
        candidate.mapping,
    ))
    assemblies = []
    capped = False

    def emit(selected):
        owner = {}
        for index, candidate in enumerate(selected):
            for atom in candidate.covered_target_atoms:
                owner[atom] = index
        formed = []
        for left, right in graph_P.edges():
            if owner[left] == owner[right]:
                continue
            if require_attachment_bonds:
                left_candidate = selected[owner[left]]
                right_candidate = selected[owner[right]]
                if (left not in left_candidate.attachment_atoms_P or
                        right not in right_candidate.attachment_atoms_P):
                    return
            formed.append(tuple(sorted((int(left), int(right)))))
        broken = tuple(sorted(
            (candidate.precursor_id, int(left), int(right))
            for candidate in selected
            for left, right in candidate.boundary_bonds
        ))
        assemblies.append(RetroAssembly(
            candidates=tuple(selected),
            formed_bonds=tuple(sorted(formed)),
            broken_bonds=broken,
        ))

    def visit(start, selected, covered, used_precursors):
        nonlocal capped
        if capped:
            return
        if covered == target_atoms:
            emit(selected)
            if len(assemblies) >= assembly_limit:
                capped = True
            return
        if len(selected) >= maximum_precursors:
            return
        for index in range(start, len(ordered)):
            candidate = ordered[index]
            coverage = frozenset(candidate.covered_target_atoms)
            if coverage & covered:
                continue
            if (not allow_repeated_precursors and
                    candidate.precursor_id in used_precursors):
                continue
            visit(
                index + 1,
                selected + [candidate],
                covered | coverage,
                used_precursors | {candidate.precursor_id},
            )

    visit(0, [], frozenset(), set())
    assemblies.sort(key=lambda assembly: (
        len(assembly.candidates),
        len(assembly.formed_bonds),
        len(assembly.broken_bonds),
        assembly.precursor_ids,
    ))
    return RetroAssemblySearchResult(
        assemblies=tuple(assemblies),
        status="capped" if capped else ("matched" if assemblies else "no_cover"),
        complete=not capped,
        assembly_limit=int(assembly_limit),
    )
