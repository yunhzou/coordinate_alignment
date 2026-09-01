"""Detect target-owned fragments from one source graph."""
from __future__ import annotations

from dataclasses import replace

from ..growth import IslandBranchLimitExceeded, grow_island
from ..matcher import (
    _atom_tuple_orbit,
    _nauty_atom_generators,
    _nauty_orbits,
)
from ..subgraph import _coerce_graph
from .augmentation import match_augmented_residuals, project_augmented_placement
from .graph_ops import partition_at_retained_fragment
from .models import (
    FragmentCandidate,
    FragmentDetectionConfig,
    FragmentDetectionResult,
    FragmentTargetContext,
)


def _initial_fragment_placements(
        source, target, config, *, target_orbits=None):
    source_orbits = _nauty_orbits(source, wbo_tol=config.iso_tolerance)
    if target_orbits is None:
        target_orbits = _nauty_orbits(
            target, wbo_tol=config.iso_tolerance)
    placements_by_identity = {}
    capped_seed_count = 0
    maximum_branch_count = 0
    candidate_capped = False

    seed_order = sorted(
        source,
        key=lambda atom: (
            -source.degree(atom),
            str(source.nodes[atom].get("element")),
            int(atom),
        ),
    )
    seed_limited = (
        config.seed_limit is not None
        and len(seed_order) > config.seed_limit)
    if config.seed_limit is not None:
        seed_order = seed_order[:config.seed_limit]
    for seed in seed_order:
        try:
            placements = grow_island(
                source,
                target,
                seed,
                {},
                graph_floor=config.graph_floor,
                iso_tol=config.iso_tolerance,
                min_lock_size=config.minimum_fragment_size,
                max_branches=config.branch_limit,
                p_orbits=target_orbits,
                r_orbits=source_orbits,
            )
        except IslandBranchLimitExceeded as exc:
            capped_seed_count += 1
            maximum_branch_count = max(maximum_branch_count, exc.count)
            continue

        maximum_branch_count = max(maximum_branch_count, len(placements))
        for placement in placements:
            retained = tuple(sorted(map(int, placement.fragment)))
            mapping = tuple(sorted(
                (int(source_atom), int(target_atom))
                for source_atom, target_atom in placement.items()
                if source_atom in placement.fragment
            ))
            placements_by_identity.setdefault(
                (retained, mapping), (retained, mapping))
            if len(placements_by_identity) >= config.candidate_limit:
                candidate_capped = True
                break
        if candidate_capped:
            break

    return (
        tuple(placements_by_identity.values()),
        capped_seed_count,
        maximum_branch_count,
        candidate_capped,
        seed_limited,
    )


def prepare_fragment_target(
        target, *, config: FragmentDetectionConfig | None = None):
    """Prepare target-owned symmetry data once for repeated source searches."""
    config = config or FragmentDetectionConfig()
    target_graph = _coerce_graph(target, config.graph_floor)
    return FragmentTargetContext(
        graph=target_graph,
        atom_orbits=_nauty_orbits(
            target_graph, wbo_tol=config.iso_tolerance),
        automorphism_generators=_nauty_atom_generators(
            target_graph, wbo_tol=config.iso_tolerance),
        graph_floor=config.graph_floor,
        iso_tolerance=config.iso_tolerance,
    )


def _candidate_identity(candidate):
    return (
        candidate.retained_atoms,
        candidate.covered_target_atoms,
        candidate.retained_fragments,
        candidate.leftover_fragments,
        candidate.boundary_bonds,
        candidate.attachment_atoms_target,
    )


def _target_automorphism_variants(candidate, generators):
    source_order = tuple(source for source, _target in candidate.mapping)
    target_images = tuple(target for _source, target in candidate.mapping)
    for images in _atom_tuple_orbit(target_images, generators):
        mapping = dict(zip(source_order, images))
        yield replace(
            candidate,
            mapping=tuple(sorted(mapping.items())),
            covered_target_atoms=tuple(sorted(mapping.values())),
            attachment_atoms_target=tuple(sorted({
                mapping[source]
                for source in candidate.attachment_atoms_source
                if source in mapping
            })),
        )


def detect_fragments(
        source, target, *, source_id="",
        config: FragmentDetectionConfig | None = None):
    """Generate augmented fragment candidates for one source-target pair."""
    config = config or FragmentDetectionConfig()
    source_graph = _coerce_graph(source, config.graph_floor)
    if isinstance(target, FragmentTargetContext):
        if (target.graph_floor != config.graph_floor
                or target.iso_tolerance != config.iso_tolerance):
            raise ValueError(
                "prepared target context and detection config disagree")
        target_context = target
    else:
        target_context = prepare_fragment_target(target, config=config)
    target_graph = target_context.graph
    (
        initial_placements,
        capped_seed_count,
        maximum_branch_count,
        candidate_capped,
        seed_limited,
    ) = _initial_fragment_placements(
        source_graph, target_graph, config,
        target_orbits=target_context.atom_orbits)

    best_initial_size = max(
        (len(retained) for retained, _mapping in initial_placements),
        default=0,
    )
    candidates = []
    seen_candidates = set()

    for retained, mapping_pairs in initial_placements:
        if len(retained) != best_initial_size:
            continue
        outside, boundary, _fragments = partition_at_retained_fragment(
            source_graph, retained)
        if (config.maximum_boundary_bonds is not None
                and len(boundary) > config.maximum_boundary_bonds):
            continue

        (
            augmented_mappings,
            augmented_capped,
            augmented_branch_count,
            augmented_atom_count,
        ) = match_augmented_residuals(
            source_graph,
            target_graph,
            dict(mapping_pairs),
            outside,
            boundary,
            graph_floor=config.graph_floor,
            iso_tolerance=config.iso_tolerance,
            branch_limit=config.branch_limit,
        )
        maximum_branch_count = max(
            maximum_branch_count, augmented_branch_count)
        if augmented_capped:
            capped_seed_count += 1
        if not augmented_mappings:
            continue

        for augmented_pairs in augmented_mappings:
            augmented_mapping = dict(augmented_pairs)
            target_mapping = {
                source_atom: target_atom
                for source_atom, target_atom in augmented_pairs
                if target_atom < len(target_graph)
            }
            retained_atoms = tuple(sorted(target_mapping))
            (
                retained_fragments,
                leftover_fragments,
                copied_residual_placements,
                attachment_atoms_source,
            ) = project_augmented_placement(
                source_graph,
                target_mapping,
                boundary,
                augmented_mapping,
                len(target_graph),
            )
            if (config.maximum_leftover_fragments is not None
                    and len(leftover_fragments)
                    > config.maximum_leftover_fragments):
                continue

            attachment_atoms_source = tuple(sorted(
                set(attachment_atoms_source) | {
                    atom for atom in retained_atoms
                    if int((source_graph.nodes[atom].get("features") or {}).get(
                        "formal_charge", 0) or 0) != 0
                }
            ))
            attachment_atoms_target = tuple(sorted({
                target_mapping[atom]
                for atom in attachment_atoms_source
                if atom in target_mapping
            }))
            candidate = FragmentCandidate(
                source_id=str(source_id),
                mapping=tuple(sorted(target_mapping.items())),
                retained_atoms=retained_atoms,
                covered_target_atoms=tuple(sorted(target_mapping.values())),
                leftover_fragments=leftover_fragments,
                boundary_bonds=boundary,
                attachment_atoms_source=attachment_atoms_source,
                attachment_atoms_target=attachment_atoms_target,
                copied_residual_placements=copied_residual_placements,
                augmented_target_atom_count=augmented_atom_count,
                retained_fragments=retained_fragments,
            )
            for variant in _target_automorphism_variants(
                    candidate, target_context.automorphism_generators):
                identity = _candidate_identity(variant)
                if identity in seen_candidates:
                    continue
                seen_candidates.add(identity)
                candidates.append(variant)
                if len(candidates) >= config.candidate_limit:
                    candidate_capped = True
                    break
            if candidate_capped:
                break
        if len(candidates) >= config.candidate_limit:
            break

    candidates.sort(key=lambda candidate: (
        candidate.covered_target_atoms,
        candidate.attachment_atoms_target,
        candidate.mapping,
    ))
    best_fragment_size = max(
        (candidate.retained_size for candidate in candidates),
        default=best_initial_size,
    )
    cap_hit = bool(capped_seed_count or candidate_capped)
    incomplete = bool(cap_hit or seed_limited)
    status = (
        "capped" if cap_hit
        else ("seed_limited" if seed_limited
              else ("matched" if candidates else "no_match"))
    )
    return FragmentDetectionResult(
        source_id=str(source_id),
        candidates=tuple(candidates),
        status=status,
        complete=not incomplete,
        branch_limit=config.branch_limit,
        maximum_branch_count=maximum_branch_count,
        capped_seed_count=capped_seed_count + int(candidate_capped),
        best_fragment_size=best_fragment_size,
    )
