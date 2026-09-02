"""Detect target-owned fragments from one source graph."""
from __future__ import annotations

from dataclasses import dataclass, replace

from ..alignment.branch import _generate_seed_orders
from ..growth import IslandBranchLimitExceeded, grow_island
from ..matcher import (
    _PartialMappingCanonicalizer,
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


@dataclass(frozen=True)
class _InitialFragmentFamily:
    """One exact endpoint-automorphism family from AAM island growth."""

    retained_atoms: tuple[int, ...]
    representative_mapping: tuple[tuple[int, int], ...]
    symmetry: dict
    encounter_count: int = 1


def _initial_fragment_placements(
        source, target, config, *, target_orbits=None,
        target_region_atoms=None):
    source_orbits = _nauty_orbits(source, wbo_tol=config.iso_tolerance)
    if target_orbits is None:
        target_orbits = _nauty_orbits(
            target, wbo_tol=config.iso_tolerance)
    placement_families = {}
    canonicalizer = _PartialMappingCanonicalizer(
        source,
        target,
        wbo_tol=config.iso_tolerance,
        target_atom_tags=(
            {int(atom): 'requested_region' for atom in target_region_atoms}
            if target_region_atoms is not None else None),
    )
    capped_seed_count = 0
    maximum_branch_count = 0
    candidate_capped = False

    if config.seed_mode == "fragment_cover":
        seed_order = _generate_seed_orders(source, n_trials=1)[0]
    else:
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
    remaining_seeds = set(seed_order)
    seed_attempt_count = 0
    rough_stop_hit = False
    for seed in seed_order:
        if seed not in remaining_seeds:
            continue
        seed_attempt_count += 1
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
            remaining_seeds.discard(seed)
            continue

        maximum_branch_count = max(maximum_branch_count, len(placements))
        for placement in placements:
            retained = tuple(sorted(map(int, placement.fragment)))
            mapping = tuple(sorted(
                (int(source_atom), int(target_atom))
                for source_atom, target_atom in placement.items()
                if source_atom in placement.fragment
            ))
            if (target_region_atoms is not None
                    and not target_region_atoms.intersection(
                        target_atom for _source_atom, target_atom in mapping)):
                continue
            certificate = canonicalizer.certificate(mapping)
            family = placement_families.get(certificate)
            if family is None:
                placement_families[certificate] = _InitialFragmentFamily(
                    retained_atoms=retained,
                    representative_mapping=mapping,
                    symmetry=dict(placement.symmetry or {}),
                )
            else:
                placement_families[certificate] = replace(
                    family, encounter_count=family.encounter_count + 1)
            if len(placement_families) >= config.candidate_limit:
                candidate_capped = True
                break
        if candidate_capped:
            break
        if config.seed_mode == "fragment_cover" and placements:
            discovered_fragment = set(map(int, placements[0].fragment))
            remaining_seeds.difference_update(discovered_fragment)
            if (len(discovered_fragment) / len(source)
                    > config.rough_retention_threshold):
                rough_stop_hit = True
                break

    return (
        tuple(placement_families.values()),
        capped_seed_count,
        maximum_branch_count,
        candidate_capped,
        seed_limited,
        seed_attempt_count,
        len(seed_order) - seed_attempt_count,
        rough_stop_hit,
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


def detect_fragments(
        source, target, *, source_id="",
        config: FragmentDetectionConfig | None = None,
        target_region_atoms=None):
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
    region = (
        frozenset(map(int, target_region_atoms))
        if target_region_atoms is not None else None
    )
    if region is not None and (
            not region or not region.issubset(set(map(int, target_graph)))):
        raise ValueError("target region must be a nonempty target-atom subset")
    (
        initial_placements,
        capped_seed_count,
        maximum_branch_count,
        candidate_capped,
        seed_limited,
        seed_attempt_count,
        seed_pruned_count,
        rough_stop_hit,
    ) = _initial_fragment_placements(
        source_graph, target_graph, config,
        target_orbits=target_context.atom_orbits,
        target_region_atoms=region)

    def placement_score(placement):
        retained = placement.retained_atoms
        mapping = placement.representative_mapping
        if region is None:
            return (len(retained),)
        images = {target for _source, target in mapping}
        return (len(images & region), len(retained))

    best_initial_score = max(
        map(placement_score, initial_placements), default=(0,))
    best_initial_size = max(
        (len(placement.retained_atoms) for placement in initial_placements
         if placement_score(placement) == best_initial_score),
        default=0)
    initial_placement_encounters = sum(
        placement.encounter_count for placement in initial_placements)
    best_initial_family_count = sum(
        placement_score(placement) == best_initial_score
        for placement in initial_placements)
    candidates = []
    seen_candidates = set()

    for placement in initial_placements:
        if placement_score(placement) != best_initial_score:
            continue
        retained = placement.retained_atoms
        mapping_pairs = placement.representative_mapping
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
            target_region_atoms=region,
            retained_symmetry=placement.symmetry,
        )
        maximum_branch_count = max(
            maximum_branch_count, augmented_branch_count)
        if augmented_capped:
            capped_seed_count += 1
        if not augmented_mappings:
            continue

        for augmented_placement in augmented_mappings:
            augmented_mapping = dict(augmented_placement.mapping)
            target_mapping = {
                source_atom: target_atom
                for source_atom, target_atom in augmented_placement.mapping
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
                aam_hierarchy=augmented_placement.hierarchy,
            )
            if (region is not None
                    and not region.intersection(
                        candidate.covered_target_atoms)):
                continue
            identity = _candidate_identity(candidate)
            if identity in seen_candidates:
                continue
            seen_candidates.add(identity)
            candidates.append(candidate)
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
    approximate = bool(
        config.seed_mode == "fragment_cover"
        and (seed_pruned_count or rough_stop_hit))
    incomplete = bool(cap_hit or seed_limited or approximate)
    status = (
        "capped" if cap_hit
        else ("seed_limited" if seed_limited
              else ("rough" if approximate
                    else ("matched" if candidates else "no_match")))
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
        initial_placement_encounters=initial_placement_encounters,
        initial_family_count=len(initial_placements),
        best_initial_family_count=best_initial_family_count,
        seed_attempt_count=seed_attempt_count,
        seed_pruned_count=seed_pruned_count,
        rough_stop_hit=rough_stop_hit,
    )
