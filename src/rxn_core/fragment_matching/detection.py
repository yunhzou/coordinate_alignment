"""Detect target-owned fragments from one source graph."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, replace

from .._native import paired_mapping_invariant as _native_paired_invariant
from ..alignment.branch import _generate_seed_orders
from ..fragment import FragmentMatchConfig, FragmentMatchContext, match_fragment
from ..matcher import (
    _PartialMappingCanonicalizer,
    _edge_wbo,
    _nauty_atom_generators,
    _nauty_orbits,
    _orbit_wbo_bucket,
)
from ..subgraph import _coerce_graph
from ..search_graph import AAMSearchGraph
from .augmentation import match_augmented_residuals, project_augmented_placement
from .graph_ops import partition_at_retained_fragment, fragment_equivalence_classes
from .models import (
    FragmentCandidate,
    FragmentDerivation,
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
    search_paths: tuple = ()


def _paired_mapping_invariant(mapping, source_orbits, target_orbits,
                              source, target):
    """Exact-invariant WL partition of a partial endpoint relation."""
    pairs = tuple(mapping)
    raw_colors = [
        (source_orbits[source_atom], target_orbits[target_atom])
        for source_atom, target_atom in pairs
    ]
    source_zero = getattr(source_orbits, "zero_bucket", None)
    target_zero = getattr(target_orbits, "zero_bucket", None)
    if source_zero is not None and target_zero is not None:
        zero_relation = (source_zero, target_zero)
        source_positions = {
            source_atom: position
            for position, (source_atom, _target_atom) in enumerate(pairs)
        }
        target_positions = {
            target_atom: position
            for position, (_source_atom, target_atom) in enumerate(pairs)
        }
        active_pairs = set()
        for position, (source_atom, target_atom) in enumerate(pairs):
            active_pairs.update(
                tuple(sorted((position, other_position)))
                for neighbor in source.neighbors(source_atom)
                if (other_position := source_positions.get(neighbor))
                is not None and other_position != position
            )
            active_pairs.update(
                tuple(sorted((position, other_position)))
                for neighbor in target.neighbors(target_atom)
                if (other_position := target_positions.get(neighbor))
                is not None and other_position != position
            )
        relations = {
            (left, right): (
                _orbit_wbo_bucket(
                    source_orbits, pairs[left][0], pairs[right][0],
                    _edge_wbo(source, pairs[left][0], pairs[right][0])),
                _orbit_wbo_bucket(
                    target_orbits, pairs[left][1], pairs[right][1],
                    _edge_wbo(target, pairs[left][1], pairs[right][1])),
            )
            for left, right in active_pairs
        }
        relations = {
            pair: relation
            for pair, relation in relations.items()
            if relation != zero_relation
        }
        return _native_paired_invariant(
            raw_colors,
            zero_relation,
            tuple(
                (left, right, relation[0], relation[1])
                for (left, right), relation in relations.items()
            ),
        )

    initial_color_counts = tuple(sorted(Counter(raw_colors).items()))

    def compact(values):
        classes = {
            value: index for index, value in enumerate(
                sorted(set(values), key=repr))
        }
        return [classes[value] for value in values]

    colors = compact(raw_colors)
    relations = {}
    for left in range(len(pairs)):
        source_left, target_left = pairs[left]
        for right in range(left + 1, len(pairs)):
            source_right, target_right = pairs[right]
            relations[(left, right)] = (
                _orbit_wbo_bucket(
                    source_orbits, source_left, source_right,
                    _edge_wbo(source, source_left, source_right)),
                _orbit_wbo_bucket(
                    target_orbits, target_left, target_right,
                    _edge_wbo(target, target_left, target_right)),
            )
    for _ in pairs:
        signatures = []
        for left in range(len(pairs)):
            neighborhood = []
            for right in range(len(pairs)):
                if left == right:
                    continue
                edge = relations[tuple(sorted((left, right)))]
                neighborhood.append((edge, colors[right]))
            signatures.append(
                (colors[left], tuple(sorted(neighborhood))))
        refined = compact(signatures)
        if len(set(refined)) == len(set(colors)):
            break
        colors = refined
    color_counts = tuple(sorted(Counter(colors).items()))
    relation_counts = tuple(sorted(Counter(
        (
            min(colors[left], colors[right]),
            max(colors[left], colors[right]),
            relation,
        )
        for (left, right), relation in relations.items()
    ).items()))
    return initial_color_counts, color_counts, relation_counts


class _InitialFamilyAccumulator:
    """Exact cross-seed quotient for initial AAM fragment families."""

    def __init__(self, source, target, config, target_region_atoms,
                 source_orbits, target_orbits):
        self.config = config
        self.target_region_atoms = target_region_atoms
        self.source_orbits = source_orbits
        self.target_orbits = target_orbits
        self.families = {}
        self.literal_families = {}
        self.coarse_buckets = {}
        self.canonicalizer = _PartialMappingCanonicalizer(
            source,
            target,
            wbo_tol=config.iso_tolerance,
            target_atom_tags=(
                {int(atom): "requested_region"
                 for atom in target_region_atoms}
                if target_region_atoms is not None else None),
        )

    def add(self, placements, graph):
        for placement, path in zip(placements, graph.paths(), strict=True):
            retained = tuple(sorted(map(int, placement.fragment)))
            mapping = tuple(sorted(
                (int(source_atom), int(target_atom))
                for source_atom, target_atom in placement.items()
                if source_atom in placement.fragment
            ))
            if (self.target_region_atoms is not None
                    and not self.target_region_atoms.intersection(
                        target_atom for _source_atom, target_atom in mapping)):
                continue
            literal_key = (retained, mapping)
            family_id = self.literal_families.get(literal_key)
            if family_id is not None:
                family = self.families[family_id]
                self.families[family_id] = replace(
                    family, encounter_count=family.encounter_count + 1,
                    search_paths=family.search_paths + (path,))
                continue
            coarse = _paired_mapping_invariant(
                mapping,
                self.source_orbits,
                self.target_orbits,
                self.canonicalizer.g_R,
                self.canonicalizer.g_P,
            )
            bucket = self.coarse_buckets.setdefault(coarse, [])
            family_id = None
            if bucket:
                for prior_id in bucket:
                    prior_mapping = self.families[
                        prior_id].representative_mapping
                    if self.canonicalizer.equivalent(
                            dict(prior_mapping), dict(mapping)):
                        family_id = prior_id
                        break
            if family_id is None:
                family_id = len(self.families)
                self.families[family_id] = _InitialFragmentFamily(
                    retained_atoms=retained,
                    representative_mapping=mapping,
                    symmetry=dict(placement.symmetry or {}),
                    search_paths=(path,),
                )
                bucket.append(family_id)
                self.literal_families[literal_key] = family_id
            else:
                family = self.families[family_id]
                self.families[family_id] = replace(
                    family, encounter_count=family.encounter_count + 1,
                    search_paths=family.search_paths + (path,))
                self.literal_families[literal_key] = family_id
            if (self.config.candidate_limit is not None
                    and len(self.families) >= self.config.candidate_limit):
                return True
        return False


def _initial_seed_order(source, config, source_orbits=None):
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
        if (config.seed_mode == "orbit_representatives"
                and source_orbits is not None):
            representatives = {}
            for atom in seed_order:
                representatives.setdefault(source_orbits[atom], atom)
            seed_order = list(representatives.values())
    symmetry_pruned = len(source) - len(seed_order)
    seed_limited = (
        config.seed_limit is not None
        and len(seed_order) > config.seed_limit)
    if config.seed_limit is not None:
        seed_order = seed_order[:config.seed_limit]
    return tuple(seed_order), seed_limited, symmetry_pruned


def _grow_initial_seed(
        source, target, seed, config, source_orbits, target_orbits):
    result = match_fragment(source, target, seed=seed,
        context=FragmentMatchContext(source_orbits=source_orbits, target_orbits=target_orbits),
        config=FragmentMatchConfig(graph_floor=config.graph_floor,
            iso_tolerance=config.iso_tolerance, minimum_size=config.minimum_fragment_size,
            branch_limit=config.branch_limit))
    graph = AAMSearchGraph.initial_fragment_search(source, target, seed, result, config)
    from ..search_symmetry import finalize_graph_symmetry
    graph, _metrics = finalize_graph_symmetry(graph, target,
                                             iso_tolerance=config.iso_tolerance)
    matches = tuple(replace(match, symmetry=edge.match['symmetry'])
                    for match, edge in zip(result.matches,
                        (e for e in graph.transitions if e.match is not None), strict=True))
    return matches, result.capped, result.branch_count, graph


def _initial_fragment_placements(
        source, target, config, *, target_orbits=None,
        target_region_atoms=None):
    source_orbits = _nauty_orbits(source, wbo_tol=config.iso_tolerance)
    if target_orbits is None:
        target_orbits = _nauty_orbits(
            target, wbo_tol=config.iso_tolerance)
    accumulator = _InitialFamilyAccumulator(
        source, target, config, target_region_atoms,
        source_orbits, target_orbits)
    capped_seed_count = 0
    maximum_branch_count = 0
    candidate_capped = False

    seed_order, seed_limited, symmetry_pruned = _initial_seed_order(
        source, config, source_orbits)
    remaining_seeds = set(seed_order)
    seed_attempt_count = 0
    rough_stop_hit = False
    search_graphs = []
    for seed in seed_order:
        if seed not in remaining_seeds:
            continue
        seed_attempt_count += 1
        placements, capped, branch_count, graph = _grow_initial_seed(
            source, target, seed, config, source_orbits, target_orbits)
        search_graphs.append(graph)
        if capped:
            capped_seed_count += 1
            maximum_branch_count = max(maximum_branch_count, branch_count)
            remaining_seeds.discard(seed)
            continue

        maximum_branch_count = max(maximum_branch_count, branch_count)
        candidate_capped = accumulator.add(placements, graph) if placements else False
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
        tuple(accumulator.families.values()),
        capped_seed_count,
        maximum_branch_count,
        candidate_capped,
        seed_limited,
        seed_attempt_count,
        symmetry_pruned + len(seed_order) - seed_attempt_count,
        rough_stop_hit,
        tuple(search_graphs),
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


def _prepare_fragment_detection(
        source, target, config, target_region_atoms=None):
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
    return source_graph, target_context, region


def _augment_initial_family(
        source_graph, target_graph, placement, config, region):
    retained = placement.retained_atoms
    mapping_pairs = placement.representative_mapping
    outside, boundary, _fragments = partition_at_retained_fragment(
        source_graph, retained)
    boundary = tuple(sorted(set(boundary) | set(placement.search_paths[0].deferred_edges)))
    if (config.maximum_boundary_bonds is not None
            and len(boundary) > config.maximum_boundary_bonds):
        return (), False, 0, ()

    augmented = match_augmented_residuals(
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
    candidates = []
    for augmented_placement in augmented.placements:
        placement_boundary = tuple(sorted(set(boundary) | {
            tuple(sorted(edge)) for fragment in augmented_placement.hierarchy.fragments
            for edge in fragment.deferred_edges}))
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
            placement_boundary,
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
            source_id="",
            mapping=tuple(sorted(target_mapping.items())),
            retained_atoms=retained_atoms,
            covered_target_atoms=tuple(sorted(target_mapping.values())),
            leftover_fragments=leftover_fragments,
            boundary_bonds=placement_boundary,
            attachment_atoms_source=attachment_atoms_source,
            attachment_atoms_target=attachment_atoms_target,
            copied_residual_placements=copied_residual_placements,
            augmented_target_atom_count=augmented.augmented_target_atom_count,
            retained_fragments=retained_fragments,
            aam_hierarchy=augmented_placement.hierarchy,
            fragment_classes=fragment_equivalence_classes(source_graph,
                placement_boundary, retained_fragments, config.iso_tolerance),
            preserved_source_bonds=tuple(sorted(tuple(sorted((a, b)))
                for a, b in source_graph.edges() if a in target_mapping and b in target_mapping
                and tuple(sorted((a, b))) not in placement_boundary)),
            derivations=(FragmentDerivation(placement.search_paths,
                                           augmented_placement.search_paths,
                                           augmented_placement.target_action, True),),
        )
        if (region is None
                or region.intersection(candidate.covered_target_atoms)):
            candidates.append(candidate)
    return (
        tuple(candidates),
        augmented.capped,
        augmented.maximum_branch_count,
        augmented.search_graphs,
    )


def _detect_fragments_from_initial(
        source_graph, target_context, initial_search, *, source_id, config,
        region, augmentation_runner=None):
    target_graph = target_context.graph
    (
        initial_placements,
        capped_seed_count,
        maximum_branch_count,
        candidate_capped,
        seed_limited,
        seed_attempt_count,
        seed_pruned_count,
        rough_stop_hit,
        search_graphs,
    ) = initial_search

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
    # Every exact AAM family is a valid occupation hypothesis.  Ranking by the
    # initial connected island before residual augmentation discards families
    # that become optimal only after the other fragments are placed.
    selected_placements = initial_placements
    if augmentation_runner is None:
        augmentation_results = (
            _augment_initial_family(
                source_graph, target_graph, placement, config, region)
            for placement in selected_placements
        )
    else:
        augmentation_results = augmentation_runner(selected_placements)

    candidates = []
    seen_candidates = {}
    search_graphs = list(search_graphs)
    for family_candidates, augmented_capped, augmented_branch_count, augmented_graphs in (
            augmentation_results):
        search_graphs.extend(augmented_graphs)
        maximum_branch_count = max(
            maximum_branch_count, augmented_branch_count)
        if augmented_capped:
            capped_seed_count += 1
        for raw_candidate in family_candidates:
            candidate = replace(raw_candidate, source_id=str(source_id))
            identity = _candidate_identity(candidate)
            if identity in seen_candidates:
                index = seen_candidates[identity]
                candidates[index] = replace(candidates[index],
                    derivations=candidates[index].derivations + candidate.derivations)
                continue
            seen_candidates[identity] = len(candidates)
            candidates.append(candidate)
            if config.candidate_limit is not None and len(candidates) >= config.candidate_limit:
                candidate_capped = True
                break
            if candidate_capped:
                break
        if config.candidate_limit is not None and len(candidates) >= config.candidate_limit:
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
        config.seed_mode in {"fragment_cover", "orbit_representatives"}
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
        search_graphs=tuple(search_graphs),
    )


def detect_fragments(
        source, target, *, source_id="",
        config: FragmentDetectionConfig | None = None,
        target_region_atoms=None):
    """Generate augmented fragment candidates for one source-target pair."""
    config = config or FragmentDetectionConfig()
    source_graph, target_context, region = _prepare_fragment_detection(
        source, target, config, target_region_atoms)
    initial_search = _initial_fragment_placements(
        source_graph,
        target_context.graph,
        config,
        target_orbits=target_context.atom_orbits,
        target_region_atoms=region,
    )
    return _detect_fragments_from_initial(
        source_graph,
        target_context,
        initial_search,
        source_id=str(source_id),
        config=config,
        region=region,
    )
