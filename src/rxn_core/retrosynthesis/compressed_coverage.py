"""Recommendation over compressed target-occupation alternatives."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import heapq
import itertools

from ..fragment_matching import materialize_target_coverage_orbit

TargetOccupation = tuple[int, ...]
CoverageSignature = tuple[TargetOccupation, ...]
CoveragePattern = tuple[CoverageSignature, ...]


def candidate_target_domains(candidate):
    """Return AAM atom domains for display and diagnostic purposes.

    These domains describe allowed atoms, but correlated AAM placements must
    not be treated as their Cartesian product during assembly.
    """
    retained = set(map(int, candidate.retained_atoms))
    witness = {int(source): int(target) for source, target in candidate.mapping}
    pools = {
        source: (target,) for source, target in witness.items()
        if source in retained
    }
    for fragment in candidate.aam_hierarchy.fragments:
        for domain in fragment.symmetry_domains:
            target_pool = tuple(map(int, domain.p_atoms))
            for source in domain.r_atoms:
                if source in retained:
                    pools[int(source)] = target_pool
    return tuple((source, pools[source]) for source in sorted(pools))


def candidate_target_occupations(
        candidate, target, *, iso_tolerance=0.5, generators=None,
        orbit_cache=None):
    """Return the distinct whole-fragment target regions allowed by AAM.

    The existing AAM candidate remains compressed internally.  We quotient
    only by target ownership/attachments, so internal atom bijections are
    never enumerated.  Each returned record is one chemically connected
    occupation spot that assembly may consume with a copy of the precursor.
    """
    if orbit_cache is None or generators is None:
        variants = materialize_target_coverage_orbit(
            candidate, target, iso_tolerance=iso_tolerance,
            generators=generators)
        return tuple({
            "covered_target_atoms": tuple(map(
                int, variant.covered_target_atoms)),
            "mapping": tuple(
                (int(source), int(target_atom))
                for source, target_atom in variant.mapping),
            "attachment_atoms_target": tuple(map(
                int, variant.attachment_atoms_target)),
        } for variant in variants)

    coverage = tuple(sorted(map(int, candidate.covered_target_atoms)))
    attachments = tuple(sorted(map(int, candidate.attachment_atoms_target)))
    cache_key = coverage, attachments
    orbit = orbit_cache.get(cache_key)
    tracked_atoms = tuple(sorted(set(coverage) | set(attachments)))
    if orbit is None:
        degree = len(target)
        generator_images = tuple(tuple(
            int(generator.get(atom, atom)) for atom in range(degree)
        ) for generator in generators)

        def atom_mask(atoms):
            return sum(1 << atom for atom in atoms)

        def transform_mask(mask, generator):
            transformed = 0
            while mask:
                bit = mask & -mask
                transformed |= 1 << generator[bit.bit_length() - 1]
                mask ^= bit
            return transformed

        initial_state = atom_mask(coverage), atom_mask(attachments)
        witnesses_by_state = {initial_state: tracked_atoms}
        queue = [initial_state]
        cursor = 0
        while cursor < len(queue):
            state = queue[cursor]
            cursor += 1
            witness = witnesses_by_state[state]
            for generator in generator_images:
                transformed_state = (
                    transform_mask(state[0], generator),
                    transform_mask(state[1], generator),
                )
                if transformed_state in witnesses_by_state:
                    continue
                witnesses_by_state[transformed_state] = tuple(
                    generator[atom] for atom in witness)
                queue.append(transformed_state)
        orbit = tuple(
            (state, witnesses_by_state[state])
            for state in sorted(witnesses_by_state))
        orbit_cache[cache_key] = orbit

    def mask_atoms(mask):
        return tuple(
            atom for atom in range(len(target)) if mask & (1 << atom))

    occupations = []
    for state, witness in orbit:
        image_by_atom = dict(zip(tracked_atoms, witness))
        occupations.append({
            "covered_target_atoms": mask_atoms(state[0]),
            "mapping": tuple(sorted(
                (int(source), image_by_atom[int(target_atom)])
                for source, target_atom in candidate.mapping)),
            "attachment_atoms_target": mask_atoms(state[1]),
        })
    return tuple(occupations)


def coverage_signature(target_occupations):
    """Source-index-free set of allowed whole-fragment occupations."""
    return tuple(sorted({
        tuple(sorted(map(int, occupation["covered_target_atoms"])))
        for occupation in target_occupations
    }))


def assign_occupation_signatures(signatures, atom_count):
    """Choose whole-fragment occupations whose union covers the target."""
    signatures = tuple(signatures)
    full = frozenset(range(atom_count))
    assigned = [None] * len(signatures)

    def visit(remaining, covered):
        if not remaining:
            return covered == full
        choices = []
        for index in remaining:
            compatible = [
                occupation for occupation in signatures[index]
                if not frozenset(occupation).issubset(covered)
            ]
            if not compatible:
                return False
            choices.append((len(compatible), index, compatible))
        _count, index, compatible = min(choices)
        rest = tuple(item for item in remaining if item != index)
        for occupation in compatible:
            assigned[index] = occupation
            if visit(rest, covered.union(occupation)):
                return True
        assigned[index] = None
        return False

    if not visit(tuple(range(len(signatures))), frozenset()):
        return None
    return tuple(assigned)


def assign_candidate_items(items, atom_count):
    """Place candidate copies so their allowed AAM regions cover the target."""
    items = tuple(items)
    full = frozenset(range(atom_count))
    selected = [None] * len(items)

    def visit(remaining, covered):
        if not remaining:
            return covered == full
        choices = []
        for index in remaining:
            compatible = [
                occupation for occupation in items[index]["target_occupations"]
                if not frozenset(occupation["covered_target_atoms"]).issubset(
                    covered)
            ]
            if not compatible:
                return False
            choices.append((len(compatible), index, compatible))
        _count, index, compatible = min(choices)
        rest = tuple(item for item in remaining if item != index)
        for occupation in compatible:
            selected[index] = occupation
            if visit(rest, covered.union(
                    occupation["covered_target_atoms"])):
                return True
        selected[index] = None
        return False

    if not visit(tuple(range(len(items))), frozenset()):
        return None
    placed = []
    for item, occupation in zip(items, selected):
        transformed = dict(item)
        transformed.update({
            "mapping": [list(pair) for pair in occupation["mapping"]],
            "covered_target_atoms": list(
                occupation["covered_target_atoms"]),
            "attachment_atoms_target": list(
                occupation["attachment_atoms_target"]),
        })
        placed.append(transformed)
    return tuple(placed)


def place_candidate_items(items, covered_target_atoms):
    """Materialize candidate mappings for an already assigned pattern.

    Items in one coverage-signature pool have the same whole-fragment target
    regions.  Assembly therefore solves region assignment once and reuses the
    selected regions while comparing precursor substitutions.
    """
    placed = []
    for item, covered_atoms in zip(items, covered_target_atoms, strict=True):
        covered_atoms = tuple(covered_atoms)
        occupation = next(
            occupation for occupation in item["target_occupations"]
            if tuple(occupation["covered_target_atoms"]) == covered_atoms
        )
        transformed = dict(item)
        transformed.update({
            "mapping": [list(pair) for pair in occupation["mapping"]],
            "covered_target_atoms": list(
                occupation["covered_target_atoms"]),
            "attachment_atoms_target": list(
                occupation["attachment_atoms_target"]),
        })
        placed.append(transformed)
    return tuple(placed)


def place_candidate_items_maximum_coverage(items, atom_count):
    """Select one correlated AAM occupation per item with maximum union.

    This is used to inspect a specified precursor set when it does not form a
    complete cover. Occupations with the same target region are equivalent for
    this objective, so only one mapping witness per region is explored.
    """
    items = tuple(items)
    choices = []
    for item in items:
        by_mask = {}
        for occupation in item["target_occupations"]:
            mask = sum(
                1 << atom for atom in occupation["covered_target_atoms"])
            by_mask.setdefault(mask, occupation)
        choices.append(tuple(by_mask.items()))
    order = tuple(sorted(range(len(items)), key=lambda index: len(choices[index])))
    possible_suffix = [0] * (len(order) + 1)
    for depth in range(len(order) - 1, -1, -1):
        possible_suffix[depth] = possible_suffix[depth + 1]
        for mask, _occupation in choices[order[depth]]:
            possible_suffix[depth] |= mask

    selected = [None] * len(items)
    best = None
    best_key = None

    def visit(depth, covered, claimed_count):
        nonlocal best, best_key
        upper_coverage = (covered | possible_suffix[depth]).bit_count()
        if best_key is not None and upper_coverage < -best_key[0]:
            return
        if depth == len(order):
            key = (-covered.bit_count(), claimed_count - covered.bit_count())
            if best_key is None or key < best_key:
                best_key = key
                best = tuple(selected)
            return
        index = order[depth]
        ranked_choices = sorted(
            choices[index],
            key=lambda pair: (-(pair[0] & ~covered).bit_count(), pair[0]),
        )
        for mask, occupation in ranked_choices:
            selected[index] = occupation
            visit(depth + 1, covered | mask, claimed_count + mask.bit_count())

    visit(0, 0, 0)
    placed = []
    for item, occupation in zip(items, best, strict=True):
        transformed = dict(item)
        transformed.update({
            "mapping": [list(pair) for pair in occupation["mapping"]],
            "covered_target_atoms": list(occupation["covered_target_atoms"]),
            "attachment_atoms_target": list(
                occupation["attachment_atoms_target"]),
        })
        placed.append(transformed)
    return tuple(placed)


@dataclass(frozen=True)
class CoverageRecommendationResult:
    patterns: tuple[CoveragePattern, ...]
    occupations: tuple[tuple[TargetOccupation, ...], ...]
    truncated: bool


@dataclass(frozen=True)
class CoverageRecommendationConfig:
    maximum_precursors: int | None = None

    def __post_init__(self):
        if (self.maximum_precursors is not None
                and self.maximum_precursors < 1):
            raise ValueError("maximum precursors must be positive")


def recommend_compressed_coverage_patterns(
        signatures, atom_count, rank_pattern, *, result_limit,
        config: CoverageRecommendationConfig | None = None):
    """Recommend full covers with a target-scaled best-first search."""
    config = config or CoverageRecommendationConfig()
    maximum_precursors = config.maximum_precursors or atom_count
    signatures = tuple(sorted(set(signatures), key=lambda signature: (
        rank_pattern((signature,), max(map(len, signature))), signature)))
    occupation_masks = tuple(tuple(
        sum(1 << atom for atom in occupation)
        for occupation in signature
    ) for signature in signatures)
    signature_sizes = tuple(len(signature[0]) for signature in signatures)
    maximum_signature_size = max(signature_sizes)
    by_target = defaultdict(list)
    for signature_index, signature in enumerate(signatures):
        for occupation in signature:
            for target in occupation:
                by_target[target].append(signature_index)
    by_target = {
        target: tuple(sorted(set(indices)))
        for target, indices in by_target.items()
    }

    completed = {}
    full = (1 << atom_count) - 1
    serial = itertools.count()

    def atoms(mask):
        return tuple(
            atom for atom in range(atom_count) if mask & (1 << atom))

    def pattern(selected):
        return tuple(signatures[index] for index in selected)

    def queue_rank(selected, covered):
        covered_count = covered.bit_count()
        overlap = sum(signature_sizes[index] for index in selected) \
            - covered_count
        return (
            -covered_count,
            len(set(selected)),
            overlap,
            next(serial),
        )

    def repeated_cover_count(masks):
        possible = 0
        for occupation in masks:
            possible |= occupation
        if full & ~possible:
            return None
        reachable = {0}
        for copy_count in range(1, maximum_precursors + 1):
            reachable = {
                covered | occupation
                for covered in reachable for occupation in masks
                if occupation & ~covered
            }
            if full in reachable:
                return copy_count
        return None

    for signature_index, masks in enumerate(occupation_masks):
        copy_count = repeated_cover_count(masks)
        if copy_count is not None:
            selected = (signature_index,) * copy_count
            selected_pattern = pattern(selected)
            witness = assign_occupation_signatures(
                selected_pattern, atom_count)
            completed[selected] = selected_pattern, witness

    # Exact complementary pairs are cheap to index and are the strongest
    # two-module constructions: they cover every target atom with no overlap.
    # Seed all of them so abundant one-module covers cannot crowd them out.
    signatures_by_mask = defaultdict(set)
    for signature_index, masks in enumerate(occupation_masks):
        for mask in masks:
            signatures_by_mask[mask].add(signature_index)
    for left_index, masks in enumerate(occupation_masks):
        for left_mask in masks:
            complement = full ^ left_mask
            for right_index in signatures_by_mask.get(complement, ()):
                selected = tuple(sorted((left_index, right_index)))
                selected_pattern = pattern(selected)
                witness = assign_occupation_signatures(
                    selected_pattern, atom_count)
                completed.setdefault(selected, (selected_pattern, witness))

    preloaded_cover_count = len(completed)
    queue = [(queue_rank((), 0), (), 0, ())]
    seen = {((), 0)}
    state_budget = result_limit * max(atom_count, 1)
    branch_budget = result_limit * maximum_precursors
    visited_states = 0
    truncated = False
    while (queue
           and len(completed) < preloaded_cover_count + result_limit
           and visited_states < state_budget):
        _priority, selected, covered, assigned_masks = heapq.heappop(queue)
        visited_states += 1
        missing = full & ~covered
        if missing == 0:
            selected_pattern = pattern(selected)
            completed.setdefault(selected, (
                selected_pattern,
                tuple(atoms(mask) for mask in assigned_masks),
            ))
            continue
        if len(selected) >= maximum_precursors:
            continue
        remaining_slots = maximum_precursors - len(selected)
        if missing.bit_count() > remaining_slots * maximum_signature_size:
            continue
        missing_atoms = tuple(
            atom for atom in range(atom_count)
            if missing & (1 << atom))
        pivot = min(missing_atoms, key=lambda target: (
            len(by_target.get(target, ())), target))
        compatible_signatures = tuple(sorted(
            by_target.get(pivot, ()),
            key=lambda signature_index: (
                -max(
                    (occupation & missing).bit_count()
                    for occupation in occupation_masks[signature_index]),
                signature_index,
            ),
        ))
        if len(compatible_signatures) > branch_budget:
            truncated = True
        for signature_index in compatible_signatures[:branch_budget]:
            for occupation in occupation_masks[signature_index]:
                if not occupation & (1 << pivot):
                    continue
                new_covered = covered | occupation
                selected_pairs = sorted(
                    zip(selected + (signature_index,),
                        assigned_masks + (occupation,)))
                new_selected = tuple(
                    index for index, _mask in selected_pairs)
                new_assigned_masks = tuple(
                    mask for _index, mask in selected_pairs)
                state_key = new_selected, new_covered
                if state_key in seen:
                    continue
                seen.add(state_key)
                heapq.heappush(queue, (
                    queue_rank(new_selected, new_covered),
                    new_selected,
                    new_covered,
                    new_assigned_masks,
                ))
    if queue:
        truncated = True
    ranked = sorted(completed.values(), key=lambda result: (
        sum(len(signature[0]) for signature in result[0]) - atom_count,
        rank_pattern(result[0], atom_count),
    ))
    by_module_count = defaultdict(list)
    for result in ranked:
        by_module_count[len(set(result[0]))].append(result)
    diverse = []
    module_counts = sorted(by_module_count)
    for offset in range(max(map(len, by_module_count.values()), default=0)):
        for module_count in module_counts:
            patterns = by_module_count[module_count]
            if offset < len(patterns):
                diverse.append(patterns[offset])
                if len(diverse) >= result_limit:
                    break
        if len(diverse) >= result_limit:
            break
    return CoverageRecommendationResult(
        patterns=tuple(result[0] for result in diverse),
        occupations=tuple(result[1] for result in diverse),
        truncated=truncated,
    )
