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
    permutations = orbit_cache.get(cache_key)
    if permutations is None:
        degree = len(target)
        identity = tuple(range(degree))
        permutations_by_state = {(coverage, attachments): identity}
        queue = [identity]
        cursor = 0
        while cursor < len(queue):
            permutation = queue[cursor]
            cursor += 1
            for generator in generators:
                transformed = tuple(
                    int(generator.get(image, image))
                    for image in permutation)
                state = (
                    tuple(sorted(transformed[atom] for atom in coverage)),
                    tuple(sorted(transformed[atom] for atom in attachments)),
                )
                if state in permutations_by_state:
                    continue
                permutations_by_state[state] = transformed
                queue.append(transformed)
        permutations = tuple(
            permutations_by_state[state]
            for state in sorted(permutations_by_state))
        orbit_cache[cache_key] = permutations
    return tuple({
        "covered_target_atoms": tuple(sorted(
            permutation[atom] for atom in coverage)),
        "mapping": tuple(sorted(
            (int(source), permutation[int(target_atom)])
            for source, target_atom in candidate.mapping)),
        "attachment_atoms_target": tuple(sorted(
            permutation[atom] for atom in attachments)),
    } for permutation in permutations)


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


@dataclass(frozen=True)
class CoverageRecommendationResult:
    patterns: tuple[CoveragePattern, ...]
    truncated: bool


@dataclass(frozen=True)
class CoverageRecommendationConfig:
    maximum_precursors: int = 3

    def __post_init__(self):
        if self.maximum_precursors < 1:
            raise ValueError("maximum precursors must be positive")


def recommend_compressed_coverage_patterns(
        signatures, atom_count, rank_pattern, *, result_limit,
        config: CoverageRecommendationConfig | None = None):
    """Recommend full covers over compressed whole-fragment occupations."""
    config = config or CoverageRecommendationConfig()
    signatures = tuple(sorted(set(signatures), key=lambda signature: (
        rank_pattern((signature,), max(map(len, signature))), signature)))
    by_target = defaultdict(list)
    for signature in signatures:
        for occupation in signature:
            for target in occupation:
                by_target[target].append(signature)

    completed = {}
    full = frozenset(range(atom_count))
    serial = itertools.count()

    def queue_rank(selected, occupations, covered):
        overlap = sum(map(len, occupations)) - len(covered)
        return (
            overlap,
            rank_pattern(selected, len(covered)),
            -len(covered),
            next(serial),
        )

    queue = [(queue_rank((), (), frozenset()), (), (), frozenset())]
    seen = {((), ())}
    while queue and len(completed) < result_limit:
        _priority, selected, occupations, covered = heapq.heappop(queue)
        missing = full - covered
        if not missing:
            if any(
                    occupation.issubset(frozenset().union(*(
                        other for position, other in enumerate(occupations)
                        if position != index)))
                    for index, occupation in enumerate(occupations)):
                continue
            key = tuple(sorted(selected))
            completed.setdefault(key, selected)
            continue
        if len(selected) >= config.maximum_precursors:
            continue
        pivot = min(missing, key=lambda target: (
            len(by_target[target]), target))
        for signature in by_target[pivot]:
            for occupation in signature:
                occupation = frozenset(occupation)
                if pivot not in occupation:
                    continue
                new_covered = covered.union(occupation)
                new_selected = selected + (signature,)
                new_occupations = occupations + (occupation,)
                state_key = (
                    tuple(sorted(new_selected)),
                    tuple(sorted(new_occupations)),
                )
                if state_key in seen:
                    continue
                seen.add(state_key)
                heapq.heappush(queue, (
                    queue_rank(new_selected, new_occupations, new_covered),
                    new_selected,
                    new_occupations,
                    new_covered,
                ))
    ranked = sorted(completed.values(), key=lambda pattern: rank_pattern(
        pattern, atom_count))
    return CoverageRecommendationResult(
        patterns=tuple(ranked[:result_limit]),
        truncated=bool(queue),
    )
