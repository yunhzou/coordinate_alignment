"""Bounded enumeration of disjoint target-coverage patterns."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache


@dataclass(frozen=True)
class CoverageEnumerationConfig:
    maximum_precursors: int = 3
    mode: str = "modular"
    beam_width: int = 200
    patterns_per_coverage: int = 4
    state_limit: int = 100_000
    allow_overlaps: bool = False

    def __post_init__(self):
        if self.maximum_precursors < 1:
            raise ValueError("maximum_precursors must be positive")
        if self.mode not in {"modular", "recommendation", "exhaustive"}:
            raise ValueError(f"unknown coverage search mode: {self.mode}")
        if min(self.beam_width, self.patterns_per_coverage,
               self.state_limit) < 1:
            raise ValueError("coverage search limits must be positive")


@dataclass(frozen=True)
class CoverageEnumerationResult:
    patterns: tuple[tuple[int, ...], ...]
    complete: bool
    truncated: bool


def _mask_index(masks, atom_count):
    by_atom = defaultdict(list)
    for mask in masks:
        for atom in range(atom_count):
            if mask & (1 << atom):
                by_atom[atom].append(mask)
    for atom in by_atom:
        by_atom[atom].sort(key=lambda mask: (-mask.bit_count(), mask))
    return by_atom


def enumerate_coverage_patterns(
        masks, atom_count, rank_pattern, *,
        config: CoverageEnumerationConfig | None = None):
    """Enumerate disjoint coverage-mask tuples that exactly cover a target.

    ``rank_pattern(pattern, covered_atom_count)`` ranks partial and complete
    paths but does not control their validity.  With ``allow_overlaps``, masks
    describe available target atoms; exact ownership is resolved afterwards.
    """
    config = config or CoverageEnumerationConfig()
    masks = tuple(sorted(set(masks), key=lambda mask: (
        -mask.bit_count(), mask)))
    full_mask = (1 << atom_count) - 1
    masks_by_atom = _mask_index(masks, atom_count)
    largest_mask_size = max(
        (mask.bit_count() for mask in masks), default=0)
    truncated = False

    @lru_cache(maxsize=None)
    def can_complete(covered, slots):
        uncovered = full_mask ^ covered
        if not uncovered:
            return True
        compatible = [
            mask for mask in masks
            if (mask & uncovered)
            and (config.allow_overlaps or not (mask & covered))
        ]
        if not compatible or slots == 0:
            return False
        possible = 0
        sizes = []
        for mask in compatible:
            contribution = mask & uncovered
            possible |= contribution
            sizes.append(contribution.bit_count())
        if uncovered & ~possible:
            return False
        sizes.sort(reverse=True)
        return uncovered.bit_count() <= sum(sizes[:slots])

    def exhaustive(covered=0, selected=()):
        if covered == full_mask:
            yield selected
            return
        if len(selected) >= config.maximum_precursors:
            return
        if not can_complete(
                covered, config.maximum_precursors - len(selected)):
            return
        uncovered = full_mask ^ covered
        pivot = (uncovered & -uncovered).bit_length() - 1
        for mask in masks_by_atom[pivot]:
            if not config.allow_overlaps and mask & covered:
                continue
            yield from exhaustive(covered | mask, selected + (mask,))

    def recommendation():
        nonlocal truncated
        states = [(0, ())]
        completed = []
        for _depth in range(config.maximum_precursors):
            next_states = []
            for covered, selected in states:
                uncovered = full_mask ^ covered
                pivot = (uncovered & -uncovered).bit_length() - 1
                for mask in masks_by_atom[pivot]:
                    if not config.allow_overlaps and mask & covered:
                        continue
                    new_covered = covered | mask
                    new_selected = selected + (mask,)
                    if new_covered == full_mask:
                        completed.append(new_selected)
                    elif ((full_mask ^ new_covered).bit_count()
                          <= (config.maximum_precursors - len(new_selected))
                          * largest_mask_size):
                        next_states.append((new_covered, new_selected))
            if len(next_states) > config.beam_width:
                truncated = True
            next_states.sort(
                key=lambda state: rank_pattern(
                    state[1], state[0].bit_count()))
            states = next_states[:config.beam_width]
            if not states:
                break
        if len(completed) > config.beam_width:
            truncated = True
        completed.sort(key=lambda pattern: rank_pattern(pattern, atom_count))
        return completed[:config.beam_width]

    def modular():
        nonlocal truncated
        states = {0: [()]}
        completed = []
        for _depth in range(config.maximum_precursors):
            next_states = defaultdict(list)
            for covered, paths in states.items():
                uncovered = full_mask ^ covered
                pivot = (uncovered & -uncovered).bit_length() - 1
                for selected in paths:
                    for mask in masks_by_atom[pivot]:
                        if not config.allow_overlaps and mask & covered:
                            continue
                        new_covered = covered | mask
                        new_selected = selected + (mask,)
                        if new_covered == full_mask:
                            completed.append(new_selected)
                            continue
                        remaining_slots = (
                            config.maximum_precursors - len(new_selected))
                        if ((full_mask ^ new_covered).bit_count()
                                > remaining_slots * largest_mask_size):
                            continue
                        bucket = next_states[new_covered]
                        bucket.append(new_selected)
                        if len(bucket) >= config.patterns_per_coverage * 2:
                            bucket.sort(key=lambda pattern: rank_pattern(
                                pattern, new_covered.bit_count()))
                            del bucket[config.patterns_per_coverage:]
            for covered, paths in next_states.items():
                paths.sort(key=lambda pattern: rank_pattern(
                    pattern, covered.bit_count()))
                del paths[config.patterns_per_coverage:]
            state_count = sum(map(len, next_states.values()))
            if state_count > config.state_limit:
                truncated = True
                ranked = sorted(
                    ((covered, pattern)
                     for covered, paths in next_states.items()
                     for pattern in paths),
                    key=lambda pair: rank_pattern(
                        pair[1], pair[0].bit_count()),
                )[:config.state_limit]
                next_states = defaultdict(list)
                for covered, pattern in ranked:
                    next_states[covered].append(pattern)
            states = dict(next_states)
            if not states:
                break
        return completed

    if config.mode == "modular":
        patterns = modular()
    elif config.mode == "recommendation":
        patterns = recommendation()
    else:
        patterns = exhaustive()
    patterns = tuple(patterns)
    return CoverageEnumerationResult(
        patterns=patterns,
        complete=not truncated,
        truncated=truncated,
    )
