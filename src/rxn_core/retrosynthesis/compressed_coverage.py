"""Correlated fragment occupations and exact target-union construction."""
from dataclasses import dataclass
from itertools import product

from ..fragment_matching import materialize_target_coverage_orbit
from .decision_graph import CoverageDecisionGraph

TargetOccupation = tuple[int, ...]
CoverageSignature = tuple[TargetOccupation, ...]
CoveragePattern = tuple[CoverageSignature, ...]


def candidate_target_domains(candidate):
    """Display-only atom domains. Never multiplied into assembly choices."""
    pools = {source: (target,) for source, target in candidate.mapping}
    target_owned = set(candidate.covered_target_atoms)
    for fragment in candidate.aam_hierarchy.fragments:
        for domain in fragment.symmetry_domains:
            for source in domain.r_atoms:
                if source in pools:
                    pools[source] = tuple(a for a in domain.p_atoms if a in target_owned)
    return tuple(sorted(pools.items()))


def candidate_target_occupations(candidate, target, *, iso_tolerance=0.5,
                                 generators=None):
    """Project exact branch actions onto joint matched-fragment regions.

    A coverage-only cache is invalid: equal atom sets can carry different
    fragment partitions and conditioned matching relations.
    """
    variants = materialize_target_coverage_orbit(candidate, target,
        iso_tolerance=iso_tolerance, generators=generators)
    return tuple({
        "covered_target_atoms": variant.covered_target_atoms,
        "mapping": variant.mapping,
        "attachment_atoms_target": variant.attachment_atoms_target,
        "retained_fragments": variant.retained_fragments,
        "target_fragment_atoms": tuple(tuple(sorted(variant.atom_mapping[a] for a in f))
                                        for f in variant.retained_fragments),
        "aam_hierarchy": variant.aam_hierarchy.to_record(),
        "target_domains": candidate_target_domains(variant),
        "derivation_actions": tuple(d.target_action for d in variant.derivations),
    } for variant in variants)


def coverage_signature(occupations):
    """Coverage index only; not an identity for matching relations."""
    return tuple(sorted({tuple(sorted(o["covered_target_atoms"])) for o in occupations}))


def place_item(item, occupation):
    """Transport the entire selected relation, including bonds and evidence."""
    placed = dict(item, **occupation)
    placed["allowed_target_regions"] = tuple(o["covered_target_atoms"]
                                               for o in item["target_occupations"])
    mapping = dict(placed["mapping"])
    placed["preserved_target_bonds"] = tuple(sorted(tuple(sorted((mapping[a], mapping[b])))
        for a, b in item.get("preserved_source_bonds", ())))
    return placed


def iter_candidate_assignments(items, atom_count):
    """All correlated choices for given precursor copies, lazily."""
    full = set(range(atom_count))
    for choices in product(*(item["target_occupations"] for item in items)):
        if set().union(*(set(o["covered_target_atoms"]) for o in choices)) == full:
            yield tuple(place_item(item, occupation)
                        for item, occupation in zip(items, choices, strict=True))


def assign_candidate_items(items, atom_count):
    return next(iter_candidate_assignments(tuple(items), atom_count), None)


def assign_occupation_signatures(signatures, atom_count):
    full = set(range(atom_count))
    return next((choices for choices in product(*signatures)
                 if set().union(*map(set, choices)) == full), None)


def iter_placed_candidate_items(items, covered_target_atoms):
    """Keep every fragment partition/attachment for assigned regions."""
    pools = [tuple(o for o in item["target_occupations"]
                   if tuple(o["covered_target_atoms"]) == tuple(atoms))
             for item, atoms in zip(items, covered_target_atoms, strict=True)]
    for choices in product(*pools):
        yield tuple(place_item(item, occupation)
                    for item, occupation in zip(items, choices, strict=True))


def place_candidate_items(items, covered_target_atoms):
    """A deterministic display witness; search uses the iterator above."""
    return next(iter_placed_candidate_items(tuple(items), covered_target_atoms))


@dataclass(frozen=True)
class CoverageRecommendationConfig:
    # Explicit caller constraint, never an inferred chemical assumption.
    maximum_precursors: int | None = None


@dataclass(frozen=True)
class CoverageRecommendationResult:
    patterns: tuple[CoveragePattern, ...]
    occupations: tuple[tuple[TargetOccupation, ...], ...]
    truncated: bool


def recommend_compressed_coverage_patterns(signatures, atom_count, rank_pattern,
        *, result_limit=None, config=None):
    config = config or CoverageRecommendationConfig()
    by_mask = {}
    for signature in sorted(set(signatures)):
        for occupation in signature:
            mask = sum(1 << atom for atom in occupation)
            by_mask.setdefault(mask, []).append(signature)
    graph = CoverageDecisionGraph.build(by_mask, atom_count)
    records = []
    for cover in graph.covers():
        if config.maximum_precursors is not None and len(cover) > config.maximum_precursors:
            continue
        atoms = tuple(tuple(a for a in range(atom_count) if mask & (1 << a)) for mask in cover)
        for pattern in product(*(by_mask[mask] for mask in cover)):
            records.append((pattern, atoms))
    records.sort(key=lambda r: (rank_pattern(r[0], atom_count), r))
    # Output limit only: the complete search is ranked before slicing.
    shown = records if result_limit is None else records[:result_limit]
    return CoverageRecommendationResult(tuple(r[0] for r in shown),
                                        tuple(r[1] for r in shown), False)
