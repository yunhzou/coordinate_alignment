"""Lossless coverage index from persisted fragment-search records."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass

from rdkit import Chem

from ..fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from ..fragment_matching.serialization import (
    FRAGMENT_DETECTION_SCHEMA,
    fragment_candidate_from_record,
    fragment_archive_from_record,
)
from ..matcher import _nauty_atom_generators
from ..subgraph import _coerce_graph
from .compressed_coverage import (
    CoverageSignature,
    candidate_target_domains,
    candidate_target_occupations,
    coverage_signature,
)
from .ranking import candidate_entry_rank


def candidate_entry(
        record, candidate, explicit_molecule, structure_key, target_domains,
        target_occupations, copy_capacity):
    leftovers = candidate["leftover_fragments"]
    retained_heavy_atoms = sum(
        explicit_molecule.GetAtomWithIdx(int(atom)).GetAtomicNum() > 1
        for atom in candidate["retained_atoms"]
    )
    retained_atom_count = len(candidate["retained_atoms"])
    total_atom_count = explicit_molecule.GetNumAtoms()
    total_heavy_atoms = explicit_molecule.GetNumHeavyAtoms()
    symmetry_copy_capacity, symmetry_retained_atom_indices = copy_capacity
    symmetry_retained_atoms = min(
        total_atom_count, retained_atom_count * symmetry_copy_capacity)
    symmetry_retained_heavy_atoms = min(
        total_heavy_atoms, retained_heavy_atoms * symmetry_copy_capacity)
    return {
        "precursor_id": record["source_id"],
        "smiles": record["representation"],
        "row_index": record["row_index"],
        "detection_reference": {
            "row_index": record["row_index"],
            "candidate_index": candidate.get("detection_candidate_index"),
        },
        "complete": record["complete"],
        "status": record["status"],
        "best_fragment_size": record["best_fragment_size"],
        "covered_target_atoms": candidate["covered_target_atoms"],
        "retained_atoms": candidate["retained_atoms"],
        "leftover_fragments": leftovers,
        "boundary_bonds": candidate["boundary_bonds"],
        "attachment_atoms_source": candidate["attachment_atoms_source"],
        "attachment_atoms_target": candidate["attachment_atoms_target"],
        "mapping": candidate["mapping"],
        "preserved_source_bonds": candidate["preserved_source_bonds"],
        "target_domains": target_domains,
        "target_occupations": target_occupations,
        "retained_fragments": candidate["retained_fragments"],
        "leftover_atom_count": sum(map(len, leftovers)),
        "structure_key": structure_key,
        "retained_heavy_atoms": retained_heavy_atoms,
        "total_heavy_atoms": total_heavy_atoms,
        "heavy_atom_retention": retained_heavy_atoms / total_heavy_atoms if total_heavy_atoms else None,
        "retained_atom_count": retained_atom_count,
        "total_atom_count": total_atom_count,
        "atom_retention": retained_atom_count / total_atom_count,
        "symmetry_copy_capacity": symmetry_copy_capacity,
        "symmetry_retained_atoms": symmetry_retained_atom_indices,
        "symmetry_retained_atom_count": symmetry_retained_atoms,
        "symmetry_retained_heavy_atoms": symmetry_retained_heavy_atoms,
        "symmetry_atom_retention": (
            symmetry_retained_atoms / total_atom_count),
        "symmetry_heavy_atom_retention": (
            symmetry_retained_heavy_atoms / total_heavy_atoms if total_heavy_atoms else None),
        "attachment_trimmed_target_atoms": candidate[
            "attachment_trimmed_target_atoms"],
        "chirality_violations": 0,
    }


def exact_source_copy_capacity(retained_atoms, generators):
    """Maximum disjoint whole-fragment copies, not independent atom capacity."""
    initial = frozenset(map(int, retained_atoms))
    regions, queue = {initial}, [initial]
    for region in queue:
        for generator in generators:
            image = frozenset(generator.get(a, a) for a in region)
            if image not in regions:
                regions.add(image)
                queue.append(image)
    masks = tuple(sum(1 << a for a in region) for region in regions)
    from functools import lru_cache
    @lru_cache(None)
    def pack(remaining):
        if not remaining:
            return 0
        if all(not (a & b) for i, a in enumerate(remaining) for b in remaining[i + 1:]):
            return len(remaining)
        first, rest = remaining[0], remaining[1:]
        return max(pack(rest), 1 + pack(tuple(m for m in rest if not m & first)))
    return pack(tuple(sorted(masks))), sorted(set().union(*regions))


@dataclass(frozen=True)
class CandidateIndexConfig:
    exclude_target_identity: bool = False
    expected_ids: tuple[str, ...] = ()
    iso_tolerance: float = 0.5

    def __post_init__(self):
        if self.iso_tolerance <= 0:
            raise ValueError("isomorphism tolerance must be positive")


@dataclass(frozen=True)
class CandidateIndex:
    groups: dict[CoverageSignature, list[dict]]
    expected: dict[str, list[dict]]
    direct_target_matches: tuple[dict, ...]
    counts: dict[str, int]


def merge_candidate_indexes(indexes):
    """Reduce independent shard indexes into one deterministic index."""
    groups = defaultdict(list)
    expected = defaultdict(list)
    direct_target_matches = []
    counts = Counter()
    for index in indexes:
        for signature, items in index.groups.items():
            groups[signature].extend(items)
        for precursor_id, items in index.expected.items():
            expected[precursor_id].extend(items)
        direct_target_matches.extend(index.direct_target_matches)
        counts.update(index.counts)
    return CandidateIndex(
        groups={
            signature: sorted(items, key=candidate_entry_rank)
            for signature, items in groups.items()
        },
        expected={
            precursor_id: sorted(items, key=candidate_entry_rank)
            for precursor_id, items in expected.items()
        },
        direct_target_matches=tuple(direct_target_matches),
        counts=dict(counts),
    )


def build_candidate_index(records, target, target_key, *, config=None):
    config = config or CandidateIndexConfig()
    groups = defaultdict(list)
    expected = defaultdict(list)
    direct_target_matches = []
    counts = Counter()
    target_graph = _coerce_graph(molecule_to_weighted_graph(target), 0.2)
    target_generators = _nauty_atom_generators(
        target_graph, wbo_tol=config.iso_tolerance)
    for record in records:
        if record["schema"] != FRAGMENT_DETECTION_SCHEMA:
            raise ValueError(
                f"unsupported fragment record schema: {record['schema']!r}")
        molecule = Chem.MolFromSmiles(record["representation"])
        molecule_key = (
            Chem.MolToSmiles(molecule, isomericSmiles=True)
            if molecule is not None else record["source_id"]
        )
        if (config.exclude_target_identity and molecule is not None
                and molecule_key == target_key):
            direct_target_matches.append({
                "precursor_id": record["source_id"],
                "smiles": record["representation"],
                "row_index": record["row_index"],
            })
            continue
        if molecule is None:
            raise ValueError(
                "invalid persisted source representation: "
                f"{record['representation']!r}")

        counts["matched_precursors"] += bool(record["candidates"])
        counts["capped_precursors"] += record["status"] == "capped"
        counts["incomplete_precursors"] += not record["complete"]
        explicit_molecule = Chem.AddHs(molecule)
        source_graph = _coerce_graph(
            molecule_to_weighted_graph(explicit_molecule), 0.2)
        source_generators = _nauty_atom_generators(source_graph, wbo_tol=config.iso_tolerance)
        source_capacity_cache = {}
        search_graphs, hierarchy_fragments = fragment_archive_from_record(record)
        candidates = []
        for candidate_index, raw_candidate in enumerate(record["candidates"]):
            typed = fragment_candidate_from_record(dict(
                raw_candidate, source_id=record["source_id"]),
                search_graphs=search_graphs, hierarchy_fragments=hierarchy_fragments)
            candidates.append((
                dict(raw_candidate, detection_candidate_index=candidate_index),
                candidate_target_domains(typed),
                candidate_target_occupations(
                    typed, target_graph,
                    iso_tolerance=config.iso_tolerance,
                    generators=target_generators),
            ))
        candidates = [
            (dict(candidate, attachment_trimmed_target_atoms=[]), domains, occupations)
            for candidate, domains, occupations in candidates
        ]
        for candidate, target_domains, target_occupations in candidates:
            retained = tuple(candidate["retained_atoms"])
            if retained not in source_capacity_cache:
                source_capacity_cache[retained] = exact_source_copy_capacity(retained, source_generators)
            item = candidate_entry(
                record,
                candidate,
                explicit_molecule,
                molecule_key,
                target_domains,
                target_occupations,
                source_capacity_cache[retained],
            )
            signature = coverage_signature(target_occupations)
            bucket = groups[signature]
            bucket.append(item)
            counts["fragment_candidates"] += 1
            if item["precursor_id"] in config.expected_ids:
                expected[item["precursor_id"]].append(item)

    ranked_groups = {
        signature: sorted(items, key=candidate_entry_rank)
        for signature, items in groups.items()
    }
    return CandidateIndex(
        groups=ranked_groups,
        expected=dict(expected),
        direct_target_matches=tuple(direct_target_matches),
        counts=dict(counts),
    )
