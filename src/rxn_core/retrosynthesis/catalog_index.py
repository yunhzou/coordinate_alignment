"""Build a bounded coverage index from persisted fragment-search records."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass

import networkx as nx
from rdkit import Chem

from ..fragment_matching.rdkit_adapter import molecule_to_weighted_graph
from ..fragment_matching.serialization import (
    FRAGMENT_DETECTION_SCHEMA,
    fragment_candidate_from_record,
    fragment_candidate_to_record,
)
from ..matcher import _nauty_orbits
from ..subgraph import _coerce_graph
from .compressed_coverage import (
    CoverageSignature,
    candidate_target_domains,
    coverage_signature,
)
from .ranking import candidate_entry_rank


def _cip(atom):
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None


def chirality_violations(candidate, source, target):
    reaction_center = {int(edge[0]) for edge in candidate["boundary_bonds"]}
    violations = 0
    for source_atom, target_atom in candidate["mapping"]:
        if int(source_atom) in reaction_center:
            continue
        source_cip = _cip(source.GetAtomWithIdx(int(source_atom)))
        target_cip = _cip(target.GetAtomWithIdx(int(target_atom)))
        if ((source_cip is not None or target_cip is not None)
                and source_cip != target_cip):
            violations += 1
    return violations


def _symmetry_copy_capacity(retained_atoms, source_orbits):
    """Return how many disjoint retained atom inventories symmetry permits.

    A source automorphism orbit is an exact pool of interchangeable atoms.  A
    repeated copy of a retained fragment consumes the same number of atoms
    from every orbit as the observed copy.  The smallest orbit capacity is
    therefore the structural upper bound on symmetric product copies.
    """
    orbit_sizes = Counter(source_orbits.values())
    retained_orbits = Counter(
        source_orbits[int(atom)] for atom in retained_atoms)
    return min(
        orbit_sizes[orbit] // count
        for orbit, count in retained_orbits.items()
    )


def candidate_entry(
        record, candidate, explicit_molecule, structure_key, target,
        chirality_ranking, source_orbits, target_domains):
    leftovers = candidate["leftover_fragments"]
    retained_heavy_atoms = sum(
        explicit_molecule.GetAtomWithIdx(int(atom)).GetAtomicNum() > 1
        for atom in candidate["retained_atoms"]
    )
    retained_atom_count = len(candidate["retained_atoms"])
    total_atom_count = explicit_molecule.GetNumAtoms()
    total_heavy_atoms = explicit_molecule.GetNumHeavyAtoms()
    symmetry_copy_capacity = _symmetry_copy_capacity(
        candidate["retained_atoms"], source_orbits)
    retained_orbits = {
        source_orbits[int(atom)] for atom in candidate["retained_atoms"]
    }
    symmetry_retained_atom_indices = sorted(
        int(atom) for atom, orbit in source_orbits.items()
        if orbit in retained_orbits
    )
    symmetry_retained_atoms = min(
        total_atom_count, retained_atom_count * symmetry_copy_capacity)
    symmetry_retained_heavy_atoms = min(
        total_heavy_atoms, retained_heavy_atoms * symmetry_copy_capacity)
    return {
        "precursor_id": record["source_id"],
        "smiles": record["representation"],
        "row_index": record["row_index"],
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
        "target_domains": target_domains,
        "retained_fragments": candidate["retained_fragments"],
        "leftover_atom_count": sum(map(len, leftovers)),
        "structure_key": structure_key,
        "retained_heavy_atoms": retained_heavy_atoms,
        "total_heavy_atoms": total_heavy_atoms,
        "heavy_atom_retention": retained_heavy_atoms / total_heavy_atoms,
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
            symmetry_retained_heavy_atoms / total_heavy_atoms),
        "attachment_trimmed_target_atoms": candidate[
            "attachment_trimmed_target_atoms"],
        "chirality_violations": (
            chirality_violations(candidate, explicit_molecule, target)
            if chirality_ranking else 0),
    }


def attachment_trim_variants(candidate, molecule):
    """Yield connected variants that relinquish one mapped attachment atom."""
    base = dict(candidate)
    base["attachment_trimmed_target_atoms"] = []
    yield base
    mapping = {
        int(source): int(target) for source, target in candidate["mapping"]
    }
    retained = {int(atom) for atom in candidate["retained_atoms"]}
    attachment_targets = set(map(
        int, candidate["attachment_atoms_target"]))
    inverse = {target: source for source, target in mapping.items()}
    seen = set()
    for target_atom in sorted(attachment_targets):
        source_atom = inverse.get(target_atom)
        if source_atom is None or source_atom not in retained:
            continue
        new_retained = retained - {source_atom}
        if not new_retained:
            continue
        graph = nx.Graph()
        graph.add_nodes_from(new_retained)
        for bond in molecule.GetBonds():
            left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if left in new_retained and right in new_retained:
                graph.add_edge(left, right)
        if not nx.is_connected(graph):
            continue
        new_mapping = {
            source: target for source, target in mapping.items()
            if source in new_retained
        }
        leftover = set(range(molecule.GetNumAtoms())) - new_retained
        leftover_graph = nx.Graph()
        leftover_graph.add_nodes_from(leftover)
        boundary = []
        attachment_r = set()
        for bond in molecule.GetBonds():
            left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if left in leftover and right in leftover:
                leftover_graph.add_edge(left, right)
            elif left in new_retained and right in leftover:
                boundary.append([left, right])
                attachment_r.add(left)
            elif right in new_retained and left in leftover:
                boundary.append([right, left])
                attachment_r.add(right)
        if not attachment_r:
            continue
        key = tuple(sorted(new_mapping.items()))
        if key in seen:
            continue
        seen.add(key)
        variant = dict(candidate)
        variant.update({
            "mapping": [list(item) for item in sorted(new_mapping.items())],
            "retained_atoms": sorted(new_retained),
            "covered_target_atoms": sorted(new_mapping.values()),
            "leftover_fragments": [
                sorted(component)
                for component in nx.connected_components(leftover_graph)
            ],
            "boundary_bonds": sorted(boundary),
            "attachment_atoms_source": sorted(attachment_r),
            "attachment_atoms_target": sorted(
                new_mapping[atom] for atom in attachment_r),
            "attachment_trimmed_target_atoms": [target_atom],
        })
        yield variant


@dataclass(frozen=True)
class CandidateIndexConfig:
    per_domain_limit: int = 200
    exclude_target_identity: bool = True
    attachment_trim_variants: bool = False
    chirality_ranking: bool = False
    expected_ids: tuple[str, ...] = ()
    iso_tolerance: float = 0.5

    def __post_init__(self):
        if self.per_domain_limit < 1:
            raise ValueError("candidate index limits must be positive")
        if self.iso_tolerance <= 0:
            raise ValueError("isomorphism tolerance must be positive")


@dataclass(frozen=True)
class CandidateIndex:
    groups: dict[CoverageSignature, list[dict]]
    expected: dict[str, list[dict]]
    direct_target_matches: tuple[dict, ...]
    counts: dict[str, int]


def build_candidate_index(records, target, target_key, *, config=None):
    config = config or CandidateIndexConfig()
    groups = defaultdict(list)
    expected = defaultdict(list)
    direct_target_matches = []
    counts = Counter()
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

        counts["matched_precursors"] += 1
        counts["capped_precursors"] += record["status"] == "capped"
        explicit_molecule = Chem.AddHs(molecule)
        source_graph = _coerce_graph(
            molecule_to_weighted_graph(explicit_molecule), 0.2)
        source_orbits = _nauty_orbits(
            source_graph, wbo_tol=config.iso_tolerance)
        if config.chirality_ranking:
            Chem.AssignStereochemistry(
                explicit_molecule, cleanIt=True, force=True)
        candidates = []
        for raw_candidate in record["candidates"]:
            typed = fragment_candidate_from_record(dict(
                raw_candidate, source_id=record["source_id"]))
            candidates.append((
                fragment_candidate_to_record(typed),
                candidate_target_domains(typed),
            ))
        if config.attachment_trim_variants:
            candidates = [
                (variant, candidate_target_domains(
                    fragment_candidate_from_record(dict(
                        variant, source_id=record["source_id"]))))
                for candidate, _domains in candidates
                for variant in attachment_trim_variants(
                    candidate, explicit_molecule)
            ]
        else:
            candidates = [
                (dict(candidate, attachment_trimmed_target_atoms=[]), domains)
                for candidate, domains in candidates
            ]
        for candidate, target_domains in candidates:
            item = candidate_entry(
                record,
                candidate,
                explicit_molecule,
                molecule_key,
                target,
                config.chirality_ranking,
                source_orbits,
                target_domains,
            )
            signature = coverage_signature(target_domains)
            bucket = groups[signature]
            bucket.append(item)
            if len(bucket) >= config.per_domain_limit * 2:
                bucket.sort(key=candidate_entry_rank)
                del bucket[config.per_domain_limit:]
            counts["fragment_candidates"] += 1
            if item["precursor_id"] in config.expected_ids:
                expected[item["precursor_id"]].append(item)

    ranked_groups = {
        signature: sorted(items, key=candidate_entry_rank)[
            :config.per_domain_limit]
        for signature, items in groups.items()
    }
    return CandidateIndex(
        groups=ranked_groups,
        expected=dict(expected),
        direct_target_matches=tuple(direct_target_matches),
        counts=dict(counts),
    )
