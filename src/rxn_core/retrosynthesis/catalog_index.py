"""Build a bounded coverage index from persisted fragment-search records."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache

import networkx as nx
from rdkit import Chem

from .ranking import candidate_entry_rank
from ..fragment_matching.serialization import FRAGMENT_DETECTION_SCHEMA


def coverage_mask(atoms):
    value = 0
    for atom in atoms:
        value |= 1 << int(atom)
    return value


def _cip(atom):
    return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None


def chirality_violations(candidate, source, target):
    return len(chirality_violation_target_atoms(
        candidate, source, target))


def chirality_violation_target_atoms(candidate, source, target):
    reaction_center = {int(edge[0]) for edge in candidate["boundary_bonds"]}
    return sorted(
        int(target_atom)
        for source_atom, target_atom in candidate["mapping"]
        if int(source_atom) not in reaction_center
        and ((_cip(source.GetAtomWithIdx(int(source_atom))) is not None
              or _cip(target.GetAtomWithIdx(int(target_atom))) is not None)
             and _cip(source.GetAtomWithIdx(int(source_atom)))
             != _cip(target.GetAtomWithIdx(int(target_atom))))
    )


def candidate_entry(
        record, candidate, explicit_molecule, structure_key, target,
        chirality_ranking):
    leftovers = candidate["leftover_fragments"]
    retained_heavy_atoms = sum(
        explicit_molecule.GetAtomWithIdx(int(atom)).GetAtomicNum() > 1
        for atom in candidate["retained_atoms"]
    )
    retained_atom_count = len(candidate["retained_atoms"])
    total_atom_count = explicit_molecule.GetNumAtoms()
    total_heavy_atoms = explicit_molecule.GetNumHeavyAtoms()
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
        "attachment_atoms_target": candidate["attachment_atoms_target"],
        "mapping": candidate["mapping"],
        "retained_fragments": candidate["retained_fragments"],
        "leftover_atom_count": sum(map(len, leftovers)),
        "structure_key": structure_key,
        "retained_heavy_atoms": retained_heavy_atoms,
        "total_heavy_atoms": total_heavy_atoms,
        "heavy_atom_retention": retained_heavy_atoms / total_heavy_atoms,
        "retained_atom_count": retained_atom_count,
        "total_atom_count": total_atom_count,
        "atom_retention": retained_atom_count / total_atom_count,
        "attachment_trimmed_target_atoms": candidate[
            "attachment_trimmed_target_atoms"],
        "chirality_violations": (
            chirality_violations(candidate, explicit_molecule, target)
            if chirality_ranking else 0),
        "chirality_violation_target_atoms": (
            chirality_violation_target_atoms(
                candidate, explicit_molecule, target)
            if chirality_ranking else []),
    }


@lru_cache(maxsize=None)
def _explicit_molecule(smiles):
    return Chem.AddHs(Chem.MolFromSmiles(smiles))


def assign_owned_target_atoms(entry, owned_target_atoms):
    """Rebuild an entry after joint assembly assigns its target atoms."""
    owned = set(map(int, owned_target_atoms))
    molecule = _explicit_molecule(entry["smiles"])
    mapping = {
        int(source): int(target) for source, target in entry["mapping"]
        if int(target) in owned
    }
    retained = set(mapping)
    outside = set(range(molecule.GetNumAtoms())) - retained
    retained_graph = nx.Graph()
    retained_graph.add_nodes_from(retained)
    outside_graph = nx.Graph()
    outside_graph.add_nodes_from(outside)
    boundary = []
    attachments = set()
    for bond in molecule.GetBonds():
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if left in retained and right in retained:
            retained_graph.add_edge(left, right)
        elif left in outside and right in outside:
            outside_graph.add_edge(left, right)
        elif left in retained:
            boundary.append([left, right])
            attachments.add(left)
        else:
            boundary.append([right, left])
            attachments.add(right)
    attachments.update(
        atom for atom in retained
        if molecule.GetAtomWithIdx(atom).GetFormalCharge() != 0)
    retained_fragments = [
        sorted(component)
        for component in nx.connected_components(retained_graph)
    ]
    leftovers = [
        sorted(component)
        for component in nx.connected_components(outside_graph)
    ]
    retained_heavy_atoms = sum(
        molecule.GetAtomWithIdx(atom).GetAtomicNum() > 1 for atom in retained)
    removed_targets = (
        set(map(int, entry["covered_target_atoms"])) - owned)
    rebuilt = dict(entry)
    rebuilt.update({
        "mapping": [list(item) for item in sorted(mapping.items())],
        "retained_atoms": sorted(retained),
        "covered_target_atoms": sorted(owned),
        "retained_fragments": retained_fragments,
        "leftover_fragments": leftovers,
        "boundary_bonds": sorted(boundary),
        "attachment_atoms_target": sorted(
            mapping[atom] for atom in attachments if atom in mapping),
        "leftover_atom_count": len(outside),
        "retained_heavy_atoms": retained_heavy_atoms,
        "retained_atom_count": len(retained),
        "heavy_atom_retention": (
            retained_heavy_atoms / entry["total_heavy_atoms"]),
        "atom_retention": len(retained) / entry["total_atom_count"],
        "attachment_trimmed_target_atoms": sorted(removed_targets),
        "chirality_violations": len(
            owned & set(entry["chirality_violation_target_atoms"])),
    })
    return rebuilt


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
    per_mask_limit: int = 200
    exclude_target_identity: bool = True
    attachment_trim_variants: bool = False
    chirality_ranking: bool = False
    expected_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class CandidateIndex:
    groups: dict[int, list[dict]]
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
        if config.chirality_ranking:
            Chem.AssignStereochemistry(
                explicit_molecule, cleanIt=True, force=True)
        candidates = record["candidates"]
        if config.attachment_trim_variants:
            candidates = [
                variant
                for candidate in candidates
                for variant in attachment_trim_variants(
                    candidate, explicit_molecule)
            ]
        else:
            candidates = [
                dict(candidate, attachment_trimmed_target_atoms=[])
                for candidate in candidates
            ]
        for candidate in candidates:
            item = candidate_entry(
                record,
                candidate,
                explicit_molecule,
                molecule_key,
                target,
                config.chirality_ranking,
            )
            mask = coverage_mask(item["covered_target_atoms"])
            bucket = groups[mask]
            bucket.append(item)
            if len(bucket) >= config.per_mask_limit * 2:
                bucket.sort(key=candidate_entry_rank)
                del bucket[config.per_mask_limit:]
            counts["fragment_candidates"] += 1
            if item["precursor_id"] in config.expected_ids:
                expected[item["precursor_id"]].append(item)

    ranked_groups = {
        mask: sorted(items, key=candidate_entry_rank)[:config.per_mask_limit]
        for mask, items in groups.items()
    }
    return CandidateIndex(
        groups=ranked_groups,
        expected=dict(expected),
        direct_target_matches=tuple(direct_target_matches),
        counts=dict(counts),
    )
