"""Pure structural ranking and assembly construction."""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction


def candidate_entry_rank(entry):
    return (
        entry["chirality_violations"],
        not entry["complete"],
        -entry["atom_retention"],
        len(entry["boundary_bonds"]),
        entry["leftover_atom_count"],
        len(entry["leftover_fragments"]),
        entry["precursor_id"],
        entry["mapping"],
    )


def validate_atom_ownership(
        precursors, target_edges, require_attachment_bonds=True):
    """Return cross-owner product bonds, or None for an invalid ownership."""
    owner = {}
    attachments = []
    for index, precursor in enumerate(precursors):
        attachments.append(set(precursor["attachment_atoms_target"]))
        for atom in precursor["covered_target_atoms"]:
            if atom in owner:
                return None
            owner[atom] = index
    formed = []
    for atom_a, atom_b in target_edges:
        owner_a, owner_b = owner.get(atom_a), owner.get(atom_b)
        if owner_a is None or owner_b is None:
            return None
        if owner_a == owner_b:
            continue
        if (require_attachment_bonds
                and (atom_a not in attachments[owner_a]
                     or atom_b not in attachments[owner_b])):
            return None
        formed.append([atom_a, atom_b])
    return formed


def build_ranked_assembly(items, formed_bonds):
    precursors = sorted(items, key=lambda item: item["precursor_id"])
    stoichiometry = Counter(item["precursor_id"] for item in precursors)
    retention_by_structure = defaultdict(list)
    for item in precursors:
        retention_by_structure[item["structure_key"]].append(Fraction(
            item["retained_heavy_atoms"], item["total_heavy_atoms"]))
    unique_retentions = [
        min(values) for values in retention_by_structure.values()
    ]
    worst_retention = min(unique_retentions)
    mean_retention = sum(
        unique_retentions, Fraction()) / len(unique_retentions)
    set_retention = Fraction(
        sum(item["retained_heavy_atoms"] for item in precursors),
        sum(item["total_heavy_atoms"] for item in precursors),
    )
    set_atom_retention = Fraction(
        sum(item["retained_atom_count"] for item in precursors),
        sum(item["total_atom_count"] for item in precursors),
    )
    return {
        "precursors": precursors,
        "precursor_stoichiometry": dict(sorted(stoichiometry.items())),
        "formed_bonds": formed_bonds,
        "score": {
            "chirality_violations": sum(
                item["chirality_violations"] for item in precursors),
            "unique_precursor_structures": len(retention_by_structure),
            "set_atom_retention": float(set_atom_retention),
            "set_heavy_atom_retention": float(set_retention),
            "worst_heavy_atom_retention": float(worst_retention),
            "mean_heavy_atom_retention": float(mean_retention),
            "capped_precursors": sum(
                not item["complete"] for item in precursors),
            "broken_bonds": sum(
                len(item["boundary_bonds"]) for item in precursors),
            "leftover_atoms": sum(
                item["leftover_atom_count"] for item in precursors),
            "formed_bonds": len(formed_bonds),
        },
    }


def assembly_rank(assembly):
    score = assembly["score"]
    return (
        score["chirality_violations"],
        score["unique_precursor_structures"],
        -score["set_atom_retention"],
        score["capped_precursors"],
        score["broken_bonds"],
        score["leftover_atoms"],
        score["formed_bonds"],
        tuple(item["precursor_id"] for item in assembly["precursors"]),
    )
