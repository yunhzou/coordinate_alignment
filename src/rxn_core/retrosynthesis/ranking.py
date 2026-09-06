"""Pure structural ranking and assembly construction."""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction


def matched_fragment_count(item):
    """Matched units in one selected R copy, including singleton H units."""
    return len(item["retained_fragments"])


def precursor_cost(item):
    """Input cost per supported product copy, using proven fragment packing."""
    return (Fraction(item["total_atom_count"] * item["retained_atom_count"],
                     item.get("symmetry_retained_atom_count", item["retained_atom_count"])),
            item["total_atom_count"])


def candidate_entry_rank(entry):
    return (
        entry["chirality_violations"],
        not entry["complete"],
        -entry["symmetry_atom_retention"],
        -entry["atom_retention"],
        len(entry["boundary_bonds"]),
        entry["leftover_atom_count"],
        len(entry["leftover_fragments"]),
        entry["precursor_id"],
        entry["mapping"],
    )


def validate_atom_ownership(
        precursors, target_edges, require_attachment_bonds=False):
    """Validate coverage without assigning an arbitrary owner to overlaps.

    Return product connections unsupported by any matched fragment. These are
    structural assembly connections, NOT an atom-balanced reaction edit count.
    """
    support = defaultdict(set)
    carried = set()
    for index, precursor in enumerate(precursors):
        for atom in precursor["covered_target_atoms"]:
            support[atom].add(index)
        carried.update(tuple(sorted(edge)) for edge in precursor.get("preserved_target_bonds", ()))
    connections = []
    for left, right in target_edges:
        if not support[left] or not support[right]:
            return None
        if tuple(sorted((left, right))) in carried:
            continue
        if require_attachment_bonds:
            if not any(left in precursors[i]["attachment_atoms_target"] for i in support[left]):
                return None
            if not any(right in precursors[i]["attachment_atoms_target"] for i in support[right]):
                return None
        connections.append([left, right])
    return sorted(connections)


def build_ranked_assembly(items, formed_bonds):
    precursors = sorted(items, key=lambda item: item["precursor_id"])
    stoichiometry = Counter(item["precursor_id"] for item in precursors)
    retention_by_structure = defaultdict(list)
    for item in precursors:
        if item["total_heavy_atoms"]:
            retention_by_structure[item["structure_key"]].append(Fraction(
                item["retained_heavy_atoms"], item["total_heavy_atoms"]))
    unique_retentions = [
        min(values) for values in retention_by_structure.values()
    ]
    worst_retention = min(unique_retentions, default=None)
    mean_retention = (sum(unique_retentions, Fraction()) / len(unique_retentions)
                      if unique_retentions else None)
    heavy_total = sum(item["total_heavy_atoms"] for item in precursors)
    set_retention = (Fraction(
        sum(item["retained_heavy_atoms"] for item in precursors),
        heavy_total) if heavy_total else None)
    set_atom_retention = Fraction(
        sum(item["retained_atom_count"] for item in precursors),
        sum(item["total_atom_count"] for item in precursors),
    )
    set_symmetry_heavy_atom_retention = (Fraction(
        sum(item.get(
            "symmetry_retained_heavy_atoms", item["retained_heavy_atoms"])
            for item in precursors),
        heavy_total) if heavy_total else None)
    claimed_target_atoms = [
        atom for item in precursors
        for atom in item.get("covered_target_atoms", ())
    ]
    overlap_count = len(claimed_target_atoms) - len(set(claimed_target_atoms))
    if claimed_target_atoms:
        set_atom_retention = Fraction(len(set(claimed_target_atoms)),
                                      sum(item["total_atom_count"] for item in precursors))
    coverage_count = (len(set(claimed_target_atoms)) if claimed_target_atoms
                      else sum(item["retained_atom_count"] for item in precursors))
    set_symmetry_atom_retention = Fraction(coverage_count) / sum(
        (precursor_cost(item)[0] for item in precursors), Fraction())
    return {
        "precursors": precursors,
        "precursor_stoichiometry": dict(sorted(stoichiometry.items())),
        "formed_bonds": formed_bonds,
        "connection_semantics": "target bonds not carried by a matched fragment; not reaction edits",
        "target_support": {atom: [i for i, item in enumerate(precursors)
                                  if atom in item.get("covered_target_atoms", ())]
                           for atom in sorted(set(claimed_target_atoms))},
        "ownership_resolved": overlap_count == 0,
        "score": {
            "matched_fragment_count": sum(matched_fragment_count(item) for item in precursors),
            "chirality_violations": sum(
                item["chirality_violations"] for item in precursors),
            "unique_precursor_structures": len({item["structure_key"] for item in precursors}),
            "set_atom_retention": float(set_atom_retention),
            "set_symmetry_atom_retention": float(
                set_symmetry_atom_retention),
            "set_symmetry_heavy_atom_retention": float(
                set_symmetry_heavy_atom_retention) if set_symmetry_heavy_atom_retention is not None else None,
            "set_heavy_atom_retention": float(set_retention) if set_retention is not None else None,
            "worst_heavy_atom_retention": float(worst_retention) if worst_retention is not None else None,
            "mean_heavy_atom_retention": float(mean_retention) if mean_retention is not None else None,
            "capped_precursors": sum(
                not item["complete"] for item in precursors),
            "broken_bonds": sum(
                len(item["boundary_bonds"]) for item in precursors),
            "leftover_atoms": sum(
                item["leftover_atom_count"] for item in precursors),
            "formed_bonds": len(formed_bonds),
            "overlapping_target_atoms": overlap_count,
        },
    }


def assembly_rank(assembly):
    score = assembly["score"]
    items = assembly["precursors"]
    covered = {a for item in items for a in item.get("covered_target_atoms", ())}
    count = len(covered) if covered else sum(item["retained_atom_count"] for item in items)
    adjusted = sum((precursor_cost(item)[0] for item in items), Fraction())
    total = sum(item["total_atom_count"] for item in items)
    return (
        sum(matched_fragment_count(item) for item in items),
        score["unique_precursor_structures"],
        -Fraction(count) / adjusted,
        -Fraction(count, total),
        tuple(item["precursor_id"] for item in assembly["precursors"]),
    )
