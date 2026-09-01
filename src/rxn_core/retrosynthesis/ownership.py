"""Joint target-atom ownership for overlapping precursor mappings."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from rdkit import Chem

from .catalog_index import assign_owned_target_atoms
from .ranking import (
    assembly_rank,
    build_ranked_assembly,
    validate_atom_ownership,
)


@dataclass(frozen=True)
class OwnershipResolution:
    assemblies: tuple[dict, ...]
    truncated: bool


def resolve_overlapping_ownership(
        items, atom_count, target_edges, *, beam_width=200,
        assembly_limit=4, require_attachment_bonds=False):
    """Partition a full union of mapped atoms and rank the exact owners."""
    available = [set(map(int, item["covered_target_atoms"])) for item in items]
    choices = [
        tuple(index for index, atoms in enumerate(available) if atom in atoms)
        for atom in range(atom_count)
    ]
    if any(not owners for owners in choices):
        return OwnershipResolution((), False)

    neighbors = [[] for _ in range(atom_count)]
    for left, right in target_edges:
        neighbors[int(left)].append(int(right))
        neighbors[int(right)].append(int(left))
    order = sorted(
        range(atom_count),
        key=lambda atom: (len(choices[atom]), -len(neighbors[atom]), atom),
    )
    states = [({}, frozenset(), 0)]
    truncated = False
    for atom in order:
        expanded = []
        for owners, used, formed in states:
            for owner in choices[atom]:
                new_formed = formed + sum(
                    owners[neighbor] != owner
                    for neighbor in neighbors[atom] if neighbor in owners)
                new_owners = dict(owners)
                new_owners[atom] = owner
                expanded.append((new_owners, used | {owner}, new_formed))
        if len(expanded) > beam_width:
            truncated = True
        expanded.sort(key=lambda state: (
            len({items[index]["structure_key"] for index in state[1]}),
            sum(items[index]["total_atom_count"] for index in state[1]),
            state[2],
            len(state[1]),
        ))
        states = expanded[:beam_width]

    source_edges = []
    inverse_mappings = []
    for item in items:
        molecule = Chem.AddHs(Chem.MolFromSmiles(item["smiles"]))
        source_edges.append(tuple(
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            for bond in molecule.GetBonds()
        ))
        inverse_mappings.append({
            int(target): int(source) for source, target in item["mapping"]
        })

    def final_rank(state):
        owners, used, formed = state
        retained_by_item = [set() for _item in items]
        for target_atom, owner in owners.items():
            retained_by_item[owner].add(inverse_mappings[owner][target_atom])
        broken = sum(
            (left in retained_by_item[index])
            != (right in retained_by_item[index])
            for index in used for left, right in source_edges[index]
        )
        total_atoms = sum(items[index]["total_atom_count"] for index in used)
        chirality = sum(
            len(set(items[index]["chirality_violation_target_atoms"])
                & {atom for atom, owner in owners.items() if owner == index})
            for index in used
        )
        return (
            chirality,
            len({items[index]["structure_key"] for index in used}),
            -Fraction(atom_count, total_atoms),
            sum(not items[index]["complete"] for index in used),
            broken,
            total_atoms - atom_count,
            formed,
        )

    states.sort(key=final_rank)
    assemblies = []
    for owners, used, _formed in states[:assembly_limit]:
        rebuilt = []
        for index in sorted(used):
            owned = [atom for atom, owner in owners.items() if owner == index]
            rebuilt.append(assign_owned_target_atoms(items[index], owned))
        formed = validate_atom_ownership(
            rebuilt, target_edges, require_attachment_bonds)
        if formed is not None:
            assemblies.append(build_ranked_assembly(rebuilt, formed))
    assemblies.sort(key=assembly_rank)
    return OwnershipResolution(tuple(assemblies[:assembly_limit]), truncated)
