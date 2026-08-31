"""RDKit conversion for fragment detection."""
from __future__ import annotations

import numpy as np

from ..frag import WeightedGraph


def molecule_to_weighted_graph(molecule):
    size = molecule.GetNumAtoms()
    matrix = np.zeros((size, size), dtype=float)
    nodes = []
    for atom in molecule.GetAtoms():
        nodes.append({
            "element": atom.GetSymbol(),
            "features": {
                "formal_charge": int(atom.GetFormalCharge()),
                "aromatic": bool(atom.GetIsAromatic()),
            },
        })
    for bond in molecule.GetBonds():
        left = int(bond.GetBeginAtomIdx())
        right = int(bond.GetEndAtomIdx())
        weight = float(bond.GetBondTypeAsDouble()) or 1.0
        matrix[left, right] = matrix[right, left] = weight
    return WeightedGraph(nodes, matrix)
