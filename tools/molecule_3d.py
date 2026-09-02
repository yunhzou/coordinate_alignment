"""Shared RDKit preparation for standalone molecular viewers."""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem


def mol_3d(smiles, *, spread_ions=False, show_hydrogens=False,
           cut_bonds=()):
    molecule = Chem.MolFromSmiles(smiles)
    with_hydrogens = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 20260827
    if AllChem.EmbedMolecule(with_hydrogens, parameters) != 0:
        raise RuntimeError(f"could not embed {smiles!r}")
    try:
        AllChem.UFFOptimizeMolecule(with_hydrogens, maxIters=300)
    except Exception:
        pass
    if cut_bonds:
        editable = Chem.RWMol(with_hydrogens)
        for left, right in cut_bonds:
            if editable.GetBondBetweenAtoms(int(left), int(right)) is not None:
                editable.RemoveBond(int(left), int(right))
        with_hydrogens = editable.GetMol()
        with_hydrogens.UpdatePropertyCache(strict=False)
    molecule = (with_hydrogens if show_hydrogens
                else Chem.RemoveHs(with_hydrogens))
    conformer = molecule.GetConformer()
    if spread_ions:
        fragments = Chem.GetMolFrags(
            molecule, asMols=False, sanitizeFrags=False)
        cursor = 0.0
        gap = 8.0
        for fragment_index, fragment in enumerate(fragments):
            x_values = [
                conformer.GetAtomPosition(atom_index).x
                for atom_index in fragment
            ]
            shift = cursor - min(x_values)
            for atom_index in fragment:
                point = conformer.GetAtomPosition(atom_index)
                point.x += shift
                point.y += 1.5 * (-1) ** fragment_index
                conformer.SetAtomPosition(atom_index, point)
            cursor += max(x_values) - min(x_values) + gap
        center_shift = (cursor - gap) / 2.0
        for atom_index in range(molecule.GetNumAtoms()):
            point = conformer.GetAtomPosition(atom_index)
            point.x -= center_shift
            conformer.SetAtomPosition(atom_index, point)
    coordinates = [
        [
            float(conformer.GetAtomPosition(index).x),
            float(conformer.GetAtomPosition(index).y),
            float(conformer.GetAtomPosition(index).z),
        ]
        for index in range(molecule.GetNumAtoms())
    ]
    return Chem.MolToMolBlock(molecule), coordinates, [
        atom.GetSymbol() for atom in molecule.GetAtoms()
    ]
