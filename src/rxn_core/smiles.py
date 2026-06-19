"""Formal-bond-order graph adapters for SMILES/CXSMILES inputs."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .frag import WeightedGraph


@dataclass
class FormalWBOEndpoint:
    """One molecule parsed from SMILES as an element + formal-WBO graph.

    The ``wbo`` values are formal bond orders from the written molecular graph,
    not quantum Wiberg bond orders.  By default, explicit hydrogens in the
    input are kept and implicit hydrogens are not expanded.  Callers may request
    hydrogen expansion to materialize atom hydrogen counts as separate H nodes.
    """

    smiles: str
    elements: list[str]
    coords: np.ndarray
    wbo: np.ndarray
    atom_maps: dict[int, int]
    nodes: list[dict[str, Any]]
    hydrogen_policy: str = "preserve_explicit_only"


def _require_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
    except Exception as exc:  # pragma: no cover - depends on local install
        raise RuntimeError(
            "SMILES support requires RDKit. Install the optional chemistry "
            "dependency, then rerun the SMILES workflow."
        ) from exc
    return Chem, rdDepictor


def _mol_from_smiles(smiles, *, sanitize=True):
    Chem, _rdDepictor = _require_rdkit()
    params = Chem.SmilesParserParams()
    params.sanitize = bool(sanitize)
    params.removeHs = False
    mol = Chem.MolFromSmiles(str(smiles), params)
    if mol is None:
        raise ValueError(f"could not parse SMILES/CXSMILES: {smiles!r}")
    return mol


def _formal_bond_order(bond):
    value = float(bond.GetBondTypeAsDouble())
    if value > 0.0:
        return value
    # Keep unusual but explicit graph bonds active in the formal-WBO graph.
    return 1.0


def _compute_2d_coords(mol, *, component_spacing=3.0):
    Chem, rdDepictor = _require_rdkit()
    mol = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    coords = np.array([
        [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, 0.0]
        for i in range(mol.GetNumAtoms())
    ], dtype=float)

    fragments = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
    if len(fragments) <= 1:
        return coords

    # RDKit may place disconnected one-atom fragments close to a molecule.
    # Spread components for display only; the AAM consumes only element + WBO.
    out = coords.copy()
    cursor = 0.0
    for frag in sorted(fragments, key=lambda f: min(f)):
        idx = np.asarray(frag, dtype=int)
        frag_coords = coords[idx]
        min_x = float(np.min(frag_coords[:, 0]))
        max_x = float(np.max(frag_coords[:, 0]))
        out[idx, 0] += cursor - min_x
        cursor += max(0.0, max_x - min_x) + float(component_spacing)
    out[:, 0] -= float(np.mean(out[:, 0]))
    return out


def smiles_to_formal_wbo(smiles, *, sanitize=True,
                         component_spacing=3.0,
                         expand_hydrogens=False):
    """Parse SMILES/CXSMILES into a formal-bond-order weighted graph.

    Returns element labels, planar display coordinates, a square formal-WBO
    matrix, atom-map labels, and node records.  Atom-map labels are metadata
    only; they are not used as AAM constraints unless the caller separately
    supplies anchors.

    If ``expand_hydrogens`` is true, RDKit materializes atom hydrogen counts
    as explicit H atoms before coordinates and formal bond orders are built.
    """
    Chem, _rdDepictor = _require_rdkit()
    mol = _mol_from_smiles(smiles, sanitize=sanitize)
    hydrogen_policy = "expand_hydrogens" if expand_hydrogens else (
        "preserve_explicit_only")
    if expand_hydrogens:
        mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    elements = []
    nodes = []
    atom_maps = {}
    for i, atom in enumerate(mol.GetAtoms()):
        symbol = atom.GetSymbol()
        amap = int(atom.GetAtomMapNum())
        if amap:
            atom_maps[i] = amap
        features = {
            "formal_charge": int(atom.GetFormalCharge()),
            "isotope": int(atom.GetIsotope()),
            "radical_electrons": int(atom.GetNumRadicalElectrons()),
            "atom_map": amap if amap else None,
            "aromatic": bool(atom.GetIsAromatic()),
        }
        elements.append(symbol)
        nodes.append({
            "element": symbol,
            "label": str(amap) if amap else None,
            "features": features,
        })

    wbo = np.zeros((n, n), dtype=float)
    for bond in mol.GetBonds():
        i = int(bond.GetBeginAtomIdx())
        j = int(bond.GetEndAtomIdx())
        w = _formal_bond_order(bond)
        wbo[i, j] = wbo[j, i] = w

    return FormalWBOEndpoint(
        smiles=str(smiles),
        elements=elements,
        coords=_compute_2d_coords(mol, component_spacing=component_spacing),
        wbo=wbo,
        atom_maps=atom_maps,
        nodes=nodes,
        hydrogen_policy=hydrogen_policy,
    )


def smiles_to_weighted_graph(smiles, *, sanitize=True,
                             component_spacing=3.0,
                             expand_hydrogens=False):
    """Return a :class:`WeightedGraph` built from SMILES formal bond orders."""
    endpoint = smiles_to_formal_wbo(
        smiles, sanitize=sanitize, component_spacing=component_spacing,
        expand_hydrogens=expand_hydrogens)
    return WeightedGraph(
        nodes=endpoint.nodes,
        weights=endpoint.wbo,
        weight_name="wbo",
        coords=endpoint.coords,
        metadata={
            "source": "smiles",
            "smiles": endpoint.smiles,
            "wbo_kind": "formal_bond_order",
            "hydrogen_policy": endpoint.hydrogen_policy,
            "atom_maps": endpoint.atom_maps,
        },
    )
