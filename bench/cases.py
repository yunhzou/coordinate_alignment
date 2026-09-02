"""Benchmark inputs for the Layer-0 replay harness.

``tempo`` is the committed 57-atom xtb example.  ``tetraphenyl`` is a synthetic
symmetry-heavy case built from formal bond orders: tetraphenylmethane (45 atoms,
four equivalent phenyl rings, 20 equivalent-looking hydrogens) undergoing a
1,2-hydrogen shift on one ring.  Coordinates are deterministic pseudo-random
and irrelevant to AAM search (they enter only the geometric post-processing).
"""
from __future__ import annotations

import random
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def tempo():
    from rxn_core.cli import _endpoint_cache
    base = ROOT / "docs/example_runs/pr1.tempo_ts3/work/endpoints"
    return _endpoint_cache(str(base / "R"), "R"), _endpoint_cache(str(base / "P"), "P")


def _tetraphenylmethane():
    elements = ["C"]  # central carbon 0
    bonds = []
    ring_atoms = []
    for _ring in range(4):
        base = len(elements)
        carbons = list(range(base, base + 6))
        elements.extend(["C"] * 6)
        for k in range(6):
            bonds.append((carbons[k], carbons[(k + 1) % 6], 1.4))
        bonds.append((0, carbons[0], 1.0))              # ipso carbon
        hydrogens = []
        for k in range(1, 6):                           # five H on non-ipso C
            h = len(elements)
            elements.append("H")
            bonds.append((carbons[k], h, 1.0))
            hydrogens.append(h)
        ring_atoms.append((carbons, hydrogens))
    return elements, bonds, ring_atoms


def tetraphenyl():
    from rxn_core.domain import MolecularEndpoint
    elements, bonds, rings = _tetraphenylmethane()
    n = len(elements)
    rng = random.Random(7)
    coords = np.array([[rng.uniform(-6, 6) for _ in range(3)] for _ in range(n)])

    def matrix(bond_list):
        w = np.zeros((n, n))
        for a, b, v in bond_list:
            w[a, b] = w[b, a] = v
        return w

    reactant = MolecularEndpoint(tuple(elements), coords, matrix(bonds), label="R")
    # product: on ring 0 move the hydrogen of carbon k=2 to carbon k=3
    carbons, hydrogens = rings[0]
    moved_h = hydrogens[1]                               # H bonded to carbons[2]
    product_bonds = [
        (a, b, v) for a, b, v in bonds
        if not (a == carbons[2] and b == moved_h)]
    product_bonds.append((carbons[3], moved_h, 1.0))
    product = MolecularEndpoint(
        tuple(elements), coords + 0.05, matrix(product_bonds), label="P")
    return reactant, product


CASES = {"tempo": tempo, "tetraphenyl": tetraphenyl}
