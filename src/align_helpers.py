"""
Helpers for loading cached xtb output and reindexing target atoms into
the reactant's atom-order (R-frame).

The PQ aligner returns a partial atom-to-atom mapping; for atoms it
couldn't pin down (typically equivalent spectator H's where the search
ran out of constraints) `fill_unmapped_greedy` matches each unmapped
R-atom to the closest still-free same-element T-atom by Euclidean
distance. `reindex_to_R_frame` then realizes the mapping into a
(elements, coords) tuple in R atom order.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np

from rxn_core_frag import parse_xyz


def load_cached_xtb(workdir):
    """Load (elements, coords, wbo, xyz_path) from a pre-computed xtb workdir.
    The workdir must contain exactly one *.xyz and a `wbo` file."""
    workdir = Path(workdir)
    xyz_files = [f for f in workdir.iterdir() if f.suffix == '.xyz']
    if len(xyz_files) != 1:
        raise RuntimeError(f"expected 1 xyz in {workdir}, found {len(xyz_files)}")
    elements, coords = parse_xyz(xyz_files[0])
    n = len(elements)
    wbo = np.zeros((n, n))
    wf = workdir / "wbo"
    if not wf.exists():
        raise RuntimeError(f"no wbo file in {workdir}")
    for ln in wf.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3:
            continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo[i, j] = v; wbo[j, i] = v
    return elements, coords, wbo, xyz_files[0]


def fill_unmapped_greedy(elR, xyzR, elT, xyzT, mapping):
    """Greedy nearest-element-match for atoms PQ left unmapped. Pairs
    each unmapped R-atom with the closest still-free same-element T-atom."""
    used_T = set(mapping.values())
    unmapped_R = [i for i in range(len(elR)) if i not in mapping]
    if not unmapped_R:
        return dict(mapping)
    free_by_el = {}
    for j in range(len(elT)):
        if j in used_T:
            continue
        free_by_el.setdefault(elT[j], []).append(j)
    out = dict(mapping)
    for i in unmapped_R:
        cands = free_by_el.get(elR[i], [])
        if not cands:
            continue
        d = [np.linalg.norm(xyzT[j] - xyzR[i]) for j in cands]
        best = cands[int(np.argmin(d))]
        out[i] = best
        cands.remove(best)
    return out


def reindex_to_R_frame(elR, xyzR, elT, xyzT, mapping):
    """Build (elements, coords, missing) of length len(elR) in R atom order."""
    n = len(elR)
    out_el = list(elR)
    out_xyz = np.array(xyzR, dtype=float).copy()
    missing = []
    for i in range(n):
        j = mapping.get(i)
        if j is None:
            missing.append(i)
            continue
        out_el[i] = elT[j]
        out_xyz[i] = xyzT[j]
    return out_el, out_xyz, missing
