"""Coordinate frame helpers."""
from __future__ import annotations

import numpy as np


def reindex_to_R_frame(elR, xyzR, elT, xyzT, mapping):
    """Build ``(elements, coords, missing)`` of length ``len(elR)`` in R order."""
    out_el = list(elR)
    out_xyz = np.array(xyzR, dtype=float).copy()
    missing = []
    for i in range(len(elR)):
        j = mapping.get(i)
        if j is None:
            missing.append(i)
            continue
        out_el[i] = elT[j]
        out_xyz[i] = xyzT[j]
    return out_el, out_xyz, missing
