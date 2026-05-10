"""
One-shot diagnostic: for one (step, IG) print every alignment branch
with its ranker score + features. Reads from the per-step containerized
output at out/ranked_views/<step>/xtb/.

Usage:
  python viewer/inspect_branches.py <step> <ig_label>
  e.g. python viewer/inspect_branches.py ru_nh3_1a iter20
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
PROJECT_ROOT = _HERE.parent

import sys
import numpy as np

from rxn_core_pq import align_from_arrays
from analyze_core_modes import (
    parse_g98_modes, core_atoms_in_R_frame, reindex_modes_to_R,
    bond_reaction_vector, bond_overlap_per_mode,
    rxn_overlap_per_mode, reaction_coord_delta,
)
from align_bgcp_coords import load_cached_xtb, fill_unmapped_greedy
from rxn_core_frag import parse_xyz

W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3


def main():
    step, label = sys.argv[1], sys.argv[2]
    cache = PROJECT_ROOT / "out" / "ranked_views" / step / "xtb"
    elR, xyzR, wboR, _ = load_cached_xtb(cache / "R")
    elP, xyzP, wboP, _ = load_cached_xtb(cache / "P")
    rp = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP)
    mapping_RP = dict(rp["mapping"])
    inv_RP = {v: k for k, v in mapping_RP.items()}
    core_R = core_atoms_in_R_frame(mapping_RP, rp["broken"], rp["formed"])
    full_RP = fill_unmapped_greedy(elR, xyzR, elP, xyzP, mapping_RP)
    delta_RP = reaction_coord_delta(np.asarray(xyzR, float),
                                     np.asarray(xyzP, float), full_RP)
    broken_R = [(int(a), int(b)) for (a, b, _, _) in rp["broken"]]
    formed_R = [(int(inv_RP[a]), int(inv_RP[b]))
                for (a, b, _, _) in rp["formed"]
                if a in inv_RP and b in inv_RP]

    hess = cache / f"hess_{label}"
    # The hess dir has both the input xyz and xtbhess.xyz; load_cached_xtb
    # expects exactly one xyz. Pick the original input (filename matches the
    # IG label).
    ig_xyz = next(p for p in hess.glob("*.xyz")
                  if "xtbhess" not in p.name)
    elT, xyzT_arr = parse_xyz(ig_xyz)
    xyzT = np.asarray(xyzT_arr, float)
    n = len(elT)
    wboT = np.zeros((n, n))
    for ln in (hess / "wbo").read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3: continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wboT[i, j] = v; wboT[j, i] = v
    freqs, modes_TS = parse_g98_modes(hess / "g98.out")
    imag_idx = list(np.where(freqs < 0)[0])
    n_imag = len(imag_idx)

    print(f"step:        {step}")
    print(f"IG:          {label}")
    print(f"n_atoms:     {len(elR)}")
    print(f"broken_R:    {broken_R}")
    print(f"formed_R:    {formed_R}")
    print(f"core_R:      {sorted(core_R)}")
    print(f"n_imag:      {n_imag}")
    print(f"imag freqs:  {[f'{freqs[k]:.1f}' for k in imag_idx]}")
    print()

    it = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT, return_all=True)
    branches = it.get("all_scored", [])

    seen = {}
    for bi, (alignment_score, br_mapping, _, _, _) in enumerate(branches):
        mkey = tuple(sorted(dict(br_mapping).items()))
        if mkey in seen:
            seen[mkey].append(bi)
            continue
        seen[mkey] = [bi]

        mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, dict(br_mapping))
        modes_R = reindex_modes_to_R(modes_TS, mapping_RT, len(elR))
        sq = (modes_R ** 2).sum(axis=2)
        total = sq.sum(axis=1)
        core_e = (sq[:, core_R].sum(axis=1) if core_R
                  else np.zeros(modes_R.shape[0]))
        kappa = np.where(total > 1e-12, core_e / total, 0.0)
        rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R)
        ts_xyz_in_R = np.zeros_like(np.asarray(xyzR))
        for r_idx, t_idx in mapping_RT.items():
            ts_xyz_in_R[r_idx] = xyzT[t_idx]
        V = bond_reaction_vector(ts_xyz_in_R, broken_R, formed_R)
        beta = bond_overlap_per_mode(modes_R, V)

        picked_k = max(imag_idx, key=lambda k: beta[k]) if imag_idx else None
        if picked_k is None:
            print(f"branch {len(seen):>2}: no imag modes")
            continue
        b = float(beta[picked_k])
        r_ = float(rho[picked_k])
        c = float(kappa[picked_k])
        score = (b * (1 + W_RXN * r_) * (1 + W_CORE * c)
                 / max(n_imag, 1) ** IMAG_PEN)
        print(f"branch {len(seen):>2}  alignment={alignment_score}  "
              f"picked_mode_k={picked_k}  freq={freqs[picked_k]:>7.1f}  "
              f"beta={b:.3f}  rho={r_:.3f}  kappa={c:.3f}  score={score:.4f}")

    print()
    print(f"unique branches: {len(seen)}  /  raw branches: {len(branches)}")


if __name__ == "__main__":
    main()
