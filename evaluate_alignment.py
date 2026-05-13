"""
Per-IG alignment + GT-similarity evaluation on the BGCP cached set,
designed to A/B different ranker formulas on the same data.

For each of 155 BGCP steps:
  1. Load cached R, P, GT, 20 IG xtb output (no xtb runs).
  2. R<->P graph alignment -> broken_R, formed_R, core_R, Delta_RP.
  3. For GT and each IG: multi-branch R<->X align (capped to keep
     runtime sane on highly symmetric systems), keep the branch that
     maximizes S = beta*(1+rho)*(1+0.2*kappa) / n_imag^0.3, take its
     picked imag mode (max-beta).
  4. Mass-weighted cosine of each IG's picked-mode displacement against
     GT's picked-mode displacement (both in R-frame).

For ranker A/B: we save per-IG (beta, rho, kappa, n_imag, freq,
mwc_to_GT) to JSON. A separate analysis re-ranks the same IGs under
different score formulas (beta-only, beta/n_imag, full formula, etc.)
and reports each formula's top-1 mwc-to-GT distribution.

Usage:
  python evaluate_alignment.py [--workers N] [--limit N] [--max-branches N]
"""
from __future__ import annotations
import argparse
import json
import multiprocessing as mp
import os
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from rxn_core import (
    align_from_arrays,
    bond_overlap_per_mode, bond_reaction_vector,
    core_atoms_in_R_frame,
    parse_g98_modes, parse_xyz, reaction_coord_delta,
    reindex_modes_to_R, rxn_overlap_per_mode,
)


PROJECT = Path(__file__).resolve().parent
WORK = PROJECT / "appendix_perparation" / "xtb_frequency_calculations"
OUT_JSON = PROJECT / "out" / "bgcp_alignment_eval.json"

# Default ranker hyperparameters for the multi-branch picker (reused
# from rxn_core.pipeline). Used only to choose which alignment branch
# to keep per IG. Downstream ranker A/B uses the captured (beta, rho,
# kappa, n_imag) so the choice of in-loop S formula doesn't pre-bias
# the ranker comparison.
W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3

ATOMIC_MASS = {
    'H': 1.008, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00,
    'F': 19.00, 'Na': 22.99, 'Mg': 24.31, 'Al': 26.98, 'Si': 28.09,
    'P': 30.97, 'S': 32.07, 'Cl': 35.45, 'K': 39.10, 'Ca': 40.08,
    'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.00, 'Mn': 54.94,
    'Fe': 55.85, 'Co': 58.93, 'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.38,
    'Ga': 69.72, 'Ge': 72.63, 'As': 74.92, 'Se': 78.97, 'Br': 79.90,
    'Mo': 95.95, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87,
    'In': 114.82, 'Sn': 118.71, 'I': 126.90, 'Pt': 195.08, 'Au': 196.97,
    'Re': 186.21, 'Hg': 200.59, 'Pb': 207.20,
}


def load_xyz_wbo(workdir: Path):
    xyz_path = next(p for p in workdir.glob("*.xyz") if "xtbhess" not in p.name)
    el, xyz = parse_xyz(xyz_path)
    xyz = np.asarray(xyz, float)
    n = len(el)
    wbo = np.zeros((n, n))
    for ln in (workdir / "wbo").read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3: continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo[i, j] = v; wbo[j, i] = v
    return el, xyz, wbo


def mass_weighted_cos(disp_a, disp_b, elements, core_atoms=None):
    """Cosine similarity on Cartesian displacement vectors.

    NOTE: name kept for compatibility but mass-weighting is OFF -- it
    over-amplified heavy-atom mismatches and suppressed the H motion
    which is usually the chemically dominant signal in TS modes.
    Restricted to core atoms when given."""
    a = np.asarray(disp_a, float); b = np.asarray(disp_b, float)
    if core_atoms is not None and len(core_atoms) > 0:
        idx = np.asarray(core_atoms, dtype=int)
        a = a[idx]; b = b[idx]
    af = a.reshape(-1); bf = b.reshape(-1)
    nA = float(np.linalg.norm(af)); nB = float(np.linalg.norm(bf))
    if nA < 1e-9 or nB < 1e-9: return 0.0
    return float(abs(af @ bf) / (nA * nB))


def best_branch_score(elR, xyzR, wboR, elT, xyzT, wboT, freqs, modes_TS,
                      broken_R, formed_R, core_R, delta_RP, max_branches):
    """Multi-branch sweep + pick highest-S branch. Critical for mwc-to-GT
    correctness: both GT and each IG must use their highest-S core-atom
    mapping so the picked imag mode vectors live in a consistent R-frame
    (i.e. atom i in IG's picked-disp == same physical atom as i in GT's
    picked-disp). Single-branch alignment can choose different swap-
    equivalent core mappings for GT vs IG and produce spurious mwc."""
    imag_idx = list(np.where(freqs < 0)[0])
    n_imag = len(imag_idx)
    if n_imag == 0:
        return None
    it = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT,
                           return_all=True, max_branches=max_branches)
    branches = it.get("all_scored", [])
    if not branches:
        branches = [(None, dict(it["mapping"]), None, None, None)]
    seen_witness = set(); seen_core = set()
    best = None
    for (_, br_mapping, _, _, _) in branches:
        br_d = dict(br_mapping)
        witness_key = tuple(sorted(br_d.items()))
        if witness_key in seen_witness: continue
        seen_witness.add(witness_key)
        # Score-equivalent core-only signature -- spectator-only forks give
        # identical S, so we only need one representative.
        core_key = tuple(sorted((c, br_d[c]) for c in core_R if c in br_d))
        if core_key in seen_core: continue
        seen_core.add(core_key)

        mapping_RT = br_d
        modes_R = reindex_modes_to_R(modes_TS, mapping_RT, len(elR))
        mode_norms = np.linalg.norm(modes_TS.reshape(modes_TS.shape[0], -1), axis=1)
        sq = (modes_R ** 2).sum(axis=2)
        total = mode_norms ** 2
        core_e = sq[:, core_R].sum(axis=1) if core_R else np.zeros(modes_R.shape[0])
        kappa = np.where(total > 1e-12, core_e / total, 0.0)
        rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R,
                                    mode_norms=mode_norms)
        ts_xyz_in_R = np.asarray(xyzR, float).copy()
        for r, t in mapping_RT.items():
            ts_xyz_in_R[r] = xyzT[t]
        V = bond_reaction_vector(ts_xyz_in_R, broken_R, formed_R)
        beta = bond_overlap_per_mode(modes_R, V, mode_norms=mode_norms)
        picked_k = max(imag_idx, key=lambda k: beta[k])
        b = float(beta[picked_k]); r = float(rho[picked_k]); c = float(kappa[picked_k])
        score = b * (1 + W_RXN * r) * (1 + W_CORE * c) / max(n_imag, 1) ** IMAG_PEN
        if best is None or score > best["S"]:
            best = {
                "S": score, "beta": b, "rho": r, "kappa": c,
                "picked_disp": modes_R[picked_k],
                "picked_freq": float(freqs[picked_k]),
                "n_imag": n_imag,
            }
    return best


def list_iter_dirs(step_modes_dir: Path):
    for hess in sorted(step_modes_dir.glob("hess_iter*")):
        m = re.match(r"hess_iter(\d+)$", hess.name)
        if not m: continue
        sp = step_modes_dir / f"sp_iter{m.group(1)}"
        if not sp.exists(): continue
        yield f"iter{m.group(1)}", hess, sp


def process_step(args_tuple):
    step_name, max_branches = args_tuple
    step_dir = WORK / step_name
    if not (step_dir / "R" / "wbo").exists() or not (step_dir / "P" / "wbo").exists():
        return {"step": step_name, "error": "missing R or P xtb cache"}
    elR, xyzR, wboR = load_xyz_wbo(step_dir / "R")
    elP, xyzP, wboP = load_xyz_wbo(step_dir / "P")
    rp = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP,
                           max_branches=max_branches)
    mapping_RP = dict(rp["mapping"])
    inv_RP = {v: k for k, v in mapping_RP.items()}
    core_R = core_atoms_in_R_frame(mapping_RP, rp["broken"], rp["formed"])
    delta_RP = reaction_coord_delta(np.asarray(xyzR, float),
                                     np.asarray(xyzP, float), mapping_RP)
    broken_R = [(int(a), int(b)) for (a, b, _, _) in rp["broken"]]
    formed_R = [(int(inv_RP[a]), int(inv_RP[b]))
                for (a, b, _, _) in rp["formed"]
                if a in inv_RP and b in inv_RP]

    gt_dir = step_dir / "hess_groundtruth"; gt_sp = step_dir / "sp_groundtruth"
    if not (gt_dir / "g98.out").exists() or not gt_sp.exists():
        return {"step": step_name, "error": "missing GT hess/sp"}
    elT, xyzT, wboT = load_xyz_wbo(gt_sp)
    freqs_gt, modes_gt = parse_g98_modes(gt_dir / "g98.out")
    gt_best = best_branch_score(
        elR, xyzR, wboR, elT, xyzT, wboT,
        freqs_gt, modes_gt, broken_R, formed_R, core_R, delta_RP, max_branches,
    )
    if gt_best is None:
        return {"step": step_name, "error": "GT has no imag mode"}

    ig_records = []
    for label, hess_dir, sp_dir in list_iter_dirs(step_dir):
        try:
            elT, xyzT, wboT = load_xyz_wbo(sp_dir)
            freqs_ig, modes_ig = parse_g98_modes(hess_dir / "g98.out")
            best = best_branch_score(
                elR, xyzR, wboR, elT, xyzT, wboT,
                freqs_ig, modes_ig, broken_R, formed_R, core_R, delta_RP, max_branches,
            )
            if best is None:
                ig_records.append({
                    "label": label, "n_imag": 0,
                    "beta": 0.0, "rho": 0.0, "kappa": 0.0,
                    "picked_freq": None, "mwc_to_GT": 0.0,
                })
                continue
            mwc = mass_weighted_cos(best["picked_disp"], gt_best["picked_disp"], elR,
                                    core_atoms=core_R)
            ig_records.append({
                "label": label,
                "beta": best["beta"], "rho": best["rho"], "kappa": best["kappa"],
                "n_imag": best["n_imag"], "picked_freq": best["picked_freq"],
                "mwc_to_GT": mwc,
            })
        except Exception as e:
            ig_records.append({"label": label, "error": f"{type(e).__name__}: {e}"})

    return {
        "step": step_name,
        "n_atoms": len(elR),
        "n_broken": len(broken_R), "n_formed": len(formed_R),
        "gt_picked_freq": gt_best["picked_freq"],
        "gt_n_imag": gt_best["n_imag"],
        "gt_beta": gt_best["beta"], "gt_rho": gt_best["rho"],
        "gt_kappa": gt_best["kappa"], "gt_score": gt_best["S"],
        "igs": ig_records,
    }


def _safe(args_tuple):
    name, _ = args_tuple
    try:
        return process_step(args_tuple)
    except Exception as e:
        return {"step": name, "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--steps", nargs="+", default=None)
    ap.add_argument("--max-branches", type=int, default=1_000_000,
                    help="cap on alignment branches per pass; default matches "
                         "the rxn_core alignment default (1_000_000). Lower "
                         "this only when debugging pathological symmetric systems.")
    args = ap.parse_args()

    all_steps = sorted(d.name for d in WORK.iterdir() if d.is_dir())
    if args.steps:
        steps = [s for s in all_steps if s in set(args.steps)]
    elif args.limit:
        steps = all_steps[:args.limit]
    else:
        steps = all_steps
    print(f"Evaluating {len(steps)} BGCP steps  workers={args.workers}  "
          f"max_branches={args.max_branches}")

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.time()
    job_args = [(s, args.max_branches) for s in steps]
    with mp.Pool(args.workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(_safe, job_args), 1):
            rows.append(rec)
            # Incremental JSON dump after every step so we have partial
            # results if the run is interrupted.
            OUT_JSON.write_text(json.dumps(rows, indent=1))
            err = rec.get("error", "")
            if err:
                print(f"  [{i:>3d}/{len(steps)}] {rec['step']:50s}  ERROR: {err[:60]}",
                      flush=True)
                continue
            igs = rec["igs"]
            if not igs:
                print(f"  [{i:>3d}/{len(steps)}] {rec['step']:50s}  no scored IGs",
                      flush=True)
                continue
            valid = [r for r in igs if 'beta' in r]
            if valid:
                top = max(valid, key=lambda r: r['beta'] * (1 + r['rho'])
                          * (1 + 0.2 * r['kappa']) / max(r['n_imag'], 1) ** 0.3)
                print(f"  [{i:>3d}/{len(steps)}] {rec['step']:50s}  "
                      f"top1={top['label']}  beta={top['beta']:.3f}  "
                      f"mwc={top['mwc_to_GT']:.3f}", flush=True)
    print(f"\ntotal: {time.time()-t0:.1f}s")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
