"""
For each step in BGCP: run xtb Hessian on the ground-truth TS and on
each initial-guess TS. For every mode, project displacements onto the
core atoms (atoms touching R->P broken/formed bonds) and pick the
mode with highest core-atom contribution. Report:
  - n_imag (count of modes with freq < 0)
  - reaction_mode_freq (frequency of the highest-core-fraction mode)
  - core_fraction (||mode_core||^2 / ||mode_all||^2)
  - imag_core_freq (highest-core-fraction freq among imaginary modes only)

Usage:
  python analyze_modes.py                # all 161 BGCP steps × all guesses
  python analyze_modes.py --limit 5      # first 5 steps
  python analyze_modes.py --steps name1  # specific step
  python analyze_modes.py --gt-only      # GT TS only, skip initial guesses
"""
from __future__ import annotations
import argparse
import csv
import json
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from rxn_core_frag import (
    run_xtb, run_xtb_hess, build_graph, find_islands, expand_mapping,
    classify_bonds, _generate_seed_orders,
)
from build_bgcp_viewer import (
    BGCP_ROOT, OUT, WORK, LOOKUP, list_step_dirs, list_initial_guesses,
    iter_num, read_xyzs,
)


WORK_HESS = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_modes")
WORK_HESS.mkdir(parents=True, exist_ok=True)


def best_mapping(g_a, g_b, wbo_a, wbo_b, n_seeds=10):
    orders = _generate_seed_orders(g_a, n_seeds)
    best = None
    for order in orders:
        m, _ = find_islands(g_a, g_b, seed_order=order)
        m = expand_mapping(m, g_a, g_b)
        br, fm, _, _ = classify_bonds(m, wbo_a, wbo_b)
        score = (len(br) + len(fm), -len(m))
        if best is None or score < best[0]:
            best = (score, m)
    return best[1]


def core_atoms_in_R_frame(map_R_to_P, broken, formed):
    inv = {v: k for k, v in map_R_to_P.items()}
    core = set()
    for (i, j, _, _) in broken:
        core.add(i); core.add(j)
    for (ip, jp, _, _) in formed:
        if ip in inv: core.add(inv[ip])
        if jp in inv: core.add(inv[jp])
    return sorted(core)


def reindex_modes(modes_ts, map_R_to_TS, n_R):
    """modes_ts is shape (n_modes, n_TS, 3); produce (n_modes, n_R, 3)
    where each R-atom takes its TS-mapped displacement (zero if unmapped)."""
    n_modes = modes_ts.shape[0]
    out = np.zeros((n_modes, n_R, 3))
    for ri, ti in map_R_to_TS.items():
        out[:, ri, :] = modes_ts[:, ti, :]
    return out


def best_core_mode(modes_R, freqs, core_indices, imag_only=False):
    """Pick the mode with highest fraction of energy on core atoms.
    Returns (mode_idx, freq, core_fraction)."""
    if not core_indices or modes_R.shape[0] == 0:
        return None, None, None
    core = np.array(core_indices, dtype=int)
    # Per-mode total energy
    sq = (modes_R ** 2).sum(axis=2)   # (n_modes, n_atoms)
    total = sq.sum(axis=1)            # (n_modes,)
    core_e = sq[:, core].sum(axis=1)  # (n_modes,)
    fraction = np.where(total > 1e-9, core_e / total, 0.0)
    candidates = np.arange(len(freqs))
    if imag_only:
        candidates = candidates[freqs < 0]
        if len(candidates) == 0:
            return None, None, None
    best_idx = candidates[np.argmax(fraction[candidates])]
    return int(best_idx), float(freqs[best_idx]), float(fraction[best_idx])


def analyze_step(step_dir, gt_only=False):
    name = step_dir.name
    chg, uhf = LOOKUP.get(name, (0, 0))
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    wd = WORK_HESS / sanitized
    wd.mkdir(parents=True, exist_ok=True)

    # R + P xyz
    rxyz = read_xyzs(step_dir / "reactants")
    pxyz = read_xyzs(step_dir / "products")
    if rxyz is None or pxyz is None:
        raise RuntimeError("missing reactant or product")
    (wd / "reactant.xyz").write_text(rxyz)
    (wd / "product.xyz").write_text(pxyz)

    elR, xyzR, wboR = run_xtb(wd / "reactant.xyz", wd / "R", charge=chg, uhf=uhf)
    elP, xyzP, wboP = run_xtb(wd / "product.xyz", wd / "P", charge=chg, uhf=uhf)
    g_R = build_graph(elR, wboR); g_P = build_graph(elP, wboP)
    map_R_to_P = best_mapping(g_R, g_P, wboR, wboP)
    br, fm, _, _ = classify_bonds(map_R_to_P, wboR, wboP)
    core = core_atoms_in_R_frame(map_R_to_P, br, fm)
    n_R = len(elR)

    # Build list of TS candidates to analyze
    ts_paths = []
    gt = sorted((step_dir / "groundtruth").glob("*.xyz"))
    if gt:
        ts_paths.append(("groundtruth", gt[0]))
    if not gt_only:
        for g in sorted(list_initial_guesses(step_dir), key=iter_num):
            ts_paths.append((f"iter{iter_num(g)}", g))

    rows = []
    for label, ts_path in ts_paths:
        sub = re.sub(r"[^A-Za-z0-9._-]", "_", ts_path.stem)[:80]
        ts_local = wd / f"{sub}.xyz"
        ts_local.write_text(ts_path.read_text())
        ts_workdir = wd / f"hess_{sub}"
        try:
            elTS, xyzTS, freqs, modes_ts = run_xtb_hess(
                ts_local, ts_workdir, charge=chg, uhf=uhf)
        except Exception as e:
            rows.append({
                "label": label, "ts_file": ts_path.name,
                "n_imag": None, "reaction_mode_freq": None,
                "core_fraction": None, "imag_core_freq": None,
                "imag_core_fraction": None,
                "error": str(e)[:120],
            })
            continue

        g_TS = build_graph(elTS, build_graph_input := __import__("numpy").asarray(np.zeros((len(elTS), len(elTS))))).copy() if False else None
        # ^ artifact; just call build_graph properly
        from rxn_core_frag import run_xtb as _run  # already imported, just being explicit
        # re-build TS WBO graph by running a regular xtb single-point
        # (caches its own dir). We'll use the TS coords/wbo from a separate
        # cached SP rather than from the hess (xtb hess also writes wbo)
        ts_sp_dir = wd / f"sp_{sub}"
        elTS2, xyzTS2, wboTS = run_xtb(ts_local, ts_sp_dir, charge=chg, uhf=uhf)
        g_TS = build_graph(elTS2, wboTS)
        map_R_to_TS = best_mapping(g_R, g_TS, wboR, wboTS)
        modes_R = reindex_modes(modes_ts, map_R_to_TS, n_R)

        n_imag = int((freqs < 0).sum())
        idx_all,  f_all,  cf_all  = best_core_mode(modes_R, freqs, core, imag_only=False)
        idx_imag, f_imag, cf_imag = best_core_mode(modes_R, freqs, core, imag_only=True)

        rows.append({
            "label": label, "ts_file": ts_path.name,
            "n_imag": n_imag,
            "reaction_mode_freq": f_all,
            "core_fraction": cf_all,
            "imag_core_freq": f_imag,
            "imag_core_fraction": cf_imag,
            "error": None,
        })
    return name, n_R, len(core), rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--steps", nargs="+", default=None)
    ap.add_argument("--gt-only", action="store_true")
    ap.add_argument("--out", default=str(OUT / "bgcp_modes.csv"))
    args = ap.parse_args()

    step_dirs = list_step_dirs()
    if args.steps:
        wanted = set(args.steps)
        step_dirs = [d for d in step_dirs if d.name in wanted]
    elif args.limit is not None:
        step_dirs = step_dirs[args.start:args.start + args.limit]
    else:
        step_dirs = step_dirs[args.start:]
    print(f"[modes] {len(step_dirs)} steps")

    out_path = Path(args.out)
    fieldnames = [
        "step_id", "n_atoms", "core_size",
        "label", "ts_file",
        "n_imag", "reaction_mode_freq", "core_fraction",
        "imag_core_freq", "imag_core_fraction",
        "error",
    ]
    with out_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for k, sd in enumerate(step_dirs, 1):
            t = time.time()
            try:
                name, n_atoms, core_size, rows = analyze_step(sd, gt_only=args.gt_only)
                for r in rows:
                    w.writerow({"step_id": name, "n_atoms": n_atoms,
                                "core_size": core_size, **r})
                f.flush()
                # Summarize first row (GT TS) as quick log
                gt_row = rows[0] if rows else None
                if gt_row:
                    rf = gt_row.get('reaction_mode_freq')
                    icf = gt_row.get('imag_core_freq')
                    cf = gt_row.get('core_fraction')
                    icf_frac = gt_row.get('imag_core_fraction')
                    f_fmt = lambda x: f"{x:7.1f}" if isinstance(x, float) else "    n/a"
                    cf_fmt = lambda x: f"{x:.2f}" if isinstance(x, float) else "n/a"
                    print(f"[{k:>3}/{len(step_dirs)}]  {time.time()-t:5.1f}s  {name:<60s}  "
                          f"core={core_size:>2}  n_imag={gt_row['n_imag']}  "
                          f"top_core: ν={f_fmt(rf)} cf={cf_fmt(cf)}  "
                          f"top_imag_core: ν={f_fmt(icf)} cf={cf_fmt(icf_frac)}")
                else:
                    print(f"[{k:>3}/{len(step_dirs)}]  {time.time()-t:5.1f}s  {name:<60s}  no TS")
            except Exception as e:
                print(f"[{k:>3}/{len(step_dirs)}]  FAIL {sd.name}: {e}")
                traceback.print_exc()
    print(f"[modes] wrote {out_path}")


if __name__ == "__main__":
    main()
