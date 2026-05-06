"""
Evaluate dwbo_overlap (per-atom WBO-environment reaction direction)
against existing rankers.

Reuses per-step HTML data (xyz_coords + mode disps already in R-frame)
so we DO NOT re-run R↔TS alignment for every initial guess. Only
R↔P alignment is recomputed once per step (~0.5s/step) to get
mapping_RP for the WBO difference matrix.

For each step:
  1. Read DATA from out/mode_viewer/<step>.html (no recompute):
     - panels/ts_list[*].xyz_coords   (TS in R-frame)
     - panels/ts_list[*].modes[*].disp (mode disp in R-frame)
     - broken_bonds, formed_bonds_R, core_atoms
  2. Load cached R, P WBO matrices (file read).
  3. Run R↔P PQ alignment ONCE per step → mapping_RP.
  4. Build V_dwbo at each TS's coords using all atom pairs weighted by
     dW_ij = WBO_R[i,j] - WBO_P[m(i), m(j)].
  5. Compute dwbo_overlap and existing metrics on every mode.
  6. Compute gt_alignment of every IG mode against GT default mode.
  7. Compare ranker top-1 picks across criteria.

Outputs:
  out/mode_analysis/dwbo_overlap_eval.csv
  prints summary table.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import csv
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

from rxn_core_pq import align_from_arrays
from align_bgcp_coords import load_cached_xtb
from analyze_core_modes import bond_reaction_vector, bond_overlap_per_mode
from bgcp_io import list_step_dirs


SRC_DIR = PROJECT_ROOT / "out" / "mode_viewer"
WORK_MODES = Path("/Users/yunhengz/empty_for_claude/rxn_core/work_modes")
OUT_CSV = PROJECT_ROOT / "out" / "mode_analysis" / "dwbo_overlap_eval.csv"


def cos_sim(a, b):
    a = np.asarray(a).reshape(-1)
    b = np.asarray(b).reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9: return 0.0
    return abs(float(a @ b)) / (na * nb)


def dwbo_reaction_vector(xyz_TS_in_R, wboR, wboP, mapping_RP, dw_threshold=0.05):
    """V_i = -sum_k dW_ik * û_ik. Only sums pairs where |dW| > threshold
    (skip noise-floor) and both atoms are mapped."""
    n_R = len(xyz_TS_in_R)
    xyz = np.asarray(xyz_TS_in_R, dtype=float)
    # Pre-compute dW matrix in R-frame, shape (n_R, n_R), zero outside
    # mapped atoms.
    mapped = list(mapping_RP.keys())
    dW = np.zeros((n_R, n_R))
    for i in mapped:
        pi = mapping_RP[i]
        for j in mapped:
            if j <= i: continue
            v = float(wboR[i, j]) - float(wboP[pi, mapping_RP[j]])
            if abs(v) >= dw_threshold:
                dW[i, j] = v
                dW[j, i] = v
    # Now build V
    V = np.zeros((n_R, 3))
    # Vectorize over j for each i:
    for i in range(n_R):
        if not np.any(dW[i]): continue
        diffs = xyz - xyz[i]                    # (n_R, 3)
        norms = np.linalg.norm(diffs, axis=1)
        mask = (norms > 1e-9) & (dW[i] != 0)
        if not mask.any(): continue
        u = diffs[mask] / norms[mask, None]     # (k, 3)
        w = dW[i][mask][:, None]                # (k, 1)
        V[i] -= (w * u).sum(axis=0)
    return V


def dwbo_overlap_per_mode(modes_R, V):
    n_modes = modes_R.shape[0]
    if n_modes == 0: return np.zeros(0)
    v_flat = np.asarray(V).reshape(-1)
    v_norm = float(np.linalg.norm(v_flat))
    if v_norm < 1e-9: return np.zeros(n_modes)
    v_unit = v_flat / v_norm
    m_flat = modes_R.reshape(n_modes, -1)
    m_norm = np.linalg.norm(m_flat, axis=1)
    dots = np.abs(m_flat @ v_unit)
    out = np.zeros(n_modes)
    valid = m_norm > 1e-9
    out[valid] = dots[valid] / m_norm[valid]
    return out


def evaluate_step(step):
    html_path = SRC_DIR / f"{step}.html"
    if not html_path.exists():
        return None
    text = html_path.read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m: return None
    data = json.loads(m.group(1))

    md = WORK_MODES / step
    if not (md / "R" / "wbo").exists() or not (md / "P" / "wbo").exists():
        return None

    elR, xyzR, wboR, _ = load_cached_xtb(md / "R")
    elP, xyzP, wboP, _ = load_cached_xtb(md / "P")
    rp = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP)
    mapping_RP = dict(rp['mapping'])

    # GT default = its default_mode_idx in payload (already chosen by bond_overlap)
    gt = next((t for t in data['ts_list'] if t['label'] == 'groundtruth'), None)
    if gt is None or not gt.get('modes'): return None
    gt_default = gt['modes'][gt['default_mode_idx']]
    gt_disp = np.asarray(gt_default['disp'])
    if gt_default['freq'] >= 0:
        return None

    ig_tss = [t for t in data['ts_list'] if t['label'] != 'groundtruth' and t.get('modes')]
    if not ig_tss: return None

    # For each TS, compute dwbo_overlap on its existing modes (no realignment).
    results_per_ts = {}
    for ts in data['ts_list']:
        ts_xyz = np.asarray(ts['xyz_coords'])
        V_dwbo = dwbo_reaction_vector(ts_xyz, wboR, wboP, mapping_RP)
        modes_R = np.asarray([mm['disp'] for mm in ts['modes']])
        dwbo_ov = dwbo_overlap_per_mode(modes_R, V_dwbo)
        results_per_ts[ts['label']] = {
            'modes': ts['modes'], 'dwbo_ov': dwbo_ov,
        }

    # Per-mode gt_alignment for IGs
    for label in [t['label'] for t in ig_tss]:
        modes = results_per_ts[label]['modes']
        for mm in modes:
            mm['_align'] = cos_sim(np.asarray(mm['disp']), gt_disp)

    def pick_mode(modes, score_fn):
        imag = [m for m in modes if m['freq'] < 0]
        if imag:
            return max(imag, key=score_fn)
        return max(modes, key=score_fn)

    rankers = {
        'bond_overlap':   lambda mode: mode.get('bond_overlap', 0.0),
        'rxn_overlap':    lambda mode: mode.get('rxn_overlap', 0.0),
        'core_fraction':  lambda mode: mode.get('core_fraction', 0.0),
    }

    out_row = {'step': step}
    for rname, sf in rankers.items():
        ig_picks = []
        for ts in ig_tss:
            picked = pick_mode(results_per_ts[ts['label']]['modes'], sf)
            ig_picks.append((sf(picked), picked['_align'], ts['label']))
        ig_picks.sort(key=lambda x: -x[0])
        out_row[rname] = ig_picks[0][1]
        out_row[rname + '_label'] = ig_picks[0][2]

    # dwbo-based rankers (need to attach dwbo_ov to each mode)
    for ts in ig_tss:
        modes = results_per_ts[ts['label']]['modes']
        dwbo_ov = results_per_ts[ts['label']]['dwbo_ov']
        for k, mm in enumerate(modes):
            mm['_dwbo'] = float(dwbo_ov[k])

    dwbo_rankers = {
        'dwbo_overlap':   lambda mode: mode['_dwbo'],
        'bond_x_dwbo':    lambda mode: mode.get('bond_overlap', 0.0) * mode['_dwbo'],
        'bond_plus_dwbo': lambda mode: mode.get('bond_overlap', 0.0) + mode['_dwbo'],
    }
    for rname, sf in dwbo_rankers.items():
        ig_picks = []
        for ts in ig_tss:
            picked = pick_mode(results_per_ts[ts['label']]['modes'], sf)
            ig_picks.append((sf(picked), picked['_align'], ts['label']))
        ig_picks.sort(key=lambda x: -x[0])
        out_row[rname] = ig_picks[0][1]
        out_row[rname + '_label'] = ig_picks[0][2]

    # Oracles
    out_row['oracle'] = max((m['_align'] for ts in ig_tss for m in results_per_ts[ts['label']]['modes']), default=0.0)
    out_row['oracle_imag'] = max((m['_align'] for ts in ig_tss for m in results_per_ts[ts['label']]['modes'] if m['freq'] < 0), default=0.0)
    return out_row


def main():
    steps = [d.name for d in list_step_dirs()]
    steps = [s for s in steps if (WORK_MODES / s / "R" / "wbo").exists()]
    print(f"Evaluating {len(steps)} steps (R↔P alignment 1×/step, no R↔TS)...")
    t0 = time.time()
    rows = []
    for i, s in enumerate(steps, 1):
        try:
            r = evaluate_step(s)
            if r is not None:
                rows.append(r)
        except Exception as e:
            print(f"  {s}: ERROR {e}")
        if i % 20 == 0:
            elapsed = time.time() - t0
            print(f"  [{i}/{len(steps)}]  ({elapsed:.0f}s, ~{elapsed/i:.1f}s/step)")
    print(f"\nDone in {time.time()-t0:.0f}s")

    if not rows: return
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with OUT_CSV.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)

    rankers = ['bond_overlap', 'rxn_overlap', 'core_fraction',
               'dwbo_overlap', 'bond_x_dwbo', 'bond_plus_dwbo']
    print(f"\n{len(rows)} steps evaluated.\n")
    print(f"{'ranker':22s}  {'mean':>7s}  {'median':>7s}  "
          f"{'≥0.7':>5s}  {'≥0.5':>5s}  {'gap':>7s}")
    oracle_vals = np.array([r['oracle'] for r in rows])
    o_imag = np.array([r['oracle_imag'] for r in rows])
    print('-' * 62)
    for name in rankers:
        v = np.array([r[name] for r in rows])
        gap = oracle_vals.mean() - v.mean()
        print(f"{name:22s}  {v.mean():7.3f}  {np.median(v):7.3f}  "
              f"{(v >= 0.7).mean()*100:4.0f}%  {(v >= 0.5).mean()*100:4.0f}%  {gap:7.3f}")
    print('-' * 62)
    print(f"{'oracle_imag':22s}  {o_imag.mean():7.3f}  {np.median(o_imag):7.3f}  "
          f"{(o_imag >= 0.7).mean()*100:4.0f}%  {(o_imag >= 0.5).mean()*100:4.0f}%")
    print(f"{'oracle':22s}  {oracle_vals.mean():7.3f}  {np.median(oracle_vals):7.3f}  "
          f"{(oracle_vals >= 0.7).mean()*100:4.0f}%  {(oracle_vals >= 0.5).mean()*100:4.0f}%")
    print(f"\nCSV: {OUT_CSV}")


if __name__ == "__main__":
    main()
