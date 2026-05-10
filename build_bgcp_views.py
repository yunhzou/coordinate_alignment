"""
For each BGCP step, emit a ranked HTML viewer at out/bgcp_views/<step>/view.html.

Layout:
  - Top row: R (left), P (middle), GT-TS (right, animated on its picked imag mode)
  - Below: all 20 IGs in a 4-col grid, sorted by ranker score, each animated
    on its picked imag mode. Each panel shows label, S, beta, rho, kappa,
    n_imag, picked_freq, AND mwc-to-GT (computed picked-imag vs GT-picked-imag,
    R-frame).

Reuses cached xtb output at appendix_perparation/xtb_frequency_calculations/.
No xtb runs; alignment graph search per step (~ tenths of a second).

Usage:
  python build_bgcp_views.py [--workers N] [--limit N] [--steps STEP ...]
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
    core_atoms_in_R_frame, fill_unmapped_greedy,
    parse_g98_modes, parse_xyz, reaction_coord_delta,
    reindex_modes_to_R, rxn_overlap_per_mode,
)


PROJECT = Path(__file__).resolve().parent
WORK = PROJECT / "appendix_perparation" / "xtb_frequency_calculations"
OUT_ROOT = PROJECT / "out" / "bgcp_views"

W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3
ATOMIC_MASS = {
    'H': 1.008, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'F': 19.00,
    'Na': 22.99, 'Mg': 24.31, 'Al': 26.98, 'Si': 28.09, 'P': 30.97, 'S': 32.07,
    'Cl': 35.45, 'K': 39.10, 'Ca': 40.08, 'Sc': 44.96, 'Ti': 47.87, 'V': 50.94,
    'Cr': 52.00, 'Mn': 54.94, 'Fe': 55.85, 'Co': 58.93, 'Ni': 58.69, 'Cu': 63.55,
    'Zn': 65.38, 'Ga': 69.72, 'Ge': 72.63, 'As': 74.92, 'Se': 78.97, 'Br': 79.90,
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
        wbo[i, j] = float(parts[2]); wbo[j, i] = wbo[i, j]
    return el, xyz, wbo


def mwc(da, db, elements, core_atoms=None):
    """Cosine of Cartesian displacement vectors, sign-blind, restricted
    to core atoms when given. Mass-weighting OFF (it over-amplified
    heavy-atom mismatches and suppressed the H motion that's usually
    the chemically dominant signal in TS modes)."""
    a = np.asarray(da, float); b = np.asarray(db, float)
    if core_atoms is not None and len(core_atoms) > 0:
        idx = np.asarray(core_atoms, dtype=int)
        a = a[idx]; b = b[idx]
    af = a.reshape(-1); bf = b.reshape(-1)
    nA = float(np.linalg.norm(af)); nB = float(np.linalg.norm(bf))
    return float(abs(af @ bf) / (nA * nB)) if (nA > 1e-9 and nB > 1e-9) else 0.0


def best_branch(elR, xyzR, wboR, elT, xyzT, wboT, freqs, modes_TS,
                broken_R, formed_R, core_R, delta_RP, max_branches=1_000_000):
    """Multi-branch sweep + pick highest-S branch. Both GT and each IG go
    through this so their picked imag modes live in a consistent R-frame
    for downstream mwc comparison."""
    imag_idx = list(np.where(freqs < 0)[0])
    n_imag = len(imag_idx)
    if n_imag == 0:
        return None
    it = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT,
                           return_all=True, max_branches=max_branches)
    branches = it.get("all_scored", [])
    if not branches:
        branches = [(None, dict(it["mapping"]), None, None, None)]
    seen_full = set(); seen_core = set()
    best = None
    for (_, br_mapping, _, _, _) in branches:
        br_d = dict(br_mapping)
        full_key = tuple(sorted(br_d.items()))
        if full_key in seen_full: continue
        seen_full.add(full_key)
        core_key = tuple(sorted((c, br_d[c]) for c in core_R if c in br_d))
        if core_key in seen_core: continue
        seen_core.add(core_key)

        mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, br_d)
        modes_R = reindex_modes_to_R(modes_TS, mapping_RT, len(elR))
        sq = (modes_R ** 2).sum(axis=2)
        total = sq.sum(axis=1)
        core_e = sq[:, core_R].sum(axis=1) if core_R else np.zeros(modes_R.shape[0])
        kappa = np.where(total > 1e-12, core_e / total, 0.0)
        rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R)
        ts_xyz_in_R = np.zeros_like(np.asarray(xyzR))
        for r, t in mapping_RT.items():
            ts_xyz_in_R[r] = xyzT[t]
        V = bond_reaction_vector(ts_xyz_in_R, broken_R, formed_R)
        beta = bond_overlap_per_mode(modes_R, V)
        picked_k = max(imag_idx, key=lambda k: beta[k])
        b = float(beta[picked_k]); r = float(rho[picked_k]); c = float(kappa[picked_k])
        score = b * (1 + W_RXN * r) * (1 + W_CORE * c) / max(n_imag, 1) ** IMAG_PEN
        if best is None or score > best["score"]:
            best = {
                "score": score, "beta": b, "rho": r, "kappa": c,
                "picked_disp": modes_R[picked_k].tolist(),
                "picked_freq": float(freqs[picked_k]),
                "n_imag": n_imag,
                "ts_xyz_in_R": ts_xyz_in_R.tolist(),
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
    full_RP = fill_unmapped_greedy(elR, xyzR, elP, xyzP, mapping_RP)
    delta_RP = reaction_coord_delta(np.asarray(xyzR, float),
                                     np.asarray(xyzP, float), full_RP)
    broken_R = [(int(a), int(b)) for (a, b, _, _) in rp["broken"]]
    formed_R = [(int(inv_RP[a]), int(inv_RP[b]))
                for (a, b, _, _) in rp["formed"]
                if a in inv_RP and b in inv_RP]

    # P reindexed to R-frame so bond pairs land on the right atoms
    xyzP_in_R = np.zeros_like(np.asarray(xyzR, float))
    for i_R, i_P in full_RP.items():
        xyzP_in_R[i_R] = xyzP[i_P]

    # GT
    gt_dir = step_dir / "hess_groundtruth"; gt_sp = step_dir / "sp_groundtruth"
    if not (gt_dir / "g98.out").exists() or not gt_sp.exists():
        return {"step": step_name, "error": "missing GT hess/sp"}
    elT, xyzT, wboT = load_xyz_wbo(gt_sp)
    freqs_gt, modes_gt = parse_g98_modes(gt_dir / "g98.out")
    gt_best = best_branch(elR, xyzR, wboR, elT, xyzT, wboT,
                          freqs_gt, modes_gt, broken_R, formed_R, core_R,
                          delta_RP, max_branches)
    if gt_best is None:
        return {"step": step_name, "error": "GT has no imag mode"}
    gt_disp = gt_best["picked_disp"]

    # IGs
    ig_records = []
    for label, hess_dir, sp_dir in list_iter_dirs(step_dir):
        try:
            elT, xyzT, wboT = load_xyz_wbo(sp_dir)
            freqs_ig, modes_ig = parse_g98_modes(hess_dir / "g98.out")
            best = best_branch(elR, xyzR, wboR, elT, xyzT, wboT,
                               freqs_ig, modes_ig, broken_R, formed_R,
                               core_R, delta_RP, max_branches)
            if best is None:
                continue
            mw = mwc(best["picked_disp"], gt_disp, elR, core_atoms=core_R)
            ig_records.append({
                "label": label,
                "score": best["score"],
                "beta": best["beta"], "rho": best["rho"], "kappa": best["kappa"],
                "n_imag": best["n_imag"], "picked_freq": best["picked_freq"],
                "mwc_to_GT": mw,
                "picked_disp": best["picked_disp"],
                "xyz_elements": list(elR),
                "xyz_coords": best["ts_xyz_in_R"],
            })
        except Exception as e:
            print(f"  {step_name} {label}: {e}", file=sys.stderr)

    if not ig_records:
        return {"step": step_name, "error": "no scored IGs"}
    # Sort by score descending
    ig_records.sort(key=lambda r: -r["score"])

    # Build payload + write HTML
    data = {
        "step": step_name,
        "n_atoms": len(elR),
        "core_atoms": list(core_R),
        "broken_bonds":   [list(b) for b in broken_R],
        "formed_bonds_R": [list(b) for b in formed_R],
        "reactant":   {"xyz_elements": elR,
                       "xyz_coords": np.asarray(xyzR, float).tolist()},
        "product":    {"xyz_elements": elR,
                       "xyz_coords": xyzP_in_R.tolist()},
        "groundtruth": {
            "xyz_elements": elR,
            "xyz_coords":   gt_best["ts_xyz_in_R"],
            "picked_disp":  gt_disp,
            "picked_freq":  gt_best["picked_freq"],
            "n_imag":       gt_best["n_imag"],
            "score":        gt_best["score"],
            "beta":         gt_best["beta"],
            "rho":          gt_best["rho"],
            "kappa":        gt_best["kappa"],
        },
        "igs": ig_records,
    }

    run_dir = OUT_ROOT / step_name
    run_dir.mkdir(parents=True, exist_ok=True)
    out_path = run_dir / "view.html"
    out_path.write_text(HTML.format(
        title=f"BGCP — {step_name}",
        n_atoms=len(elR),
        n_broken=len(broken_R),
        n_formed=len(formed_R),
        n_core=len(core_R),
        gt_freq_str=f"{gt_best['picked_freq']:.0f}i",
        gt_n_imag=gt_best["n_imag"],
        gt_score=f"{gt_best['score']:.3f}",
        gt_beta=f"{gt_best['beta']:.3f}",
        gt_rho=f"{gt_best['rho']:.3f}",
        gt_kappa=f"{gt_best['kappa']:.3f}",
        data_json=json.dumps(data),
    ))
    return {
        "step": step_name,
        "out": str(out_path),
        "n_ig": len(ig_records),
        "top1_label": ig_records[0]["label"],
        "top1_score": ig_records[0]["score"],
        "top1_mwc": ig_records[0]["mwc_to_GT"],
    }


def _safe(args_tuple):
    name, _ = args_tuple
    try:
        return process_step(args_tuple)
    except Exception as e:
        return {"step": name, "error": f"{type(e).__name__}: {e}",
                "trace": traceback.format_exc()}


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>{title}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 html, body {{ margin:0; padding:0; }}
 body {{ font-family:-apple-system,sans-serif; background:#fafafa; padding:14px; box-sizing:border-box; }}
 h2 {{ margin:0 0 4px; font-size:18px; }}
 .sub {{ color:#444; font-size:13px; margin-bottom:14px; }}
 .legend span {{ display:inline-block; padding:2px 8px; border-radius:4px; font-size:12px; margin-right:6px; }}
 .leg-broken {{ background:#fcd3d3; color:#a00; }}
 .leg-formed {{ background:#cdebd0; color:#070; }}
 .leg-mode   {{ background:#d6e7ff; color:#024; }}
 .leg-static {{ background:#eee;    color:#666; }}
 .ref-row {{ display:grid; grid-template-columns:repeat(3,1fr); gap:12px; margin-bottom:18px; }}
 .ig-grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; }}
 .panel  {{ background:white; border:1px solid #ddd; border-radius:6px; padding:6px 8px 8px; }}
 .panel.no-imag {{ background:#f7f7f7; }}
 .ph    {{ display:flex; justify-content:space-between; align-items:baseline; font-size:12px; margin-bottom:4px; }}
 .ph .lbl {{ font-weight:600; font-size:13px; }}
 .ph .rk  {{ font-family:ui-monospace,monospace; color:#024; }}
 .vw    {{ position:relative; width:100%; height:230px; }}
 .ref-row .vw {{ height:300px; }}
 .vwbox {{ position:absolute; inset:0; }}
 .meta  {{ font-family:ui-monospace,monospace; font-size:11px; color:#444; padding:3px 0 0; line-height:1.4; }}
 .meta b {{ color:#024; }}
 .badge-static {{ display:inline-block; background:#eee; color:#666; padding:1px 6px; border-radius:3px; font-size:11px; margin-left:6px; }}
</style></head><body>
<h2>{title}</h2>
<div class="sub">
  N atoms: <b>{n_atoms}</b> &nbsp;|&nbsp;
  broken: <b>{n_broken}</b> &nbsp;|&nbsp;
  formed: <b>{n_formed}</b> &nbsp;|&nbsp;
  core atoms: <b>{n_core}</b><br>
  <span class="legend">
    <span class="leg-broken">breaking bond</span>
    <span class="leg-formed">forming bond</span>
    <span class="leg-mode">picked imag mode</span>
    <span class="leg-static">no imag (static)</span>
  </span>
  <br>
  IGs sorted by ranker score
  <code>S = &beta;(1 + &rho;)(1 + 0.2&kappa;) / n_imag<sup>0.3</sup></code>;
  per panel we show S, &beta;, &rho;, &kappa;, n_imag, picked freq, and mwc-to-GT.
</div>
<div class="ref-row">
  <div class="panel"><div class="ph"><span class="lbl">Reactant</span><span class="rk">static</span></div>
    <div class="vw"><div id="vw_R" class="vwbox"></div></div></div>
  <div class="panel"><div class="ph"><span class="lbl">Product</span><span class="rk">static</span></div>
    <div class="vw"><div id="vw_P" class="vwbox"></div></div></div>
  <div class="panel"><div class="ph"><span class="lbl">Ground-truth TS</span>
    <span class="rk">S = {gt_score}</span></div>
    <div class="vw"><div id="vw_GT" class="vwbox"></div></div>
    <div class="meta">
      <b>&beta;</b>={gt_beta} &nbsp;
      <b>&rho;</b>={gt_rho} &nbsp;
      <b>&kappa;</b>={gt_kappa} &nbsp;
      <b>n_imag</b>={gt_n_imag} &nbsp;
      <b>freq</b>={gt_freq_str} cm⁻¹
    </div></div>
</div>
<div class="ig-grid" id="grid"></div>
<script>
const DATA = {data_json};
function buildBody(elements, xyz) {{
  let s = `${{xyz.length}}\nframe\n`;
  for (let i=0;i<xyz.length;i++) s += `${{elements[i]}}  ${{xyz[i][0].toFixed(6)}}  ${{xyz[i][1].toFixed(6)}}  ${{xyz[i][2].toFixed(6)}}\n`;
  return s;
}}
function buildBodyAt(elements, xyz, disp, scale) {{
  let s = `${{xyz.length}}\nframe\n`;
  for (let i=0;i<xyz.length;i++) {{
    const x=xyz[i][0]+scale*disp[i][0], y=xyz[i][1]+scale*disp[i][1], z=xyz[i][2]+scale*disp[i][2];
    s += `${{elements[i]}}  ${{x.toFixed(6)}}  ${{y.toFixed(6)}}  ${{z.toFixed(6)}}\n`;
  }}
  return s;
}}
function xyzAt(xyz, disp, scale) {{
  const out = new Array(xyz.length);
  for (let i=0;i<xyz.length;i++)
    out[i]=[xyz[i][0]+scale*disp[i][0], xyz[i][1]+scale*disp[i][1], xyz[i][2]+scale*disp[i][2]];
  return out;
}}
function decorateEvent(viewer, xyz, pair, color) {{
  const [i,j] = pair;
  if (i>=xyz.length || j>=xyz.length) return;
  const a = xyz[i], b = xyz[j];
  viewer.addCylinder({{start:{{x:a[0],y:a[1],z:a[2]}}, end:{{x:b[0],y:b[1],z:b[2]}},
                      color:color, radius:0.20, dashed:true}});
}}
function decorateBonds(viewer, which, xyz) {{
  if (which !== 'formed') for (const p of DATA.broken_bonds) decorateEvent(viewer,xyz,p,'red');
  if (which !== 'broken') for (const p of DATA.formed_bonds_R) decorateEvent(viewer,xyz,p,'green');
}}
function drawArrows(viewer, xyz, disp) {{
  for (const i of DATA.core_atoms) {{
    if (i>=xyz.length || !disp || !disp[i]) continue;
    const d = disp[i]; const len = Math.hypot(d[0],d[1],d[2]); if (len<1e-3) continue;
    const ax=xyz[i][0], ay=xyz[i][1], az=xyz[i][2];
    viewer.addArrow({{start:{{x:ax,y:ay,z:az}},
                     end:{{x:ax+d[0]*1.5,y:ay+d[1]*1.5,z:az+d[2]*1.5}},
                     color:'#0066cc', radius:0.06}});
  }}
}}
function makeStatic(divId, ts, which) {{
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(buildBody(ts.xyz_elements, ts.xyz_coords), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  decorateBonds(v, which||'both', ts.xyz_coords);
  v.zoomTo(); v.render(); return v;
}}
function makeAnimated(divId, ts, disp, which) {{
  const w = which||'both';
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(buildBody(ts.xyz_elements, ts.xyz_coords), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  decorateBonds(v, w, ts.xyz_coords); drawArrows(v, ts.xyz_coords, disp);
  v.zoomTo(); v.render();
  let t=0; const period=30; const amp=0.6;
  setInterval(()=>{{
    t=(t+1)%period;
    const scale = amp*Math.sin(2*Math.PI*t/period);
    const cur = xyzAt(ts.xyz_coords, disp, scale);
    v.removeAllModels(); v.removeAllShapes();
    v.addModel(buildBodyAt(ts.xyz_elements, ts.xyz_coords, disp, scale), 'xyz');
    v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
    decorateBonds(v, w, cur); drawArrows(v, cur, disp);
    v.render();
  }}, 60);
  return v;
}}
window.addEventListener('load', () => {{
  makeStatic('vw_R', DATA.reactant, 'broken');
  makeStatic('vw_P', DATA.product,  'formed');
  if (DATA.groundtruth.picked_disp) makeAnimated('vw_GT', DATA.groundtruth, DATA.groundtruth.picked_disp, 'both');
  else                               makeStatic('vw_GT', DATA.groundtruth, 'both');
  const grid = document.getElementById('grid');
  const S_GT = DATA.groundtruth.score || 0;
  DATA.igs.forEach((ig, i) => {{
    const div = document.createElement('div');
    const hasMode = !!ig.picked_disp;
    div.className = 'panel' + (hasMode?'':' no-imag');
    const staticBadge = hasMode ? '' : `<span class="badge-static">n_imag=0</span>`;
    const freqStr = ig.picked_freq != null ? ig.picked_freq.toFixed(0)+'i' : '-';
    const sRatio = (S_GT > 1e-9) ? (ig.score / S_GT) : null;
    const sRatioStr = (sRatio != null) ? sRatio.toFixed(3) : '-';
    div.innerHTML = `
      <div class="ph">
        <span class="lbl">${{ig.label}}${{staticBadge}}</span>
        <span class="rk">S = ${{ig.score.toFixed(3)}} &nbsp; (S/S<sub>GT</sub> = ${{sRatioStr}})</span>
      </div>
      <div class="vw"><div id="vw_ig${{i}}" class="vwbox"></div></div>
      <div class="meta">
        <b>&beta;</b>=${{ig.beta.toFixed(3)}} &nbsp;
        <b>&rho;</b>=${{ig.rho.toFixed(3)}} &nbsp;
        <b>&kappa;</b>=${{ig.kappa.toFixed(3)}} &nbsp;
        <b>n_imag</b>=${{ig.n_imag}} &nbsp;
        <b>freq</b>=${{freqStr}} cm⁻¹ &nbsp;
        <b>mwc</b>=${{ig.mwc_to_GT.toFixed(3)}}
      </div>`;
    grid.appendChild(div);
    if (ig.picked_disp) makeAnimated(`vw_ig${{i}}`, ig, ig.picked_disp, 'both');
    else                makeStatic(`vw_ig${{i}}`, ig, 'both');
  }});
}});
</script>
</body></html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--steps", nargs="+", default=None)
    ap.add_argument("--max-branches", type=int, default=1_000_000)
    args = ap.parse_args()

    all_steps = sorted(d.name for d in WORK.iterdir() if d.is_dir())
    if args.steps: steps = [s for s in all_steps if s in set(args.steps)]
    elif args.limit: steps = all_steps[:args.limit]
    else: steps = all_steps
    print(f"Building {len(steps)} BGCP views  workers={args.workers}  max_branches={args.max_branches}")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_ok = 0; n_err = 0
    job_args = [(s, args.max_branches) for s in steps]
    # Build a top-level index
    index_rows = []
    with mp.Pool(args.workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(_safe, job_args), 1):
            err = rec.get("error", "")
            if err:
                print(f"  [{i:>3d}/{len(steps)}] {rec['step']:50s}  ERROR: {err[:60]}", flush=True)
                n_err += 1
            else:
                print(f"  [{i:>3d}/{len(steps)}] {rec['step']:50s}  "
                      f"top1={rec['top1_label']}  S={rec['top1_score']:.3f}  "
                      f"mwc={rec['top1_mwc']:.3f}", flush=True)
                index_rows.append(rec)
                n_ok += 1
    print(f"\n{n_ok} ok, {n_err} errors in {time.time()-t0:.0f}s")

    # Write a top-level index.html
    rows_html = "".join(
        f"<tr><td><a href='{r['step']}/view.html'>{r['step']}</a></td>"
        f"<td>{r['n_ig']}</td>"
        f"<td>{r['top1_label']}</td>"
        f"<td>{r['top1_score']:.3f}</td>"
        f"<td>{r['top1_mwc']:.3f}</td></tr>"
        for r in sorted(index_rows, key=lambda r: r['step'])
    )
    (OUT_ROOT / "index.html").write_text(f"""<!doctype html><html><head><meta charset="utf-8">
<title>BGCP ranked views</title>
<style>body{{font-family:-apple-system,sans-serif;margin:20px;max-width:1200px}}
table{{border-collapse:collapse}} th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}</style>
</head><body>
<h2>BGCP ranked views ({n_ok} steps)</h2>
<table>
<tr><th>step</th><th>n_ig</th><th>top1 label</th><th>top1 S</th><th>top1 mwc</th></tr>
{rows_html}
</table></body></html>""")
    print(f"index: {OUT_ROOT / 'index.html'}")


if __name__ == "__main__":
    main()
