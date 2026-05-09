"""
External-dataset version of build_ranked_view_one_step.py.

Runs the full pipeline from scratch (no cached payload, no groundtruth
required) on a directory shaped like the El Agente Pathways "plain"
mode output:

    <step_dir>/
      source/reactant*.xyz       # one (or more) reactant xyz
      source/product*.xyz        # product xyz
      initial_ts_guesses/        # 20 LLM-generated IG-TS xyz files
      generation_report.json     # used to read charge, multiplicity

Pipeline:
  1. xtb GFN2 single-point on R, P -> WBO matrices
  2. PQ alignment R <-> P -> atom mapping + broken / formed bonds
  3. core_atoms in R-frame; bond-reaction vector V at TS coords
  4. xtb hess on each IG (parallel) -> g98.out -> normal modes
  5. align each IG <-> R, reindex modes to R-frame
  6. per-mode features: bond_overlap (beta), rxn_overlap (rho),
     core_fraction (kappa); pick max-bond_overlap imag mode per IG
  7. score = beta * (1 + w_r * rho) * (1 + w_c * kappa) / n_imag^p,
     sort descending. Every IG with >=1 imaginary mode gets animated
     on its picked mode; IGs with n_imag == 0 render as static.
     No n_imag-count or rho filter on the viewer side.

Caching: all xtb output lands in
  work_modes/<workflow_name>/{R, P, hess_iter<i>}/
so re-runs are cache hits (xtb is the dominant cost).

Usage:
  python viewer/build_ranked_view_external.py <step_dir> [step_dir ...]

  e.g.
  python viewer/build_ranked_view_external.py \\
      /Users/.../run_20260508_202120_plain_parallel/{1a,1b,1c}

Output:
  out/ranked_views/<workflow_name>.html
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent

import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from rxn_core_pq import align_from_arrays
from rxn_core_frag import run_xtb, parse_xyz
from analyze_core_modes import (
    parse_g98_modes, core_atoms_in_R_frame, reindex_modes_to_R,
    bond_reaction_vector, bond_overlap_per_mode,
    rxn_overlap_per_mode, reaction_coord_delta,
)
from align_bgcp_coords import fill_unmapped_greedy


WORK_MODES = PROJECT_ROOT / "work_modes"
OUT_DIR = PROJECT_ROOT / "out" / "ranked_views"

# Score-formula hyperparameters (same shape as rk_clean_v2 score, but
# applied without the n_imag<=2 / rho>=0.10 filters — the viewer shows
# every IG that has any imaginary mode, ranked by score).
W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3

# xtb parallelism: each worker uses OMP_NUM_THREADS, and we run N_WORKERS
# in parallel. 4 workers x 4 threads = 16, fine on a 14-core box.
N_WORKERS = 4
OMP_THREADS = 4


def parse_xyz_file(path: Path):
    el, xyz = parse_xyz(path)
    return el, xyz


def read_first_xyz(d: Path):
    files = sorted(d.glob("*.xyz"))
    if not files:
        raise FileNotFoundError(f"no xyz in {d}")
    return parse_xyz_file(files[0])


def run_xtb_hess(xyz_path: Path, workdir: Path, charge: int = 0, uhf: int = 0):
    """Run `xtb input.xyz --gfn 2 --hess`. Cached: skip if g98.out + wbo
    are already there with a matching xyz copy. Returns (elements,
    coords, wbo, freqs, modes_TS) where freqs is (n_modes,) and
    modes_TS is (n_modes, n_atoms, 3)."""
    workdir.mkdir(parents=True, exist_ok=True)
    local = workdir / Path(xyz_path).name
    src_text = Path(xyz_path).read_text()
    cached_text = local.read_text() if local.exists() else None
    g98 = workdir / "g98.out"
    wbo = workdir / "wbo"
    cached = (cached_text == src_text) and g98.exists() and wbo.exists()
    if not cached:
        shutil.copy(xyz_path, local)
        cmd = ["xtb", local.name, "--gfn", "2", "--hess"]
        if charge: cmd += ["--chrg", str(charge)]
        if uhf:    cmd += ["--uhf",  str(uhf)]
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(OMP_THREADS)
        res = subprocess.run(cmd, cwd=workdir, capture_output=True,
                             text=True, env=env)
        if res.returncode != 0:
            raise RuntimeError(f"xtb hess failed in {workdir}: "
                               f"{res.stderr[-500:]}")
        if not g98.exists() or not wbo.exists():
            raise RuntimeError(f"missing g98.out or wbo in {workdir}")
    elements, coords = parse_xyz(local)
    n = len(elements)
    wbo_arr = np.zeros((n, n))
    for ln in wbo.read_text().splitlines():
        parts = ln.split()
        if len(parts) < 3: continue
        i, j = int(parts[0]) - 1, int(parts[1]) - 1
        v = float(parts[2])
        wbo_arr[i, j] = v; wbo_arr[j, i] = v
    freqs, modes_TS = parse_g98_modes(g98)
    return elements, coords, wbo_arr, freqs, modes_TS


def _hess_worker(args):
    """Standalone for ProcessPoolExecutor."""
    xyz_path, workdir, charge, uhf = args
    t0 = time.time()
    try:
        elT, xyzT, wboT, freqs, modes_TS = run_xtb_hess(
            Path(xyz_path), Path(workdir), charge, uhf)
        return {"ok": True, "xyz_path": str(xyz_path),
                "secs": time.time() - t0,
                "elT": elT, "xyzT": xyzT.tolist() if hasattr(xyzT, 'tolist') else xyzT,
                "wboT": wboT.tolist(), "freqs": freqs.tolist(),
                "modes_TS": modes_TS.tolist()}
    except Exception as e:
        return {"ok": False, "xyz_path": str(xyz_path),
                "secs": time.time() - t0, "error": str(e)}


def label_for_ig(path: Path):
    m = re.search(r"_iter(\d+)_", path.name)
    return f"iter{m.group(1)}" if m else path.stem


def process_step(step_dir: Path):
    rep = json.loads((step_dir / "generation_report.json").read_text())
    spec = rep["generation_spec"]
    workflow_name = spec["workflow_name"]
    charge = int(spec.get("charge", 0))
    mult = int(spec.get("multiplicity", 1))
    uhf = max(0, mult - 1)
    print(f"\n=== {workflow_name}  charge={charge} mult={mult} ===", flush=True)

    cache = WORK_MODES / workflow_name
    cache.mkdir(parents=True, exist_ok=True)

    # 1. R, P single-points
    rxyz_path = sorted((step_dir / "source").glob("reactant*.xyz"))[0]
    pxyz_path = sorted((step_dir / "source").glob("product*.xyz"))[0]
    print(f"  xtb sp on R, P ...", flush=True)
    t0 = time.time()
    elR, xyzR, wboR = run_xtb(rxyz_path, cache / "R", charge=charge, uhf=uhf)
    elP, xyzP, wboP = run_xtb(pxyz_path, cache / "P", charge=charge, uhf=uhf)
    print(f"    sp done in {time.time()-t0:.1f}s ({len(elR)} atoms each)", flush=True)

    # 2. R<->P alignment + core atoms + bond-reaction vector at R-coords
    print(f"  PQ alignment R<->P ...", flush=True)
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
    print(f"    mapping={len(mapping_RP)}/{len(elR)}, broken={len(broken_R)}, "
          f"formed={len(formed_R)}, core_atoms={len(core_R)}", flush=True)

    # 3. xtb hess on each IG, in parallel
    ig_paths = sorted((step_dir / "initial_ts_guesses").glob("*.xyz"))
    print(f"  xtb hess on {len(ig_paths)} IGs (workers={N_WORKERS}, omp={OMP_THREADS}) ...",
          flush=True)
    t0 = time.time()
    args = [(str(p), str(cache / f"hess_{label_for_ig(p)}"), charge, uhf)
            for p in ig_paths]
    results = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futures = {ex.submit(_hess_worker, a): a for a in args}
        for done in as_completed(futures):
            r = done.result()
            results.append(r)
            tag = "ok" if r["ok"] else f"FAIL({r.get('error','?')[:60]})"
            print(f"    [{r['secs']:>5.1f}s]  {Path(r['xyz_path']).name}  {tag}",
                  flush=True)
    print(f"  all hess done in {time.time()-t0:.1f}s", flush=True)

    # 4. Per-IG: align IG<->R, reindex modes, compute features, score
    n_R = len(elR)
    bond_V = None  # built per-IG using IG coords in R-frame
    ig_records = []
    for r in sorted(results, key=lambda r: r["xyz_path"]):
        label = label_for_ig(Path(r["xyz_path"]))
        if not r["ok"]:
            ig_records.append({
                "label": label, "score": 0.0,
                "beta": 0.0, "rho": 0.0, "kappa": 0.0, "n_imag": 0,
                "picked_freq": None, "picked_disp": None,
                "xyz_elements": elR, "xyz_coords": xyzR.tolist() if hasattr(xyzR,'tolist') else xyzR,
                "error": r.get("error", ""),
            })
            continue
        elT = r["elT"]
        xyzT = np.asarray(r["xyzT"], float)
        wboT = np.asarray(r["wboT"], float)
        freqs = np.asarray(r["freqs"], float)
        modes_TS = np.asarray(r["modes_TS"], float)

        # Align IG<->R
        try:
            it = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT)
            mapping_RT = dict(it["mapping"])
            mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, mapping_RT)
        except Exception as e:
            print(f"    align fail for {label}: {e}")
            ig_records.append({
                "label": label, "score": 0.0,
                "beta": 0.0, "rho": 0.0, "kappa": 0.0, "n_imag": 0,
                "picked_freq": None, "picked_disp": None,
                "xyz_elements": elT, "xyz_coords": xyzT.tolist(),
                "error": f"align: {e}",
            })
            continue

        modes_R = reindex_modes_to_R(modes_TS, mapping_RT, n_R)
        sq = (modes_R ** 2).sum(axis=2)
        total = sq.sum(axis=1)
        core_e = sq[:, core_R].sum(axis=1) if core_R else np.zeros(modes_R.shape[0])
        kappa = np.where(total > 1e-12, core_e / total, 0.0)
        rho = rxn_overlap_per_mode(modes_R, delta_RP, core_R)
        # IG xyz expressed in R-frame (for bond-reaction vector)
        ts_xyz_in_R = np.zeros_like(np.asarray(xyzR))
        for r_idx, t_idx in mapping_RT.items():
            ts_xyz_in_R[r_idx] = xyzT[t_idx]
        V = bond_reaction_vector(ts_xyz_in_R, broken_R, formed_R)
        beta = bond_overlap_per_mode(modes_R, V)

        imag_idx = list(np.where(freqs < 0)[0])
        n_imag = len(imag_idx)
        if n_imag == 0:
            ig_records.append({
                "label": label, "score": 0.0,
                "beta": 0.0, "rho": 0.0, "kappa": 0.0, "n_imag": 0,
                "picked_freq": None, "picked_disp": None,
                "xyz_elements": elT, "xyz_coords": xyzT.tolist(),
            })
            continue
        # Pick max-bond_overlap imag mode (in R-frame indexing).
        # NO n_imag-count filter: every IG with >=1 imag mode gets ranked
        # and animated; the score formula already down-weights large n_imag
        # via the / n_imag^p term.
        picked_k = max(imag_idx, key=lambda k: beta[k])
        b = float(beta[picked_k]); r_ = float(rho[picked_k]); c = float(kappa[picked_k])
        score = b * (1 + W_RXN * r_) * (1 + W_CORE * c) / max(n_imag, 1) ** IMAG_PEN
        # Render the picked mode in R-frame on R coords (so animation
        # uses consistent atom indexing across panels).
        ig_records.append({
            "label": label,
            "score": float(score),
            "beta": b, "rho": r_, "kappa": c,
            "n_imag": int(n_imag),
            "picked_freq": float(freqs[picked_k]),
            "picked_disp": modes_R[picked_k].tolist(),
            "xyz_elements": elR,
            "xyz_coords":   ts_xyz_in_R.tolist(),
        })

    # 5. Sort by score descending. IGs with no imag mode (score=0) sink
    # to the bottom naturally.
    ig_records.sort(key=lambda x: -x["score"])

    # 6. Build payload + render HTML
    data = {
        "step": workflow_name,
        "n_atoms": n_R,
        "core_atoms": list(core_R),
        "broken_bonds":  [list(b) for b in broken_R],
        "formed_bonds_R":[list(b) for b in formed_R],
        "reactant":  {"xyz_elements": elR,
                      "xyz_coords":   xyzR.tolist() if hasattr(xyzR,'tolist') else xyzR},
        "product":   {"xyz_elements": elP,
                      "xyz_coords":   xyzP.tolist() if hasattr(xyzP,'tolist') else xyzP},
        "igs": ig_records,
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{workflow_name}.html"
    out_path.write_text(HTML.format(
        title=f"Ranked view (external) — {workflow_name}",
        n_atoms=n_R,
        n_broken=len(broken_R),
        n_formed=len(formed_R),
        n_core=len(core_R),
        data_json=json.dumps(data),
    ))
    has_mode = sum(1 for ig in ig_records if ig["n_imag"] > 0)
    print(f"  IGs with at least one imag mode: {has_mode}/{len(ig_records)}", flush=True)
    print(f"  top-3 by score:")
    for ig in ig_records[:3]:
        print(f"    {ig['label']:>8s}  S={ig['score']:.3f}  "
              f"beta={ig['beta']:.3f}  rho={ig['rho']:.3f}  "
              f"kappa={ig['kappa']:.3f}  n_imag={ig['n_imag']}", flush=True)
    print(f"  wrote {out_path}", flush=True)
    return out_path


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>{title}</title>
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
 .ref-row {{ display:grid; grid-template-columns:repeat(2,1fr); gap:12px; margin-bottom:18px; }}
 .ig-grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; }}
 .panel  {{ background:white; border:1px solid #ddd; border-radius:6px; padding:6px 8px 8px; }}
 .panel.no-imag {{ background:#f7f7f7; }}
 .ph    {{ display:flex; justify-content:space-between; align-items:baseline; font-size:12px; margin-bottom:4px; }}
 .ph .lbl {{ font-weight:600; font-size:13px; }}
 .ph .rk  {{ font-family:ui-monospace,monospace; color:#024; }}
 .panel.no-imag .ph .rk {{ color:#888; }}
 .vw    {{ position:relative; width:100%; height:230px; }}
 .ref-row .vw {{ height:300px; }}
 .vwbox {{ position:absolute; inset:0; }}
 .meta  {{ font-family:ui-monospace,monospace; font-size:11px; color:#444; padding:3px 0 0; line-height:1.4; }}
 .meta b {{ color:#024; }}
 .badge-static {{ display:inline-block; background:#eee; color:#666; padding:1px 6px; border-radius:3px; font-size:11px; margin-left:6px; }}
</style>
</head>
<body>
<h2>{title}</h2>
<div class="sub">
  N atoms: <b>{n_atoms}</b> &nbsp;|&nbsp;
  broken: <b>{n_broken}</b> &nbsp;|&nbsp;
  formed: <b>{n_formed}</b> &nbsp;|&nbsp;
  core atoms: <b>{n_core}</b><br>
  <span class="legend">
    <span class="leg-broken">broken bond</span>
    <span class="leg-formed">formed bond</span>
    <span class="leg-mode">mode arrow (core atoms)</span>
    <span class="leg-static">n_imag = 0 (rendered static)</span>
  </span>
  <br>
  IGs are sorted by score
  <code>S = &beta; (1 + w_r &rho;) (1 + w_c &kappa;) / n_imag^p</code>
  with <code>w_r=1.0, w_c=0.2, p=0.3</code>, descending. Every IG with
  &ge;1 imaginary mode is animated on its picked mode (max-&beta; imag);
  IGs with <code>n_imag = 0</code> render as static structures.
</div>

<div class="ref-row">
  <div class="panel"><div class="ph"><span class="lbl">Reactant</span><span class="rk">static</span></div>
    <div class="vw"><div id="vw_R" class="vwbox"></div></div></div>
  <div class="panel"><div class="ph"><span class="lbl">Product</span><span class="rk">static</span></div>
    <div class="vw"><div id="vw_P" class="vwbox"></div></div></div>
</div>

<div class="ig-grid" id="grid"></div>

<script>
const DATA = {data_json};

function buildBody(elements, xyz) {{
  const n = xyz.length;
  let s = `${{n}}\nframe\n`;
  for (let i = 0; i < n; i++) {{
    s += `${{elements[i]}}  ${{xyz[i][0].toFixed(6)}}  ${{xyz[i][1].toFixed(6)}}  ${{xyz[i][2].toFixed(6)}}\n`;
  }}
  return s;
}}
function buildBodyAt(elements, xyz, disp, scale) {{
  const n = xyz.length;
  let s = `${{n}}\nframe\n`;
  for (let i = 0; i < n; i++) {{
    const x = xyz[i][0] + scale * disp[i][0];
    const y = xyz[i][1] + scale * disp[i][1];
    const z = xyz[i][2] + scale * disp[i][2];
    s += `${{elements[i]}}  ${{x.toFixed(6)}}  ${{y.toFixed(6)}}  ${{z.toFixed(6)}}\n`;
  }}
  return s;
}}
function decorateBonds(viewer) {{
  const atoms = viewer.selectedAtoms({{}});
  for (const [i, j] of DATA.broken_bonds) {{
    if (i < atoms.length && j < atoms.length) {{
      const a = atoms[i], b = atoms[j];
      viewer.addCylinder({{start:{{x:a.x,y:a.y,z:a.z}}, end:{{x:b.x,y:b.y,z:b.z}},
                          color:'red', radius:0.08, dashed:true}});
    }}
  }}
  for (const [i, j] of DATA.formed_bonds_R) {{
    if (i < atoms.length && j < atoms.length) {{
      const a = atoms[i], b = atoms[j];
      viewer.addCylinder({{start:{{x:a.x,y:a.y,z:a.z}}, end:{{x:b.x,y:b.y,z:b.z}},
                          color:'green', radius:0.08, dashed:true}});
    }}
  }}
}}
function drawArrows(viewer, xyz, disp) {{
  const atoms = viewer.selectedAtoms({{}});
  for (const i of DATA.core_atoms) {{
    if (i >= atoms.length || !disp || !disp[i]) continue;
    const d = disp[i];
    const len = Math.hypot(d[0], d[1], d[2]);
    if (len < 1e-3) continue;
    const a = atoms[i];
    viewer.addArrow({{
      start:{{x:a.x, y:a.y, z:a.z}},
      end:  {{x:a.x + d[0]*1.5, y:a.y + d[1]*1.5, z:a.z + d[2]*1.5}},
      color:'#0066cc', radius:0.06,
    }});
  }}
}}
function makeStatic(divId, ts) {{
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(buildBody(ts.xyz_elements, ts.xyz_coords), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  decorateBonds(v);
  v.zoomTo();
  v.render();
  return v;
}}
function makeAnimated(divId, ts, disp) {{
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(buildBody(ts.xyz_elements, ts.xyz_coords), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  decorateBonds(v); drawArrows(v, ts.xyz_coords, disp);
  v.zoomTo(); v.render();
  let t = 0; const period = 30; const amp = 0.6;
  setInterval(() => {{
    t = (t + 1) % period;
    const scale = amp * Math.sin(2 * Math.PI * t / period);
    v.removeAllModels(); v.removeAllShapes();
    v.addModel(buildBodyAt(ts.xyz_elements, ts.xyz_coords, disp, scale), 'xyz');
    v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
    decorateBonds(v); drawArrows(v, ts.xyz_coords, disp);
    v.render();
  }}, 60);
  return v;
}}

window.addEventListener('load', () => {{
  makeStatic('vw_R', DATA.reactant);
  makeStatic('vw_P', DATA.product);
  const grid = document.getElementById('grid');
  DATA.igs.forEach((ig, i) => {{
    const div = document.createElement('div');
    const hasMode = !!ig.picked_disp;
    div.className = 'panel' + (hasMode ? '' : ' no-imag');
    const staticBadge = hasMode ? '' :
      `<span class="badge-static">n_imag = 0</span>`;
    const freqStr = ig.picked_freq != null ? ig.picked_freq.toFixed(0) + 'i' : '—';
    div.innerHTML = `
      <div class="ph">
        <span class="lbl">${{ig.label}}${{staticBadge}}</span>
        <span class="rk">S = ${{ig.score.toFixed(3)}}</span>
      </div>
      <div class="vw"><div id="vw_ig${{i}}" class="vwbox"></div></div>
      <div class="meta">
        <b>&beta;</b>=${{ig.beta.toFixed(3)}} &nbsp;
        <b>&rho;</b>=${{ig.rho.toFixed(3)}} &nbsp;
        <b>&kappa;</b>=${{ig.kappa.toFixed(3)}} &nbsp;
        <b>n_imag</b>=${{ig.n_imag}} &nbsp;
        <b>freq</b>=${{freqStr}} cm⁻¹
      </div>`;
    grid.appendChild(div);
    if (ig.picked_disp) {{
      makeAnimated(`vw_ig${{i}}`, ig, ig.picked_disp);
    }} else {{
      makeStatic(`vw_ig${{i}}`, ig);
    }}
  }});
}});
</script>
</body></html>
"""


def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    for arg in sys.argv[1:]:
        step_dir = Path(arg)
        if not step_dir.is_dir():
            print(f"  WARN: skip {step_dir} (not a directory)")
            continue
        try:
            process_step(step_dir)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  ERROR processing {step_dir}: {e}")


if __name__ == "__main__":
    main()
