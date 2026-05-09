"""
Single-step end-to-end view for El Agente Disco.

Pipeline (one step):
  1. atom alignment (R <-> P)              [rxn_core_pq, cached]
  2. core-atom identification               [classify_bonds]
  3. xtb Hessian + normal modes per IG     [cached parses]
  4. picked imaginary mode per IG           [max bond_overlap, n_imag<=2 filter]
  5. rk_clean_v2 score per IG               [b * (1+wr*r) * (1+wc*c) / n_imag^p]
  6. one HTML page: R | P | GT (top row) +
     20 IG panels in a grid, sorted by score descending. Every IG with
     >=1 imaginary mode is animated on its picked mode; IGs with
     n_imag == 0 render as static structures. No n_imag-count or rho
     filter. Each panel shows score, beta, rho, kappa, n_imag.

Reuses the per-step mode_viewer payload at
  appendix_perparation/viewer/mode_viewer/<step>.html
which already contains all 20 IGs + groundtruth with their parsed modes
and reactive-bond lists.

Usage:
  python viewer/build_ranked_view_one_step.py [STEP_NAME]
  default STEP_NAME = pr16.carbocation_ts11  (18 atoms, clean 1/1 bonds)

Output:
  out/ranked_views/<step>.html
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent

import json
import re
import sys
from pathlib import Path

import numpy as np

from rxn_core_pq import align_from_arrays
from align_bgcp_coords import load_cached_xtb, fill_unmapped_greedy


MODE_VIEWER_DIR = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
BGCP_ROOT = PROJECT_ROOT / 'appendix_perparation' / 'Pure_Geometries_Elementary_Step' / 'Benchmark_Guesses_Collective_Package'
WORK_MODES = PROJECT_ROOT / 'work_modes'
OUT_DIR = PROJECT_ROOT / 'out' / 'ranked_views'

# Score-formula hyperparameters (same shape as rk_clean_v2 score, but
# applied without the n_imag<=2 / rho>=0.10 filters — the viewer shows
# every IG with any imag mode, ranked by score).
W_RXN, W_CORE, IMAG_PEN = 1.0, 0.2, 0.3


def read_xyz(path: Path):
    lines = path.read_text().strip().splitlines()
    n = int(lines[0])
    el, xyz = [], []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        el.append(parts[0])
        xyz.append([float(x) for x in parts[1:4]])
    return el, xyz


def load_step_payload(step):
    html = MODE_VIEWER_DIR / f'{step}.html'
    if not html.exists():
        raise SystemExit(f"per-step payload not found: {html}\n"
                         f"Run viewer/build_mode_viewer.py first to populate "
                         f"appendix_perparation/viewer/mode_viewer/.")
    text = html.read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m:
        raise SystemExit(f"could not find DATA block in {html}")
    return json.loads(m.group(1))


def imag_modes(ts):
    return [m for m in ts.get('modes', []) if m.get('freq', 0) < 0]


def pick_and_score(ts):
    """Return (picked_mode, score, b, r, c, n_imag).
    picked_mode is the max-bond_overlap imag mode; None if no imag modes."""
    imag = imag_modes(ts)
    n_imag = len(imag)
    if not imag:
        return None, 0.0, 0.0, 0.0, 0.0, 0
    picked = max(imag, key=lambda m: m.get('bond_overlap', 0.0))
    b = picked.get('bond_overlap', 0.0)
    r = picked.get('rxn_overlap',  0.0)
    c = picked.get('core_fraction',0.0)
    score = b * (1 + W_RXN * r) * (1 + W_CORE * c) / max(n_imag, 1) ** IMAG_PEN
    return picked, score, b, r, c, n_imag


def build_view_data(step):
    payload = load_step_payload(step)
    n_atoms = payload['n_atoms']

    gt = next((t for t in payload['ts_list']
               if t['label'] == 'groundtruth' and t.get('modes')), None)
    if gt is None:
        raise SystemExit(f"no groundtruth with modes in {step}")
    igs = [t for t in payload['ts_list']
           if t['label'] != 'groundtruth']

    # Score every IG
    ig_records = []
    for ig in igs:
        picked, score, b, r, c, n_imag = pick_and_score(ig)
        ig_records.append({
            'label': ig['label'],
            'xyz_elements': ig['xyz_elements'],
            'xyz_coords': ig['xyz_coords'],
            'picked_disp': picked['disp'] if picked else None,
            'picked_freq': float(picked['freq']) if picked else None,
            'score': score,
            'beta': b,
            'rho':  r,
            'kappa': c,
            'n_imag': n_imag,
        })
    # Sort by score descending. IGs with no imag mode (score=0) sink
    # to the bottom naturally.
    ig_records.sort(key=lambda x: -x['score'])

    # GT picked mode (for the GT panel)
    gt_picked, _, _, _, _, gt_n_imag = pick_and_score(gt)

    # R, P xyz
    # Load cached xtb output for R and P (created by build_mode_viewer.py)
    # so we can re-derive mapping_RP and reindex P into R-frame for display.
    # Without this, the P panel would draw bond cylinders at wrong atom
    # positions (formed_bonds_R is in R-frame indices; raw P xyz is in
    # P-frame).
    cache = WORK_MODES / step
    elR, xyzR_arr, wboR, _ = load_cached_xtb(cache / 'R')
    elP, xyzP_arr, wboP, _ = load_cached_xtb(cache / 'P')
    rp = align_from_arrays(elR, xyzR_arr, wboR, elP, xyzP_arr, wboP)
    full_RP = fill_unmapped_greedy(elR, xyzR_arr, elP, xyzP_arr, dict(rp['mapping']))
    xyzP_in_R = np.zeros_like(np.asarray(xyzR_arr, float))
    for i_R, i_P in full_RP.items():
        xyzP_in_R[i_R] = xyzP_arr[i_P]

    return {
        'step': step,
        'n_atoms': n_atoms,
        'core_atoms': payload.get('core_atoms', []),
        'broken_bonds': payload.get('broken_bonds', []),
        'formed_bonds_R': payload.get('formed_bonds_R', []),
        'reactant':   {'xyz_elements': elR,
                       'xyz_coords':   xyzR_arr.tolist() if hasattr(xyzR_arr,'tolist') else xyzR_arr},
        # P shown in R-frame so atom indices line up with broken/formed pairs
        'product':    {'xyz_elements': elR,
                       'xyz_coords':   xyzP_in_R.tolist()},
        'groundtruth': {
            'xyz_elements': gt['xyz_elements'],
            'xyz_coords':   gt['xyz_coords'],
            'picked_disp':  gt_picked['disp'] if gt_picked else None,
            'picked_freq':  float(gt_picked['freq']) if gt_picked else None,
            'n_imag':       gt_n_imag,
        },
        'igs': ig_records,
    }


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
 .ref-row {{ display:grid; grid-template-columns:repeat(3,1fr); gap:12px; margin-bottom:18px; }}
 .ig-grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; }}
 .panel  {{ background:white; border:1px solid #ddd; border-radius:6px; padding:6px 8px 8px; }}
 .panel.no-imag {{ background:#f7f7f7; }}
 .ph    {{ display:flex; justify-content:space-between; align-items:baseline; font-size:12px; margin-bottom:4px; }}
 .ph .lbl {{ font-weight:600; font-size:13px; }}
 .ph .rk  {{ font-family:ui-monospace,monospace; color:#024; }}
 .panel.no-imag .ph .rk {{ color:#888; }}
 .vw    {{ position:relative; width:100%; height:230px; }}
 .ref-row .vw {{ height:280px; }}
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
  <div class="panel"><div class="ph"><span class="lbl">Ground-truth TS</span>
    <span class="rk">freq {gt_freq_str} cm&#x207B;&#xB9;</span></div>
    <div class="vw"><div id="vw_GT" class="vwbox"></div></div>
    <div class="meta">n_imag = <b>{gt_n_imag}</b></div></div>
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
function xyzAt(xyz, disp, scale) {{
  const out = new Array(xyz.length);
  for (let i = 0; i < xyz.length; i++) {{
    out[i] = [
      xyz[i][0] + scale * disp[i][0],
      xyz[i][1] + scale * disp[i][1],
      xyz[i][2] + scale * disp[i][2],
    ];
  }}
  return out;
}}

// which: 'broken' (R panel only) | 'formed' (P panel only) | 'both' (TS-like)
// Uses xyz coords directly rather than viewer.selectedAtoms({{}}) so the
// cylinders can be added before the first render() (selectedAtoms can
// return empty if called before the viewer is initialized).
function decorateBonds(viewer, which, xyz) {{
  if (which !== 'formed') {{
    for (const [i, j] of DATA.broken_bonds) {{
      if (i < xyz.length && j < xyz.length) {{
        viewer.addCylinder({{
          start:{{x:xyz[i][0], y:xyz[i][1], z:xyz[i][2]}},
          end:  {{x:xyz[j][0], y:xyz[j][1], z:xyz[j][2]}},
          color:'red', radius:0.08, dashed:true,
        }});
      }}
    }}
  }}
  if (which !== 'broken') {{
    for (const [i, j] of DATA.formed_bonds_R) {{
      if (i < xyz.length && j < xyz.length) {{
        viewer.addCylinder({{
          start:{{x:xyz[i][0], y:xyz[i][1], z:xyz[i][2]}},
          end:  {{x:xyz[j][0], y:xyz[j][1], z:xyz[j][2]}},
          color:'green', radius:0.08, dashed:true,
        }});
      }}
    }}
  }}
}}
function drawArrows(viewer, xyz, disp) {{
  for (const i of DATA.core_atoms) {{
    if (i >= xyz.length || !disp || !disp[i]) continue;
    const d = disp[i];
    const len = Math.hypot(d[0], d[1], d[2]);
    if (len < 1e-3) continue;
    const ax = xyz[i][0], ay = xyz[i][1], az = xyz[i][2];
    viewer.addArrow({{
      start:{{x:ax, y:ay, z:az}},
      end:  {{x:ax + d[0]*1.5, y:ay + d[1]*1.5, z:az + d[2]*1.5}},
      color:'#0066cc', radius:0.06,
    }});
  }}
}}

// Static (R, P).  `which` selects 'broken' | 'formed' | 'both'.
function makeStatic(divId, ts, which) {{
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(buildBody(ts.xyz_elements, ts.xyz_coords), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  decorateBonds(v, which || 'both', ts.xyz_coords);
  v.zoomTo();
  v.render();
  return v;
}}

// Animated panel (GT, IG-with-picked-mode)
function makeAnimated(divId, ts, disp, which) {{
  const w = which || 'both';
  const v = $3Dmol.createViewer(divId, {{backgroundColor:'white'}});
  v.addModel(buildBody(ts.xyz_elements, ts.xyz_coords), 'xyz');
  v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
  decorateBonds(v, w, ts.xyz_coords);
  drawArrows(v, ts.xyz_coords, disp);
  v.zoomTo();
  v.render();

  // animation loop
  let t = 0;
  const period = 30;        // frames per cycle
  const amp = 0.6;
  setInterval(() => {{
    t = (t + 1) % period;
    const scale = amp * Math.sin(2 * Math.PI * t / period);
    const cur = xyzAt(ts.xyz_coords, disp, scale);
    v.removeAllModels();
    v.removeAllShapes();
    v.addModel(buildBodyAt(ts.xyz_elements, ts.xyz_coords, disp, scale), 'xyz');
    v.setStyle({{}}, {{stick:{{radius:0.10}}, sphere:{{scale:0.20}}}});
    decorateBonds(v, w, cur);
    drawArrows(v, cur, disp);
    v.render();
  }}, 60);
  return v;
}}

window.addEventListener('load', () => {{
  // Reference row.
  // R panel: only the bonds that are BREAKING (red dashed) -- exist
  //          in R now, about to disappear.
  // P panel: only the bonds that just FORMED (green dashed) -- exist
  //          in P. P xyz is reindexed to R-frame above so the
  //          formed_bonds_R pairs land on the right atoms.
  makeStatic('vw_R', DATA.reactant, 'broken');
  makeStatic('vw_P', DATA.product,  'formed');
  if (DATA.groundtruth.picked_disp) {{
    makeAnimated('vw_GT', DATA.groundtruth, DATA.groundtruth.picked_disp);
  }} else {{
    makeStatic('vw_GT', DATA.groundtruth);
  }}

  // IG grid
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
    step = sys.argv[1] if len(sys.argv) > 1 else 'pr16.carbocation_ts11'
    print(f"Building one-step ranked view for: {step}")
    data = build_view_data(step)
    print(f"  n_atoms={data['n_atoms']}, "
          f"broken={len(data['broken_bonds'])}, formed={len(data['formed_bonds_R'])}, "
          f"core_atoms={len(data['core_atoms'])}, "
          f"IGs={len(data['igs'])}")

    has_mode = sum(1 for ig in data['igs'] if ig['n_imag'] > 0)
    print(f"  IGs with at least one imag mode: {has_mode}/{len(data['igs'])}")
    print(f"  top-3 by score:")
    for ig in data['igs'][:3]:
        print(f"    {ig['label']:>8s}  S={ig['score']:.3f}  "
              f"beta={ig['beta']:.3f}  rho={ig['rho']:.3f}  "
              f"kappa={ig['kappa']:.3f}  n_imag={ig['n_imag']}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / f"{step}.html"
    gt_freq_str = (f"{data['groundtruth']['picked_freq']:.0f}i"
                   if data['groundtruth']['picked_freq'] else "—")
    html = HTML.format(
        title=f"Ranked view — {step}",
        n_atoms=data['n_atoms'],
        n_broken=len(data['broken_bonds']),
        n_formed=len(data['formed_bonds_R']),
        n_core=len(data['core_atoms']),
        gt_freq_str=gt_freq_str,
        gt_n_imag=data['groundtruth']['n_imag'],
        data_json=json.dumps(data),
    )
    out_path.write_text(html)
    print(f"\nwrote {out_path}  ({out_path.stat().st_size/1e6:.2f} MB)")


if __name__ == '__main__':
    main()
