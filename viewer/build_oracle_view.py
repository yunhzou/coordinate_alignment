"""
Oracle viewer: GT + top-5 initial-guess TS picked by gt_alignment
(oracle ranking), each IG showing the mode within it that aligns
best with the GT mode.

Per panel label includes:
  - the gt_alignment value (cosine similarity in 3N-Cartesian space)
  - the mode's freq, bond_overlap, rxn_overlap, core_fraction
  - whether the mode is imaginary

Output: out/mode_viewer/oracle_view.html

Reuses cached payload from out/mode_viewer/<step>.html (no recompute).
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np


def kabsch(P, Q):
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    Pc = P.mean(0); Qc = Q.mean(0)
    H = (Q - Qc).T @ (P - Pc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Pc - R @ Qc
    return R, t


def cos_sim_disp(a, b):
    a = np.asarray(a).reshape(-1)
    b = np.asarray(b).reshape(-1)
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-9 or nb < 1e-9: return 0.0
    return abs(float(a @ b)) / (na * nb)


SRC_DIR = PROJECT_ROOT / "out" / "mode_viewer"


def load_step_payload(html_path):
    text = html_path.read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m:
        raise ValueError(f"could not find DATA= in {html_path}")
    return json.loads(m.group(1))


def best_aligned_mode(ts, gt_disp):
    """The mode within `ts` whose displacement aligns best with gt_disp.
    Considered: ALL modes (not restricted to imaginary)."""
    best = None; best_a = -1
    for m in ts['modes']:
        a = cos_sim_disp(np.asarray(m['disp']), gt_disp)
        if a > best_a:
            best_a = a; best = m
    return best, best_a


def build_oracle_payload(step_payload, top_k=5):
    ts_list = step_payload['ts_list']
    by_label = {ts['label']: ts for ts in ts_list}
    gt = by_label.get('groundtruth')
    if gt is None:
        return None
    gt_default_idx = gt.get('default_mode_idx', 0)
    gt_mode = gt['modes'][gt_default_idx]
    gt_disp = np.asarray(gt_mode['disp'])

    # Score every IG by best gt_alignment over its modes
    igs = []
    for ts in ts_list:
        if ts['label'] == 'groundtruth' or not ts.get('modes'):
            continue
        m, a = best_aligned_mode(ts, gt_disp)
        if m is None: continue
        igs.append((a, ts, m))
    if len(igs) < 1: return None
    igs.sort(key=lambda t: -t[0])
    top_igs = igs[:top_k]

    chosen = [(gt, gt_mode, 1.0)]  # GT vs itself = 1.0
    for a, ts, m in top_igs:
        chosen.append((ts, m, a))

    def to_panel(ts, m, gt_align):
        return {
            'label': ts['label'],
            'mode_idx': m['idx'],
            'freq': m['freq'],
            'bond_overlap': m.get('bond_overlap', 0.0),
            'rxn_overlap': m.get('rxn_overlap', 0.0),
            'core_fraction': m.get('core_fraction', 0.0),
            'gt_align': gt_align,
            'is_imag': m['freq'] < 0,
            'n_imag': ts['n_imag'],
            'n_modes_total': ts['n_modes_total'],
            'xyz_elements': ts['xyz_elements'],
            'xyz_coords': ts['xyz_coords'],
            'disp': m['disp'],
        }

    panels = [to_panel(ts, m, a) for (ts, m, a) in chosen]

    # Kabsch-align IG panels to GT panel (operates on R-frame coords).
    # Coords already share atom indexing. Apply rotation to disp too.
    target_xyz = np.asarray(panels[0]['xyz_coords'])
    for p in panels[1:]:
        Q = np.asarray(p['xyz_coords'])
        R, t = kabsch(target_xyz, Q)
        new_xyz = (Q @ R.T) + t
        new_disp = np.asarray(p['disp']) @ R.T
        p['xyz_coords'] = [[round(float(x), 4) for x in v] for v in new_xyz]
        p['disp']       = [[round(float(x), 4) for x in v] for v in new_disp]

    return {
        'step': step_payload['step'],
        'n_atoms': step_payload['n_atoms'],
        'core_atoms': step_payload['core_atoms'],
        'broken_bonds': step_payload['broken_bonds'],
        'formed_bonds_R': step_payload['formed_bonds_R'],
        'panels': panels,
    }


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Oracle view — GT + top-5 IGs by gt_alignment</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select { padding: 4px 6px; font-size: 13px; min-width: 600px; }
 input { padding: 4px 6px; font-size: 13px; }
 .grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 12px; }
 .pane { background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; min-width: 0; }
 .pane.gt { border: 2px solid #333; background: #fffbf0; }
 .viewer { width: 100%; height: 360px; position: relative; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 .lb { background: #ffd6d6; color: #800; }
 .lf { background: #d6f0d6; color: #060; }
 .imag { color: #c00; font-weight: 700; }
 .real { color: #888; }
 h2, h3 { margin: 4px 0 6px; }
 h3 { font-size: 14px; }
 .stats { color: #444; font-size: 11px; line-height: 1.5; }
 .align { font-size: 16px; font-weight: 700; padding: 2px 6px; border-radius: 4px; }
 .align.high  { background: #c8f0c8; color: #060; }
 .align.mid   { background: #fff0c8; color: #860; }
 .align.low   { background: #ffc8c8; color: #800; }
 .input-range { width: 160px; }
 button { padding: 4px 10px; }
</style></head><body>

<h2>BGCP — oracle view: GT + top-5 IGs ranked by gt_alignment</h2>
<div class="ctl">
  <label><b>Step:</b><select id="stepSel"></select></label>
  <input id="filter" placeholder="filter steps" oninput="rebuildOptions()">
  <button onclick="prevStep()">◀</button>
  <button onclick="nextStep()">▶</button>
  <label>Amplitude
    <input type="range" id="amp" min="0.05" max="1.5" step="0.05" value="0.5" class="input-range">
    <span id="ampVal">0.50</span>
  </label>
  <label>Speed
    <input type="range" id="speed" min="50" max="600" step="50" value="200" class="input-range">
    <span id="speedVal">200</span> ms
  </label>
  <button onclick="togglePlay()" id="playBtn">⏸ Pause</button>
  <span class="legend" style="margin-left:auto"><span class="lb">broken</span><span class="lf">formed</span></span>
</div>

<div class="grid" id="panes"></div>

<script>
const DATA = __DATA__;
const stepNames = Object.keys(DATA);
const sel = document.getElementById('stepSel');
const panesDiv = document.getElementById('panes');
let curStep = null;
let viewers = [];
let amp = 0.5;
let speed = 200;
let playing = true;
let timer = null;
let phase = 0;

function alignClass(a) {
  if (a >= 0.7) return 'high';
  if (a >= 0.3) return 'mid';
  return 'low';
}

function rebuildOptions() {
  const f = document.getElementById('filter').value.toLowerCase();
  sel.innerHTML = '';
  for (const n of stepNames) {
    if (!n.toLowerCase().includes(f)) continue;
    const d = DATA[n];
    const opt = document.createElement('option');
    opt.value = n;
    const aligns = d.panels.slice(1).map(p => p.gt_align.toFixed(2)).join('/');
    opt.textContent = `${n}  top5_align=${aligns}`;
    sel.appendChild(opt);
  }
  if (sel.options.length) render(sel.value);
}
sel.addEventListener('change', () => render(sel.value));
function prevStep() { if (sel.selectedIndex > 0) { sel.selectedIndex--; render(sel.value); } }
function nextStep() { if (sel.selectedIndex < sel.options.length - 1) { sel.selectedIndex++; render(sel.value); } }

document.getElementById('amp').addEventListener('input', e => {
  amp = +e.target.value; document.getElementById('ampVal').textContent = amp.toFixed(2);
});
document.getElementById('speed').addEventListener('input', e => {
  speed = +e.target.value; document.getElementById('speedVal').textContent = speed;
  restartTimer();
});
function togglePlay() {
  playing = !playing;
  document.getElementById('playBtn').textContent = playing ? '⏸ Pause' : '▶ Play';
  if (playing) restartTimer(); else { if (timer) clearInterval(timer); timer = null; }
}
function restartTimer() {
  if (timer) clearInterval(timer);
  if (!playing) return;
  timer = setInterval(() => { phase += 0.18; drawAll(Math.sin(phase) * amp); }, speed);
}

function buildPanes(d) {
  panesDiv.innerHTML = '';
  for (let i = 0; i < d.panels.length; i++) {
    const p = d.panels[i];
    const tag = i === 0 ? 'GT' : `IG #${i}`;
    const imagTag = p.is_imag ? '<span class="imag">imag</span>' : '<span class="real">real</span>';
    const ac = alignClass(p.gt_align);
    const div = document.createElement('div');
    div.className = 'pane' + (i === 0 ? ' gt' : '');
    div.innerHTML = `<h3>${tag}: ${p.label}
      <span class="align ${ac}">align=${p.gt_align.toFixed(3)}</span></h3>
      <div class="stats">
        mode #${p.mode_idx} · ${imagTag} freq=${p.freq.toFixed(1)} ·
        bond_ovlp=${p.bond_overlap.toFixed(2)} · rxn_ovlp=${p.rxn_overlap.toFixed(2)} ·
        core_frac=${p.core_fraction.toFixed(2)} · ${p.n_imag} imag/${p.n_modes_total}
      </div>
      <div id="v${i}" class="viewer"></div>`;
    panesDiv.appendChild(div);
  }
}

function render(name) {
  curStep = name;
  const d = DATA[name];
  buildPanes(d);
  viewers = [];
  for (let i = 0; i < d.panels.length; i++) {
    viewers.push($3Dmol.createViewer('v' + i, {backgroundColor: 'white'}));
    drawPane(i, 0);
    viewers[i].zoomTo();
  }
}

function drawPane(i, scale) {
  const d = DATA[curStep];
  const p = d.panels[i];
  const v = viewers[i];
  if (!v) return;
  const xyz = p.xyz_coords;
  const els = p.xyz_elements;
  const disp = p.disp;
  const n = xyz.length;
  let body = `${n}\nframe\n`;
  for (let k = 0; k < n; k++) {
    const x = xyz[k][0] + scale * disp[k][0];
    const y = xyz[k][1] + scale * disp[k][1];
    const z = xyz[k][2] + scale * disp[k][2];
    body += `${els[k]}  ${x.toFixed(6)}  ${y.toFixed(6)}  ${z.toFixed(6)}\n`;
  }
  v.removeAllModels();
  v.removeAllShapes();
  v.addModel(body, 'xyz');
  v.setStyle({}, {stick: {radius: 0.10}, sphere: {scale: 0.20}});
  const atoms = v.selectedAtoms({});
  for (const [a, b] of d.broken_bonds) {
    if (a < atoms.length && b < atoms.length) {
      v.addCylinder({start:{x:atoms[a].x,y:atoms[a].y,z:atoms[a].z},
                     end:{x:atoms[b].x,y:atoms[b].y,z:atoms[b].z},
                     color:'red', radius:0.10, dashed:true});
    }
  }
  for (const [a, b] of d.formed_bonds_R) {
    if (a < atoms.length && b < atoms.length) {
      v.addCylinder({start:{x:atoms[a].x,y:atoms[a].y,z:atoms[a].z},
                     end:{x:atoms[b].x,y:atoms[b].y,z:atoms[b].z},
                     color:'green', radius:0.10, dashed:true});
    }
  }
  for (const ai of d.core_atoms) {
    const a = atoms[ai]; if (!a) continue;
    const dd = p.disp[ai];
    const len = Math.hypot(dd[0], dd[1], dd[2]);
    if (len < 1e-3) continue;
    v.addArrow({start:{x:a.x, y:a.y, z:a.z},
                end:{x:a.x+dd[0]*1.5, y:a.y+dd[1]*1.5, z:a.z+dd[2]*1.5},
                color:'#0066cc', radius:0.08});
  }
  v.render();
}
function drawAll(scale) {
  const d = DATA[curStep];
  for (let i = 0; i < d.panels.length; i++) drawPane(i, scale);
}

window.addEventListener('load', () => {
  rebuildOptions();
  restartTimer();
});
</script>
</body></html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--top-k', type=int, default=5)
    ap.add_argument('--out', type=str, default=None,
                    help='Output filename (default: oracle_view_top<k>.html)')
    args = ap.parse_args()
    out_html = SRC_DIR / (args.out or f'oracle_view_top{args.top_k}.html')
    t0 = time.time()
    files = sorted(SRC_DIR.glob("*.html"))
    files = [f for f in files
             if not (f.name in ('index.html', 'flat_view.html', 'guess_quality.html')
                     or f.name.startswith('oracle_view'))]
    print(f"Reading {len(files)} per-step HTML payloads from {SRC_DIR}")

    out = {}
    n_skip = 0
    for i, hp in enumerate(files, 1):
        try:
            payload = load_step_payload(hp)
            flat = build_oracle_payload(payload, top_k=args.top_k)
            if flat is None:
                n_skip += 1
                continue
            out[flat['step']] = flat
        except Exception as e:
            n_skip += 1
            print(f"  skip {hp.name}: {e}", file=sys.stderr)
        if i % 20 == 0:
            print(f"  [{i}/{len(files)}]")

    ordered = {k: out[k] for k in sorted(out.keys())}
    html = HTML.replace('__DATA__', json.dumps(ordered))
    out_html.write_text(html)
    size_mb = out_html.stat().st_size / 1e6
    print(f"\nDone in {time.time()-t0:.1f}s. {len(ordered)} steps, {n_skip} skipped.")
    print(f"Output: {out_html}  ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
