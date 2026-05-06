"""
Single-page flat viewer: GT + top-2 initial-guess TS for every BGCP step,
3 animated 3Dmol viewers side by side per step, all playing each TS's
highest-core-fraction imaginary mode in sync.

Reuses data from the per-step HTMLs at out/mode_viewer/<step>.html
(no xtb, no alignment recomputation — those HTMLs already contain
the full payload as embedded JSON; we just pick the right subset and
re-bundle into one combined HTML).

Output: out/mode_viewer/flat_view.html
"""
from __future__ import annotations
import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np


def kabsch(P, Q):
    """Optimal rotation+translation aligning Q to P.
    Returns (R, t) such that (R @ Q.T).T + t ≈ P (least-squares)."""
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


SRC_DIR = Path(__file__).parent / "out" / "mode_viewer"
OUT_HTML = SRC_DIR / "flat_view.html"


def load_step_payload(html_path):
    """Parse the embedded JSON payload from a per-step HTML."""
    text = html_path.read_text()
    m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
    if not m:
        raise ValueError(f"could not find DATA= in {html_path}")
    return json.loads(m.group(1))


def best_imag_mode(ts):
    """Best mode = highest bond_overlap among imag modes (else among all).
    bond_overlap projects the mode onto the bond-stretching/-compressing
    direction at TS coords; this rewards true reaction modes."""
    imag = [m for m in ts['modes'] if m['freq'] < 0]
    key = lambda m: m.get('bond_overlap', m.get('rxn_overlap',
                                                m.get('core_fraction', 0)))
    if imag:
        return max(imag, key=key)
    return max(ts['modes'], key=key)


def build_flat_payload(step_payload):
    """Pick GT + top-2 IG (by their best imag-mode core_fraction) and
    minimize each panel to just the one mode being shown."""
    ts_list = step_payload['ts_list']
    by_label = {ts['label']: ts for ts in ts_list}
    gt = by_label.get('groundtruth')
    if gt is None:
        return None
    igs = [(ts, best_imag_mode(ts)) for ts in ts_list if ts['label'] != 'groundtruth']
    if len(igs) < 2:
        return None
    rank_key = lambda t: -t[1].get('bond_overlap',
                                    t[1].get('rxn_overlap',
                                             t[1].get('core_fraction', 0)))
    igs.sort(key=rank_key)
    chosen = [(gt, best_imag_mode(gt)), igs[0], igs[1]]

    def to_panel(ts, m):
        return {
            'label': ts['label'],
            'mode_idx': m['idx'],
            'freq': m['freq'],
            'bond_overlap': m.get('bond_overlap', 0.0),
            'rxn_overlap': m.get('rxn_overlap', 0.0),
            'core_fraction': m['core_fraction'],
            'is_imag': m['freq'] < 0,
            'n_imag': ts['n_imag'],
            'n_modes_total': ts['n_modes_total'],
            'xyz_elements': ts['xyz_elements'],
            'xyz_coords': ts['xyz_coords'],
            'disp': m['disp'],
        }

    panels = [to_panel(ts, m) for (ts, m) in chosen]

    # Kabsch-align panels 1 and 2 to panel 0 (GT). Operates on the R-frame
    # aligned coordinates so atom indices already correspond. The same
    # rotation is applied to displacement vectors (rotation only — no
    # translation, since displacements are differences not positions).
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
<title>Flat view — GT + top-2 initial-guess TS vibration modes</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select { padding: 4px 6px; font-size: 13px; min-width: 480px; }
 input { padding: 4px 6px; font-size: 13px; }
 .row { display: flex; gap: 12px; }
 .pane { flex: 1; background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; min-width: 0; }
 .viewer { width: 100%; height: 540px; position: relative; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 .lb { background: #ffd6d6; color: #800; }
 .lf { background: #d6f0d6; color: #060; }
 .imag { color: #c00; font-weight: 700; }
 h2, h3 { margin: 4px 0 6px; }
 h3 { font-size: 14px; }
 .stats { color: #444; font-size: 12px; }
 .input-range { width: 160px; }
 button { padding: 4px 10px; }
</style></head><body>

<h2>BGCP — flat view: groundtruth + top-2 initial-guess TS vibration modes</h2>
<div class="ctl">
  <label><b>Step:</b><select id="stepSel"></select></label>
  <input id="filter" placeholder="filter" oninput="rebuildOptions()">
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

<div class="row" id="panes"></div>

<script>
const DATA = __DATA__;
const stepNames = Object.keys(DATA);
const sel = document.getElementById('stepSel');
const panesDiv = document.getElementById('panes');
let curStep = null;
let viewers = [null, null, null];
let amp = 0.5;
let speed = 200;
let playing = true;
let timer = null;
let phase = 0;

function rebuildOptions() {
  const f = document.getElementById('filter').value.toLowerCase();
  sel.innerHTML = '';
  for (const n of stepNames) {
    if (!n.toLowerCase().includes(f)) continue;
    const d = DATA[n];
    const opt = document.createElement('option');
    opt.value = n;
    const gt = d.panels[0], a = d.panels[1], b = d.panels[2];
    opt.textContent = `${n}  GT=${gt.bond_overlap.toFixed(2)}/${gt.label}  A=${a.bond_overlap.toFixed(2)}/${a.label}  B=${b.bond_overlap.toFixed(2)}/${b.label}`;
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
  for (let i = 0; i < 3; i++) {
    const p = d.panels[i];
    const tag = i === 0 ? 'GT' : (i === 1 ? 'IG #1' : 'IG #2');
    const imagTag = p.is_imag ? '<span class="imag">imag</span>' : 'real';
    const div = document.createElement('div');
    div.className = 'pane';
    div.innerHTML = `<h3>${tag}: ${p.label}</h3>
      <div class="stats">mode #${p.mode_idx} · ${imagTag} freq=${p.freq.toFixed(2)} cm⁻¹ · bond_ovlp=${p.bond_overlap.toFixed(3)} · rxn_ovlp=${p.rxn_overlap.toFixed(3)} · core_frac=${p.core_fraction.toFixed(3)} · ${p.n_imag} imag modes</div>
      <div id="v${i}" class="viewer"></div>`;
    panesDiv.appendChild(div);
  }
}

function render(name) {
  curStep = name;
  const d = DATA[name];
  buildPanes(d);
  for (let i = 0; i < 3; i++) {
    viewers[i] = $3Dmol.createViewer('v' + i, {backgroundColor: 'white'});
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
function drawAll(scale) { for (let i = 0; i < 3; i++) drawPane(i, scale); }

window.addEventListener('load', () => {
  rebuildOptions();
  restartTimer();
});
</script>
</body></html>
"""


def main():
    ap = argparse.ArgumentParser()
    args = ap.parse_args()
    t0 = time.time()
    files = sorted(SRC_DIR.glob("*.html"))
    files = [f for f in files if f.name not in ('index.html', 'flat_view.html')]
    print(f"Reading {len(files)} per-step HTML payloads from {SRC_DIR}")

    out = {}
    n_skip = 0
    for i, hp in enumerate(files, 1):
        try:
            payload = load_step_payload(hp)
            flat = build_flat_payload(payload)
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
    OUT_HTML.write_text(html)
    size_mb = OUT_HTML.stat().st_size / 1e6
    print(f"\nDone in {time.time()-t0:.1f}s. {len(ordered)} steps, {n_skip} skipped.")
    print(f"Output: {OUT_HTML}  ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
