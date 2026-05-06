"""
Combined HTML viewer: vibration modes for every BGCP TS, ranked by
core-atom contribution.

For each step, generates `out/mode_viewer/<step>.html` with:
  - 3D viewer showing the TS structure (R-frame atom indices)
  - Broken bonds (red dashed) and formed bonds (green dashed)
  - Animated vibration mode (Cartesian displacement, sin(t)·δ)
  - Dropdown to switch TS (groundtruth, iter1, …, iter20)
  - Mode list ranked by core_fraction; click to switch
  - Default: highest-core-fraction imaginary mode of groundtruth TS

Top-level `out/mode_viewer/index.html` lists every step.

Multiprocessing parallelizes per-step HTML generation.
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path
_HERE = _Path(__file__).resolve().parent
_sys.path.insert(0, str(_HERE.parent / "src"))
_sys.path.insert(0, str(_HERE))
PROJECT_ROOT = _HERE.parent  # _RXN_CORE_PATH_SETUP

import argparse
import csv
import json
import multiprocessing as mp
import os
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from rxn_core_pq import align_from_arrays
from rxn_core_frag import write_xyz_str
from bgcp_io import BGCP_ROOT, LOOKUP, list_step_dirs
from align_bgcp_coords import load_cached_xtb, fill_unmapped_greedy, reindex_to_R_frame
from analyze_core_modes import (
    parse_g98_modes, core_atoms_in_R_frame, reindex_modes_to_R, list_ts_caches,
    reaction_coord_delta, rxn_overlap_per_mode,
    bond_reaction_vector, bond_overlap_per_mode,
    WORK_MODES,
)


OUT_DIR = PROJECT_ROOT / "out" / "mode_viewer"
OUT_DIR.mkdir(parents=True, exist_ok=True)


HTML_TEMPLATE = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Vibration modes — __STEP__</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select { padding: 4px 6px; font-size: 13px; min-width: 200px; }
 .row { display: flex; gap: 12px; }
 .pane { background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
 .v { flex: 2; }
 .l { flex: 1; min-width: 320px; max-width: 420px; }
 #viewer { width: 100%; height: 600px; position: relative; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 .lb { background: #ffd6d6; color: #800; }
 .lf { background: #d6f0d6; color: #060; }
 table { border-collapse: collapse; font-size: 13px; width: 100%; }
 td, th { border: 1px solid #ddd; padding: 4px 8px; text-align: left; }
 tr.imag { color: #c00; font-weight: 600; }
 tr.selected { background: #ffd700; }
 tr.clickable:hover { background: #fffacc; cursor: pointer; }
 .stats { color: #444; font-size: 13px; }
 h2, h3 { margin: 4px 0 6px; }
 .core-atoms { font-family: ui-monospace, monospace; font-size: 12px; color: #666; }
 .input-range { width: 200px; }
</style></head><body>

<h2>__STEP__</h2>
<div class="ctl">
  <label>TS:
    <select id="tsSel"></select>
  </label>
  <label>Amplitude
    <input type="range" id="amp" min="0.05" max="1.5" step="0.05" value="0.5" class="input-range">
    <span id="ampVal">0.5</span>
  </label>
  <label>Speed
    <input type="range" id="speed" min="50" max="600" step="50" value="200" class="input-range">
    <span id="speedVal">200</span> ms
  </label>
  <button onclick="togglePlay()" id="playBtn">⏸ Pause</button>
  <span class="stats" id="info"></span>
  <span class="legend" style="margin-left:auto"><span class="lb">broken</span><span class="lf">formed</span></span>
</div>

<div class="row">
  <div class="pane v">
    <h3>Animated mode</h3>
    <div id="viewer"></div>
  </div>
  <div class="pane l">
    <h3>Modes (click to switch)</h3>
    <div class="core-atoms" id="coreAtomsInfo"></div>
    <div style="overflow-y:auto;max-height:540px;margin-top:6px">
      <table>
        <thead><tr><th>rank</th><th>idx</th><th>freq cm⁻¹</th><th>bond ovlp</th><th>rxn ovlp</th><th>core frac</th></tr></thead>
        <tbody id="modeTab"></tbody>
      </table>
    </div>
  </div>
</div>

<script>
const DATA = __DATA__;
let curTSIdx = 0;
let curModeIdx = 0;
let amp = 0.5;
let speed = 200;
let playing = true;
let timer = null;
let phase = 0;
let viewer = null;

const tsSel = document.getElementById('tsSel');
DATA.ts_list.forEach((t, i) => {
  const o = document.createElement('option');
  o.value = i; o.textContent = `${t.label}  (${t.n_imag} imag, ${t.n_modes_total} modes)`;
  tsSel.appendChild(o);
});
tsSel.value = DATA.default_ts_idx;
tsSel.addEventListener('change', () => { curTSIdx = +tsSel.value; rebuildModeList(); selectMode(DATA.ts_list[curTSIdx].default_mode_idx); });

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
  timer = setInterval(() => {
    phase += 0.18;
    drawFrame(Math.sin(phase) * amp);
  }, speed);
}

function selectMode(idx) {
  curModeIdx = idx;
  // Highlight in mode list
  document.querySelectorAll('#modeTab tr').forEach(tr => tr.classList.remove('selected'));
  const row = document.getElementById('mode-row-' + idx);
  if (row) row.classList.add('selected');
  updateInfo();
  drawFrame(0);
}

function updateInfo() {
  const ts = DATA.ts_list[curTSIdx];
  const m = ts.modes[curModeIdx];
  const f = m.freq.toFixed(2);
  const tag = m.freq < 0 ? '<span style="color:#c00;font-weight:600">imag</span>' : 'real';
  document.getElementById('info').innerHTML =
    `mode #${m.idx} · ${tag} freq=${f} cm⁻¹ · bond_ovlp=${m.bond_overlap.toFixed(3)} · rxn_ovlp=${m.rxn_overlap.toFixed(3)} · core_frac=${m.core_fraction.toFixed(3)}`;
  document.getElementById('coreAtomsInfo').innerHTML =
    `<b>Core atoms (R-frame, ${DATA.core_atoms.length}):</b> ${DATA.core_atoms.join(', ')}`;
}

function rebuildModeList() {
  const ts = DATA.ts_list[curTSIdx];
  const tbody = document.getElementById('modeTab');
  tbody.innerHTML = '';
  // Show imag modes first (already sorted by rank), then real top-15
  const imags = ts.modes.filter((m, i) => m.freq < 0);
  const reals = ts.modes.filter((m, i) => m.freq >= 0).slice(0, 30);
  const showOrder = [...imags.map((_, i) => ts.modes.indexOf(imags[i])),
                     ...reals.map((_, i) => ts.modes.indexOf(reals[i]))];
  showOrder.forEach((mIdx, rank) => {
    const m = ts.modes[mIdx];
    const tr = document.createElement('tr');
    tr.id = 'mode-row-' + mIdx;
    tr.classList.add('clickable');
    if (m.freq < 0) tr.classList.add('imag');
    tr.innerHTML = `<td>${rank}</td><td>${m.idx}</td><td>${m.freq.toFixed(2)}</td><td>${m.bond_overlap.toFixed(3)}</td><td>${m.rxn_overlap.toFixed(3)}</td><td>${m.core_fraction.toFixed(3)}</td>`;
    tr.onclick = () => selectMode(mIdx);
    tbody.appendChild(tr);
  });
}

function drawFrame(scale) {
  if (!viewer) return;
  const ts = DATA.ts_list[curTSIdx];
  const m = ts.modes[curModeIdx];
  const xyz = ts.xyz_coords;  // [n][3]
  const elements = ts.xyz_elements;
  const disp = m.disp;       // [n][3]
  const n = xyz.length;
  // Build xyz body string
  let body = `${n}\nframe\n`;
  for (let i = 0; i < n; i++) {
    const x = xyz[i][0] + scale * disp[i][0];
    const y = xyz[i][1] + scale * disp[i][1];
    const z = xyz[i][2] + scale * disp[i][2];
    body += `${elements[i]}  ${x.toFixed(6)}  ${y.toFixed(6)}  ${z.toFixed(6)}\n`;
  }
  viewer.removeAllModels();
  viewer.removeAllLabels();
  viewer.removeAllShapes();
  viewer.addModel(body, 'xyz');
  viewer.setStyle({}, {stick: {radius: 0.10}, sphere: {scale: 0.20}});
  // bond cylinders
  const atoms = viewer.selectedAtoms({});
  for (const [i, j] of DATA.broken_bonds) {
    if (i < atoms.length && j < atoms.length) {
      const a = atoms[i], b = atoms[j];
      viewer.addCylinder({start:{x:a.x,y:a.y,z:a.z}, end:{x:b.x,y:b.y,z:b.z},
                          color:'red', radius:0.10, dashed:true});
    }
  }
  for (const [i, j] of DATA.formed_bonds_R) {
    if (i < atoms.length && j < atoms.length) {
      const a = atoms[i], b = atoms[j];
      viewer.addCylinder({start:{x:a.x,y:a.y,z:a.z}, end:{x:b.x,y:b.y,z:b.z},
                          color:'green', radius:0.10, dashed:true});
    }
  }
  // small displacement arrows for current mode (only on core atoms)
  for (const i of DATA.core_atoms) {
    const a = atoms[i]; if (!a) continue;
    const d = m.disp[i];
    const len = Math.hypot(d[0], d[1], d[2]);
    if (len < 1e-3) continue;
    viewer.addArrow({
      start: {x:a.x, y:a.y, z:a.z},
      end:   {x:a.x + d[0]*1.5, y:a.y + d[1]*1.5, z:a.z + d[2]*1.5},
      color: '#0066cc', radius: 0.08,
    });
  }
  viewer.render();
}

window.addEventListener('load', () => {
  viewer = $3Dmol.createViewer('viewer', {backgroundColor: 'white'});
  curTSIdx = DATA.default_ts_idx;
  curModeIdx = DATA.ts_list[curTSIdx].default_mode_idx;
  rebuildModeList();
  drawFrame(0);
  viewer.zoomTo();
  selectMode(curModeIdx);
  restartTimer();
});
</script>
</body></html>
"""


def build_step_data(name):
    step_modes_dir = WORK_MODES / name
    elR, xyzR, wboR, _ = load_cached_xtb(step_modes_dir / "R")
    elP, xyzP, wboP, _ = load_cached_xtb(step_modes_dir / "P")
    rp_res = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP)
    mapping_RP = dict(rp_res['mapping'])
    inv_RP = {v: k for k, v in mapping_RP.items()}
    core_R = core_atoms_in_R_frame(mapping_RP, rp_res['broken'], rp_res['formed'])
    full_RP = fill_unmapped_greedy(elR, xyzR, elP, xyzP, mapping_RP)
    delta_RP = reaction_coord_delta(xyzR, xyzP, full_RP)
    n_R = len(elR)
    broken_R = [(int(a), int(b)) for (a, b, _, _) in rp_res['broken']]
    formed_R = [(int(inv_RP[a]), int(inv_RP[b]))
                for (a, b, _, _) in rp_res['formed']
                if a in inv_RP and b in inv_RP]
    broken_bonds = [[int(a), int(b)] for (a, b, _, _) in rp_res['broken']]
    formed_bonds_R = []
    for (a, b, _, _) in rp_res['formed']:
        ra = inv_RP.get(a); rb = inv_RP.get(b)
        if ra is not None and rb is not None:
            formed_bonds_R.append([int(ra), int(rb)])

    ts_list = []
    default_ts_idx = 0
    for ti, (label, hess_dir, sp_dir) in enumerate(list_ts_caches(step_modes_dir)):
        g98 = hess_dir / "g98.out"
        if not g98.exists():
            continue
        try:
            freqs, modes_TS = parse_g98_modes(g98)
        except Exception:
            continue
        if modes_TS.shape[0] == 0:
            continue
        elT, xyzT, wboT, _ = load_cached_xtb(sp_dir)
        ts_res = align_from_arrays(elR, xyzR, wboR, elT, xyzT, wboT)
        mapping_RT = fill_unmapped_greedy(elR, xyzR, elT, xyzT, dict(ts_res['mapping']))
        # Reindex TS coords to R-frame for display
        aligned_el, aligned_xyz, _ = reindex_to_R_frame(
            elR, xyzR, elT, xyzT, mapping_RT)
        modes_R = reindex_modes_to_R(modes_TS, mapping_RT, n_R)
        sq = (modes_R ** 2).sum(axis=2)
        total = sq.sum(axis=1)
        core_e = sq[:, core_R].sum(axis=1) if core_R else np.zeros(modes_R.shape[0])
        fraction = np.where(total > 1e-12, core_e / total, 0.0)
        rxn_ov = rxn_overlap_per_mode(modes_R, delta_RP, core_R)
        # bond_overlap: project onto bond-stretching / -compressing direction
        ts_xyz_in_R = np.zeros_like(np.asarray(xyzR))
        for r, t in mapping_RT.items():
            ts_xyz_in_R[r] = xyzT[t]
        V = bond_reaction_vector(ts_xyz_in_R, broken_R, formed_R)
        bond_ov = bond_overlap_per_mode(modes_R, V)
        imag_mask = freqs < 0
        n_imag = int(imag_mask.sum())

        # Pick default mode: imag mode with highest bond_overlap; else mode 0.
        if n_imag > 0:
            imag_idx = np.where(imag_mask)[0]
            default_mode_idx = int(imag_idx[np.argmax(bond_ov[imag_idx])])
        else:
            default_mode_idx = int(np.argmax(bond_ov))

        # Keep all imag, plus top-30 real by bond_overlap, for the table.
        keep_idx = list(np.where(imag_mask)[0])
        real_idx = [i for i in range(len(freqs)) if not imag_mask[i]]
        real_idx.sort(key=lambda i: -bond_ov[i])
        keep_idx += real_idx[:30]
        keep_idx = sorted(set(keep_idx))
        modes_payload = []
        for i in keep_idx:
            modes_payload.append({
                'idx': int(i),
                'freq': round(float(freqs[i]), 4),
                'bond_overlap': round(float(bond_ov[i]), 4),
                'rxn_overlap': round(float(rxn_ov[i]), 4),
                'core_fraction': round(float(fraction[i]), 4),
                'disp': [[round(float(x), 4) for x in v] for v in modes_R[i]],
            })

        # Update default_mode_idx to be index within the kept list
        try:
            default_mode_position = keep_idx.index(default_mode_idx)
        except ValueError:
            default_mode_position = 0

        ts_list.append({
            'label': label,
            'n_imag': n_imag,
            'n_modes_total': int(modes_R.shape[0]),
            'xyz_elements': aligned_el,
            'xyz_coords': [[round(float(x), 4) for x in v] for v in aligned_xyz],
            'modes': modes_payload,
            'default_mode_idx': default_mode_position,
        })
        if label == 'groundtruth':
            default_ts_idx = ti

    payload = {
        'step': name,
        'n_atoms': n_R,
        'core_atoms': core_R,
        'broken_bonds': broken_bonds,
        'formed_bonds_R': formed_bonds_R,
        'ts_list': ts_list,
        'default_ts_idx': default_ts_idx,
    }
    return payload


def render_one(name):
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    out_path = OUT_DIR / f"{sanitized}.html"
    payload = build_step_data(name)
    if not payload['ts_list']:
        return None
    html = (HTML_TEMPLATE
            .replace('__STEP__', name)
            .replace('__DATA__', json.dumps(payload)))
    out_path.write_text(html)
    # Top imag bond_overlap across all TS, for index summary.
    best_imag_ov = 0.0
    best_imag_label = '—'
    n_ts_with_imag = 0
    for ts in payload['ts_list']:
        for m in ts['modes']:
            if m['freq'] < 0 and m['bond_overlap'] > best_imag_ov:
                best_imag_ov = m['bond_overlap']
                best_imag_label = ts['label']
        if ts['n_imag'] > 0:
            n_ts_with_imag += 1
    return {
        'name': name,
        'file': out_path.name,
        'n_atoms': payload['n_atoms'],
        'n_core': len(payload['core_atoms']),
        'n_ts': len(payload['ts_list']),
        'n_ts_with_imag': n_ts_with_imag,
        'best_imag_ov': best_imag_ov,
        'best_imag_label': best_imag_label,
    }


def _safe(name):
    try:
        return (name, True, render_one(name))
    except Exception as e:
        return (name, False, f"{type(e).__name__}: {e}\n{traceback.format_exc()}")


INDEX_HTML = r"""<!doctype html>
<html><head><meta charset="utf-8"><title>BGCP vibration mode viewer</title>
<style>
 body { font-family: -apple-system, sans-serif; margin: 20px; max-width: 1100px; }
 input { padding: 4px 6px; font-size: 13px; width: 320px; }
 table { border-collapse: collapse; font-size: 13px; width: 100%; margin-top: 12px; }
 td, th { border: 1px solid #ccc; padding: 4px 8px; text-align: left; }
 a { color: #06c; text-decoration: none; }
 a:hover { text-decoration: underline; }
</style></head><body>
<h2>BGCP — vibration modes ranked by core-atom contribution</h2>
<p>For each step the per-step page shows the TS structure with broken (red dashed)
and formed (green dashed) bonds plus the animated vibration of the selected
mode. The default selection is the imaginary mode whose displacement is most
concentrated on the reaction-core atoms (the atoms touching broken/formed bonds).</p>

<input id="filter" placeholder="filter (substring)" oninput="filt()">

<table id="t"><thead><tr>
<th>step</th><th>N</th><th>core</th><th>n_TS</th><th>n_TS_with_imag</th>
<th>best imag bond_overlap</th><th>best label</th>
</tr></thead><tbody>
__ROWS__
</tbody></table>

<script>
function filt() {
  const f = document.getElementById('filter').value.toLowerCase();
  document.querySelectorAll('#t tbody tr').forEach(tr => {
    tr.style.display = tr.dataset.name.toLowerCase().includes(f) ? '' : 'none';
  });
}
</script>
</body></html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    ap.add_argument('--steps', nargs='+', default=None)
    ap.add_argument('--limit', type=int, default=None)
    args = ap.parse_args()

    all_steps = [d.name for d in list_step_dirs()]
    if args.steps:
        steps = [s for s in all_steps if s in set(args.steps)]
    elif args.limit:
        steps = all_steps[:args.limit]
    else:
        steps = all_steps
    steps = [s for s in steps if (WORK_MODES / s / "R" / "wbo").exists()]

    print(f"Building viewer HTML for {len(steps)} steps using {args.workers} workers")
    print(f"  output: {OUT_DIR}")
    t0 = time.time()
    results = []
    with mp.Pool(args.workers) as pool:
        for i, (name, ok, payload) in enumerate(pool.imap_unordered(_safe, steps), 1):
            if ok and payload:
                results.append(payload)
                print(f"[{i:3d}/{len(steps)}] {name[:55]:55s}  "
                      f"core={payload['n_core']}  ts={payload['n_ts']}  "
                      f"best_imag={payload['best_imag_ov']:.3f}/{payload['best_imag_label']}")
            else:
                print(f"[{i:3d}/{len(steps)}] {name[:55]:55s}  ERROR/skipped: {payload}")
            sys.stdout.flush()

    # Index page
    rows = ""
    for r in sorted(results, key=lambda x: -x['best_imag_ov']):
        rows += (f"<tr data-name='{r['name']}'>"
                 f"<td><a href='{r['file']}' target='_blank'>{r['name']}</a></td>"
                 f"<td>{r['n_atoms']}</td><td>{r['n_core']}</td>"
                 f"<td>{r['n_ts']}</td><td>{r['n_ts_with_imag']}</td>"
                 f"<td>{r['best_imag_ov']:.3f}</td><td>{r['best_imag_label']}</td>"
                 f"</tr>")
    (OUT_DIR / "index.html").write_text(INDEX_HTML.replace('__ROWS__', rows))
    print(f"\nDone in {time.time()-t0:.1f}s")
    print(f"Index: {OUT_DIR / 'index.html'}")


if __name__ == "__main__":
    main()
