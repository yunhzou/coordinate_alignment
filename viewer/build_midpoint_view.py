"""
Mid-point alignment-inspection viewer.

For each elementary step, render three structures side-by-side in a
single 3Dmol scene:

   R   |   midpoint = (R + Kabsch(P→R)) / 2   |   P

If the atom mapping (R[i] ↔ P[i]) is correct, the midpoint should
look like a plausible TS-ish geometry — bonds at intermediate
lengths, atoms in reasonable neighbourhoods. If the mapping is wrong,
the midpoint will be visually garbage: atoms inside other atoms,
disconnected fragments, super-long bonds.

Reads:
  appendix_perparation/Pure_Geometries_Elementary_Step/
    Benchmark_Guesses_Coordinate_Aligned_Version/<step>/{reactants,products}/
  appendix_perparation/viewer/mode_viewer/<step>.html  (broken/formed bonds)

Output:
  appendix_perparation/viewer/midpoint_view.html
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
import time
from pathlib import Path

import numpy as np


ALIGNED_DIR = (PROJECT_ROOT / 'appendix_perparation' / 'Pure_Geometries_Elementary_Step'
               / 'Benchmark_Guesses_Coordinate_Aligned_Version')
MODE_DIR    = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
OUT_HTML    = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'midpoint_view.html'


def parse_xyz(path):
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines[0])
    elements, coords = [], []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        elements.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return elements, np.asarray(coords, dtype=float)


def kabsch(P, Q):
    P = np.asarray(P, dtype=float); Q = np.asarray(Q, dtype=float)
    Pc = P.mean(0); Qc = Q.mean(0)
    H = (Q - Qc).T @ (P - Pc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Pc - R @ Qc
    return R, t


def parse_alignment_qc(comment_line):
    out = {}
    for k, default in (('pq_mapped', '0/0'), ('fallback', '0'), ('missing', '0')):
        m = re.search(rf'{k}=([\d/]+)', comment_line)
        out[k] = m.group(1) if m else default
    return out


def process_step(step):
    rdir = ALIGNED_DIR / step / 'reactants'
    pdir = ALIGNED_DIR / step / 'products'
    rfiles = list(rdir.glob('*.xyz')) if rdir.exists() else []
    pfiles = list(pdir.glob('*.xyz')) if pdir.exists() else []
    if not rfiles or not pfiles: return None
    elR, xyzR = parse_xyz(rfiles[0])
    elP, xyzP = parse_xyz(pfiles[0])
    if len(elR) != len(elP): return None
    p_comment = Path(pfiles[0]).read_text().splitlines()[1] if len(pfiles) else ''
    qc = parse_alignment_qc(p_comment)

    # Bonds + core from per-step HTML
    bro_pairs = []; fmd_pairs_R = []; core_R = []
    html = MODE_DIR / f"{step}.html"
    if html.exists():
        m = re.search(r"const DATA = (\{.*?\});\n", html.read_text(), re.DOTALL)
        if m:
            data = json.loads(m.group(1))
            bro_pairs   = [(int(i), int(j)) for (i, j) in data.get('broken_bonds', [])]
            fmd_pairs_R = [(int(i), int(j)) for (i, j) in data.get('formed_bonds_R', [])]
            core_R      = sorted({int(x) for x in data.get('core_atoms', [])})

    # Kabsch P→R using ALL atoms (the most charitable alignment so that
    # any midpoint anomaly really is the chemistry, not a bad fit).
    Rmat, t = kabsch(xyzR, xyzP)
    xyzP_kabsch = xyzP @ Rmat.T + t

    # Midpoint: simple linear interpolation. If mapping is correct,
    # bonds are at ~half-formed lengths and atoms are in plausible
    # neighbourhoods. If mapping is wrong, atoms collide / bonds explode.
    xyzM = 0.5 * (xyzR + xyzP_kabsch)

    # Side-by-side layout: R at x=0, M offset right, P offset further right
    span = float(np.ptp(xyzR, axis=0).max() if len(xyzR) else 5.0)
    offset = max(span * 1.4, 8.0)
    xyzM_view = xyzM         + np.array([offset,     0, 0])
    xyzP_view = xyzP_kabsch  + np.array([2 * offset, 0, 0])

    # Per-atom Δ at midpoint = ||P-R||/2 = half of the displacement.
    # Big midpoint Δ on a non-reactive atom = misalignment signature.
    deltas = np.linalg.norm(xyzP_kabsch - xyzR, axis=1)

    # Header positions (a label above each panel)
    header_R = [float(xyzR[:,0].min()) , float(xyzR[:,1].max() + 2.0), float(xyzR[:,2].max())]
    header_M = [float(xyzM_view[:,0].mean()), header_R[1], header_R[2]]
    header_P = [float(xyzP_view[:,0].mean()), header_R[1], header_R[2]]

    return {
        'step': step,
        'n_atoms': len(elR),
        'elements': list(elR),
        'coords_R': [[round(float(x), 4) for x in v] for v in xyzR],
        'coords_M': [[round(float(x), 4) for x in v] for v in xyzM_view],
        'coords_P': [[round(float(x), 4) for x in v] for v in xyzP_view],
        'deltas':   [round(float(x), 4) for x in deltas],
        'broken':   bro_pairs,
        'formed_R': fmd_pairs_R,
        'core_R':   core_R,
        'pq_mapped': qc['pq_mapped'],
        'fallback':  qc['fallback'],
        'missing':   qc['missing'],
        'offset_x':  offset,
        'header_R': header_R, 'header_M': header_M, 'header_P': header_P,
    }


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Mid-point alignment inspection — R | midpoint | P</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select { padding: 4px 6px; font-size: 13px; min-width: 540px; }
 .pane { background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
 .viewer { width: 100%; height: 720px; position: relative; }
 .stats { color: #444; font-size: 12px; margin-bottom: 6px; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 button { padding: 4px 10px; }
</style></head><body>

<h2>Mid-point alignment inspection — R | (R + P)/2 | P</h2>
<div class="ctl">
  <label><b>Step:</b><select id="stepSel"></select></label>
  <input id="filter" placeholder="filter steps" oninput="rebuildOptions()">
  <button onclick="prev()">◀</button>
  <button onclick="next()">▶</button>
  <label><input type="checkbox" id="showBonds" checked> broken/formed bond cylinders</label>
  <span class="legend" style="margin-left:auto">
    <span style="background:#888; color:white; padding:2px 8px; border-radius:4px;">R / midpoint / P</span>
    <span style="background:#cc3333; color:white; padding:2px 8px; border-radius:4px;">broken</span>
    <span style="background:#2a8a2a; color:white; padding:2px 8px; border-radius:4px;">formed</span>
  </span>
</div>

<div class="pane">
  <div class="stats" id="info"></div>
  <div id="v" class="viewer"></div>
</div>

<script>
const DATA = __DATA__;
const stepNames = Object.keys(DATA);
const sel = document.getElementById('stepSel');
let viewer = null;
let curStep = null;

function rebuildOptions() {
  const f = document.getElementById('filter').value.toLowerCase();
  sel.innerHTML = '';
  for (const n of stepNames) {
    if (!n.toLowerCase().includes(f)) continue;
    const d = DATA[n];
    const opt = document.createElement('option');
    opt.value = n;
    opt.textContent = `${n}   N=${d.n_atoms}   pq=${d.pq_mapped}   broken=${d.broken.length}   formed=${d.formed_R.length}`;
    sel.appendChild(opt);
  }
  if (sel.options.length) render(sel.value);
}
sel.addEventListener('change', () => render(sel.value));
function prev() { if (sel.selectedIndex > 0) { sel.selectedIndex--; render(sel.value); } }
function next() { if (sel.selectedIndex < sel.options.length - 1) { sel.selectedIndex++; render(sel.value); } }
document.getElementById('showBonds').addEventListener('change', () => render(curStep));

function buildXyz(elements, coords) {
  let body = `${elements.length}\nframe\n`;
  for (let k = 0; k < elements.length; k++) {
    const c = coords[k];
    body += `${elements[k]}  ${c[0].toFixed(4)}  ${c[1].toFixed(4)}  ${c[2].toFixed(4)}\n`;
  }
  return body;
}

function render(name) {
  curStep = name;
  const d = DATA[name];
  document.getElementById('info').innerHTML =
    `<b>${name}</b> · ${d.n_atoms} atoms · `
    + `PQ-mapped ${d.pq_mapped} · fallback ${d.fallback} · missing ${d.missing} · `
    + `broken = ${d.broken.length} · formed = ${d.formed_R.length} · `
    + `<i>If the midpoint looks nonsensical (atoms inside each other, broken connectivity), the atom mapping is wrong.</i>`;

  if (viewer) { viewer.removeAllModels(); viewer.removeAllShapes(); viewer.removeAllLabels(); }
  viewer = $3Dmol.createViewer('v', {backgroundColor: 'white'});

  // Three structures, all in one viewer — translated along x.
  // Models are added in the order: R (0), midpoint (1), P (2).
  viewer.addModel(buildXyz(d.elements, d.coords_R), 'xyz');
  viewer.addModel(buildXyz(d.elements, d.coords_M), 'xyz');
  viewer.addModel(buildXyz(d.elements, d.coords_P), 'xyz');
  // Distinct styling so the three panels read at a glance:
  //  R: full CPK
  //  M: a touch lighter sticks (this is the diagnostic structure)
  //  P: full CPK
  viewer.setStyle({model: 0}, {stick: {radius: 0.13}, sphere: {scale: 0.22}});
  viewer.setStyle({model: 1}, {stick: {radius: 0.10, opacity: 0.95}, sphere: {scale: 0.20, opacity: 0.95}});
  viewer.setStyle({model: 2}, {stick: {radius: 0.13}, sphere: {scale: 0.22}});

  // Header labels above each panel
  viewer.addLabel('R', {position:{x:d.header_R[0], y:d.header_R[1], z:d.header_R[2]},
      fontSize:18, fontColor:'white', backgroundColor:'#3a6dbf',
      backgroundOpacity:0.95, showBackground:true, inFront:true});
  viewer.addLabel('midpoint  =  (R + P)/2',
      {position:{x:d.header_M[0], y:d.header_M[1], z:d.header_M[2]},
      fontSize:18, fontColor:'white', backgroundColor:'#aa6600',
      backgroundOpacity:0.95, showBackground:true, inFront:true});
  viewer.addLabel('P', {position:{x:d.header_P[0], y:d.header_P[1], z:d.header_P[2]},
      fontSize:18, fontColor:'white', backgroundColor:'#cc3366',
      backgroundOpacity:0.95, showBackground:true, inFront:true});

  // Optional broken/formed bond cylinders.
  // - broken: dashed red inside R (bond was here, now gone)
  // - formed: dashed green inside P (bond is here, was absent)
  // - both: dashed grey at the midpoint (visible at half-bond)
  if (document.getElementById('showBonds').checked) {
    for (const [a, b] of d.broken) {
      const ra = d.coords_R[a], rb = d.coords_R[b];
      viewer.addCylinder({start:{x:ra[0],y:ra[1],z:ra[2]},
                           end:{x:rb[0],y:rb[1],z:rb[2]},
                           color:'red', radius: 0.10, dashed: true});
      const ma = d.coords_M[a], mb = d.coords_M[b];
      viewer.addCylinder({start:{x:ma[0],y:ma[1],z:ma[2]},
                           end:{x:mb[0],y:mb[1],z:mb[2]},
                           color:'red', radius: 0.06, dashed: true, opacity: 0.6});
    }
    for (const [a, b] of d.formed_R) {
      const pa = d.coords_P[a], pb = d.coords_P[b];
      viewer.addCylinder({start:{x:pa[0],y:pa[1],z:pa[2]},
                           end:{x:pb[0],y:pb[1],z:pb[2]},
                           color:'green', radius: 0.10, dashed: true});
      const ma = d.coords_M[a], mb = d.coords_M[b];
      viewer.addCylinder({start:{x:ma[0],y:ma[1],z:ma[2]},
                           end:{x:mb[0],y:mb[1],z:mb[2]},
                           color:'green', radius: 0.06, dashed: true, opacity: 0.6});
    }
  }

  // Hover support: hovering an atom in any of the 3 panels highlights
  // the same atom index in the other two panels.
  let hoverShapes = [];
  function clearHover() { for (const s of hoverShapes) viewer.removeShape(s); hoverShapes = []; }
  viewer.setHoverable({}, true,
    function(atom) {
      const k = (atom.serial !== undefined) ? (atom.serial - 1) : atom.index;
      if (k === null || k < 0 || k >= d.coords_R.length) return;
      clearHover();
      for (const c of [d.coords_R[k], d.coords_M[k], d.coords_P[k]]) {
        hoverShapes.push(viewer.addSphere(
          {center:{x:c[0],y:c[1],z:c[2]}, radius: 0.5,
           color: '#ffaa00', opacity: 0.65}));
      }
      viewer.addLabel(`#${k} ${d.elements[k]}`,
        {position:{x:d.coords_M[k][0], y:d.coords_M[k][1], z:d.coords_M[k][2] + 0.7},
         fontSize:14, fontColor:'white', backgroundColor:'#553300',
         backgroundOpacity:0.9, showBackground:true, inFront:true});
      viewer.render();
    },
    function(_) { clearHover(); viewer.render(); }
  );

  viewer.zoomTo();
  viewer.render();
}

window.addEventListener('load', rebuildOptions);
</script>
</body></html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--steps', nargs='+', default=None)
    args = ap.parse_args()

    steps = sorted(d.name for d in ALIGNED_DIR.iterdir() if d.is_dir())
    if args.steps: steps = [s for s in steps if s in set(args.steps)]
    if args.limit: steps = steps[:args.limit]
    print(f"Processing {len(steps)} steps from {ALIGNED_DIR}")

    out = {}; t0 = time.time(); n_skip = 0
    for i, s in enumerate(steps, 1):
        try:
            r = process_step(s)
            if r is not None: out[s] = r
            else:             n_skip += 1
        except Exception as e:
            n_skip += 1
            print(f"  skip {s}: {e}")
        if i % 30 == 0: print(f"  [{i}/{len(steps)}]  ({time.time()-t0:.0f}s)")

    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    html = HTML.replace('__DATA__', json.dumps(out))
    OUT_HTML.write_text(html)
    sz = OUT_HTML.stat().st_size / 1e6
    print(f"\nDone in {time.time()-t0:.0f}s. {len(out)} steps, {n_skip} skipped.")
    print(f"Output: {OUT_HTML}  ({sz:.1f} MB)")


if __name__ == '__main__':
    main()
