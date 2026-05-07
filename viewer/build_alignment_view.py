"""
Atom-alignment inspection viewer.

For each elementary step, render reactant and product side-by-side in
a single 3Dmol scene, with thin grey cylinders connecting each atom
i in R to atom i in P. Atoms touching broken or formed bonds are
highlighted (broken-bond endpoints in red, formed-bond endpoints in
green) so the user can sanity-check the PQ alignment + bond
classification.

Reads cached xtb output:
  appendix_perparation/xtb_frequency_calculations/<step>/R/{reactant.xyz, wbo}
  appendix_perparation/xtb_frequency_calculations/<step>/P/{product.xyz, wbo}

Runs `align_from_arrays` per step to produce the mapping (~0.5 s/step).
P is then Kabsch-aligned to R using the mapped atoms, indexed into
R-frame, and translated rightward so it sits to the right of R in the
viewer.

Output:
  appendix_perparation/viewer/alignment_view.html
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
import time
from pathlib import Path

import numpy as np

from rxn_core_pq import align_from_arrays
from rxn_core_frag import classify_bonds, expand_mapping, build_graph
from align_bgcp_coords import load_cached_xtb, fill_unmapped_greedy


XTB_DIR = PROJECT_ROOT / 'appendix_perparation' / 'xtb_frequency_calculations'
OUT_HTML = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'alignment_view.html'


def kabsch(P, Q):
    """Rotation+translation aligning Q to P."""
    P = np.asarray(P, dtype=float); Q = np.asarray(Q, dtype=float)
    Pc = P.mean(0); Qc = Q.mean(0)
    H = (Q - Qc).T @ (P - Pc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Pc - R @ Qc
    return R, t


def reindex_to_R_frame(elR, elP, xyzP, mapping):
    """Return (elements_in_R_order, coords_in_R_order). For unmapped
    R atoms we fill from the nearest free P atom (already done by
    fill_unmapped_greedy upstream)."""
    n = len(elR)
    out_el = list(elR); out_xyz = np.zeros((n, 3))
    for i in range(n):
        j = mapping.get(i)
        if j is None:
            out_xyz[i] = np.zeros(3)  # gap — should be rare after fallback
            continue
        out_el[i] = elP[j]
        out_xyz[i] = xyzP[j]
    return out_el, out_xyz


def process_step(step):
    rdir = XTB_DIR / step / 'R'
    pdir = XTB_DIR / step / 'P'
    if not (rdir / 'wbo').exists() or not (pdir / 'wbo').exists():
        return None
    elR, xyzR, wboR, _ = load_cached_xtb(rdir)
    elP, xyzP, wboP, _ = load_cached_xtb(pdir)
    res = align_from_arrays(elR, xyzR, wboR, elP, xyzP, wboP)
    pq_mapping = dict(res['mapping'])
    full = fill_unmapped_greedy(elR, xyzR, elP, xyzP, pq_mapping)
    elP_R, xyzP_R = reindex_to_R_frame(elR, elP, xyzP, full)
    # Kabsch P (in R-frame index) onto R using mapped atoms
    mapped_idx = sorted(pq_mapping.keys())
    if mapped_idx:
        Rmat, t = kabsch(np.asarray(xyzR)[mapped_idx],
                         np.asarray(xyzP_R)[mapped_idx])
        xyzP_aligned = (np.asarray(xyzP_R) @ Rmat.T) + t
    else:
        xyzP_aligned = np.asarray(xyzP_R)
    broken = res['broken']; formed = res['formed']
    core_R = sorted({int(i) for (i,j,_,_) in broken} |
                    {int(j) for (i,j,_,_) in broken})
    inv = {v: k for k, v in full.items()}
    core_P = sorted({int(inv[i]) for (i,j,_,_) in formed if i in inv} |
                    {int(inv[j]) for (i,j,_,_) in formed if j in inv})
    bro_pairs = [(int(i), int(j)) for (i,j,_,_) in broken]
    fmd_pairs_R = []
    for (ip, jp, _, _) in formed:
        if ip in inv and jp in inv:
            fmd_pairs_R.append((int(inv[ip]), int(inv[jp])))
    # Translate P so it sits to the right of R in the viewer
    span = float(np.ptp(xyzR, axis=0).max() if len(xyzR) else 5.0)
    offset = max(span * 1.4, 8.0)
    xyzP_view = xyzP_aligned + np.array([offset, 0, 0])
    return {
        'step': step,
        'n_atoms': len(elR),
        'elements_R': list(elR),
        'coords_R':   [[round(float(x), 4) for x in v] for v in xyzR],
        'elements_P': elP_R,
        'coords_P':   [[round(float(x), 4) for x in v] for v in xyzP_view],
        'mapping':    {int(k): int(v) for k, v in full.items()},
        'pq_mapped':  len(pq_mapping),
        'core_R':     core_R,
        'core_P':     core_P,    # already R-frame indices for core
        'broken':     bro_pairs,
        'formed_R':   fmd_pairs_R,
        'offset_x':   offset,
    }


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Atom-alignment inspection — R | P with mapping lines</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select { padding: 4px 6px; font-size: 13px; min-width: 540px; }
 input { padding: 4px 6px; font-size: 13px; }
 .pane { background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
 .viewer { width: 100%; height: 720px; position: relative; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 .lb { background: #ffd6d6; color: #800; }
 .lf { background: #d6f0d6; color: #060; }
 .ll { background: #ddd; color: #444; }
 .stats { color: #444; font-size: 12px; margin-bottom: 6px; }
 button { padding: 4px 10px; }
 input[type=range] { width: 160px; }
</style></head><body>

<h2>BGCP — atom-alignment inspection (R left  ↔  P right)</h2>
<div class="ctl">
  <label><b>Step:</b><select id="stepSel"></select></label>
  <input id="filter" placeholder="filter steps" oninput="rebuildOptions()">
  <button onclick="prev()">◀</button>
  <button onclick="next()">▶</button>
  <label><input type="checkbox" id="lines" checked> show mapping lines</label>
  <label>line opacity <input type="range" id="alpha" min="0.1" max="1.0" step="0.05" value="0.45"></label>
  <label>line radius <input type="range" id="rad" min="0.005" max="0.05" step="0.005" value="0.018"></label>
  <span class="legend" style="margin-left:auto">
    <span class="lb">broken-bond atoms</span>
    <span class="lf">formed-bond atoms</span>
    <span class="ll">mapping lines</span>
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
    opt.textContent = `${n}   N=${d.n_atoms}  pq=${d.pq_mapped}/${d.n_atoms}  broken=${d.broken.length} formed=${d.formed_R.length}`;
    sel.appendChild(opt);
  }
  if (sel.options.length) render(sel.value);
}
sel.addEventListener('change', () => render(sel.value));
function prev() { if (sel.selectedIndex > 0) { sel.selectedIndex--; render(sel.value); } }
function next() { if (sel.selectedIndex < sel.options.length - 1) { sel.selectedIndex++; render(sel.value); } }
document.getElementById('lines').addEventListener('change', () => render(curStep));
document.getElementById('alpha').addEventListener('input', () => render(curStep));
document.getElementById('rad').addEventListener('input', () => render(curStep));

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
    `<b>${name}</b> · ${d.n_atoms} atoms · PQ-mapped ${d.pq_mapped}/${d.n_atoms} · `
    + `broken bonds = ${d.broken.length} · formed bonds = ${d.formed_R.length} · `
    + `core atoms (R) = ${d.core_R.length}`;
  if (viewer) { viewer.removeAllModels(); viewer.removeAllShapes(); }
  viewer = $3Dmol.createViewer('v', {backgroundColor: 'white'});
  // Two models: R first, then P (offset to right by d.offset_x already applied)
  viewer.addModel(buildXyz(d.elements_R, d.coords_R), 'xyz');
  viewer.addModel(buildXyz(d.elements_P, d.coords_P), 'xyz');
  viewer.setStyle({}, {stick: {radius: 0.10}, sphere: {scale: 0.20}});

  // Build a set of broken/formed atom indices
  const broken_atoms = new Set();
  for (const [i, j] of d.broken) { broken_atoms.add(i); broken_atoms.add(j); }
  const formed_atoms = new Set();
  for (const [i, j] of d.formed_R) { formed_atoms.add(i); formed_atoms.add(j); }

  // Highlight broken/formed atoms as colored spheres in BOTH R and P
  for (const i of broken_atoms) {
    const r = d.coords_R[i], p = d.coords_P[i];
    viewer.addSphere({center:{x:r[0], y:r[1], z:r[2]}, radius: 0.35,
                       color: '#cc3333', opacity: 0.55});
    viewer.addSphere({center:{x:p[0], y:p[1], z:p[2]}, radius: 0.35,
                       color: '#cc3333', opacity: 0.55});
  }
  for (const i of formed_atoms) {
    const r = d.coords_R[i], p = d.coords_P[i];
    viewer.addSphere({center:{x:r[0], y:r[1], z:r[2]}, radius: 0.35,
                       color: '#2a8a2a', opacity: 0.55});
    viewer.addSphere({center:{x:p[0], y:p[1], z:p[2]}, radius: 0.35,
                       color: '#2a8a2a', opacity: 0.55});
  }

  // Draw broken bond cylinders inside R, formed bond cylinders inside P
  for (const [a, b] of d.broken) {
    const ra = d.coords_R[a], rb = d.coords_R[b];
    viewer.addCylinder({start:{x:ra[0],y:ra[1],z:ra[2]},
                         end:{x:rb[0],y:rb[1],z:rb[2]},
                         color:'red', radius: 0.10, dashed: true});
  }
  for (const [a, b] of d.formed_R) {
    const pa = d.coords_P[a], pb = d.coords_P[b];
    viewer.addCylinder({start:{x:pa[0],y:pa[1],z:pa[2]},
                         end:{x:pb[0],y:pb[1],z:pb[2]},
                         color:'green', radius: 0.10, dashed: true});
  }

  // Mapping lines: thin cylinders connecting R[i] to P[i]
  if (document.getElementById('lines').checked) {
    const alpha = +document.getElementById('alpha').value;
    const radius = +document.getElementById('rad').value;
    for (let i = 0; i < d.coords_R.length; i++) {
      const r = d.coords_R[i], p = d.coords_P[i];
      // Color: highlight if reactive, else light grey
      let color = '#888';
      if (broken_atoms.has(i)) color = '#cc3333';
      else if (formed_atoms.has(i)) color = '#2a8a2a';
      viewer.addCylinder({start:{x:r[0], y:r[1], z:r[2]},
                          end:  {x:p[0], y:p[1], z:p[2]},
                          color: color, radius: radius, opacity: alpha});
    }
  }

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

    steps = sorted(d.name for d in XTB_DIR.iterdir() if d.is_dir())
    if args.steps:
        steps = [s for s in steps if s in set(args.steps)]
    if args.limit:
        steps = steps[:args.limit]
    print(f"Processing {len(steps)} steps from {XTB_DIR}")

    out = {}
    t0 = time.time(); n_skip = 0
    for i, s in enumerate(steps, 1):
        try:
            r = process_step(s)
            if r is not None:
                out[s] = r
            else:
                n_skip += 1
        except Exception as e:
            n_skip += 1
            print(f"  skip {s}: {e}")
        if i % 20 == 0:
            print(f"  [{i}/{len(steps)}]  ({time.time()-t0:.0f}s)")

    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    html = HTML.replace('__DATA__', json.dumps(out))
    OUT_HTML.write_text(html)
    sz = OUT_HTML.stat().st_size / 1e6
    print(f"\nDone in {time.time()-t0:.0f}s. {len(out)} steps, {n_skip} skipped.")
    print(f"Output: {OUT_HTML}  ({sz:.1f} MB)")


if __name__ == '__main__':
    main()
