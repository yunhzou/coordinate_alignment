"""
Atom-alignment inspection viewer.

For each elementary step, render reactant and product side-by-side in
a single 3Dmol scene, with thin grey cylinders connecting each atom
i in R to atom i in P. Atoms touching broken or formed bonds are
highlighted (broken-bond endpoints in red, formed-bond endpoints in
green) so the user can sanity-check the alignment + bond
classification.

Reads ALREADY-ALIGNED geometries (no recomputation):
  appendix_perparation/Pure_Geometries_Elementary_Step/
    Benchmark_Guesses_Coordinate_Aligned_Version/<step>/
      reactants/reactant_aligned.xyz   — R in its own indexing
      products/product_aligned.xyz     — P reindexed into R-frame

Bonds (broken / formed / core) are pulled from the per-step viewer
HTMLs which already contain the classified output:
  appendix_perparation/viewer/mode_viewer/<step>.html

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
import re
import time
from pathlib import Path

import numpy as np


ALIGNED_DIR = (PROJECT_ROOT / 'appendix_perparation' / 'Pure_Geometries_Elementary_Step'
               / 'Benchmark_Guesses_Coordinate_Aligned_Version')
MODE_DIR    = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'mode_viewer'
OUT_HTML    = PROJECT_ROOT / 'appendix_perparation' / 'viewer' / 'alignment_view.html'


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
    """Rotation+translation aligning Q to P."""
    P = np.asarray(P, dtype=float); Q = np.asarray(Q, dtype=float)
    Pc = P.mean(0); Qc = Q.mean(0)
    H = (Q - Qc).T @ (P - Pc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Pc - R @ Qc
    return R, t


def parse_alignment_qc(comment_line):
    """Read 'pq_mapped=N/M  fallback=K  missing=L' from xyz comment."""
    out = {}
    for k, default in (('pq_mapped', '0/0'), ('fallback', '0'), ('missing', '0')):
        m = re.search(rf'{k}=([\d/]+)', comment_line)
        if m: out[k] = m.group(1)
        else: out[k] = default
    return out


def pick_spectator_anchors(xyzR, core_R):
    """Use ALL spectator atoms (those NOT touching a broken/formed
    bond) as the Kabsch reference. If too few spectators (< 5), fall
    back to all atoms — small molecules need every atom to define a
    rigid frame."""
    n_atoms = len(xyzR)
    core_set = set(int(x) for x in core_R)
    spectators = [i for i in range(n_atoms) if i not in core_set]
    if len(spectators) >= 5:
        return sorted(spectators), 'spectators'
    return list(range(n_atoms)), 'all_atoms_fallback'


def process_step(step):
    rdir = ALIGNED_DIR / step / 'reactants'
    pdir = ALIGNED_DIR / step / 'products'
    rfiles = list(rdir.glob('*.xyz')) if rdir.exists() else []
    pfiles = list(pdir.glob('*.xyz')) if pdir.exists() else []
    if not rfiles or not pfiles:
        return None
    elR, xyzR = parse_xyz(rfiles[0])
    elP, xyzP = parse_xyz(pfiles[0])
    if len(elR) != len(elP):
        return None  # shouldn't happen for aligned files

    # P comment line carries pq_mapped / fallback / missing fields
    p_comment = Path(pfiles[0]).read_text().splitlines()[1] if len(pfiles) else ''
    qc = parse_alignment_qc(p_comment)

    # Bonds + core from the per-step viewer HTML (already R-frame)
    bro_pairs = []; fmd_pairs_R = []; core_R = []
    html = MODE_DIR / f"{step}.html"
    if html.exists():
        text = html.read_text()
        m = re.search(r"const DATA = (\{.*?\});\n", text, re.DOTALL)
        if m:
            data = json.loads(m.group(1))
            bro_pairs   = [(int(i), int(j)) for (i, j) in data.get('broken_bonds', [])]
            fmd_pairs_R = [(int(i), int(j)) for (i, j) in data.get('formed_bonds_R', [])]
            core_R      = sorted({int(x) for x in data.get('core_atoms', [])})

    # Spectator-only Kabsch (used to compute "should-be-tiny" deltas).
    anchors, anchor_kind = pick_spectator_anchors(xyzR, core_R)
    Rmat, t = kabsch(xyzR[anchors], xyzP[anchors])
    xyzP_kabsch = xyzP @ Rmat.T + t
    # Per-atom Δ (the "spectator misalignment" signal)
    deltas = np.linalg.norm(xyzP_kabsch - xyzR, axis=1)
    spec_disp = np.array([deltas[i] for i in anchors]) if anchor_kind == 'spectators' else np.array([])

    # SIDE-BY-SIDE layout: translate P along x by ~1.4× molecular span
    # so it sits to the right of R in the viewer. Lines drawn between
    # corresponding atoms span the gap.
    span = float(np.ptp(xyzR, axis=0).max() if len(xyzR) else 5.0)
    offset = max(span * 1.4, 8.0)
    xyzP_view = xyzP_kabsch + np.array([offset, 0, 0])
    return {
        'step': step,
        'n_atoms': len(elR),
        'elements_R': list(elR),
        'coords_R':   [[round(float(x), 4) for x in v] for v in xyzR],
        'elements_P': list(elP),
        'coords_P':   [[round(float(x), 4) for x in v] for v in xyzP_view],
        'anchors':    [int(a) for a in anchors],
        'anchor_kind': anchor_kind,
        'deltas':     [round(float(x), 4) for x in deltas],
        'spec_disp_mean': round(float(spec_disp.mean()), 4) if spec_disp.size else 0.0,
        'spec_disp_max':  round(float(spec_disp.max()),  4) if spec_disp.size else 0.0,
        'pq_mapped':  qc['pq_mapped'],
        'fallback':   qc['fallback'],
        'missing':    qc['missing'],
        'core_R':     core_R,
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
 .viewer { width: 100%; height: 820px; position: relative; }
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
  <label><input type="checkbox" id="numbers"> atom #</label>
  <button onclick="clearPins()">clear pins</button>
  <span style="margin-left:8px; padding:2px 6px; border:1px solid #ddd; border-radius:4px;">
    Mapping batch:
    <button onclick="prevBatch()">◀ prev 5</button>
    <span id="batchInfo">1–5 / N</span>
    <button onclick="nextBatch()">next 5 ▶</button>
    <label style="margin-left:6px;">sort by
      <select id="sortBy">
        <option value="delta_desc">Δ ↓ (largest first)</option>
        <option value="delta_asc">Δ ↑ (smallest first)</option>
        <option value="index">atom # ↑</option>
      </select>
    </label>
    <label style="margin-left:6px;"><input type="checkbox" id="alwaysReactive" checked> always show reactive</label>
  </span>
  <span class="legend" style="margin-left:auto">
    <span class="lb">broken</span>
    <span class="lf">formed</span>
    <span style="background:#aa0033; color:white; padding:2px 8px; border-radius:4px; font-size:12px;">suspect spectator (Δ&gt;thr)</span>
    <span style="background:#ffaa00; color:#553300; padding:2px 8px; border-radius:4px; font-size:12px;">hover</span>
    <span style="background:#ff6699; color:#552233; padding:2px 8px; border-radius:4px; font-size:12px;">click = pin</span>
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
    opt.textContent = `${n}   N=${d.n_atoms}  pq=${d.pq_mapped}  fallback=${d.fallback}  missing=${d.missing}  broken=${d.broken.length} formed=${d.formed_R.length}`;
    sel.appendChild(opt);
  }
  if (sel.options.length) render(sel.value);
}
sel.addEventListener('change', () => render(sel.value));
function prev() { if (sel.selectedIndex > 0) { sel.selectedIndex--; render(sel.value); } }
function next() { if (sel.selectedIndex < sel.options.length - 1) { sel.selectedIndex++; render(sel.value); } }
// Parameter changes rebuild the scene but PRESERVE the camera (no zoomTo).
for (const id of ['numbers', 'sortBy', 'alwaysReactive']) {
  const el = document.getElementById(id);
  el.addEventListener('change', () => { batchStart = 0; render(curStep, true); });
}

// Mapping-line pagination state
const BATCH_SIZE = 5;
let batchStart = 0;     // index into the *sorted* atom list
function prevBatch() { batchStart = Math.max(0, batchStart - BATCH_SIZE); render(curStep, true); }
function nextBatch() {
  const d = DATA[curStep]; if (!d) return;
  if (batchStart + BATCH_SIZE < d.n_atoms) batchStart += BATCH_SIZE;
  render(curStep, true);
}

// Pin state — survives across renders within the same step.
// pinHandles tracks the 3Dmol shape/label objects per pinned atom so we
// can remove them incrementally on click without rebuilding the scene.
let pinned = new Set();
let pinHandles = new Map();   // atomIdx -> {spheres:[2], labels:[2]}

function clearPins() {
  for (const k of Array.from(pinned)) unpinAtom(k);
}

function pinAtom(k) {
  if (pinned.has(k)) return;
  const d = DATA[curStep]; if (!d) return;
  const handles = {spheres: [], labels: []};
  for (const side of [0, 1]) {
    const c = (side === 0) ? d.coords_R[k] : d.coords_P[k];
    if (!c) continue;
    handles.spheres.push(viewer.addSphere(
      {center:{x:c[0], y:c[1], z:c[2]}, radius: 0.50,
       color: '#ff6699', opacity: 0.78}));
    handles.labels.push(viewer.addLabel(`#${k}`,
      {position:{x:c[0], y:c[1], z:c[2] + 0.7},
       fontSize: 13, fontColor: 'white',
       backgroundColor: '#993355', backgroundOpacity: 0.92,
       showBackground: true, inFront: true}));
  }
  pinHandles.set(k, handles);
  pinned.add(k);
  viewer.render();   // does NOT reset camera
}

function unpinAtom(k) {
  const h = pinHandles.get(k);
  if (h) {
    for (const s of h.spheres) viewer.removeShape(s);
    for (const l of h.labels)  viewer.removeLabel(l);
  }
  pinHandles.delete(k);
  pinned.delete(k);
  viewer.render();
}

function buildXyz(elements, coords) {
  let body = `${elements.length}\nframe\n`;
  for (let k = 0; k < elements.length; k++) {
    const c = coords[k];
    body += `${elements[k]}  ${c[0].toFixed(4)}  ${c[1].toFixed(4)}  ${c[2].toFixed(4)}\n`;
  }
  return body;
}

function render(name, preserveView=false) {
  const isStepChange = (name !== curStep);
  if (isStepChange) {
    pinned = new Set();
    pinHandles = new Map();
    batchStart = 0;
  }
  curStep = name;
  const d = DATA[name];
  document.getElementById('info').innerHTML =
    `<b>${name}</b> · ${d.n_atoms} atoms · `
    + `PQ-mapped ${d.pq_mapped} · fallback ${d.fallback} · missing ${d.missing} · `
    + `broken = ${d.broken.length} · formed = ${d.formed_R.length} · `
    + `core (R) = ${d.core_R.length} · `
    + `spectator disp: mean=${d.spec_disp_mean.toFixed(3)}Å, max=${d.spec_disp_max.toFixed(3)}Å`;
  if (viewer) { viewer.removeAllModels(); viewer.removeAllShapes(); viewer.removeAllLabels(); }
  viewer = $3Dmol.createViewer('v', {backgroundColor: 'white'});

  // Side-by-side: R (CPK) on the left, P (CPK) translated to the right.
  viewer.addModel(buildXyz(d.elements_R, d.coords_R), 'xyz');
  viewer.addModel(buildXyz(d.elements_P, d.coords_P), 'xyz');
  viewer.setStyle({}, {stick: {radius: 0.12}, sphere: {scale: 0.22}});

  // Build a set of broken/formed atom indices
  const broken_atoms = new Set();
  for (const [i, j] of d.broken) { broken_atoms.add(i); broken_atoms.add(j); }
  const formed_atoms = new Set();
  for (const [i, j] of d.formed_R) { formed_atoms.add(i); formed_atoms.add(j); }

  // Halo reactive atoms in BOTH R and P
  for (const i of broken_atoms) {
    for (const c of [d.coords_R[i], d.coords_P[i]]) {
      viewer.addSphere({center:{x:c[0],y:c[1],z:c[2]}, radius: 0.35,
                         color: '#cc3333', opacity: 0.45});
    }
  }
  for (const i of formed_atoms) {
    for (const c of [d.coords_R[i], d.coords_P[i]]) {
      viewer.addSphere({center:{x:c[0],y:c[1],z:c[2]}, radius: 0.35,
                         color: '#2a8a2a', opacity: 0.45});
    }
  }
  // Broken bonds in R (left side) and formed bonds in P (right side)
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

  // Build the *atom list to inspect*, sorted per the dropdown.
  // batchStart picks a window of BATCH_SIZE=5 from this list.
  const sortMode = document.getElementById('sortBy').value;
  let order = Array.from({length: d.n_atoms}, (_, i) => i);
  if (sortMode === 'delta_desc') order.sort((a, b) => d.deltas[b] - d.deltas[a]);
  else if (sortMode === 'delta_asc') order.sort((a, b) => d.deltas[a] - d.deltas[b]);
  // 'index' → already in 0..N order

  if (batchStart >= d.n_atoms) batchStart = 0;
  const showFrom = batchStart;
  const showTo = Math.min(batchStart + BATCH_SIZE, d.n_atoms);
  const inBatch = new Set(order.slice(showFrom, showTo));

  // Reactive atoms can be force-shown (toggle) so they're always visible
  // alongside whatever batch the user is looking at.
  const alwaysReactive = document.getElementById('alwaysReactive').checked;
  if (alwaysReactive) {
    for (const i of broken_atoms) inBatch.add(i);
    for (const i of formed_atoms) inBatch.add(i);
  }

  // Update batch info text
  document.getElementById('batchInfo').textContent =
      `${showFrom + 1}–${showTo} / ${d.n_atoms}  (sorted: ${sortMode.replace('_',' ')})`;

  // Mapping cylinders: only for the current batch + (optionally) reactive.
  // Color: red/green for reactive, dark red for "long Δ" suspects, grey
  // for short displacements.
  function lengthColor(L) {
    if (L < 0.10) return '#888';
    if (L < 0.50) return '#dd8800';
    return '#aa0033';
  }
  for (const i of inBatch) {
    const r = d.coords_R[i], p = d.coords_P[i];
    const L = d.deltas[i];
    const isBroken = broken_atoms.has(i);
    const isFormed = formed_atoms.has(i);
    let color, radius, opacity;
    if (isBroken) { color = '#cc3333'; radius = 0.06; opacity = 0.85; }
    else if (isFormed) { color = '#2a8a2a'; radius = 0.06; opacity = 0.85; }
    else { color = lengthColor(L); radius = 0.05; opacity = 0.85; }
    viewer.addCylinder({start:{x:r[0],y:r[1],z:r[2]},
                        end:  {x:p[0],y:p[1],z:p[2]},
                        color: color, radius: radius, opacity: opacity});
    // Atom-# labels at BOTH ends of the line, regardless of the global
    // "atom #" toggle, so the user can read the correspondence.
    for (const c of [r, p]) {
      viewer.addLabel(`#${i}`,
        {position:{x:c[0], y:c[1], z:c[2] + 0.5},
         fontSize: 12, fontColor: 'white',
         backgroundColor: isBroken ? '#cc3333' : (isFormed ? '#2a8a2a' : '#444'),
         backgroundOpacity: 0.9, showBackground: true, inFront: true});
    }
  }

  // Atom-number labels (one per atom, placed at R position only since
  // R and P now overlap in the same coord space). The mapping line
  // shows where P sits relative to that label.
  if (document.getElementById('numbers').checked) {
    for (let i = 0; i < d.coords_R.length; i++) {
      const r = d.coords_R[i];
      viewer.addLabel(String(i),
        {position:{x:r[0], y:r[1], z:r[2] + 0.35},
         fontSize: 9, fontColor: '#444',
         backgroundColor: 'white', backgroundOpacity: 0.6,
         showBackground: true, inFront: true});
    }
  }

  // Re-add pin highlights for atoms still pinned (only happens on
  // parameter change, since step change resets `pinned`).
  pinHandles = new Map();
  const stillPinned = Array.from(pinned);
  pinned = new Set();
  for (const k of stillPinned) pinAtom(k);

  // Hover state — transient, cleared on mouseout. We track shapes/labels
  // in arrays so we can remove them precisely.
  let hoverShapes = [];
  let hoverLabels = [];
  function clearHover() {
    for (const s of hoverShapes) viewer.removeShape(s);
    for (const l of hoverLabels) viewer.removeLabel(l);
    hoverShapes = []; hoverLabels = [];
  }
  function addHoverHighlight(k) {
    for (const side of [0, 1]) {
      const c = (side === 0) ? d.coords_R[k] : d.coords_P[k];
      const e = (side === 0) ? d.elements_R[k] : d.elements_P[k];
      if (!c) continue;
      hoverShapes.push(viewer.addSphere(
        {center:{x:c[0], y:c[1], z:c[2]}, radius: 0.55,
         color: '#ffaa00', opacity: 0.65}));
      hoverLabels.push(viewer.addLabel(`#${k} ${e}` + (side===0 ? ' (R)' : ' (P)'),
        {position:{x:c[0], y:c[1], z:c[2] + 0.7},
         fontSize: 14, fontColor: 'white',
         backgroundColor: '#553300', backgroundOpacity: 0.9,
         showBackground: true, inFront: true}));
    }
  }
  // Resolve a 3Dmol atom-event object to the shared atom index k.
  // The xyz parser assigns atom.serial (1-based, per model) and the
  // model index sits in atom.model. We use .serial-1 as the in-model
  // index k; since R and P share the R-frame indexing, the same k
  // identifies the counterpart on the other side.
  function atomToK(atom) {
    if (atom.serial !== undefined) return atom.serial - 1;
    if (atom.index  !== undefined) return atom.index;
    return null;
  }

  viewer.setHoverable({}, true,
    function(atom) {
      const k = atomToK(atom);
      if (k === null || k < 0 || k >= d.coords_R.length) return;
      clearHover();
      if (!pinned.has(k)) addHoverHighlight(k);
      viewer.render();
    },
    function(_) {
      clearHover();
      viewer.render();
    }
  );

  // Click to pin / unpin. Incremental: only adds or removes the
  // hovered atom's pin shapes; never rebuilds the whole scene, never
  // calls zoomTo. The camera stays exactly where the user left it.
  viewer.setClickable({}, true, function(atom, _viewer) {
    const k = atomToK(atom);
    if (k === null || k < 0 || k >= d.coords_R.length) return;
    if (pinned.has(k)) unpinAtom(k);
    else                pinAtom(k);
  });

  if (isStepChange && !preserveView) viewer.zoomTo();
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
    if args.steps:
        steps = [s for s in steps if s in set(args.steps)]
    if args.limit:
        steps = steps[:args.limit]
    print(f"Processing {len(steps)} steps from {ALIGNED_DIR}")

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
