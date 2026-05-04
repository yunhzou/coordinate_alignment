"""
Run analyze on every elementary step in the TSDisco benchmark
(tsdisco_benchmark_visualization_plain_portable/index.html) and emit
a single combined viewer with a dropdown for switching steps.

For steps with multiple reactants or products, the per-fragment xyz
files are concatenated into single R / P xyz files (xtb single-point
handles disconnected fragments naturally; the WBO graph will have
multiple connected components).
"""

from __future__ import annotations
import argparse
import json
import re
import time
import traceback
from pathlib import Path

from rxn_core_frag import analyze, write_xyz_str


TSDISCO_INDEX = Path(
    "/Users/yunhengz/empty_for_claude/tsdisco_benchmark_visualization_plain_portable/index.html"
)
ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core")
OUT = ROOT / "out"
WORK = ROOT / "work_tsdisco"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)


def load_tsdisco_data():
    html = TSDISCO_INDEX.read_text()
    m = re.search(r"const TSDISCO_DATA = (\{.*?\});\s*\n", html, re.DOTALL)
    if m is None:
        raise SystemExit("Could not locate TSDISCO_DATA in index.html")
    return json.loads(m.group(1))


def _parse_xyz_text(txt):
    import numpy as np
    lines = txt.strip().splitlines()
    n = int(lines[0].split()[0])
    elements = []
    coords = []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        elements.append(parts[0])
        coords.append([float(x) for x in parts[1:4]])
    return elements, np.array(coords)


def concat_xyz(xyz_texts, gap=10.0):
    """Concatenate multiple xyz texts. Multi-fragment inputs are
    translated apart along the x-axis with `gap` Å between bounding
    boxes, so xtb sees them as well-separated species (no spurious
    inter-fragment WBOs and no atom overlap from shared TS-frame
    pre-alignment)."""
    if len(xyz_texts) == 1:
        return xyz_texts[0]
    import numpy as np
    parsed = [_parse_xyz_text(t) for t in xyz_texts]
    elements = []
    coords_list = []
    cursor = 0.0
    for els, coords in parsed:
        if coords.size == 0:
            continue
        bb_min = coords.min(axis=0)
        bb_max = coords.max(axis=0)
        # Place fragment so its bounding box starts at x=cursor; center y/z at 0
        shift = np.array([
            cursor - bb_min[0],
            -(bb_max[1] + bb_min[1]) / 2.0,
            -(bb_max[2] + bb_min[2]) / 2.0,
        ])
        coords_list.append(coords + shift)
        elements.extend(els)
        cursor += (bb_max[0] - bb_min[0]) + gap
    all_coords = np.vstack(coords_list)
    n = len(elements)
    body = "\n".join(f"{el}  {c[0]:.6f}  {c[1]:.6f}  {c[2]:.6f}"
                     for el, c in zip(elements, all_coords))
    return f"{n}\nconcatenated_separated\n{body}\n"


def write_xyz_file(text, path):
    Path(path).write_text(text)


def step_inputs(step):
    """Return (reactant_xyz_text, product_xyz_text, charge, uhf)."""
    reacts = [s for s in step["input_structures"] if s["role"] == "reactant"]
    prods = [s for s in step["input_structures"] if s["role"] == "product"]
    rxyz = concat_xyz([s["xyz"] for s in reacts])
    pxyz = concat_xyz([s["xyz"] for s in prods])
    chg = step.get("charge", 0) or 0
    mult = step.get("multiplicity", 1) or 1
    uhf = max(0, mult - 1)
    return rxyz, pxyz, chg, uhf


def bond_table_html(rows):
    if not rows:
        return "<i>none</i>"
    head = "<tr><th>i</th><th>j</th><th>WBO_R</th><th>WBO_P</th></tr>"
    body = []
    for i, j, wR, wP in rows:
        wR_s = "—" if wR is None else f"{wR:.2f}"
        wP_s = "—" if wP is None else f"{wP:.2f}"
        body.append(f"<tr><td>{i}</td><td>{j}</td><td>{wR_s}</td><td>{wP_s}</td></tr>")
    return f"<table class='bondtab'>{head}{''.join(body)}</table>"


def analyze_step(step, work_dir):
    rxyz_text, pxyz_text, chg, uhf = step_inputs(step)
    work_dir.mkdir(parents=True, exist_ok=True)
    r_path = work_dir / "reactant.xyz"
    p_path = work_dir / "product.xyz"
    write_xyz_file(rxyz_text, r_path)
    write_xyz_file(pxyz_text, p_path)
    res = analyze(r_path, p_path, work_dir, charge=chg, uhf=uhf)
    elR = res["elements_R"]; elP = res["elements_P"]
    return {
        "name": f"{step['dataset']}/{step['step_id']}",
        "xyzR": write_xyz_str(elR, res["coords_R"], comment="reactant"),
        "xyzP": write_xyz_str(elP, res["coords_P"], comment="product"),
        "broken_idx": [[i, j] for (i, j, _, _) in res["broken"]],
        "formed_idx": [[i, j] for (i, j, _, _) in res["formed"]],
        "natoms": len(elR),
        "n_broken": len(res["broken"]),
        "n_formed": len(res["formed"]),
        "n_mapped": len(res["mapping"]),
        "broken_table": bond_table_html(res["broken"]),
        "formed_table": bond_table_html(res["formed"]),
        "charge": chg,
        "uhf": uhf,
        "mechanism": step.get("mechanism", ""),
        "step": step.get("step", ""),
    }


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>TSDisco bond-change viewer</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select, input { padding: 4px 6px; font-size: 13px; }
 select { min-width: 460px; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 .b { background: #ffd6d6; color: #800; }
 .f { background: #d6f0d6; color: #060; }
 .stats { color: #444; font-size: 13px; }
 .row { display: flex; gap: 12px; }
 .pane { flex: 1; background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
 .viewer { width: 100%; height: 520px; position: relative; }
 .bondtab { border-collapse: collapse; font-size: 12px; }
 .bondtab td, .bondtab th { border: 1px solid #ccc; padding: 2px 6px; }
 h3 { margin: 6px 0; }
</style></head><body>

<div class="ctl">
  <label><b>Step:</b><select id="stepSel"></select></label>
  <input id="filter" placeholder="filter (substring)">
  <button onclick="prevStep()">◀ Prev</button>
  <button onclick="nextStep()">Next ▶</button>
  <span class="stats" id="stats"></span>
  <span class="legend" style="margin-left:auto">
    <span class="b">broken</span><span class="f">formed</span>
  </span>
</div>

<div class="row">
  <div class="pane"><h3>Reactant</h3><div id="vR" class="viewer"></div></div>
  <div class="pane"><h3>Product</h3><div id="vP" class="viewer"></div></div>
</div>
<div class="row" style="margin-top:12px">
  <div class="pane"><h3>Broken bonds</h3><div id="brokenTab"></div></div>
  <div class="pane"><h3>Formed bonds</h3><div id="formedTab"></div></div>
</div>

<script>
const STEPS = __DATA__;
const stepNames = Object.keys(STEPS);
const sel = document.getElementById('stepSel');
function rebuildOptions(filter='') {
  sel.innerHTML = '';
  const f = filter.toLowerCase();
  stepNames.filter(n => n.toLowerCase().includes(f)).forEach(n => {
    const d = STEPS[n];
    const o = document.createElement('option');
    o.value = n;
    o.textContent = `${n}  (N=${d.natoms} br=${d.n_broken} fm=${d.n_formed})`;
    sel.appendChild(o);
  });
}
rebuildOptions();
document.getElementById('filter').addEventListener('input', e => {
  rebuildOptions(e.target.value);
  if (sel.options.length) loadStep(sel.options[0].value);
});

function parseXYZCoords(xyz) {
  const lines = xyz.trim().split('\n');
  const n = parseInt(lines[0]);
  const out = [];
  for (let i = 0; i < n; i++) {
    const parts = lines[2+i].trim().split(/\s+/);
    out.push([+parts[1], +parts[2], +parts[3]]);
  }
  return out;
}

function dashedCylinder(viewer, p1, p2, color, radius, nDashes) {
  const dx = p2[0]-p1[0], dy = p2[1]-p1[1], dz = p2[2]-p1[2];
  const onFrac = 0.55;
  for (let k = 0; k < nDashes; k++) {
    const t1 = (k + 0.0) / nDashes;
    const t2 = (k + onFrac) / nDashes;
    viewer.addCylinder({
      start: {x: p1[0]+dx*t1, y: p1[1]+dy*t1, z: p1[2]+dz*t1},
      end:   {x: p1[0]+dx*t2, y: p1[1]+dy*t2, z: p1[2]+dz*t2},
      radius: radius, fromCap: 2, toCap: 2, color: color,
    });
  }
}

let vR = null, vP = null;
function setupOne(divId, prev, xyz, changedBonds, color, label) {
  let v;
  if (prev) { prev.clear(); v = prev; }
  else { v = $3Dmol.createViewer(divId, {backgroundColor: 'white'}); }
  v.addModel(xyz, 'xyz');
  v.setStyle({}, {stick: {radius: 0.12}, sphere: {scale: 0.22}});
  const coords = parseXYZCoords(xyz);
  changedBonds.forEach(([i, j]) => {
    const p1 = coords[i], p2 = coords[j];
    const len = Math.hypot(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]);
    const nDashes = Math.max(4, Math.round(len * 4));
    dashedCylinder(v, p1, p2, color, 0.16, nDashes);
    v.addLabel(label, {
      position: {x: (p1[0]+p2[0])/2, y: (p1[1]+p2[1])/2, z: (p1[2]+p2[2])/2},
      backgroundColor: color, fontColor: 'white', fontSize: 11,
      borderThickness: 0, padding: 2, inFront: true,
    });
  });
  v.zoomTo();
  v.render();
  return v;
}

function loadStep(name) {
  const d = STEPS[name];
  if (!d) return;
  vR = setupOne('vR', vR, d.xyzR, d.broken_idx, 'red',   'BROKEN');
  vP = setupOne('vP', vP, d.xyzP, d.formed_idx, 'green', 'FORMED');
  document.getElementById('brokenTab').innerHTML = d.broken_table;
  document.getElementById('formedTab').innerHTML = d.formed_table;
  document.getElementById('stats').textContent =
      `N=${d.natoms} | mapped=${d.n_mapped} | broken=${d.n_broken} formed=${d.n_formed} | chg/uhf=${d.charge}/${d.uhf}`;
  if (sel.value !== name) sel.value = name;
}
sel.addEventListener('change', e => loadStep(e.target.value));
function prevStep() {
  const opts = [...sel.options].map(o => o.value);
  const i = opts.indexOf(sel.value);
  if (i > 0) loadStep(opts[i-1]);
}
function nextStep() {
  const opts = [...sel.options].map(o => o.value);
  const i = opts.indexOf(sel.value);
  if (i >= 0 && i < opts.length - 1) loadStep(opts[i+1]);
}
window.addEventListener('keydown', e => {
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
  if (e.key === 'ArrowRight') nextStep();
  if (e.key === 'ArrowLeft') prevStep();
});

if (sel.options.length) loadStep(sel.options[0].value);
</script>
</body></html>
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None,
                    help="run only first N steps (debugging)")
    ap.add_argument("--out", default=str(OUT / "tsdisco_viewer.html"))
    ap.add_argument("--start", type=int, default=0)
    args = ap.parse_args()

    tsdisco = load_tsdisco_data()
    steps = tsdisco["steps"]
    if args.limit is not None:
        steps = steps[args.start:args.start + args.limit]
    else:
        steps = steps[args.start:]
    print(f"[tsdisco] {len(steps)} steps to process")

    data = {}
    for k, step in enumerate(steps, 1):
        name = f"{step['dataset']}/{step['step_id']}"
        sanitized = re.sub(r'[^A-Za-z0-9._-]', '_', name)
        wd = WORK / sanitized
        t = time.time()
        try:
            d = analyze_step(step, wd)
            data[name] = d
            print(f"[{k:>3}/{len(steps)}]  {time.time()-t:5.1f}s  OK   "
                  f"{name:<70s}  N={d['natoms']:>3}  br/fm={d['n_broken']}/{d['n_formed']}")
        except Exception as e:
            print(f"[{k:>3}/{len(steps)}]  FAIL {name}: {e}")
            traceback.print_exc()

    if not data:
        raise SystemExit("no successful steps")
    html = HTML.replace("__DATA__", json.dumps(data))
    Path(args.out).write_text(html)
    print(f"[tsdisco] wrote {args.out}  ({len(data)} steps)")


if __name__ == "__main__":
    main()
