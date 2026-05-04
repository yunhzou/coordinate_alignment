"""
Run analysis on a chosen subset of benchmark steps and build a single
combined viewer.html with a dropdown to switch between them.

Usage:
  python build_combined.py            # 10 random steps
  python build_combined.py 20         # 20 random steps
  python build_combined.py --steps a b c   # specific steps
"""

from __future__ import annotations
import argparse
import json
import random
import re
import time
import traceback
from pathlib import Path

from rxn_core_frag import analyze, write_xyz_str


BENCH = Path("/Users/yunhengz/empty_for_claude/Benchmark")
ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core")
OUT = ROOT / "out"
WORK = ROOT / "work"
OUT.mkdir(parents=True, exist_ok=True)


def parse_charge_uhf(xyz_path):
    title = Path(xyz_path).read_text().splitlines()[1] if Path(xyz_path).exists() else ""
    chg, uhf = 0, 0
    m = re.search(r"charge\s*=\s*(-?\d+)", title)
    if m:
        chg = int(m.group(1))
    m = re.search(r"multiplicity\s*=\s*(\d+)", title)
    if m:
        uhf = max(0, int(m.group(1)) - 1)
    return chg, uhf


def list_all_steps():
    out = []
    for d in sorted(BENCH.iterdir()):
        if not d.is_dir():
            continue
        r = d / "plain" / "stage0" / "reactant.xyz"
        p = d / "plain" / "stage0" / "product.xyz"
        if r.exists() and p.exists():
            out.append((d.name, r, p))
    return out


def bond_table_html(rows, side):
    if not rows:
        return "<i>none</i>"
    head = "<tr><th>i</th><th>j</th><th>WBO_R</th><th>WBO_P</th></tr>"
    body = []
    for i, j, wR, wP in rows:
        wR_s = "—" if wR is None else f"{wR:.2f}"
        wP_s = "—" if wP is None else f"{wP:.2f}"
        body.append(f"<tr><td>{i}</td><td>{j}</td><td>{wR_s}</td><td>{wP_s}</td></tr>")
    return f"<table class='bondtab'>{head}{''.join(body)}</table>"


def analyze_step(name, r, p):
    chg, uhf = parse_charge_uhf(r)
    res = analyze(r, p, WORK / name, charge=chg, uhf=uhf)
    elR = res["elements_R"]; elP = res["elements_P"]
    mapping = res["mapping"]
    inv = {v: k for k, v in mapping.items()}
    map_lines = []
    for i in sorted(mapping):
        map_lines.append(f"R[{i:>3}]({elR[i]}) -> P[{mapping[i]:>3}]({elP[mapping[i]]})")
    unmapped_R = [i for i in range(len(elR)) if i not in mapping]
    unmapped_P = [j for j in range(len(elP)) if j not in inv]
    if unmapped_R or unmapped_P:
        map_lines.append("")
        map_lines.append(f"unmapped R: {unmapped_R}")
        map_lines.append(f"unmapped P: {unmapped_P}")

    return {
        "name": name,
        "xyzR": write_xyz_str(elR, res["coords_R"], comment="reactant"),
        "xyzP": write_xyz_str(elP, res["coords_P"], comment="product"),
        "broken_idx": [[i, j] for (i, j, _, _) in res["broken"]],
        "formed_idx": [[i, j] for (i, j, _, _) in res["formed"]],
        "coreR": list(res["core_R"]),
        "coreP": list(res["core_P"]),
        "natoms": len(elR),
        "n_anchors": res["n_anchors"],
        "n_spectator": res["n_after_merge"],
        "n_broken": len(res["broken"]),
        "n_formed": len(res["formed"]),
        "broken_table": bond_table_html(res["broken"], "R"),
        "formed_table": bond_table_html(res["formed"], "P"),
        "mapping_text": "\n".join(map_lines),
        "charge": chg, "uhf": uhf,
    }


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Reaction-core viewer</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select, input { padding: 4px 6px; font-size: 13px; }
 select { min-width: 360px; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 .b { background: #ffd6d6; color: #800; }
 .f { background: #d6f0d6; color: #060; }
 .c { background: #ffe0a8; color: #804000; }
 .stats { color: #444; font-size: 13px; }
 .row { display: flex; gap: 12px; }
 .pane { flex: 1; background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
 .viewer { width: 100%; height: 520px; position: relative; }
 .bondtab { border-collapse: collapse; font-size: 12px; }
 .bondtab td, .bondtab th { border: 1px solid #ccc; padding: 2px 6px; }
 pre { font-size: 11px; background: #f4f4f4; padding: 6px; max-height: 220px; overflow: auto; }
 h3 { margin: 6px 0; }
</style></head><body>

<div class="ctl">
  <label><b>Step:</b>
    <select id="stepSel"></select>
  </label>
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
<div class="pane" style="margin-top:12px">
  <h3>Atom mapping</h3>
  <pre id="mapPre"></pre>
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
function setupOne(divId, prev, xyz, coreSet, changedBonds, color, label) {
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
  vR = setupOne('vR', vR, d.xyzR, d.coreR, d.broken_idx, 'red',   'BROKEN');
  vP = setupOne('vP', vP, d.xyzP, d.coreP, d.formed_idx, 'green', 'FORMED');
  document.getElementById('brokenTab').innerHTML = d.broken_table;
  document.getElementById('formedTab').innerHTML = d.formed_table;
  document.getElementById('mapPre').textContent = d.mapping_text;
  document.getElementById('stats').textContent =
      `N=${d.natoms} | anchors=${d.n_anchors} | spectator=${d.n_spectator} | ` +
      `broken=${d.n_broken} formed=${d.n_formed} | chg/uhf=${d.charge}/${d.uhf}`;
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
    ap.add_argument("n", type=int, nargs="?", default=10)
    ap.add_argument("--steps", nargs="+", default=None)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=str(OUT / "viewer.html"))
    args = ap.parse_args()

    all_steps = list_all_steps()
    if args.steps:
        wanted = set(args.steps)
        chosen = [s for s in all_steps if s[0] in wanted]
    else:
        random.seed(args.seed)
        chosen = random.sample(all_steps, min(args.n, len(all_steps)))

    print(f"[combined] {len(chosen)} steps")
    data = {}
    for k, (name, r, p) in enumerate(chosen, 1):
        t = time.time()
        try:
            d = analyze_step(name, r, p)
            data[name] = d
            print(f"[{k:>2}/{len(chosen)}] {time.time()-t:5.1f}s  OK   {name:<60s} "
                  f"N={d['natoms']:>3}  br/fm={d['n_broken']}/{d['n_formed']}")
        except Exception as e:
            print(f"[{k:>2}/{len(chosen)}] FAIL {name}: {e}")
            traceback.print_exc()

    if not data:
        raise SystemExit("no successful steps")
    html = HTML.replace("__DATA__", json.dumps(data))
    Path(args.out).write_text(html)
    print(f"[combined] wrote {args.out}  ({len(data)} steps)")


if __name__ == "__main__":
    main()
