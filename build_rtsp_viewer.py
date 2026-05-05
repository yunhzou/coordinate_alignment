"""
Build a 3-panel R -> TS -> P viewer for the tsdisco benchmark.

For each step, we align the rank-0 candidate TS to R using the same
multi-seed algorithm. The TS is then re-indexed into R-atom space, so
all three structures share one atom-numbering convention. Bonds are
classified at three transitions:

  R-bonds: WBO_R >= 0.5 (existing in reactant)
  R->TS:   bonds where |WBO_R - WBO_TS| >= 0.5
  TS->P:   bonds where |WBO_TS - WBO_P| >= 0.5
  R->P:    bonds where |WBO_R - WBO_P| >= 0.5  (overall chemistry)

Each bond gets categorized as broken (lost between two consecutive
panels) or formed (gained). The TS panel highlights "in-between" bonds
whose WBO is partway between the broken/formed extremes.

Usage:
  python build_rtsp_viewer.py            # all 160 steps
  python build_rtsp_viewer.py --limit 10
  python build_rtsp_viewer.py --steps name1 name2 ...
"""
from __future__ import annotations
import argparse
import json
import re
import time
import traceback
from pathlib import Path

import numpy as np

from rxn_core_frag import (
    run_xtb, build_graph, find_islands, expand_mapping,
    classify_bonds, write_xyz_str,
    _generate_seed_orders,
)
from build_tsdisco_viewer import step_inputs, _parse_xyz_text


TSDISCO_INDEX = Path(
    "/Users/yunhengz/empty_for_claude/tsdisco_benchmark_visualization_plain_portable/index.html"
)
ROOT = Path("/Users/yunhengz/empty_for_claude/rxn_core")
OUT = ROOT / "out"
WORK = ROOT / "work_rtsp"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)


def load_data():
    html = TSDISCO_INDEX.read_text()
    m = re.search(r"const TSDISCO_DATA = (\{.*?\});\s*\n", html, re.DOTALL)
    return json.loads(m.group(1))


def best_mapping(g_R, g_P, wboR, wboP, n_seeds=10):
    """Multi-seed mapping search; same scoring as analyze."""
    orders = _generate_seed_orders(g_R, n_seeds)
    best = None
    for order in orders:
        m, _ = find_islands(g_R, g_P, seed_order=order)
        m = expand_mapping(m, g_R, g_P)
        br, fm, _, _ = classify_bonds(m, wboR, wboP)
        score = (len(br) + len(fm), -len(m))
        if best is None or score < best[0]:
            best = (score, m)
    return best[1]


def reindex_xyz_to_target(target_idx_to_src_idx, src_elements, src_coords,
                          n_target, fallback_elements, fallback_coords):
    """Build (elements, coords) in target-atom-index order. Unmapped
    target indices fall back to the corresponding R-frame atom (so the
    viewer renders the R geometry for atoms that didn't align)."""
    out_elements = list(fallback_elements)
    out_coords = np.array(fallback_coords, dtype=float).copy()
    for tgt, src in target_idx_to_src_idx.items():
        out_elements[tgt] = src_elements[src]
        out_coords[tgt] = src_coords[src]
    return out_elements, out_coords


def changed_bonds(wbo_a, wbo_b, threshold=0.5, scan=0.5):
    """Return list of (i, j, w_a, w_b, delta) for atom pairs whose
    WBO changes by >= threshold and at least one side has wbo >= scan."""
    out = []
    n = wbo_a.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            w_a = float(wbo_a[i, j]); w_b = float(wbo_b[i, j])
            if max(w_a, w_b) < scan:
                continue
            d = w_b - w_a
            if abs(d) >= threshold:
                out.append((i, j, w_a, w_b, d))
    return out


def all_bonds(wbo, threshold=0.5):
    """Return list of (i, j, wbo) for all bonds with wbo >= threshold."""
    n = wbo.shape[0]
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            w = float(wbo[i, j])
            if w >= threshold:
                out.append((i, j, w))
    return out


def analyze_step(step):
    """Run xtb on R, TS, P; align all three to a common (R) index space.
    Return everything the viewer needs."""
    name = f"{step['dataset']}/{step['step_id']}"
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)
    chg = step.get("charge", 0) or 0
    uhf = max(0, (step.get("multiplicity", 1) or 1) - 1)

    rxyz_text, pxyz_text, _, _ = step_inputs(step)
    cands = step.get("candidates", [])
    if not cands or not cands[0].get("xyz"):
        raise RuntimeError("no TS candidate")
    ts_xyz_text = cands[0]["xyz"]

    (wd / "reactant.xyz").write_text(rxyz_text)
    (wd / "product.xyz").write_text(pxyz_text)
    (wd / "ts.xyz").write_text(ts_xyz_text)

    # xtb on all three (cached)
    elR, xyzR, wboR = run_xtb(wd / "reactant.xyz", wd / "R", charge=chg, uhf=uhf)
    elP, xyzP, wboP = run_xtb(wd / "product.xyz", wd / "P", charge=chg, uhf=uhf)
    elTS, xyzTS, wboTS = run_xtb(wd / "ts.xyz", wd / "TS", charge=chg, uhf=uhf)

    g_R = build_graph(elR, wboR)
    g_P = build_graph(elP, wboP)
    g_TS = build_graph(elTS, wboTS)

    # Align P to R, and TS to R, both using multi-seed
    map_R_to_P = best_mapping(g_R, g_P, wboR, wboP)
    map_R_to_TS = best_mapping(g_R, g_TS, wboR, wboTS)

    # Reindex P and TS into R-atom-index space (fall back to R for any
    # unmapped target indices so the viewer can still render the geometry)
    n_R = len(elR)
    elP_r, xyzP_r = reindex_xyz_to_target(map_R_to_P, elP, xyzP, n_R, elR, xyzR)
    wboP_r = np.zeros_like(wboR)
    for ri, pi in map_R_to_P.items():
        for rj, pj in map_R_to_P.items():
            if ri < rj:
                wboP_r[ri, rj] = wboP[pi, pj]
                wboP_r[rj, ri] = wboP[pi, pj]

    elTS_r, xyzTS_r = reindex_xyz_to_target(map_R_to_TS, elTS, xyzTS, n_R, elR, xyzR)
    wboTS_r = np.zeros_like(wboR)
    for ri, ti in map_R_to_TS.items():
        for rj, tj in map_R_to_TS.items():
            if ri < rj:
                wboTS_r[ri, rj] = wboTS[ti, tj]
                wboTS_r[rj, ri] = wboTS[ti, tj]

    # Classify bonds at each transition
    rt = changed_bonds(wboR, wboTS_r)   # R -> TS
    tp = changed_bonds(wboTS_r, wboP_r) # TS -> P
    rp = changed_bonds(wboR, wboP_r)    # R -> P (overall)

    # Categorize: bond is "in transition" at TS if its WBO is partway
    # between extreme. e.g., R has WBO 0.9 (bond), P has WBO 0.0 (no
    # bond), TS has 0.4 -- the bond is half-broken at TS.
    inflight = []
    for (i, j, wR_, wP_, d) in rp:
        wT = float(wboTS_r[i, j])
        # in-flight iff TS wbo is not at either extreme (within scan_threshold)
        # i.e., midway: wT is between min(wR,wP) and max(wR,wP)
        lo = min(wR_, wP_); hi = max(wR_, wP_)
        if lo + 0.15 < wT < hi - 0.15:
            inflight.append((i, j, wR_, wT, wP_))

    return {
        "name": name,
        "natoms": len(elR),
        "elements": elR,
        "xyzR": write_xyz_str(elR, xyzR, comment="reactant (R-frame)"),
        "xyzTS": write_xyz_str(elTS_r, xyzTS_r, comment="TS (re-indexed to R)"),
        "xyzP": write_xyz_str(elP_r, xyzP_r, comment="product (re-indexed to R)"),
        "bonds_R": all_bonds(wboR),
        "bonds_TS": all_bonds(wboTS_r),
        "bonds_P": all_bonds(wboP_r),
        "rt_changes": rt,   # R -> TS
        "tp_changes": tp,   # TS -> P
        "rp_changes": rp,   # R -> P overall
        "inflight": inflight,
        "n_R_to_TS": len(rt),
        "n_TS_to_P": len(tp),
        "n_R_to_P": len(rp),
        "n_inflight": len(inflight),
        "charge": chg, "uhf": uhf,
    }


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>R -> TS -> P viewer</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select { min-width: 460px; padding: 4px 6px; font-size: 13px; }
 input { padding: 4px 6px; font-size: 13px; }
 .row { display: flex; gap: 12px; }
 .pane { flex: 1; background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
 .viewer { width: 100%; height: 480px; position: relative; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 .b { background: #ffd6d6; color: #800; }
 .f { background: #d6f0d6; color: #060; }
 .i { background: #ffe9b3; color: #804000; }
 .stats { color: #444; font-size: 13px; }
 .bondtab { border-collapse: collapse; font-size: 11px; width: 100%; }
 .bondtab td, .bondtab th { border: 1px solid #ccc; padding: 2px 6px; }
 h3 { margin: 6px 0; font-size: 14px; }
</style></head><body>

<div class="ctl">
  <label><b>Step:</b><select id="stepSel"></select></label>
  <input id="filter" placeholder="filter (substring)">
  <button onclick="prevStep()">◀ Prev</button>
  <button onclick="nextStep()">Next ▶</button>
  <span class="stats" id="stats"></span>
  <span class="legend" style="margin-left:auto">
    <span class="b">broken</span><span class="f">formed</span><span class="i">in-flight</span>
  </span>
</div>

<div class="row">
  <div class="pane"><h3>Reactant</h3><div id="vR" class="viewer"></div></div>
  <div class="pane"><h3>Transition state</h3><div id="vT" class="viewer"></div></div>
  <div class="pane"><h3>Product</h3><div id="vP" class="viewer"></div></div>
</div>

<div class="row" style="margin-top:12px">
  <div class="pane"><h3>R → TS changes</h3><div id="rtTab"></div></div>
  <div class="pane"><h3>TS → P changes</h3><div id="tpTab"></div></div>
  <div class="pane"><h3>R → P (overall)</h3><div id="rpTab"></div></div>
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
    o.textContent = `${n}  (R→TS=${d.n_R_to_TS} TS→P=${d.n_TS_to_P} flight=${d.n_inflight})`;
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

let vR = null, vT = null, vP = null;
function setupViewer(divId, prev, xyz, dashedBonds) {
  let v;
  if (prev) { prev.clear(); v = prev; }
  else { v = $3Dmol.createViewer(divId, {backgroundColor: 'white'}); }
  v.addModel(xyz, 'xyz');
  v.setStyle({}, {stick: {radius: 0.10}, sphere: {scale: 0.20}});
  const coords = parseXYZCoords(xyz);
  dashedBonds.forEach(([i, j, color, label]) => {
    const p1 = coords[i], p2 = coords[j];
    const len = Math.hypot(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]);
    const nDashes = Math.max(4, Math.round(len * 4));
    dashedCylinder(v, p1, p2, color, 0.16, nDashes);
    v.addLabel(label, {
      position: {x: (p1[0]+p2[0])/2, y: (p1[1]+p2[1])/2, z: (p1[2]+p2[2])/2},
      backgroundColor: color, fontColor: 'white', fontSize: 10,
      borderThickness: 0, padding: 2, inFront: true,
    });
  });
  v.zoomTo();
  v.render();
  return v;
}

function bondTable(rows, columns) {
  if (!rows || !rows.length) return "<i>none</i>";
  const head = "<tr>" + columns.map(c => `<th>${c}</th>`).join("") + "</tr>";
  const body = rows.map(r => "<tr>" + r.map(c => `<td>${c}</td>`).join("") + "</tr>").join("");
  return `<table class="bondtab">${head}${body}</table>`;
}

function loadStep(name) {
  const d = STEPS[name];
  if (!d) return;
  const els = d.elements;

  // Per-panel bond cylinders use that panel's native atom indices.
  // R panel:  rp_changes (R-indexed). delta<0 = BROKEN.
  // P panel:  rp_changes_P (P-indexed). delta>0 = FORMED.
  // TS panel: rp_changes_TS + inflight_TS (TS-indexed).
  const rDashes = [];
  const pDashes = [];
  const tDashes = [];
  d.rp_changes.forEach(([i, j, wR, wP, delta]) => {
    if (delta < 0) {
      rDashes.push([i, j, 'red', 'BROKEN']);
    }
  });
  (d.rp_changes_P || []).forEach(([i, j, wR, wP, delta]) => {
    if (delta > 0) {
      pDashes.push([i, j, 'green', 'FORMED']);
    }
  });
  (d.rp_changes_TS || []).forEach(([i, j, wR, wP, delta]) => {
    const c = delta < 0 ? 'red' : 'green';
    tDashes.push([i, j, c, delta < 0 ? 'breaking' : 'forming']);
  });
  (d.inflight_TS || []).forEach(([i, j, wR, wT, wP]) => {
    tDashes.push([i, j, '#ffa500', `flight WBO=${wT.toFixed(2)}`]);
  });

  vR = setupViewer('vR', vR, d.xyzR, rDashes);
  vT = setupViewer('vT', vT, d.xyzTS, tDashes);
  vP = setupViewer('vP', vP, d.xyzP, pDashes);

  document.getElementById('rtTab').innerHTML = bondTable(
    d.rt_changes.map(([i,j,wR,wT,delta]) => [
      `${i}(${els[i]})`, `${j}(${els[j]})`,
      wR.toFixed(2), wT.toFixed(2), delta.toFixed(2)
    ]),
    ['i','j','WBO_R','WBO_TS','ΔWBO']
  );
  document.getElementById('tpTab').innerHTML = bondTable(
    d.tp_changes.map(([i,j,wT,wP,delta]) => [
      `${i}(${els[i]})`, `${j}(${els[j]})`,
      wT.toFixed(2), wP.toFixed(2), delta.toFixed(2)
    ]),
    ['i','j','WBO_TS','WBO_P','ΔWBO']
  );
  document.getElementById('rpTab').innerHTML = bondTable(
    d.rp_changes.map(([i,j,wR,wP,delta]) => [
      `${i}(${els[i]})`, `${j}(${els[j]})`,
      wR.toFixed(2), wP.toFixed(2), delta.toFixed(2)
    ]),
    ['i','j','WBO_R','WBO_P','ΔWBO']
  );

  document.getElementById('stats').textContent =
      `N=${d.natoms} | R→TS=${d.n_R_to_TS} TS→P=${d.n_TS_to_P} R→P=${d.n_R_to_P} | in-flight=${d.n_inflight} | chg/uhf=${d.charge}/${d.uhf}`;
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
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--steps", nargs="+", default=None)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--out", default=str(OUT / "rtsp_viewer.html"))
    args = ap.parse_args()

    data = load_data()
    steps = data["steps"]
    if args.steps:
        wanted = set(args.steps)
        steps = [s for s in steps if f"{s['dataset']}/{s['step_id']}" in wanted]
    elif args.limit is not None:
        steps = steps[args.start:args.start + args.limit]
    else:
        steps = steps[args.start:]

    print(f"[rtsp] {len(steps)} steps to process")
    out_data = {}
    for k, step in enumerate(steps, 1):
        name = f"{step['dataset']}/{step['step_id']}"
        t = time.time()
        try:
            d = analyze_step(step)
            out_data[name] = d
            print(f"[{k:>3}/{len(steps)}]  {time.time()-t:5.1f}s  OK   "
                  f"{name:<70s}  R→TS={d['n_R_to_TS']} TS→P={d['n_TS_to_P']} flight={d['n_inflight']}")
        except Exception as e:
            print(f"[{k:>3}/{len(steps)}]  FAIL {name}: {e}")
            traceback.print_exc()

    if not out_data:
        raise SystemExit("no successful steps")
    html = HTML.replace("__DATA__", json.dumps(out_data))
    Path(args.out).write_text(html)
    print(f"[rtsp] wrote {args.out}  ({len(out_data)} steps)")


if __name__ == "__main__":
    main()
