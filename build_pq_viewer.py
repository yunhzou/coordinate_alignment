"""
General BGCP viewer using the NEW PQ algorithm (rxn_core_pq.analyze_pq).

Renders all 160 BGCP steps as a single HTML with a dropdown, R/P 3D viewers
side by side, broken/formed bond cylinders, and a regression marker per step
(tie / win / loose-regression / strict-regression vs OLD rxn_core_frag).

Output:
  out/bgcp_pq_viewer.html        — single combined viewer
  out/bgcp_pq_regressions.txt    — text list of regressed steps with details

Usage: python build_pq_viewer.py
"""
from __future__ import annotations
import csv
import json
import re
import sys
import time
import traceback
from pathlib import Path

from rxn_core_pq import analyze_pq
from rxn_core_frag import write_xyz_str
from build_bgcp_viewer import (
    BGCP_ROOT, WORK, LOOKUP, list_step_dirs, read_xyzs,
)


ROOT = Path(__file__).parent
OUT_HTML = ROOT / "out" / "bgcp_pq_viewer.html"
OUT_REG = ROOT / "out" / "bgcp_pq_regressions.txt"
OLD_CSV = ROOT / "out" / "bgcp_old_bonds.csv"


def classify(o, x):
    """Compare OLD csv row vs NEW (live PQ) row. Return one of:
    tie / win / loose_regress / strict_regress."""
    if not o or not o.get('n_atoms'):
        return 'unknown'
    ot = int(o['n_broken']) + int(o['n_formed'])
    nt = x['n_broken'] + x['n_formed']
    om = int(o['n_mapped']); nm = x['n_mapped']
    if nt < ot:
        return 'win'
    elif nt > ot:
        return 'strict_regress' if nm < om else 'loose_regress'
    else:
        return 'tie'


def bond_table_rows(bonds, elements):
    if not bonds: return "<tr><td colspan=4><i>none</i></td></tr>"
    out = []
    for i, j, wR, wP in bonds:
        wR_s = '—' if wR is None else f'{wR:.2f}'
        wP_s = '—' if wP is None else f'{wP:.2f}'
        out.append(f"<tr><td>{i}({elements[i]})</td><td>{j}({elements[j]})</td>"
                   f"<td>{wR_s}</td><td>{wP_s}</td></tr>")
    return ''.join(out)


def main():
    # Load OLD baseline for classification.
    old_rows = {r['step_id']: r for r in csv.DictReader(open(OLD_CSV))}

    print(f"Running PQ analyze on all BGCP steps...")
    step_data = {}
    by_class = {'tie': [], 'win': [], 'loose_regress': [], 'strict_regress': [],
                'unknown': []}

    for i, sd in enumerate(list_step_dirs(), 1):
        name = sd.name
        chg, uhf = LOOKUP.get(name, (0, 0))
        sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", name)
        wd = WORK / sanitized
        wd.mkdir(parents=True, exist_ok=True)
        rxyz = read_xyzs(sd / "reactants")
        pxyz = read_xyzs(sd / "products")
        if rxyz is None or pxyz is None:
            continue
        (wd / "reactant.xyz").write_text(rxyz)
        (wd / "product.xyz").write_text(pxyz)
        try:
            t0 = time.time()
            out = analyze_pq(wd / "reactant.xyz", wd / "product.xyz", wd,
                             charge=chg, uhf=uhf)
            elapsed = time.time() - t0
        except Exception as e:
            print(f"[{i}/160] {name}: ERROR {e}")
            traceback.print_exc()
            continue

        elR = out['elements_R']; elP = out['elements_P']
        broken = out['broken']; formed = out['formed']
        mapping = out['mapping']
        inv = {v: k for k, v in mapping.items()}

        # Bond indices for cylinder drawing.
        # broken bonds are R-side indices (drawn on R viewer in red dashed)
        # formed bonds are P-side indices (drawn on P viewer in green dashed)
        broken_idx = [[int(a), int(b)] for (a, b, _, _) in broken]
        formed_idx_P = [[int(a), int(b)] for (a, b, _, _) in formed]

        cls = classify(old_rows.get(name), {
            'n_broken': out['n_broken'],
            'n_formed': out['n_formed'],
            'n_mapped': out['n_mapped'],
        })
        by_class[cls].append((name, out, old_rows.get(name)))

        step_data[name] = {
            'name': name,
            'class': cls,
            'n_atoms': len(elR),
            'n_mapped': out['n_mapped'],
            'n_broken': out['n_broken'],
            'n_formed': out['n_formed'],
            'chir': out['chirality_violations'],
            'charge': chg, 'uhf': uhf,
            'old_broken': int(old_rows.get(name, {}).get('n_broken', 0) or 0),
            'old_formed': int(old_rows.get(name, {}).get('n_formed', 0) or 0),
            'old_mapped': int(old_rows.get(name, {}).get('n_mapped', 0) or 0),
            'xyzR': write_xyz_str(elR, out['coords_R'], comment='reactant'),
            'xyzP': write_xyz_str(elP, out['coords_P'], comment='product'),
            'broken_idx': broken_idx,
            'formed_idx_P': formed_idx_P,
            'broken_table': bond_table_rows(broken, elR),
            'formed_table': bond_table_rows(formed, elP),
        }
        print(f"[{i:3d}/160] {name[:55]:55s} {cls:14s} "
              f"old={old_rows.get(name, {}).get('n_broken','?')}/"
              f"{old_rows.get(name, {}).get('n_formed','?')}  "
              f"new={out['n_broken']}/{out['n_formed']}  t={elapsed:.1f}s")

    # ---- Write regression text file ----
    with OUT_REG.open('w') as f:
        f.write("PQ vs OLD comparison summary across BGCP\n")
        f.write("=" * 70 + "\n\n")
        for cat in ['strict_regress', 'loose_regress', 'win', 'tie']:
            lst = by_class[cat]
            f.write(f"## {cat.upper()} ({len(lst)} steps)\n")
            if cat in ('strict_regress', 'loose_regress'):
                for name, out, old in sorted(
                        lst, key=lambda t: -((t[1]['n_broken']+t[1]['n_formed'])
                                             - (int(t[2]['n_broken'])+int(t[2]['n_formed'])))):
                    ot = int(old['n_broken']) + int(old['n_formed'])
                    nt = out['n_broken'] + out['n_formed']
                    f.write(f"  +{nt-ot:<2}  OLD {old['n_broken']}/{old['n_formed']}/"
                            f"{old['n_mapped']:>3}  PQ {out['n_broken']}/"
                            f"{out['n_formed']}/{out['n_mapped']:>3}  {name}\n")
                    elR = out['elements_R']; elP = out['elements_P']
                    br_str = ','.join(f'{i}({elR[i]})-{j}({elR[j]})' for i, j, _, _ in out['broken'])
                    fm_str = ','.join(f'{i}({elP[i]})-{j}({elP[j]})' for i, j, _, _ in out['formed'])
                    f.write(f"      PQ broken: {br_str}\n")
                    f.write(f"      PQ formed: {fm_str}\n")
            elif cat == 'win':
                for name, out, old in sorted(
                        lst, key=lambda t: -((int(t[2]['n_broken'])+int(t[2]['n_formed']))
                                             - (t[1]['n_broken']+t[1]['n_formed']))):
                    ot = int(old['n_broken']) + int(old['n_formed'])
                    nt = out['n_broken'] + out['n_formed']
                    f.write(f"  -{ot-nt:<2}  OLD {old['n_broken']}/{old['n_formed']}/"
                            f"{old['n_mapped']:>3}  PQ {out['n_broken']}/"
                            f"{out['n_formed']}/{out['n_mapped']:>3}  {name}\n")
            f.write("\n")
    print(f"\nRegression text: {OUT_REG}")

    # ---- Write HTML viewer ----
    html = HTML.replace('__DATA__', json.dumps(step_data))
    OUT_HTML.write_text(html)
    print(f"HTML viewer:     {OUT_HTML}")
    print()
    for cat in ['strict_regress', 'loose_regress', 'win', 'tie', 'unknown']:
        print(f"  {cat:14s}: {len(by_class[cat])}")


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>BGCP PQ alignment viewer</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
 body { font-family: -apple-system, sans-serif; margin: 12px; background: #fafafa; }
 .ctl { background: white; padding: 8px 12px; border: 1px solid #ddd; border-radius: 6px; margin-bottom: 12px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
 select { min-width: 480px; padding: 4px 6px; font-size: 13px; }
 input { padding: 4px 6px; font-size: 13px; }
 .filter-btns button { margin-right: 4px; padding: 3px 8px; font-size: 12px; cursor: pointer; }
 .filter-btns button.active { background: #444; color: white; }
 .badge { display: inline-block; padding: 2px 8px; border-radius: 4px; font-size: 12px; font-weight: 600; }
 .b-tie { background: #eee; color: #444; }
 .b-win { background: #d6f0d6; color: #060; }
 .b-loose_regress { background: #ffe6c2; color: #803; }
 .b-strict_regress { background: #ffd6d6; color: #800; }
 .b-unknown { background: #ddd; color: #666; }
 .stats { color: #444; font-size: 13px; }
 .row { display: flex; gap: 12px; }
 .pane { flex: 1; background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
 .viewer { width: 100%; height: 480px; position: relative; }
 .legend span { display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }
 .lb { background: #ffd6d6; color: #800; }
 .lf { background: #d6f0d6; color: #060; }
 table { border-collapse: collapse; font-size: 12px; }
 td, th { border: 1px solid #ccc; padding: 2px 6px; }
 h3 { margin: 6px 0; }
</style></head><body>

<h2>BGCP — NEW PQ algorithm (rxn_core_pq) vs OLD (rxn_core_frag)</h2>
<div class="ctl">
  <label><b>Step:</b><select id="stepSel"></select></label>
  <input id="filter" placeholder="filter (substring)" oninput="rebuildOptions()">
  <span class="filter-btns">
    <button id="b-all" onclick="setCat('all')" class="active">all</button>
    <button id="b-strict_regress" onclick="setCat('strict_regress')">strict reg</button>
    <button id="b-loose_regress" onclick="setCat('loose_regress')">loose reg</button>
    <button id="b-win" onclick="setCat('win')">wins</button>
    <button id="b-tie" onclick="setCat('tie')">ties</button>
  </span>
  <button onclick="prevStep()">◀</button>
  <button onclick="nextStep()">▶</button>
  <span class="stats" id="stats"></span>
  <span class="legend" style="margin-left:auto"><span class="lb">broken</span><span class="lf">formed</span></span>
</div>

<div class="row">
  <div class="pane"><h3>Reactant <span id="sR" class="stats"></span></h3><div id="vR" class="viewer"></div></div>
  <div class="pane"><h3>Product <span id="sP" class="stats"></span></h3><div id="vP" class="viewer"></div></div>
</div>
<div class="row" style="margin-top:12px">
  <div class="pane"><h3>Broken bonds (R indices)</h3>
    <table><thead><tr><th>i</th><th>j</th><th>WBO_R</th><th>WBO_P</th></tr></thead>
    <tbody id="brokenTab"></tbody></table>
  </div>
  <div class="pane"><h3>Formed bonds (P indices)</h3>
    <table><thead><tr><th>i</th><th>j</th><th>WBO_R</th><th>WBO_P</th></tr></thead>
    <tbody id="formedTab"></tbody></table>
  </div>
</div>

<script>
const STEPS = __DATA__;
const stepNames = Object.keys(STEPS);
const sel = document.getElementById('stepSel');
let curCat = 'all';

function rebuildOptions() {
  const f = document.getElementById('filter').value.toLowerCase();
  sel.innerHTML = '';
  let names = stepNames.filter(n => n.toLowerCase().includes(f));
  if (curCat !== 'all') names = names.filter(n => STEPS[n].class === curCat);
  for (const n of names) {
    const d = STEPS[n];
    const opt = document.createElement('option');
    opt.value = n;
    const cls = d.class.replace('_regress', '_reg');
    opt.textContent = `[${cls}] ${n}  (br/fm OLD ${d.old_broken}/${d.old_formed} → PQ ${d.n_broken}/${d.n_formed})`;
    sel.appendChild(opt);
  }
  if (sel.options.length) render(sel.value);
}

function setCat(c) {
  curCat = c;
  document.querySelectorAll('.filter-btns button').forEach(b => b.classList.remove('active'));
  document.getElementById('b-' + c).classList.add('active');
  rebuildOptions();
}
function prevStep() { if (sel.selectedIndex > 0) { sel.selectedIndex--; render(sel.value); } }
function nextStep() { if (sel.selectedIndex < sel.options.length - 1) { sel.selectedIndex++; render(sel.value); } }
sel.addEventListener('change', () => render(sel.value));

let vR, vP;
function render(name) {
  const d = STEPS[name];
  if (!d) return;
  document.getElementById('stats').innerHTML =
    `<span class='badge b-${d.class}'>${d.class}</span> ` +
    `N=${d.n_atoms} · OLD ${d.old_broken}/${d.old_formed} mapped ${d.old_mapped} · ` +
    `PQ ${d.n_broken}/${d.n_formed} mapped ${d.n_mapped} · chir ${d.chir}`;
  document.getElementById('sR').textContent = `(charge=${d.charge}, uhf=${d.uhf})`;
  document.getElementById('sP').textContent = '';
  document.getElementById('vR').innerHTML = '';
  document.getElementById('vP').innerHTML = '';
  vR = $3Dmol.createViewer('vR', {backgroundColor: 'white'});
  vR.addModel(d.xyzR, 'xyz');
  vR.setStyle({}, {stick: {radius: 0.10}, sphere: {scale: 0.20}});
  // broken bonds drawn on R
  const atomsR = vR.selectedAtoms({});
  for (const [i, j] of d.broken_idx) {
    if (i < atomsR.length && j < atomsR.length) {
      vR.addCylinder({start: {x:atomsR[i].x,y:atomsR[i].y,z:atomsR[i].z},
                      end:   {x:atomsR[j].x,y:atomsR[j].y,z:atomsR[j].z},
                      color: 'red', radius: 0.10, dashed: true});
    }
  }
  vR.zoomTo(); vR.render();

  vP = $3Dmol.createViewer('vP', {backgroundColor: 'white'});
  vP.addModel(d.xyzP, 'xyz');
  vP.setStyle({}, {stick: {radius: 0.10}, sphere: {scale: 0.20}});
  const atomsP = vP.selectedAtoms({});
  for (const [i, j] of d.formed_idx_P) {
    if (i < atomsP.length && j < atomsP.length) {
      vP.addCylinder({start: {x:atomsP[i].x,y:atomsP[i].y,z:atomsP[i].z},
                      end:   {x:atomsP[j].x,y:atomsP[j].y,z:atomsP[j].z},
                      color: 'green', radius: 0.10, dashed: true});
    }
  }
  vP.zoomTo(); vP.render();

  document.getElementById('brokenTab').innerHTML = d.broken_table;
  document.getElementById('formedTab').innerHTML = d.formed_table;
}

rebuildOptions();
</script>
</body></html>
"""


if __name__ == "__main__":
    main()
