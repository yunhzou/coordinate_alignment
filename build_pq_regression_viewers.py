"""
For each of the 7 PQ-vs-OLD regression cases, run analyze_pq and write
a static 3D R/P viewer to out/regressions/<step>/pq_result.html.

The OLD-algorithm 10-seed traces already live alongside as seed_*.html.
This adds the NEW-algorithm result so they sit side-by-side.

Usage: python build_pq_regression_viewers.py
"""
from __future__ import annotations
import json
import re
from pathlib import Path

from rxn_core_pq import analyze_pq
from rxn_core_frag import write_xyz_str
from build_bgcp_viewer import BGCP_ROOT, LOOKUP, read_xyzs


REGRESSIONS = [
    "pr7.V.dodh_ts56-triplet",
    "pr7.V.dodh_ts56-singlet",
    "pr7.V.dodh_ts71",
    "pr7.V.dodh_ts1314",
    "Jackie_TS_10",
    "pr19.heck_ts1",
    "pr16.carbocation_ts5",
]

OUT_ROOT = Path(__file__).parent / "out" / "regressions"
WORK = Path(__file__).parent / "work_bgcp"


HTML_TEMPLATE = r"""<!doctype html>
<html><head><meta charset="utf-8">
<title>__TITLE__ — PQ result</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
body{font-family:-apple-system,sans-serif;margin:12px;background:#fafafa}
.row{display:flex;gap:12px}
.pane{flex:1;background:white;border:1px solid #ddd;border-radius:6px;padding:8px}
.viewer{width:100%;height:560px;position:relative}
.legend span{display:inline-block;padding:2px 8px;margin-right:6px;border-radius:4px;font-size:12px}
.b{background:#ffd6d6;color:#800}
.f{background:#d6f0d6;color:#060}
.u{background:#f8e3a3;color:#663300}
table{border-collapse:collapse;font-size:12px;margin-top:6px}
td,th{border:1px solid #ccc;padding:2px 6px}
h3{margin:6px 0}
.stats{font-size:13px;color:#444;margin-bottom:8px}
</style></head><body>
<h2>__TITLE__ — PQ algorithm result (rxn_core_pq)</h2>
<div class="stats">
  N atoms: __N__ | mapped: __MAPPED__/__N__ | broken: __NBR__ | formed: __NFM__ |
  chirality violations: __CHIR__ | charge: __CHG__ | uhf: __UHF__
</div>
<div><span class="legend"><span class="b">broken bonds (in R)</span>
<span class="f">formed bonds (in P)</span>
<span class="u">unmapped atoms (highlighted on both)</span></span></div>
<p><a href="index.html">↑ this step's seed comparison</a> |
<a href="../index.html">↑↑ all regressions</a></p>

<div class="row">
  <div class="pane"><h3>Reactant</h3><div id="vR" class="viewer"></div></div>
  <div class="pane"><h3>Product (R-frame indices)</h3><div id="vP" class="viewer"></div></div>
</div>

<div class="row" style="margin-top:12px">
  <div class="pane"><h3>Broken bonds (R indices)</h3>__BROKEN_TABLE__</div>
  <div class="pane"><h3>Formed bonds (P indices, mapped to R where possible)</h3>__FORMED_TABLE__</div>
</div>

<script>
const xyzR = __XYZR__;
const xyzP = __XYZP__;
const broken = __BROKEN_IDX__;     // [[i,j], ...] in R indices
const formed_R = __FORMED_R_IDX__; // [[i,j], ...] using mapping^-1; null if unmapped
const formed_P = __FORMED_P_IDX__; // [[ip,jp], ...] in P indices
const unmapped_R = __UNMAPPED_R__;
const unmapped_P = __UNMAPPED_P__;

function setup(div, xyz, broken_pairs, formed_pairs, unmapped) {
  const v = $3Dmol.createViewer(div, {backgroundColor: 'white'});
  v.addModel(xyz, 'xyz');
  v.setStyle({}, {stick: {radius: 0.12}, sphere: {radius: 0.30}});
  // highlight unmapped atoms in orange
  for (const i of unmapped) {
    v.setStyle({serial: i + 1}, {stick: {radius: 0.12}, sphere: {radius: 0.50, color: 'orange'}});
  }
  // dashed bond cylinders
  const atoms = v.selectedAtoms({});
  function pos(i) {
    const a = atoms[i];
    return {x: a.x, y: a.y, z: a.z};
  }
  for (const [i, j] of broken_pairs) {
    if (i == null || j == null || i >= atoms.length || j >= atoms.length) continue;
    v.addCylinder({start: pos(i), end: pos(j), color: 'red',
                   radius: 0.10, dashed: true});
  }
  for (const [i, j] of formed_pairs) {
    if (i == null || j == null || i >= atoms.length || j >= atoms.length) continue;
    v.addCylinder({start: pos(i), end: pos(j), color: 'green',
                   radius: 0.10, dashed: true});
  }
  v.zoomTo();
  v.render();
}

setup('vR', xyzR, broken, formed_R.filter(p => p[0] != null && p[1] != null),
      unmapped_R);
setup('vP', xyzP, broken.filter(p => false /* draw broken on R only */),
      formed_P, unmapped_P);
</script>
</body></html>
"""


def bond_table_html(rows):
    if not rows: return "<i>none</i>"
    head = "<tr><th>i</th><th>j</th><th>WBO_R</th><th>WBO_P</th></tr>"
    body = []
    for i, j, wR, wP in rows:
        wR_s = "—" if wR is None else f"{wR:.2f}"
        wP_s = "—" if wP is None else f"{wP:.2f}"
        body.append(f"<tr><td>{i}</td><td>{j}</td><td>{wR_s}</td><td>{wP_s}</td></tr>")
    return f"<table>{head}{''.join(body)}</table>"


def run_step(step):
    chg, uhf = LOOKUP.get(step, (0, 0))
    sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step)
    wd = WORK / sanitized
    wd.mkdir(parents=True, exist_ok=True)
    rxyz = read_xyzs(BGCP_ROOT / step / "reactants")
    pxyz = read_xyzs(BGCP_ROOT / step / "products")
    if rxyz is None or pxyz is None:
        raise RuntimeError(f"missing R or P for {step}")
    (wd / "reactant.xyz").write_text(rxyz)
    (wd / "product.xyz").write_text(pxyz)

    out = analyze_pq(wd / "reactant.xyz", wd / "product.xyz", wd,
                     charge=chg, uhf=uhf)

    elR = out["elements_R"]; elP = out["elements_P"]
    xyzR_arr = out["coords_R"]; xyzP_arr = out["coords_P"]
    mapping = out["mapping"]
    inv = {v: k for k, v in mapping.items()}

    broken_idx = [[i, j] for (i, j, _, _) in out["broken"]]
    formed_idx_P = [[ip, jp] for (ip, jp, _, _) in out["formed"]]
    formed_idx_R = [[inv.get(ip), inv.get(jp)] for (ip, jp, _, _) in out["formed"]]

    unmapped_R = [i for i in range(len(elR)) if i not in mapping]
    unmapped_P = [i for i in range(len(elP)) if i not in inv]

    xyzR_str = write_xyz_str(elR, xyzR_arr, comment="R")
    xyzP_str = write_xyz_str(elP, xyzP_arr, comment="P")

    out_dir = OUT_ROOT / sanitized
    out_dir.mkdir(parents=True, exist_ok=True)
    html = HTML_TEMPLATE
    repls = {
        "__TITLE__": step,
        "__N__": str(len(elR)),
        "__MAPPED__": str(out["n_mapped"]),
        "__NBR__": str(out["n_broken"]),
        "__NFM__": str(out["n_formed"]),
        "__CHIR__": str(out["chirality_violations"]),
        "__CHG__": str(chg),
        "__UHF__": str(uhf),
        "__BROKEN_TABLE__": bond_table_html(out["broken"]),
        "__FORMED_TABLE__": bond_table_html(out["formed"]),
        "__XYZR__": json.dumps(xyzR_str),
        "__XYZP__": json.dumps(xyzP_str),
        "__BROKEN_IDX__": json.dumps(broken_idx),
        "__FORMED_R_IDX__": json.dumps(formed_idx_R),
        "__FORMED_P_IDX__": json.dumps(formed_idx_P),
        "__UNMAPPED_R__": json.dumps(unmapped_R),
        "__UNMAPPED_P__": json.dumps(unmapped_P),
    }
    for k, v in repls.items():
        html = html.replace(k, v)
    out_path = out_dir / "pq_result.html"
    out_path.write_text(html)
    return out, out_path


def main():
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary = []
    for step in REGRESSIONS:
        try:
            out, path = run_step(step)
            print(f"{step:50s} br/fm={out['n_broken']}/{out['n_formed']:<2} "
                  f"mapped={out['n_mapped']}/{len(out['elements_R'])} "
                  f"chir={out['chirality_violations']:<2} -> {path}")
            summary.append((step, out, path))
        except Exception as e:
            print(f"{step}: ERROR {e}")
            import traceback; traceback.print_exc()

    # update top-level index to add a "PQ result" link per step
    rows = ""
    for step, out, path in summary:
        sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", step)
        rows += (
            f"<tr><td>{step}</td>"
            f"<td>{len(out['elements_R'])}</td>"
            f"<td>{out['n_broken']}/{out['n_formed']} (mapped {out['n_mapped']})</td>"
            f"<td><a href='{sanitized}/pq_result.html'>PQ viewer</a></td>"
            f"<td><a href='{sanitized}/index.html'>OLD 10-seed traces</a></td></tr>"
        )
    top = f"""<!doctype html><html><head><meta charset="utf-8">
<title>PQ regression diagnostic</title>
<style>body{{font-family:-apple-system,sans-serif;margin:20px;max-width:1100px}}
table{{border-collapse:collapse}}
th,td{{border:1px solid #ccc;padding:6px 10px;font-size:13px}}</style>
</head><body>
<h2>PQ-vs-OLD regressions: side-by-side</h2>
<p>Per step: <b>PQ viewer</b> shows the new algorithm's static R/P result with
broken (red dashed) / formed (green dashed) bonds and unmapped atoms in orange.
<b>OLD 10-seed traces</b> shows the old algorithm with 10 random seed orderings
as slider-driven animations.</p>
<table>
<tr><th>step</th><th>N</th><th>PQ br/fm (mapped)</th><th>NEW</th><th>OLD</th></tr>
{rows}
</table>
</body></html>"""
    (OUT_ROOT / "index.html").write_text(top)
    print(f"\ntop index: {OUT_ROOT / 'index.html'}")


if __name__ == "__main__":
    main()
