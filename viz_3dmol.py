"""
3Dmol.js HTML visualization for reaction-core analysis.

Outputs a single self-contained HTML file with two side-by-side viewers
(reactant left, product right). Broken bonds are drawn as red dashed
cylinders on the reactant; formed bonds as green dashed cylinders on the
product. Reaction-core atoms are highlighted with an orange sphere outline.
"""

from __future__ import annotations
import json
from pathlib import Path

from rxn_core_wbo import write_xyz_str


HTML_TEMPLATE = """<!doctype html>
<html><head>
<meta charset="utf-8">
<title>Reaction core: {title}</title>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  body {{ font-family: -apple-system, sans-serif; margin: 0; padding: 12px; background: #fafafa; }}
  h2 {{ margin: 6px 0; }}
  .row {{ display: flex; gap: 12px; }}
  .pane {{ flex: 1; background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }}
  .viewer {{ width: 100%; height: 520px; position: relative; }}
  .legend span {{ display: inline-block; padding: 2px 8px; margin-right: 6px; border-radius: 4px; font-size: 12px; }}
  .b {{ background: #ffd6d6; color: #800; }}
  .f {{ background: #d6f0d6; color: #060; }}
  .c {{ background: #ffe0a8; color: #804000; }}
  pre {{ font-size: 11px; background: #f4f4f4; padding: 6px; max-height: 220px; overflow: auto; }}
  table {{ font-size: 12px; border-collapse: collapse; }}
  td, th {{ border: 1px solid #ccc; padding: 2px 6px; }}
</style>
</head><body>
<h2>{title}</h2>
<div class="legend">
  <span class="b">broken bond</span>
  <span class="f">formed bond</span>
  <span class="c">core atom</span>
  <span style="color:#666">spectator atoms in default CPK</span>
</div>
<div class="row">
  <div class="pane">
    <h3>Reactant</h3>
    <div id="vR" class="viewer"></div>
  </div>
  <div class="pane">
    <h3>Product</h3>
    <div id="vP" class="viewer"></div>
  </div>
</div>
<div class="row">
  <div class="pane"><h3>Broken bonds</h3>{table_broken}</div>
  <div class="pane"><h3>Formed bonds</h3>{table_formed}</div>
</div>
<div class="pane" style="margin-top:12px"><h3>Atom mapping (R -> P, 0-indexed)</h3><pre>{mapping_text}</pre></div>

<script>
const xyzR = {xyzR_json};
const xyzP = {xyzP_json};
const brokenBonds = {broken_json};   // [[i,j], ...] indices into reactant
const formedBonds = {formed_json};   // [[i,j], ...] indices into product
const coreR = {coreR_json};          // 0-indexed atom indices in reactant
const coreP = {coreP_json};          // 0-indexed atom indices in product
const coordsR = {coordsR_json};
const coordsP = {coordsP_json};

function dashedCylinder(viewer, p1, p2, color, radius, nDashes) {{
  // Build a manual dashed line out of short solid cylinders (3Dmol's
  // built-in dashed: true on addCylinder is unreliable across versions).
  const dx = p2[0] - p1[0], dy = p2[1] - p1[1], dz = p2[2] - p1[2];
  const onFrac = 0.55;  // duty cycle: 55% on, 45% off
  for (let k = 0; k < nDashes; k++) {{
    const t1 = (k + 0.0) / nDashes;
    const t2 = (k + onFrac) / nDashes;
    viewer.addCylinder({{
      start: {{x: p1[0]+dx*t1, y: p1[1]+dy*t1, z: p1[2]+dz*t1}},
      end:   {{x: p1[0]+dx*t2, y: p1[1]+dy*t2, z: p1[2]+dz*t2}},
      radius: radius, fromCap: 2, toCap: 2, color: color,
    }});
  }}
}}

function setupViewer(id, xyz, coreSet, changedBonds, coords, color, label) {{
  let v = $3Dmol.createViewer(id, {{backgroundColor: 'white'}});
  v.addModel(xyz, 'xyz');
  v.setStyle({{}}, {{stick: {{radius: 0.12}}, sphere: {{scale: 0.22}}}});
  if (coreSet.length) {{
    v.setStyle({{serial: coreSet.map(i => i + 1)}}, {{
      stick: {{radius: 0.18, color: 'orange'}},
      sphere: {{scale: 0.34, color: 'orange'}},
    }});
  }}
  changedBonds.forEach(([i, j]) => {{
    const p1 = coords[i], p2 = coords[j];
    const len = Math.hypot(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]);
    const nDashes = Math.max(4, Math.round(len * 4));
    dashedCylinder(v, p1, p2, color, 0.16, nDashes);
    // Midpoint label so the change is unmistakable
    v.addLabel(label, {{
      position: {{x: (p1[0]+p2[0])/2, y: (p1[1]+p2[1])/2, z: (p1[2]+p2[2])/2}},
      backgroundColor: color, fontColor: 'white', fontSize: 11,
      borderThickness: 0, padding: 2, inFront: true,
    }});
  }});
  v.zoomTo();
  v.render();
  return v;
}}

setupViewer('vR', xyzR, coreR, brokenBonds, coordsR, 'red',   'BROKEN');
setupViewer('vP', xyzP, coreP, formedBonds, coordsP, 'green', 'FORMED');
</script>
</body></html>
"""


def _bond_table(rows, side):
    if not rows:
        return "<i>none</i>"
    head = "<tr><th>i</th><th>j</th><th>WBO_R</th><th>WBO_P</th></tr>"
    body = []
    for i, j, wR, wP in rows:
        wR_s = "—" if wR is None else f"{wR:.2f}"
        wP_s = "—" if wP is None else f"{wP:.2f}"
        body.append(f"<tr><td>{i}</td><td>{j}</td><td>{wR_s}</td><td>{wP_s}</td></tr>")
    return f"<table>{head}{''.join(body)}</table>"


def render_html(result, title, out_path):
    elR, xyzR = result["elements_R"], result["coords_R"]
    elP, xyzP = result["elements_P"], result["coords_P"]
    xyzR_str = write_xyz_str(elR, xyzR, comment="reactant")
    xyzP_str = write_xyz_str(elP, xyzP, comment="product")

    broken_idx = [[i, j] for (i, j, _, _) in result["broken"]]
    formed_idx = [[i, j] for (i, j, _, _) in result["formed"]]

    mapping = result["mapping"]
    map_lines = []
    for i in sorted(mapping):
        map_lines.append(f"  R[{i:>3}]({elR[i]}) -> P[{mapping[i]:>3}]({elP[mapping[i]]})")
    unmapped_R = [i for i in range(len(elR)) if i not in mapping]
    if unmapped_R:
        inv = {v: k for k, v in mapping.items()}
        unmapped_P = [j for j in range(len(elP)) if j not in inv]
        map_lines.append("")
        map_lines.append(f"  Unmapped R atoms: {unmapped_R}")
        map_lines.append(f"  Unmapped P atoms: {unmapped_P}")

    html = HTML_TEMPLATE.format(
        title=title,
        xyzR_json=json.dumps(xyzR_str),
        xyzP_json=json.dumps(xyzP_str),
        broken_json=json.dumps(broken_idx),
        formed_json=json.dumps(formed_idx),
        coreR_json=json.dumps(list(result["core_R"])),
        coreP_json=json.dumps(list(result["core_P"])),
        coordsR_json=json.dumps(xyzR.tolist()),
        coordsP_json=json.dumps(xyzP.tolist()),
        table_broken=_bond_table(result["broken"], "R"),
        table_formed=_bond_table(result["formed"], "P"),
        mapping_text="\n".join(map_lines),
    )
    Path(out_path).write_text(html)
    return out_path
