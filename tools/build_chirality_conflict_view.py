#!/usr/bin/env python3
"""Build a self-contained viewer for rejected index-chirality mechanisms."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]


def read_xyz(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[0].strip())
    elements, coords = [], []
    for line in lines[2:2 + count]:
        fields = line.split()
        elements.append(fields[0])
        coords.append([float(value) for value in fields[1:4]])
    if len(elements) != count:
        raise ValueError(f"truncated XYZ: {path}")
    return elements, coords


def rejected_records(case_root: Path):
    summary = json.loads((case_root / "summary.json").read_text())
    if summary.get("status") == "ok":
        stage = json.loads((case_root / "rp_stage.json").read_text())
        return summary, stage.get("rejected_index_chirality") or []
    diagnostics = summary.get("diagnostics") or {}
    return summary, diagnostics.get("rejected_mechanisms") or []


def build(case_root: Path, selection_manifest: Path, output: Path):
    summary, rejected = rejected_records(case_root)
    if not rejected:
        raise ValueError(f"{case_root.name}: no rejected mechanisms")
    selection = json.loads(selection_manifest.read_text())
    source = next(
        case for case in selection["cases"]
        if case["step_id"] == case_root.name)
    elements_R, coords_R = read_xyz(Path(source["reactant_xyz"]))
    elements_P, coords_P = read_xyz(Path(source["product_xyz"]))
    mechanisms = []
    for rejected_mechanism in rejected:
        mapping = {
            int(r): int(p) for r, p in
            rejected_mechanism["source_mapping_RP"].items()
        }
        product_in_R = [coords_P[mapping[r]] for r in range(len(elements_R))]
        product_elements_in_R = [elements_P[mapping[r]]
                                 for r in range(len(elements_R))]
        if product_elements_in_R != elements_R:
            raise ValueError("rejected mapping changes atom elements")
        diagnostics = rejected_mechanism.get("diagnostics") or {}
        mechanisms.append({
            "id": rejected_mechanism["source_mechanism_id"],
            "reason": rejected_mechanism["reason"],
            "mapping_RP": [mapping[r] for r in range(len(elements_R))],
            "product_in_R": product_in_R,
            "switchable_r_atoms": diagnostics.get("switchable_r_atoms") or [],
            "frames": diagnostics.get("source_mismatch_frames") or [],
        })
    data = {
        "step_id": case_root.name,
        "status": summary.get("status"),
        "elements": elements_R,
        "reactant": coords_R,
        "rejected_mechanisms": mechanisms,
    }
    three_dmol = (
        REPOSITORY / "src/rxn_core/static/3Dmol-min.js").read_text()
    html = HTML.replace("__THREEDMOL__", three_dmol).replace(
        "__DATA__", json.dumps(data))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    return output


HTML = r'''<!doctype html><html><head><meta charset="utf-8">
<title>Index-chirality conflict</title><script>__THREEDMOL__</script>
<style>
body{font-family:-apple-system,sans-serif;margin:14px;background:#f6f8fb;color:#182235}
h2{margin:0 0 4px}.sub{font:12px ui-monospace,monospace;color:#53657a;margin-bottom:10px}
.controls,.panel,.details{background:white;border:1px solid #d5dce5;border-radius:7px;padding:9px}
.controls{display:flex;gap:12px;align-items:center;margin-bottom:10px}.controls select{max-width:550px}
select,label{font:12px ui-monospace,monospace}.grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}
.viewer{height:520px;position:relative}.viewer>div{position:absolute;inset:0}
.ph{font-weight:650}.details{margin-top:10px;font:12px/1.55 ui-monospace,monospace;white-space:pre-wrap}
.legend{font-size:12px;margin:5px 0}.center{color:#d7191c}.neighbor{color:#f28e2b}.mutable{color:#0099b8}
</style></head><body>
<h2 id="title"></h2><div class="sub">Rejected AAM mechanism diagnostic · source mapping shown without replacement</div>
<div class="controls"><label>Mechanism <select id="mechanism"></select></label>
<label>Conflicting frame <select id="frame"></select></label>
<label><input id="allLabels" type="checkbox"> all atom labels</label></div>
<div class="legend"><span class="center">red = frame center</span> · <span class="neighbor">orange = ordered neighbor triple</span> · <span class="mutable">cyan halo = AAM-authorized mutable atom</span></div>
<div class="grid"><div class="panel"><div class="ph">Reactant (R indices)</div><div class="viewer"><div id="viewR"></div></div></div>
<div class="panel"><div class="ph">Product under rejected source mapping (R indices)</div><div class="viewer"><div id="viewP"></div></div></div></div>
<div class="details" id="details"></div>
<script>
const DATA=__DATA__;let vr=null,vp=null;
const $=id=>document.getElementById(id);
function xyz(coords){return DATA.elements.length+'\nconflict\n'+DATA.elements.map((e,i)=>e+' '+coords[i].join(' ')).join('\n')+'\n'}
function currentMechanism(){return DATA.rejected_mechanisms[Number($('mechanism').value)]}
function currentFrame(){return currentMechanism().frames[Number($('frame').value)]}
function draw(div,coords){div.innerHTML='';const v=$3Dmol.createViewer(div,{backgroundColor:'white'});v.addModel(xyz(coords),'xyz');v.setStyle({},{stick:{radius:.10},sphere:{scale:.20}});const m=currentMechanism(),f=currentFrame();for(const a of m.switchable_r_atoms)v.addSphere({center:{x:coords[a][0],y:coords[a][1],z:coords[a][2]},radius:.29,color:'#00bcd4',alpha:.22});v.setStyle({serial:f.center_R},{stick:{radius:.16,color:'#d7191c'},sphere:{scale:.48,color:'#d7191c'}});f.neighbors_R_index_order.forEach((a,k)=>{v.setStyle({serial:a},{stick:{radius:.14,color:'#f28e2b'},sphere:{scale:.40,color:'#f28e2b'}});v.addLabel(String(k+1),{position:{x:coords[a][0],y:coords[a][1],z:coords[a][2]},fontSize:12,fontColor:'#8a3f00',backgroundColor:'white',backgroundOpacity:.8,inFront:true})});const labels=$('allLabels').checked?[...coords.keys()]:[f.center_R,...f.neighbors_R_index_order];for(const a of labels)v.addLabel('R'+a,{position:{x:coords[a][0],y:coords[a][1],z:coords[a][2]},fontSize:9,fontColor:'black',backgroundColor:'white',backgroundOpacity:.65,inFront:true});v.zoomTo();v.render();return v}
function rebuildFrames(){const m=currentMechanism();$('frame').innerHTML='';m.frames.forEach((f,i)=>{const o=document.createElement('option');o.value=i;o.textContent=f.id;$('frame').appendChild(o)});render()}
function render(){const m=currentMechanism(),f=currentFrame();vr=draw($('viewR'),DATA.reactant);vp=draw($('viewP'),m.product_in_R);$('details').textContent='mechanism #'+m.id+'\n'+m.reason+'\n\nframe: '+f.id+'\ncenter: R'+f.center_R+'\nordered neighbors: '+f.neighbors_R_index_order.map(x=>'R'+x).join(' → ')+'\nR sign / normalized volume: '+f.reactant_orientation_sign+' / '+f.reactant_normalized_orientation.toFixed(6)+'\nP sign / normalized volume: '+f.source_product_orientation_sign+' / '+f.source_product_normalized_orientation.toFixed(6)}
$('title').textContent=DATA.step_id+' index-chirality conflict';DATA.rejected_mechanisms.forEach((m,i)=>{const o=document.createElement('option');o.value=i;o.textContent='#'+m.id+' · '+m.frames.length+' source mismatch frame(s)';$('mechanism').appendChild(o)});$('mechanism').onchange=rebuildFrames;$('frame').onchange=render;$('allLabels').onchange=render;rebuildFrames();
</script></body></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, required=True)
    parser.add_argument("--selection-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(build(
        args.case_root.resolve(), args.selection_manifest.resolve(),
        args.output.resolve()))


if __name__ == "__main__":
    main()
