#!/usr/bin/env python3
"""Build a standalone viewer for ranked per-precursor AAM mappings."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from rdkit import Chem

from molecule_3d import mol_3d


COLORS = {"mapped": "#1565c0", "leaving": "#c62828",
          "missing": "#ef6c00"}


def _model(smiles, mapped, *, target=False):
    block, coords, elements = mol_3d(smiles, show_hydrogens=True)
    size = len(elements)
    mapped = set(map(int, mapped))
    classes = {
        "mapped": sorted(mapped),
        ("missing" if target else "leaving"):
            sorted(set(range(size)) - mapped),
    }
    return {
        "mol": block, "coords": coords, "elements": elements,
        "styles": classes,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parts", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--json-output")
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--title", default="Precursor AAM candidates")
    args = parser.parse_args()

    target = Chem.AddHs(Chem.MolFromSmiles(args.target_smiles))
    candidates = []
    for path in sorted(Path(args.parts).glob("part_*.jsonl.gz")):
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            for line in stream:
                record = json.loads(line)
                source = Chem.AddHs(Chem.MolFromSmiles(record["smiles"]))
                for candidate in record["candidates"]:
                    retained = len(candidate["retained_atoms"])
                    covered = len(candidate["covered_target_atoms"])
                    candidates.append({
                        "id": record["precursor_id"],
                        "smiles": record["smiles"],
                        "complete": record["complete"],
                        "retained": retained,
                        "source_atoms": source.GetNumAtoms(),
                        "coverage": covered,
                        "target_atoms": target.GetNumAtoms(),
                        "retention_fraction": retained / source.GetNumAtoms(),
                        "coverage_fraction": covered / target.GetNumAtoms(),
                        "boundary_bonds": candidate["boundary_bonds"],
                        "leftover_fragments": candidate["leftover_fragments"],
                        "mapping": candidate["mapping"],
                    })
    candidates.sort(key=lambda item: (
        -item["retention_fraction"], -item["coverage_fraction"],
        len(item["boundary_bonds"]),
        sum(map(len, item["leftover_fragments"])), item["id"],
    ))
    candidates = candidates[:args.top]
    for item in candidates:
        source_atoms = [source for source, _ in item["mapping"]]
        target_atoms = [target_atom for _, target_atom in item["mapping"]]
        item["models"] = [
            _model(item["smiles"], source_atoms),
            _model(args.target_smiles, target_atoms, target=True),
        ]
        item["source_labels"] = [
            {"atom": int(source), "text": f"P{int(target_atom)}"}
            for source, target_atom in item["mapping"]
        ]
        item["target_labels"] = [
            {"atom": atom, "text": f"P{atom}"}
            for atom in range(target.GetNumAtoms())
        ]
    payload = {
        "title": args.title,
        "target_smiles": args.target_smiles,
        "target_atoms": target.GetNumAtoms(),
        "ranking": "all-atom retention, then target coverage",
        "results": candidates,
    }
    if args.json_output:
        json_output = Path(args.json_output)
        json_output.parent.mkdir(parents=True, exist_ok=True)
        json_output.write_text(json.dumps(payload, indent=2) + "\n")

    library = (Path(__file__).parents[1] / "src" / "rxn_core" / "static" /
               "3Dmol-min.js").read_text()
    data = json.dumps(payload, separators=(",", ":"))
    colors = json.dumps(COLORS)
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>{args.title}</title><style>
body{{margin:0;height:100vh;overflow:hidden;font:13px system-ui;background:#eef2f6;color:#172033}}header{{height:78px;padding:12px 18px;background:#101828;color:white}}h1{{margin:0 0 5px;font-size:19px}}.legend{{color:#b7c4d8}}.dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin:0 5px 0 12px}}#layout{{height:calc(100vh - 78px);display:grid;grid-template-columns:330px 1fr}}aside{{overflow:auto;background:white;border-right:1px solid #d9e1ea}}button{{display:block;width:100%;padding:11px;text-align:left;border:0;border-bottom:1px solid #d9e1ea;background:white;cursor:pointer}}button.active{{background:#eaf3ff;box-shadow:inset 4px 0 #2684ff}}button b,button small{{display:block}}button small{{color:#64748b;margin-top:3px}}main{{display:grid;grid-template-columns:1fr 1fr;gap:8px;padding:8px}}.panel{{position:relative;background:white;border:1px solid #d9e1ea;border-radius:9px;overflow:hidden}}.view{{position:absolute;inset:0}}.label{{position:absolute;z-index:4;left:10px;top:9px;right:10px;background:#ffffffdf;border:1px solid #d9e1ea;border-radius:7px;padding:7px 9px}}.label b,.label small{{display:block}}.label small{{color:#64748b;margin-top:2px}}.controls{{position:absolute;z-index:10;right:18px;top:90px;background:#ffffffed;padding:7px 9px;border-radius:7px}}</style><script>{library}</script></head><body><header><h1>{args.title}</h1><div class="legend"><b>Colors and P# labels come from AAM.</b><span class="dot" style="background:#1565c0"></span>retained/mapped <span class="dot" style="background:#c62828"></span>R leaving <span class="dot" style="background:#ef6c00"></span>P unresolved · same P# means the same mapped atom</div></header><div id="layout"><aside id="list"></aside><main><div class="controls"><label><input id="labels" type="checkbox"> show P# atom identities</label></div><section class="panel"><div class="label" id="LR"></div><div class="view" id="R"></div></section><section class="panel"><div class="label" id="LP"></div><div class="view" id="P"></div></section></main></div><script>
const payload={data},mappingColors={colors},mappingViewers={{}},mappingCurrent={{}};
function pt(m,i){{return{{x:m.coords[i][0],y:m.coords[i][1],z:m.coords[i][2]}}}}
function show(id,m,labels){{let v=mappingViewers[id];if(!v){{v=$3Dmol.createViewer(id,{{backgroundColor:'white'}});mappingViewers[id]=v}}else{{v.removeAllModels();v.removeAllLabels()}}v.addModel(m.mol,'sdf');v.setStyle({{}},{{stick:{{radius:.12}},sphere:{{scale:.23}}}});Object.entries(m.styles).forEach(([name,atoms])=>v.addStyle({{index:atoms}},{{stick:{{color:mappingColors[name],radius:.19}},sphere:{{color:mappingColors[name],scale:.34}}}}));if(document.getElementById('labels').checked)labels.filter(l=>m.elements[l.atom]!=='H').forEach(l=>v.addLabel(l.text,{{position:pt(m,l.atom),fontSize:9,fontColor:'#111',backgroundColor:'white',backgroundOpacity:.72,inFront:true}}));v.zoomTo();v.render()}}
function select(i){{mappingCurrent.i=i;const x=payload.results[i];document.querySelectorAll('button').forEach((b,j)=>b.classList.toggle('active',i===j));document.getElementById('LR').innerHTML='<b>R · '+x.id+'</b><small>'+x.smiles+'</small><small>retained '+x.retained+'/'+x.source_atoms+' atoms</small>';document.getElementById('LP').innerHTML='<b>P · phenylacetic acid</b><small>covered '+x.coverage+'/'+x.target_atoms+' atoms; orange atoms require the next route expansion</small>';show('R',x.models[0],x.source_labels);show('P',x.models[1],x.target_labels)}}
const list=document.getElementById('list');payload.results.forEach((x,i)=>{{const b=document.createElement('button');b.innerHTML='<b>#'+(i+1)+' · '+x.id+'</b><small>retention '+(100*x.retention_fraction).toFixed(1)+'% · target coverage '+(100*x.coverage_fraction).toFixed(1)+'%</small><small>'+x.smiles+'</small>';b.onclick=()=>select(i);list.appendChild(b)}});document.getElementById('labels').onchange=()=>select(mappingCurrent.i||0);if(payload.results.length)select(0);
</script></body></html>"""
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html)
    print(output.resolve())


if __name__ == "__main__":
    main()
