#!/usr/bin/env python3
"""Build a selectable 3D viewer for one-precursor rearrangement results."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from rdkit import Chem

from build_retro_demo_viewer import _mol_3d


FRAGMENT_COLORS = [
    "#1565c0", "#ef6c00", "#7b1fa2", "#2e7d32", "#c2185b",
    "#f9a825", "#c62828", "#00838f", "#3949ab", "#558b2f",
    "#6d4c41", "#546e7a", "#ad1457", "#00695c", "#4527a0",
]


def _model(smiles, *, fragments, broken=(), formed=(), spread=False,
           cut_bonds=(), labels=()):
    block, coords, elements = _mol_3d(
        smiles, spread_ions=spread, show_hydrogens=True,
        cut_bonds=cut_bonds)
    fragment_labels = []
    for index, fragment in enumerate(fragments):
        if not fragment:
            continue
        fragment_labels.append({
            "text": f"F{index + 1}",
            "color": FRAGMENT_COLORS[index % len(FRAGMENT_COLORS)],
            "position": [
                sum(coords[atom][axis] for atom in fragment) / len(fragment)
                for axis in range(3)
            ],
        })
    return {
        "mol": block, "coords": coords, "elements": elements,
        "fragments": [sorted(fragment) for fragment in fragments],
        "fragment_labels": fragment_labels,
        "labels": list(labels),
        "broken": [event[:2] for event in broken],
        "formed": [event[:2] for event in formed],
    }


def _retained_fragments(molecule, mapping, broken, formed, target_atoms):
    """Connected R components after every changed R bond is cut."""
    cut_edges = {tuple(sorted(map(int, event[:2]))) for event in broken}
    inverse = {int(image): int(source) for source, image in mapping}
    for event in formed:
        left, right = int(event[0]), int(event[1])
        if left in inverse and right in inverse:
            source_edge = tuple(sorted((inverse[left], inverse[right])))
            if molecule.GetBondBetweenAtoms(*source_edge) is not None:
                cut_edges.add(source_edge)
    graph = Chem.RWMol(molecule)
    for left, right in cut_edges:
        if graph.GetBondBetweenAtoms(left, right) is not None:
            graph.RemoveBond(left, right)
    components = sorted(
        (set(fragment) for fragment in Chem.GetMolFrags(graph.GetMol())),
        key=lambda fragment: (-len(fragment), min(fragment)),
    )
    mapping_dict = {int(source): int(image) for source, image in mapping}
    r_fragments, p_fragments = [], []
    for component in components:
        retained_r = {atom for atom in component
                      if mapping_dict.get(atom, target_atoms) < target_atoms}
        retained_p = {mapping_dict[atom] for atom in retained_r}
        if retained_r:
            r_fragments.append(retained_r)
            p_fragments.append(retained_p)
    return r_fragments, p_fragments, tuple(sorted(cut_edges))


def _payload(report, target_smiles, top_count):
    target_atoms = Chem.AddHs(Chem.MolFromSmiles(target_smiles)).GetNumAtoms()
    known_rank = report["expected_rank"]
    ranks = list(range(1, min(top_count, report["result_count"]) + 1))
    if known_rank not in ranks:
        ranks.append(known_rank)
    output = []
    for rank in ranks:
        item = report["results"][rank - 1]
        molecule = Chem.AddHs(Chem.MolFromSmiles(item["smiles"]))
        r_fragments, p_fragments, cut_edges = _retained_fragments(
            molecule, item["mapping"], item["broken"], item["formed"],
            target_atoms)
        target_formed = [
            event for event in item["formed"]
            if int(event[0]) < target_atoms and int(event[1]) < target_atoms
        ]
        output.append({
            "rank": rank,
            "known": item["precursor_id"] == report["expected_id"],
            "id": item["precursor_id"],
            "smiles": item["smiles"],
            "score": item["score"],
            "complete": item["complete"],
            "cap_hits": item["cap_hits"],
            "excess": item["excess_elements"],
            "models": [
                _model(item["smiles"], fragments=r_fragments,
                       broken=item["broken"], labels=[
                           {"atom": int(source), "text": f"P{int(target)}"}
                           for source, target in item["mapping"]
                           if int(target) < target_atoms
                       ]),
                _model(target_smiles, fragments=p_fragments,
                       formed=target_formed, labels=[
                           {"atom": atom, "text": f"P{atom}"}
                           for atom in range(target_atoms)
                       ]),
            ],
            "fragment_count": len(r_fragments),
            "pattern_key": tuple(sorted(
                tuple(sorted(fragment)) for fragment in p_fragments)),
        })
    pattern_ids = {}
    patterns = []
    for item in output:
        key = item.pop("pattern_key")
        if key not in pattern_ids:
            pattern_id = len(pattern_ids) + 1
            pattern_ids[key] = pattern_id
            patterns.append({
                "pattern": pattern_id,
                "coverage_atom_sets": [list(fragment) for fragment in key],
                "fragment_sizes": sorted(
                    (len(fragment) for fragment in key), reverse=True),
            })
        item["pattern"] = pattern_ids[key]
    output.sort(key=lambda item: (item["pattern"], item["rank"]))
    return {
        "summary": {
            "rows": report["scan_counts"]["rows"],
            "searched": report["scan_counts"]["searched"],
            "results": report["result_count"],
            "cap_hits": report["scan_counts"]["cap_hits"],
            "known_rank": known_rank,
            "pattern_count": len(patterns),
        },
        "patterns": patterns,
        "results": output,
    }


def _html(payload, *, title, known_label):
    library = (Path(__file__).parents[1] / "src" / "rxn_core" / "static" /
               "3Dmol-min.js").read_text()
    data = json.dumps(payload, separators=(",", ":"))
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title><style>
:root{{--bg:#f3f5f8;--card:#fff;--ink:#172033;--muted:#64748b;--line:#dbe2ea}}
*{{box-sizing:border-box}}body{{margin:0;height:100vh;overflow:hidden;font:13px system-ui;background:var(--bg);color:var(--ink)}}
header{{height:72px;padding:11px 18px;background:#101828;color:white;display:flex;align-items:center;justify-content:space-between}}h1{{font-size:18px;margin:0 0 4px}}.muted{{color:#9fb0c8}}.metrics{{display:flex;gap:20px}}.metric b{{display:block;font-size:17px}}.metric small{{color:#9fb0c8}}
#layout{{display:grid;grid-template-columns:330px 1fr;height:calc(100vh - 72px)}}aside{{background:white;border-right:1px solid var(--line);overflow:auto}}.intro{{padding:12px;border-bottom:1px solid var(--line);line-height:1.5}}
.result{{padding:10px 12px;border:0;border-bottom:1px solid var(--line);width:100%;text-align:left;background:white;cursor:pointer}}.result:hover{{background:#f5f8fc}}.result.active{{background:#eaf3ff;box-shadow:inset 4px 0 #2684ff}}.rank{{font-weight:750;font-size:14px}}.badge{{background:#0f9d66;color:white;border-radius:10px;padding:2px 7px;margin-left:7px;font-size:10px}}.ids{{color:#475569;margin-top:4px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}.score{{color:#64748b;font-size:11px;margin-top:4px}}
.patternhead{{padding:9px 12px;background:#e2e8f0;border-top:2px solid #94a3b8;border-bottom:1px solid var(--line);font-weight:800}}.patternhead small{{display:block;color:#475569;font-weight:500;margin-top:2px}}
main{{display:grid;grid-template-columns:1fr 1fr;gap:8px;padding:8px;min-width:0}}.panel{{position:relative;background:white;border:1px solid var(--line);border-radius:9px;overflow:hidden}}.view{{position:absolute;inset:0}}.label{{position:absolute;z-index:3;left:10px;top:9px;max-width:85%;background:#ffffffdf;border:1px solid var(--line);border-radius:7px;padding:7px 9px;pointer-events:none}}.label b,.label small{{display:block}}.label small{{color:#64748b;margin-top:2px;overflow:hidden;text-overflow:ellipsis}}.controls{{position:absolute;z-index:5;right:18px;top:82px;background:#ffffffdc;padding:6px 9px;border:1px solid var(--line);border-radius:7px;line-height:1.7}}.controls label{{display:block}}.dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px}}
</style><script>{library}</script></head><body><header><div><h1>{title}</h1><div class="muted">Explicit-H AAM · select a precursor to inspect its event mapping</div></div><div class="metrics"><span class="metric"><b>{payload['summary']['rows']:,}</b><small>catalog R</small></span><span class="metric"><b>{payload['summary']['searched']:,}</b><small>AAM searched</small></span><span class="metric"><b>{payload['summary']['results']:,}</b><small>low-event results</small></span></div></header>
<div id="layout"><aside><div class="intro"><b>Distinct colors = connected retained fragments.</b><br>The same fragment has the same color in R and P.<br><span style="color:#d33">Red = broken</span>; <span style="color:#159447">green = formed</span>.<br>Turn fragment colors off to restore original element colors.<br><br>Known {known_label} recovered at <b>rank {payload['summary']['known_rank']}</b>.<br>All cap events across scan: {payload['summary']['cap_hits']:,}.</div><div id="list"></div></aside><main><div class="controls"><label><input id="fragments" type="checkbox" checked> color retained fragments</label><label><input id="labels" type="checkbox"> mapped P# identities</label></div><section class="panel"><div class="label" id="LR"></div><div class="view" id="R"></div></section><section class="panel"><div class="label" id="LP"></div><div class="view" id="P"></div></section></main></div>
<script>const data={data},fragmentColors={json.dumps(FRAGMENT_COLORS)},viewers={{}};function pt(m,i){{return{{x:m.coords[i][0],y:m.coords[i][1],z:m.coords[i][2]}}}}function point(p){{return{{x:p[0],y:p[1],z:p[2]}}}}function show(id,m){{let v=viewers[id];if(!v){{v=$3Dmol.createViewer(id,{{backgroundColor:'white'}});viewers[id]=v}}else{{v.removeAllModels();v.removeAllShapes();v.removeAllLabels()}}v.addModel(m.mol,'sdf');v.setStyle({{}},{{stick:{{radius:.15}},sphere:{{scale:.28}}}});if(document.getElementById('fragments').checked){{m.fragments.forEach((atoms,i)=>{{const color=fragmentColors[i%fragmentColors.length];v.addStyle({{index:atoms}},{{stick:{{color:color,radius:.21}},sphere:{{color:color,scale:.39}}}})}});m.fragment_labels.forEach(l=>v.addLabel(l.text,{{position:point(l.position),fontSize:12,fontColor:'white',backgroundColor:l.color,backgroundOpacity:.9,inFront:true}}))}}m.broken.forEach(b=>v.addCylinder({{start:pt(m,b[0]),end:pt(m,b[1]),radius:.10,color:'#e5484d'}}));m.formed.forEach(b=>v.addCylinder({{start:pt(m,b[0]),end:pt(m,b[1]),radius:.11,color:'#16a34a'}}));if(document.getElementById('labels').checked)(m.labels||[]).forEach(l=>v.addLabel(l.text,{{position:pt(m,l.atom),fontSize:10,fontColor:'#111',backgroundColor:'white',backgroundOpacity:.72,inFront:true}}));v.zoomTo();v.render()}}function patternInfo(id){{return data.patterns.find(x=>x.pattern===id)}}function patternText(id){{const p=patternInfo(id);return p?(p.fragment_sizes.length+' modules · atom sizes '+p.fragment_sizes.join(' + ')):''}}
function select(i){{const x=data.results[i];document.querySelectorAll('.result').forEach(e=>e.classList.toggle('active',Number(e.dataset.index)===i));document.getElementById('LR').innerHTML='<b>R · '+x.id+(x.known?' · known {known_label}':'')+'</b><small>'+x.smiles+'</small><small>Pattern '+x.pattern+' · '+x.fragment_count+' retained fragment(s) · '+x.score[0]+' bond edits · '+x.score[1]+' chirality violations · cap hits '+x.cap_hits+'</small>';document.getElementById('LP').innerHTML='<b>P target only · Pattern '+x.pattern+'</b><small>'+patternText(x.pattern)+' · colored regions define the construction pattern</small>';show('R',x.models[0]);show('P',x.models[1])}}const list=document.getElementById('list');let lastPattern=null;data.results.forEach((x,i)=>{{if(x.pattern!==lastPattern){{const h=document.createElement('div');h.className='patternhead';h.innerHTML='Pattern '+x.pattern+'<small>'+patternText(x.pattern)+' · colored regions on P define this construction</small>';list.appendChild(h);lastPattern=x.pattern}}const b=document.createElement('button');b.className='result';b.dataset.index=i;b.innerHTML='<span class="rank">#'+x.rank+'</span>'+(x.known?'<span class="badge">KNOWN {known_label.upper()}</span>':'')+'<div class="ids">'+x.id+'</div><div class="score">fragments '+x.fragment_count+' · edits '+x.score[0]+' · chirality '+x.score[1]+' · '+(x.complete?'complete':'CAPPED')+'</div>';b.onclick=()=>select(i);list.appendChild(b)}});function redraw(){{const active=[...document.querySelectorAll('.result')].find(x=>x.classList.contains('active'));select(active?Number(active.dataset.index):0)}}['labels','fragments'].forEach(id=>document.getElementById(id).onchange=redraw);const knownIndex=data.results.findIndex(x=>x.known);select(Math.max(0,knownIndex));window.onresize=()=>Object.values(viewers).forEach(v=>v.resize());</script></body></html>"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top", type=int, default=24)
    parser.add_argument("--title", default="Blind rearrangement search")
    parser.add_argument("--known-label", default="precursor")
    parser.add_argument(
        "--replacement-result",
        help="Full mapping JSON used to replace the known precursor result.")
    args = parser.parse_args()
    report = json.loads(Path(args.results).read_text())
    if args.replacement_result:
        replacement = json.loads(Path(args.replacement_result).read_text())
        known = report["results"][report["expected_rank"] - 1]
        for key in ("mapping", "broken", "formed", "score",
                    "chirality_violations"):
            known[key] = replacement[key]
        known["raw_mapping"] = replacement["mapping"]
        known["complete"] = False
        known["cap_hits"] = max(1, int(known.get("cap_hits", 0)))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_html(_payload(report, args.target_smiles, args.top),
                            title=args.title,
                            known_label=args.known_label))
    print(output.resolve())


if __name__ == "__main__":
    main()
