#!/usr/bin/env python3
"""Build a detection-only viewer from persisted fragment mappings."""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
from collections import Counter
from pathlib import Path

from rdkit import Chem

from build_retro_demo_viewer import _mol_3d
from rxn_core.retrosynthesis.catalog_index import chirality_violations


COLORS = (
    "#2684ff", "#ff8b00", "#8b5cf6", "#00a896",
    "#e83e8c", "#795548", "#00acc1",
)


def _combine_source_models(sources):
    combined = None
    for index, source in enumerate(sources):
        molecule = Chem.MolFromMolBlock(
            source["model"]["mol"], removeHs=False)
        conformer = molecule.GetConformer()
        positions = [conformer.GetAtomPosition(atom)
                     for atom in range(molecule.GetNumAtoms())]
        center = (
            sum(point.x for point in positions) / len(positions),
            sum(point.y for point in positions) / len(positions),
            sum(point.z for point in positions) / len(positions),
        )
        destination = ((index % 2) * 32.0, -(index // 2) * 24.0, 0.0)
        for atom, point in enumerate(positions):
            conformer.SetAtomPosition(atom, (
                point.x + destination[0] - center[0],
                point.y + destination[1] - center[1],
                point.z + destination[2] - center[2],
            ))
        source["combined_offset"] = (
            0 if combined is None else combined.GetNumAtoms())
        combined = molecule if combined is None else Chem.CombineMols(
            combined, molecule)

    conformer = combined.GetConformer()
    return {
        "mol": Chem.MolToMolBlock(combined),
        "coords": [
            list(conformer.GetAtomPosition(atom))
            for atom in range(combined.GetNumAtoms())
        ],
        "elements": [atom.GetSymbol() for atom in combined.GetAtoms()],
    }


def _validates_mapping(source, target, candidate, tolerance):
    mapping = {int(left): int(right) for left, right in candidate["mapping"]}
    retained = set(map(int, candidate["retained_atoms"]))
    elements_match = all(
        source.GetAtomWithIdx(left).GetAtomicNum()
        == target.GetAtomWithIdx(right).GetAtomicNum()
        for left, right in mapping.items()
    )
    bonds_match = True
    for bond in source.GetBonds():
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if left not in retained or right not in retained:
            continue
        target_bond = target.GetBondBetweenAtoms(mapping[left], mapping[right])
        if (target_bond is None
                or abs(bond.GetBondTypeAsDouble()
                       - target_bond.GetBondTypeAsDouble()) > tolerance):
            bonds_match = False
            break
    return elements_match, bonds_match


def _payload(records_path, bank_path, target_smiles, tolerance):
    names = {}
    with open(bank_path, newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            names[str(row["Literature ID"])] = row["Name"]

    target = Chem.AddHs(Chem.MolFromSmiles(target_smiles))
    Chem.AssignStereochemistry(target, cleanIt=True, force=True)
    target_block, target_coords, target_elements = _mol_3d(
        target_smiles, show_hydrogens=True)
    with gzip.open(records_path, "rt", encoding="utf-8") as stream:
        records = sorted(
            (json.loads(line) for line in stream),
            key=lambda record: int(record["source_id"]),
        )

    sources = []
    for source_index, record in enumerate(records):
        source = Chem.AddHs(Chem.MolFromSmiles(record["representation"]))
        Chem.AssignStereochemistry(source, cleanIt=True, force=True)
        source_block, source_coords, source_elements = _mol_3d(
            record["representation"], show_hydrogens=True)
        candidates = []
        for candidate_index, candidate in enumerate(record["candidates"]):
            elements_match, bonds_match = _validates_mapping(
                source, target, candidate, tolerance)
            chiral = chirality_violations(candidate, source, target)
            candidates.append({
                "index": candidate_index,
                "mapping": candidate["mapping"],
                "retained_atoms": candidate["retained_atoms"],
                "covered_target_atoms": candidate["covered_target_atoms"],
                "retained_fragments": candidate["retained_fragments"],
                "boundary_bonds": candidate["boundary_bonds"],
                "element_match": elements_match,
                "bond_match": bonds_match,
                "chirality_violations": chiral,
            })
        default = min(
            range(len(candidates)),
            key=lambda index: (
                candidates[index]["chirality_violations"],
                -len(candidates[index]["retained_atoms"]),
                candidates[index]["covered_target_atoms"],
            ),
        )
        sources.append({
            "id": record["source_id"],
            "name": names.get(record["source_id"], record["source_id"]),
            "smiles": record["representation"],
            "color": COLORS[source_index % len(COLORS)],
            "complete": record["complete"],
            "maximum_branch_count": record["maximum_branch_count"],
            "model": {
                "mol": source_block,
                "coords": source_coords,
                "elements": source_elements,
            },
            "candidates": candidates,
            "default_candidate": default,
        })

    def placement_set_rank(choices):
        counts = Counter(
            atom
            for source, choice in zip(sources, choices)
            for atom in source["candidates"][choice]["covered_target_atoms"]
        )
        return (
            len(counts),
            -sum(count - 1 for count in counts.values() if count > 1),
            -sum(
                source["candidates"][choice]["chirality_violations"]
                for source, choice in zip(sources, choices)
            ),
            tuple(-choice for choice in choices),
        )

    default_choices = max(
        itertools.product(*(
            range(len(source["candidates"])) for source in sources
        )),
        key=placement_set_rank,
    )
    for source, choice in zip(sources, default_choices):
        source["default_candidate"] = choice
    combined_sources = _combine_source_models(sources)
    defaults = [
        (source_index, source["candidates"][source["default_candidate"]])
        for source_index, source in enumerate(sources)
    ]

    occupancy = Counter()
    heavy_occupancy = set()
    for _source_index, candidate in defaults:
        occupancy.update(map(int, candidate["covered_target_atoms"]))
    for atom in occupancy:
        if target.GetAtomWithIdx(atom).GetAtomicNum() > 1:
            heavy_occupancy.add(atom)
    heavy_count = target.GetNumHeavyAtoms()
    return {
        "target": {
            "smiles": target_smiles,
            "model": {
                "mol": target_block,
                "coords": target_coords,
                "elements": target_elements,
            },
            "atom_count": target.GetNumAtoms(),
            "heavy_atom_count": heavy_count,
        },
        "sources": sources,
        "combined_sources": combined_sources,
        "summary": {
            "source_count": len(sources),
            "candidate_count": sum(len(source["candidates"])
                                   for source in sources),
            "cap_hits": sum(not source["complete"] for source in sources),
            "default_union_atoms": len(occupancy),
            "default_union_heavy_atoms": len(heavy_occupancy),
            "default_chirality_violations": sum(
                candidate["chirality_violations"]
                for _source_index, candidate in defaults
            ),
            "overlap_atoms": sorted(
                atom for atom, count in occupancy.items() if count > 1),
            "uncovered_atoms": sorted(
                set(range(target.GetNumAtoms())) - set(occupancy)),
            "tolerance": tolerance,
        },
    }


def _html(payload):
    library = (Path(__file__).parents[1] / "src" / "rxn_core" / "static" /
               "3Dmol-min.js").read_text()
    data = json.dumps(payload, separators=(",", ":"))
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>Literature component occupations</title><style>
*{{box-sizing:border-box}}body{{margin:0;height:100vh;overflow:hidden;font:13px system-ui;color:#172033;background:#f3f5f8}}
header{{height:88px;background:#101828;color:white;padding:12px 18px;display:flex;justify-content:space-between;align-items:center}}
h1{{font-size:19px;margin:0 0 5px}}.sub{{color:#b8c4d8}}.warning{{color:#ffd38a;font-weight:700;max-width:590px;text-align:right}}
#layout{{height:calc(100vh - 88px);display:grid;grid-template-columns:360px 1fr}}aside{{background:white;border-right:1px solid #dbe2ea;overflow:auto}}
.summary{{padding:12px;line-height:1.45;border-bottom:1px solid #dbe2ea}}button{{font:inherit}}.all,.source,.candidate{{width:100%;border:0;border-bottom:1px solid #e4e9ef;background:white;text-align:left;padding:10px 12px;cursor:pointer}}
.all:hover,.source:hover,.candidate:hover{{background:#f3f7fc}}.source.active,.all.active{{background:#eaf3ff;box-shadow:inset 4px 0 #2684ff}}.source b{{font-size:14px}}.small{{font-size:11px;color:#64748b;margin-top:3px}}
#candidates{{padding:8px 12px;background:#f8fafc;border-bottom:1px solid #dbe2ea}}.candidate{{display:inline-block;width:auto;border:1px solid #cbd5e1;border-radius:6px;padding:5px 8px;margin:3px}}.candidate.active{{background:#172033;color:white}}
main{{min-width:0;display:grid;grid-template-rows:auto 1fr}}#info{{background:white;padding:9px 14px;border-bottom:1px solid #dbe2ea;line-height:1.4}}
#views{{display:grid;grid-template-columns:1fr 1.5fr;gap:8px;padding:8px;min-height:0}}.panel{{position:relative;background:white;border:1px solid #dbe2ea;border-radius:9px;overflow:hidden}}.view{{position:absolute;inset:0}}
.label{{position:absolute;z-index:3;left:10px;top:10px;background:#ffffffdf;border:1px solid #dbe2ea;border-radius:7px;padding:6px 9px;max-width:78%}}.label small{{display:block;color:#64748b;margin-top:2px}}
.controls{{position:absolute;z-index:5;right:18px;top:99px;background:#ffffffdf;border:1px solid #dbe2ea;border-radius:7px;padding:5px 8px}}
</style><script>{library}</script></head><body>
<header><div><h1>Seven literature components → final vancomycin aglycon</h1><div class="sub">Saved blind explicit-H AAM results · tolerance 0.5 · intact detected fragments</div></div><div class="warning">DETECTION-ONLY VIEW: occupations may overlap and do not constitute a one-step assembly.</div></header>
<div id="layout"><aside><div class="summary"><b>All seven produced structural matches.</b><br><span id="summary"></span><br>Choose an R to inspect every equal-size placement. The combined view maximizes union coverage; it is not a ground-truth assignment.</div><button class="all" id="all">Coverage-maximizing occupations</button><div id="sources"></div></aside>
<main><div id="info"></div><div class="controls"><label><input id="labels" type="checkbox"> atom-map labels</label></div><div id="views"><section class="panel"><div class="label" id="RL"></div><div class="view" id="R"></div></section><section class="panel"><div class="label" id="PL"></div><div class="view" id="P"></div></section></div></main></div>
<script>const data={data};let selected=-1,candidate=0;const viewers={{}};
function pt(m,i){{return {{x:m.coords[i][0],y:m.coords[i][1],z:m.coords[i][2]}}}}
function draw(id,m,styles,labels,broken=[],titles=[]){{let v=viewers[id];if(!v){{v=$3Dmol.createViewer(id,{{backgroundColor:'white'}});viewers[id]=v}}else{{v.removeAllModels();v.removeAllShapes();v.removeAllLabels()}}if(m.mol){{v.addModel(m.mol,'sdf');v.setStyle({{}},{{stick:{{radius:.12}},sphere:{{scale:.23}}}});styles.forEach(s=>v.addStyle({{index:s.indices}},{{stick:{{color:s.color,radius:.20}},sphere:{{color:s.color,scale:.35}}}}));broken.forEach(b=>v.addCylinder({{start:pt(m,b[0]),end:pt(m,b[1]),radius:.10,color:'#e5484d'}}));titles.forEach(l=>v.addLabel(l.text,{{position:pt(m,l.atom),fontSize:14,fontColor:l.color,backgroundColor:'white',backgroundOpacity:.85,inFront:true}}));if(document.getElementById('labels').checked)labels.forEach(l=>v.addLabel(l.text,{{position:pt(m,l.atom),fontSize:10,fontColor:'#111',backgroundColor:'white',backgroundOpacity:.75,inFront:true}}));v.zoomTo()}}v.render()}}
function showSource(index,choice){{selected=index;candidate=choice;document.querySelectorAll('.source').forEach((b,i)=>b.classList.toggle('active',i===index));document.getElementById('all').classList.remove('active');const s=data.sources[index],c=s.candidates[choice],color=s.color;document.getElementById('candidates')?.remove();const box=document.createElement('div');box.id='candidates';s.candidates.forEach((x,i)=>{{const b=document.createElement('button');b.className='candidate'+(i===choice?' active':'');b.textContent='placement '+(i+1)+(x.chirality_violations?' · stereo ⚠':' · stereo ✓');b.onclick=()=>showSource(index,i);box.appendChild(b)}});document.querySelectorAll('.source')[index].after(box);document.getElementById('info').innerHTML='<b>'+s.name+'</b> · placement '+(choice+1)+'/'+s.candidates.length+' · '+c.retained_atoms.length+' explicit atoms mapped · '+c.retained_fragments.length+' intact connected fragment · element '+(c.element_match?'✓':'✗')+' · bonds '+(c.bond_match?'✓':'✗')+' · chirality violations '+c.chirality_violations;document.getElementById('RL').innerHTML='<b>R'+s.id+' · '+s.name+'</b><small>colored = intact matched fragment; original element colors = unmatched/protecting groups</small>';document.getElementById('PL').innerHTML='<b>P_target occupation by R'+s.id+'</b><small>'+c.covered_target_atoms.length+' of '+data.target.atom_count+' explicit target atoms</small>';draw('R',s.model,[{{indices:c.retained_atoms,color}}],c.mapping.map(x=>({{atom:x[0],text:'P'+x[1]}})),c.boundary_bonds);draw('P',data.target.model,[{{indices:c.covered_target_atoms,color}}],c.mapping.map(x=>({{atom:x[1],text:'R'+x[0]}})))}}
function showAll(){{selected=-1;document.querySelectorAll('.source').forEach(b=>b.classList.remove('active'));document.getElementById('all').classList.add('active');document.getElementById('candidates')?.remove();const owners={{}};data.sources.forEach((s,i)=>s.candidates[s.default_candidate].covered_target_atoms.forEach(a=>(owners[a]??=[]).push(i)));const styles=[];data.sources.forEach((s,i)=>styles.push({{indices:Object.keys(owners).filter(a=>owners[a].length===1&&owners[a][0]===i).map(Number),color:s.color}}));styles.push({{indices:Object.keys(owners).filter(a=>owners[a].length>1).map(Number),color:'#d62828'}});const rstyles=[],rlabels=[],rbroken=[],rtitles=[];data.sources.forEach(s=>{{const c=s.candidates[s.default_candidate],o=s.combined_offset;rstyles.push({{indices:c.retained_atoms.map(a=>a+o),color:s.color}});rtitles.push({{atom:o,text:'R'+s.id,color:s.color}});c.mapping.forEach(x=>rlabels.push({{atom:x[0]+o,text:'R'+s.id+'→P'+x[1]}}));c.boundary_bonds.forEach(b=>rbroken.push([b[0]+o,b[1]+o]))}});document.getElementById('info').innerHTML='<b>Coverage-maximizing combination, one placement per R.</b> The left panel contains every intact R and the right panel shows their occupations on P. Red P atoms are claimed by more than one R; uncolored P atoms are not occupied. This is not a verified literature assignment.';document.getElementById('RL').innerHTML='<b>All seven R structures</b><small>spatially separated · colored = selected intact mapped fragment · original colors = unmatched groups</small>';document.getElementById('PL').innerHTML='<b>P_target combined occupation</b><small>'+data.summary.default_union_atoms+'/'+data.target.atom_count+' explicit atoms · '+data.summary.default_union_heavy_atoms+'/'+data.target.heavy_atom_count+' heavy atoms · '+data.summary.overlap_atoms.length+' overlap atoms · '+data.summary.default_chirality_violations+' chirality conflicts</small>';draw('R',data.combined_sources,rstyles,rlabels,rbroken,rtitles);draw('P',data.target.model,styles,[])}}
const list=document.getElementById('sources');data.sources.forEach((s,i)=>{{const b=document.createElement('button');b.className='source';b.innerHTML='<b><span style="color:'+s.color+'">R'+s.id+'</span> · '+s.name+'</b><div class="small">'+s.candidates.length+' placement'+(s.candidates.length===1?'':'s')+' · best '+s.candidates[s.default_candidate].retained_atoms.length+' atoms · branch max '+s.maximum_branch_count+'</div>';b.onclick=()=>showSource(i,s.default_candidate);list.appendChild(b)}});document.getElementById('summary').textContent=data.summary.candidate_count+' saved placements; '+data.summary.cap_hits+' cap hits; default union '+data.summary.default_union_atoms+'/'+data.target.atom_count+' explicit atoms.';document.getElementById('all').onclick=showAll;document.getElementById('labels').onchange=()=>selected<0?showAll():showSource(selected,candidate);showAll();window.onresize=()=>Object.values(viewers).forEach(v=>v.resize());</script></body></html>"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--records", required=True)
    parser.add_argument("--bank", required=True)
    parser.add_argument("--target-smiles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tolerance", type=float, default=0.5)
    args = parser.parse_args()
    payload = _payload(
        args.records, args.bank, args.target_smiles, args.tolerance)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_html(payload), encoding="utf-8")
    print(json.dumps({
        "output": str(output.resolve()),
        **payload["summary"],
    }, indent=2))


if __name__ == "__main__":
    main()
