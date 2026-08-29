#!/usr/bin/env python3
"""Build a self-contained 3D viewer for a Suzuki-style coupling demo."""
from __future__ import annotations

import html
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from rxn_core import (
    RetroFragmentSearchConfig,
    assemble_fragment_cover,
    discover_retained_fragments,
)
from rxn_core.smiles import smiles_to_weighted_graph


COLORS = {"R1": "#2684ff", "R2": "#ff8b00"}


def _mol_3d(smiles, *, spread_ions=False, show_hydrogens=False,
            cut_bonds=()):
    molecule = Chem.MolFromSmiles(smiles)
    with_hydrogens = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 20260827
    if AllChem.EmbedMolecule(with_hydrogens, parameters) != 0:
        raise RuntimeError(f"could not embed {smiles!r}")
    try:
        AllChem.UFFOptimizeMolecule(with_hydrogens, maxIters=300)
    except Exception:
        pass
    if cut_bonds:
        editable = Chem.RWMol(with_hydrogens)
        for left, right in cut_bonds:
            if editable.GetBondBetweenAtoms(int(left), int(right)) is not None:
                editable.RemoveBond(int(left), int(right))
        with_hydrogens = editable.GetMol()
        with_hydrogens.UpdatePropertyCache(strict=False)
    molecule = (with_hydrogens if show_hydrogens
                else Chem.RemoveHs(with_hydrogens))
    conformer = molecule.GetConformer()
    if spread_ions:
        fragments = Chem.GetMolFrags(molecule, asMols=False, sanitizeFrags=False)
        cursor = 0.0
        gap = 8.0
        for fragment_index, fragment in enumerate(fragments):
            x_values = [conformer.GetAtomPosition(atom_index).x
                        for atom_index in fragment]
            shift = cursor - min(x_values)
            for atom_index in fragment:
                point = conformer.GetAtomPosition(atom_index)
                point.x += shift
                point.y += 1.5 * (-1) ** fragment_index
                conformer.SetAtomPosition(atom_index, point)
            cursor += max(x_values) - min(x_values) + gap
        center_shift = (cursor - gap) / 2.0
        for atom_index in range(molecule.GetNumAtoms()):
            point = conformer.GetAtomPosition(atom_index)
            point.x -= center_shift
            conformer.SetAtomPosition(atom_index, point)
    coordinates = [
        [
            float(conformer.GetAtomPosition(index).x),
            float(conformer.GetAtomPosition(index).y),
            float(conformer.GetAtomPosition(index).z),
        ]
        for index in range(molecule.GetNumAtoms())
    ]
    return Chem.MolToMolBlock(molecule), coordinates, [
        atom.GetSymbol() for atom in molecule.GetAtoms()
    ]


def _build_payload():
    target_smiles = "Clc1ccccc1-c1ccccc1"
    target_graph = smiles_to_weighted_graph(
        target_smiles, expand_hydrogens=False)
    config = RetroFragmentSearchConfig(
        minimum_fragment_size=3,
        branch_limit=100,
        candidate_limit=100,
        maximum_boundary_bonds=1,
        maximum_leftover_fragments=1,
    )
    definitions = (
        ("R1", "Bromobenzene · MCULE-5539191636", "BrC1=CC=CC=C1"),
        ("R2", "2-Chlorophenylboronic acid · MCULE-6011753091",
         "C1(=CC=CC=C1Cl)B(O)O"),
    )
    results = []
    for key, name, smiles in definitions:
        result = discover_retained_fragments(
            smiles_to_weighted_graph(smiles, expand_hydrogens=False),
            target_graph,
            precursor_id=key,
            config=config,
        )
        results.append(result)
    assembly_result = assemble_fragment_cover(
        target_graph,
        tuple(candidate for result in results for candidate in result.candidates),
        maximum_precursors=2,
        require_attachment_bonds=True,
    )
    if not assembly_result.assemblies:
        raise RuntimeError("Suzuki demo did not produce an assembly")
    assembly = assembly_result.assemblies[0]
    selected = {candidate.precursor_id: candidate for candidate in assembly.candidates}

    models = {}
    for key, name, smiles in definitions:
        block, coords, elements = _mol_3d(smiles, spread_ions="." in smiles)
        candidate = selected[key]
        models[key] = {
            "title": f"{key}: {name}",
            "subtitle": (
                f"matched atoms {list(candidate.retained_atoms)}; "
                f"leftover {list(candidate.leftover_fragments)}"),
            "mol": block,
            "coords": coords,
            "elements": elements,
            "styles": [{
                "indices": list(candidate.retained_atoms),
                "color": COLORS[key],
            }],
            "broken": [list(bond) for bond in candidate.boundary_bonds],
            "formed": [],
        }

    product_block, product_coords, product_elements = _mol_3d(target_smiles)
    product_styles = []
    for candidate in assembly.candidates:
        product_styles.append({
            "indices": list(candidate.covered_target_atoms),
            "color": COLORS[candidate.precursor_id],
        })
    models["P"] = {
        "title": "P target: 2-chlorobiphenyl",
        "subtitle": "complete fragment-unit coverage",
        "mol": product_block,
        "coords": product_coords,
        "elements": product_elements,
        "styles": product_styles,
        "broken": [],
        "formed": [list(bond) for bond in assembly.formed_bonds],
    }
    side_block, side_coords, side_elements = _mol_3d("OB(O)Br")
    models["P2"] = {
        "title": "Reconstructed unmatched pool: B(OH)2Br",
        "subtitle": "reconstructed from unmatched R fragments",
        "mol": side_block,
        "coords": side_coords,
        "elements": side_elements,
        "styles": [],
        "broken": [],
        "formed": [],
    }
    return {
        "models": models,
        "summary": {
            "status": assembly_result.status,
            "branch_cap": 100,
            "formed_bonds": [list(bond) for bond in assembly.formed_bonds],
            "broken_bonds": [list(bond) for bond in assembly.broken_bonds],
            "coverage": {
                candidate.precursor_id: list(candidate.covered_target_atoms)
                for candidate in assembly.candidates
            },
        },
    }


def _html(payload):
    library = (Path(__file__).parents[1] / "src" / "rxn_core" / "static" /
               "3Dmol-min.js").read_text()
    data = json.dumps(payload, separators=(",", ":"))
    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Fragment-unit retrosynthesis demo</title>
<style>
:root{{--bg:#f3f5f8;--card:#fff;--ink:#172033;--muted:#64748b;--line:#dbe2ea}}
*{{box-sizing:border-box}} body{{margin:0;font:14px system-ui;background:var(--bg);color:var(--ink)}}
header{{height:72px;padding:12px 20px;background:#101828;color:white;display:flex;align-items:center;justify-content:space-between}}
h1{{font-size:19px;margin:0 0 4px}} .sub{{color:#b8c4d6}} .controls{{display:flex;gap:16px;align-items:center}}
#grid{{display:grid;grid-template-columns:1fr 1fr;gap:10px;padding:10px;height:calc(100vh - 144px)}}
.panel{{position:relative;min-height:260px;background:var(--card);border:1px solid var(--line);border-radius:10px;overflow:hidden;box-shadow:0 2px 8px #0000000c}}
.view{{position:absolute;inset:0}} .label{{position:absolute;z-index:3;left:12px;top:10px;background:#ffffffdc;border:1px solid var(--line);border-radius:7px;padding:7px 10px;pointer-events:none}}
.label b{{display:block;font-size:14px}} .label small{{color:var(--muted)}}
footer{{height:72px;background:white;border-top:1px solid var(--line);display:flex;align-items:center;gap:22px;padding:10px 20px;overflow:auto}}
.legend{{display:flex;align-items:center;gap:7px;white-space:nowrap}} .dot{{width:13px;height:13px;border-radius:50%}}
.bond{{width:24px;height:5px;border-radius:4px}} code{{background:#eef2f6;padding:3px 6px;border-radius:4px}}
</style><script>{library}</script></head>
<body><header><div><h1>Suzuki-style coupling — fragment-unit reconstruction</h1>
<div class="sub">R1 + R2 → P target + hidden P2 · no cut sweep · branch cap 100</div></div>
<div class="controls"><label><input id="labels" type="checkbox" checked> atom indices</label>
<button id="reset">reset views</button></div></header>
<main id="grid">
<section class="panel"><div class="label"><b>R1: bromobenzene</b><small>blue phenyl retained; Br element-colored</small></div><div id="R1" class="view"></div></section>
<section class="panel"><div class="label"><b>R2: 2-chlorophenylboronic acid</b><small>orange chlorophenyl retained; B(OH)2 element-colored</small></div><div id="R2" class="view"></div></section>
<section class="panel"><div class="label"><b>P target: 2-chlorobiphenyl</b><small>blue + orange coverage; green C–C bond formed</small></div><div id="P" class="view"></div></section>
<section class="panel"><div class="label"><b>Reconstructed unmatched pool</b><small>B(OH)2 + Br shown in element colors</small></div><div id="P2" class="view"></div></section>
</main>
<footer><span class="legend"><i class="dot" style="background:{COLORS['R1']}"></i>R1 retained fragment</span>
<span class="legend"><i class="dot" style="background:{COLORS['R2']}"></i>R2 retained fragment</span>
<span class="legend"><i class="dot" style="background:linear-gradient(90deg,#8c510a,#a855f7)"></i>unmatched: element colors</span>
<span class="legend"><i class="bond" style="background:#e5484d"></i>broken R bond</span>
<span class="legend"><i class="bond" style="background:#16a34a"></i>formed P bond</span>
<span>Result: <code>1 complete assembly</code></span></footer>
<script>
const data={data}; const viewers={{}};
function point(coords,i){{return {{x:coords[i][0],y:coords[i][1],z:coords[i][2]}}}}
function draw(key){{const d=data.models[key],v=$3Dmol.createViewer(key,{{backgroundColor:'#ffffff'}});viewers[key]=v;
 v.addModel(d.mol,'sdf');v.setStyle({{}},{{stick:{{radius:.16}},sphere:{{scale:.30}}}});
 d.styles.forEach(s=>v.addStyle({{index:s.indices}},{{stick:{{color:s.color,radius:.22}},sphere:{{color:s.color,scale:.42}}}}));
 d.broken.forEach(b=>v.addCylinder({{start:point(d.coords,b[0]),end:point(d.coords,b[1]),radius:.10,color:'#e5484d',alpha:.92}}));
 d.formed.forEach(b=>v.addCylinder({{start:point(d.coords,b[0]),end:point(d.coords,b[1]),radius:.11,color:'#16a34a',alpha:.95}}));
 addLabels(key);v.zoomTo();v.render();}}
function addLabels(key){{const v=viewers[key],d=data.models[key];v.removeAllLabels();if(!document.getElementById('labels').checked)return;
 d.coords.forEach((p,i)=>v.addLabel(d.elements[i]+i,{{position:point(d.coords,i),fontSize:10,fontColor:'#111827',backgroundColor:'white',backgroundOpacity:.72,inFront:true}}));}}
['R1','R2','P','P2'].forEach(draw);
document.getElementById('labels').onchange=()=>Object.keys(viewers).forEach(k=>{{addLabels(k);viewers[k].render()}});
document.getElementById('reset').onclick=()=>Object.values(viewers).forEach(v=>{{v.zoomTo();v.render()}});
window.onresize=()=>Object.values(viewers).forEach(v=>v.resize());
</script></body></html>"""


def main():
    output = Path("data/retro_demo/view.html")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(_html(_build_payload()))
    print(output.resolve())


if __name__ == "__main__":
    main()
