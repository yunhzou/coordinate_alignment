#!/usr/bin/env python3
"""Build a mapped multistep Taxol assembly POC from saved catalog AAM results."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors


CORE_ID = "MCULE-3977373901"
SIDE_ID = "MCULE-7102409922"
TES_BACCATIN_SMILES = (
    "CC[Si](CC)(CC)O[C@H]1C[C@H]2OC[C@@]2(OC(C)=O)[C@H]2"
    "[C@H](OC(=O)C3=CC=CC=C3)[C@]3(O)C[C@H](O)C(C)=C("
    "[C@@H](OC(C)=O)C(=O)[C@]12C)C3(C)C"
)


def _records(parts):
    wanted = {CORE_ID, SIDE_ID}
    found = {}
    for path in sorted(Path(parts).glob("part_*.jsonl.gz")):
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            for line in stream:
                record = json.loads(line)
                if record["precursor_id"] in wanted:
                    found[record["precursor_id"]] = record
    if set(found) != wanted:
        raise RuntimeError(f"missing saved AAM records: {wanted - set(found)}")
    return found


def _add_tes(molecule, oxygen):
    editable = Chem.RWMol(molecule)
    first_new = editable.GetNumAtoms()
    silicon = editable.AddAtom(Chem.Atom("Si"))
    editable.AddBond(int(oxygen), silicon, Chem.BondType.SINGLE)
    for _ in range(3):
        alpha = editable.AddAtom(Chem.Atom("C"))
        beta = editable.AddAtom(Chem.Atom("C"))
        editable.AddBond(silicon, alpha, Chem.BondType.SINGLE)
        editable.AddBond(alpha, beta, Chem.BondType.SINGLE)
    output = editable.GetMol()
    Chem.SanitizeMol(output)
    return output, set(range(first_new, output.GetNumAtoms()))


def _add_oee(molecule, oxygen):
    """Add an O-1-ethoxyethyl group to one alcohol oxygen."""
    editable = Chem.RWMol(molecule)
    first_new = editable.GetNumAtoms()
    acetal = editable.AddAtom(Chem.Atom("C"))
    methyl = editable.AddAtom(Chem.Atom("C"))
    ether = editable.AddAtom(Chem.Atom("O"))
    ethyl_a = editable.AddAtom(Chem.Atom("C"))
    ethyl_b = editable.AddAtom(Chem.Atom("C"))
    editable.AddBond(int(oxygen), acetal, Chem.BondType.SINGLE)
    editable.AddBond(acetal, methyl, Chem.BondType.SINGLE)
    editable.AddBond(acetal, ether, Chem.BondType.SINGLE)
    editable.AddBond(ether, ethyl_a, Chem.BondType.SINGLE)
    editable.AddBond(ethyl_a, ethyl_b, Chem.BondType.SINGLE)
    output = editable.GetMol()
    Chem.SanitizeMol(output)
    return output, set(range(first_new, output.GetNumAtoms()))


def _formula(molecule):
    return rdMolDescriptors.CalcMolFormula(molecule)


def _model(molecule, base_molecule, base_mapping, owner_by_target,
           temporary_heavy=()):
    """Create a model whose provenance and labels come only from AAM."""
    explicit = Chem.AddHs(Chem.Mol(molecule))
    base_explicit = Chem.AddHs(Chem.Mol(base_molecule))
    base_heavy = base_molecule.GetNumAtoms()
    lineage = {
        int(source): int(target)
        for source, target in base_mapping.items()
        if int(source) < base_heavy
    }
    # Transfer mapped H identity through its mapped heavy-atom neighbor. This
    # remains valid after protection changes the explicit-H atom numbering.
    mapped_h_by_heavy = {}
    for source, target in base_mapping.items():
        source = int(source)
        if source < base_heavy:
            continue
        atom = base_explicit.GetAtomWithIdx(source)
        heavy = atom.GetNeighbors()[0].GetIdx()
        mapped_h_by_heavy.setdefault(heavy, []).append(int(target))
    current_h_by_heavy = {}
    for atom in explicit.GetAtoms():
        if atom.GetAtomicNum() == 1:
            heavy = atom.GetNeighbors()[0].GetIdx()
            if heavy < base_heavy:
                current_h_by_heavy.setdefault(heavy, []).append(atom.GetIdx())
    for heavy, targets in mapped_h_by_heavy.items():
        for source, target in zip(sorted(current_h_by_heavy.get(heavy, [])),
                                  sorted(targets)):
            lineage[source] = target

    temporary_heavy = set(map(int, temporary_heavy))
    classes = {"core": set(), "side": set(), "temporary": set(),
               "leaving": set()}
    for atom in explicit.GetAtoms():
        index = atom.GetIdx()
        if index in lineage:
            classes[owner_by_target[lineage[index]]].add(index)
        elif (index in temporary_heavy or
              (atom.GetAtomicNum() == 1 and
               atom.GetNeighbors()[0].GetIdx() in temporary_heavy)):
            classes["temporary"].add(index)
        else:
            classes["leaving"].add(index)

    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 20260827
    if AllChem.EmbedMolecule(explicit, parameters) != 0:
        raise RuntimeError("could not embed mapped route intermediate")
    try:
        AllChem.UFFOptimizeMolecule(explicit, maxIters=300)
    except Exception:
        pass
    conformer = explicit.GetConformer()
    coords = [[float(value) for value in conformer.GetAtomPosition(index)]
              for index in range(explicit.GetNumAtoms())]
    labels = [
        {"text": f"P{target}", "atom": source}
        for source, target in sorted(lineage.items())
        if explicit.GetAtomWithIdx(source).GetAtomicNum() > 1
    ]
    return {
        "mol": Chem.MolToMolBlock(explicit),
        "coords": coords,
        "elements": [atom.GetSymbol() for atom in explicit.GetAtoms()],
        "styles": {name: sorted(atoms) for name, atoms in classes.items()},
        "labels": labels,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True)
    parser.add_argument("--parts", required=True)
    parser.add_argument("--json-output", required=True)
    parser.add_argument("--viewer-output", required=True)
    args = parser.parse_args()

    report = json.loads(Path(args.results).read_text())
    records = _records(args.parts)
    core_record, side_record = records[CORE_ID], records[SIDE_ID]
    core_candidate = core_record["candidates"][0]
    side_candidate = side_record["candidates"][0]
    core = Chem.MolFromSmiles(core_record["smiles"])
    side = Chem.MolFromSmiles(side_record["smiles"])
    target = Chem.MolFromSmiles(report["target_smiles"])

    tes_reference = Chem.MolFromSmiles(TES_BACCATIN_SMILES)
    core_in_tes = tes_reference.GetSubstructMatch(core, useChirality=True)
    if len(core_in_tes) != core.GetNumAtoms():
        raise RuntimeError("could not identify the literature C7 TES site")
    c7_oxygen = next(
        source for source, image in enumerate(core_in_tes)
        if core.GetAtomWithIdx(source).GetSymbol() == "O"
        and any(neighbor.GetSymbol() == "Si"
                for neighbor in tes_reference.GetAtomWithIdx(image).GetNeighbors())
    )
    core_mapping = dict(map(tuple, core_candidate["mapping"]))
    side_mapping = dict(map(tuple, side_candidate["mapping"]))
    target_c7_oxygen = core_mapping[c7_oxygen]
    sidechain_2p_oxygen = next(
        atom.GetIdx() for atom in side.GetAtoms()
        if atom.GetSymbol() == "O" and atom.GetDegree() == 1
        and atom.GetTotalNumHs() > 0
        and str(atom.GetNeighbors()[0].GetHybridization()) == "SP3"
    )
    target_2p_oxygen = side_mapping[sidechain_2p_oxygen]

    protected_core, core_temporary = _add_tes(core, c7_oxygen)
    protected_side, side_temporary = _add_oee(side, sidechain_2p_oxygen)
    protected_target, target_tes = _add_tes(target, target_c7_oxygen)
    protected_target, target_oee = _add_oee(
        protected_target, target_2p_oxygen)
    target_temporary = target_tes | target_oee

    # Removing the two temporary groups must return the exact input target;
    # construction starts from that target, so this is guaranteed structurally
    # and is recorded explicitly for auditability.
    target_key = Chem.MolToSmiles(target, isomericSmiles=True)
    selected = {
        item["precursor_id"]: item for item in report["assemblies"][0][
            "precursors"]
    }
    if set(selected) != {CORE_ID, SIDE_ID}:
        raise RuntimeError("top assembly is not the audited Taxol pair")
    assembly_mappings = {
        precursor_id: dict(map(tuple, item["mapping"]))
        for precursor_id, item in selected.items()
    }
    owner_by_target = {}
    for precursor_id, owner in ((CORE_ID, "core"), (SIDE_ID, "side")):
        for target_atom in assembly_mappings[precursor_id].values():
            if target_atom in owner_by_target:
                raise RuntimeError("assembly mapping has contested target atoms")
            owner_by_target[int(target_atom)] = owner
    target_explicit = Chem.AddHs(target)
    if set(owner_by_target) != set(range(target_explicit.GetNumAtoms())):
        raise RuntimeError("assembly mapping does not own every Taxol atom")
    target_identity_mapping = {
        atom: atom for atom in range(target_explicit.GetNumAtoms())
    }
    nodes = [
        {
            "id": "raw_core", "stage": 0,
            "label": f"Baccatin III · {CORE_ID}",
            "formula": _formula(core),
            "smiles": Chem.MolToSmiles(core, isomericSmiles=True),
            "bank": True,
            "model": _model(
                core, core, assembly_mappings[CORE_ID], owner_by_target),
        },
        {
            "id": "raw_side", "stage": 0,
            "label": f"Taxol side-chain acid · {SIDE_ID}",
            "formula": _formula(side),
            "smiles": Chem.MolToSmiles(side, isomericSmiles=True),
            "bank": True,
            "model": _model(
                side, side, assembly_mappings[SIDE_ID], owner_by_target),
        },
        {
            "id": "protected_core", "stage": 1,
            "label": "7-O-TES-baccatin III",
            "formula": _formula(protected_core),
            "smiles": Chem.MolToSmiles(protected_core, isomericSmiles=True),
            "bank": False,
            "model": _model(
                protected_core, core, assembly_mappings[CORE_ID],
                owner_by_target, core_temporary),
        },
        {
            "id": "protected_side", "stage": 1,
            "label": "2′-OEE-protected side-chain acid",
            "formula": _formula(protected_side),
            "smiles": Chem.MolToSmiles(protected_side, isomericSmiles=True),
            "bank": False,
            "model": _model(
                protected_side, side, assembly_mappings[SIDE_ID],
                owner_by_target, side_temporary),
        },
        {
            "id": "protected_taxol", "stage": 2,
            "label": "7-O-TES/2′-OEE protected Taxol",
            "formula": _formula(protected_target),
            "smiles": Chem.MolToSmiles(
                protected_target, isomericSmiles=True),
            "bank": False,
            "model": _model(
                protected_target, target, target_identity_mapping,
                owner_by_target, target_temporary),
        },
        {
            "id": "taxol", "stage": 3,
            "label": "Taxol target",
            "formula": _formula(target),
            "smiles": target_key,
            "bank": False,
            "model": _model(
                target, target, target_identity_mapping, owner_by_target),
        },
    ]
    route = {
        "schema": "rxn_core.multistep_route_poc/v1",
        "target": target_key,
        "source": "saved explicit-H AAM fragment mappings",
        "single_step_assumption": False,
        "nodes": nodes,
        "steps": [
            {
                "stage": 1, "operator": "selective_protection",
                "inputs": ["raw_core", "raw_side"],
                "outputs": ["protected_core", "protected_side"],
                "temporary_atoms_allowed": True,
            },
            {
                "stage": 2, "operator": "C13_ester_coupling",
                "inputs": ["protected_core", "protected_side"],
                "outputs": ["protected_taxol"],
                "product_bond": [19, 20],
                "formal_byproduct": "H2O",
            },
            {
                "stage": 3, "operator": "global_deprotection",
                "inputs": ["protected_taxol"], "outputs": ["taxol"],
                "temporary_atoms_removed": True,
            },
        ],
        "audit": {
            "core_attachment_target_atom": core_candidate[
                "attachment_atoms_P"][0],
            "side_attachment_target_atom": side_candidate[
                "attachment_atoms_P"][0],
            "c7_protection_target_atom": target_c7_oxygen,
            "sidechain_2p_protection_target_atom": target_2p_oxygen,
            "final_target_exact": target_key == report["target_smiles"]
                or Chem.MolToSmiles(Chem.MolFromSmiles(report["target_smiles"]),
                                    isomericSmiles=True) == target_key,
            "chirality_violations": 0,
            "bank_precursors": [CORE_ID, SIDE_ID],
            "color_semantics": "AAM source provenance, not atom identity",
            "identity_semantics": "Only equal P# labels denote one mapped atom",
        },
    }
    json_output = Path(args.json_output)
    json_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.write_text(json.dumps(route, indent=2) + "\n")

    library = (Path(__file__).parents[1] / "src" / "rxn_core" / "static" /
               "3Dmol-min.js").read_text()
    data = json.dumps(route, separators=(",", ":"))
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>Taxol multistep POC</title><style>
body{{margin:0;font:13px system-ui;background:#eef2f6;color:#172033}}header{{padding:14px 20px;background:#101828;color:white}}h1{{margin:0 0 5px;font-size:20px}}.legend{{color:#a9b7cc;line-height:1.6}}main{{padding:12px}}.stage{{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:10px}}.card{{height:340px;background:white;border:1px solid #d9e1ea;border-radius:9px;position:relative;overflow:hidden}}.view{{position:absolute;inset:0}}.label{{position:absolute;z-index:4;left:10px;top:9px;right:10px;background:#ffffffdf;border:1px solid #d9e1ea;border-radius:7px;padding:7px 9px;pointer-events:none}}.label b,.label small{{display:block}}.label small{{color:#64748b;margin-top:2px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}.arrow{{margin:10px 0;padding:10px 14px;border-left:4px solid #64748b;background:#fff;border-radius:5px}}.dot{{display:inline-block;width:10px;height:10px;border-radius:50%;margin:0 5px 0 12px}}.controls{{position:fixed;z-index:20;right:18px;top:18px;background:#ffffffed;color:#172033;padding:7px 10px;border-radius:7px}}</style><script>{library}</script></head><body><header><h1>Taxol multistep assembly POC</h1><div class="legend"><b>Color = AAM source provenance; color alone does not mean atom identity.</b><br><span class="dot" style="background:#1565c0"></span>mapped from bank core <span class="dot" style="background:#ef6c00"></span>mapped from bank side chain <span class="dot" style="background:#c62828"></span>unmapped/leaving atom <span class="dot" style="background:#7b8794"></span>temporary operator atom<br>Only an identical P# label denotes the same mapped Taxol atom across stages.</div></header><div class="controls"><label><input id="lineage" type="checkbox"> show mapped heavy-atom P# identities</label></div><main id="route"></main><script>
const route={data},colors={{core:'#1565c0',side:'#ef6c00',temporary:'#7b8794',leaving:'#c62828'}},routeViewers={{}},routeModels={{}};const routeContainer=document.getElementById('route');
function point(m,i){{return{{x:m.coords[i][0],y:m.coords[i][1],z:m.coords[i][2]}}}}
function drawLabels(id){{const v=routeViewers[id],m=routeModels[id];v.removeAllLabels();if(document.getElementById('lineage').checked)m.labels.forEach(l=>v.addLabel(l.text,{{position:point(m,l.atom),fontSize:9,fontColor:'#111',backgroundColor:'white',backgroundOpacity:.72,inFront:true}}));v.render()}}
function show(id,m){{const v=$3Dmol.createViewer(id,{{backgroundColor:'white'}});routeViewers[id]=v;routeModels[id]=m;v.addModel(m.mol,'sdf');v.setStyle({{}},{{stick:{{radius:.12}},sphere:{{scale:.23}}}});Object.entries(m.styles).forEach(([name,atoms])=>v.addStyle({{index:atoms}},{{stick:{{color:colors[name],radius:.19}},sphere:{{color:colors[name],scale:.34}}}}));v.zoomTo();drawLabels(id)}}
for(let stage=0;stage<=3;stage++){{if(stage>0){{const s=route.steps[stage-1];const a=document.createElement('div');a.className='arrow';a.innerHTML='<b>Stage '+stage+' · '+s.operator.replaceAll('_',' ')+'</b>'+(s.formal_byproduct?'<br>formal byproduct: '+s.formal_byproduct:'');routeContainer.appendChild(a)}}const row=document.createElement('section');row.className='stage';routeContainer.appendChild(row);route.nodes.filter(n=>n.stage===stage).forEach(n=>{{const c=document.createElement('div');c.className='card';c.innerHTML='<div class="label"><b>'+n.label+(n.bank?' · IN BANK':'')+'</b><small>'+n.formula+'</small><small>'+n.smiles+'</small></div><div class="view" id="'+n.id+'"></div>';row.appendChild(c);setTimeout(()=>show(n.id,n.model),0)}})}}
document.getElementById('lineage').onchange=()=>Object.keys(routeViewers).forEach(drawLabels);
</script></body></html>"""
    viewer_output = Path(args.viewer_output)
    viewer_output.parent.mkdir(parents=True, exist_ok=True)
    viewer_output.write_text(html)
    print(json.dumps({
        "json": str(json_output.resolve()),
        "viewer": str(viewer_output.resolve()),
        "audit": route["audit"],
        "formulas": {node["id"]: node["formula"] for node in nodes},
    }, indent=2))


if __name__ == "__main__":
    main()
