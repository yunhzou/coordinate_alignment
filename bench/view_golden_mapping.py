"""Offline 2D reference/actual AAM comparison from an existing checkpoint."""
import argparse
import ast
from collections import Counter
from dataclasses import asdict
import html
import json
from pathlib import Path

import numpy as np
import pynauty
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from golden_evaluation import colored_graph, project
from rxn_core.artifacts import read_aam_checkpoint

PALETTE = ['#009e73', '#e69f00', '#0072b2', '#cc79a7', '#d55e00', '#56b4e9', '#8c6bb1', '#a6761d', '#737373']


def molecules(reaction, problem):
    result = []
    for smiles, endpoint in zip(reaction.split('>>'), (problem.reactant, problem.product)):
        mol = Chem.MolFromSmiles(smiles)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.MolToSmiles(mol, canonical=True)
        mol = Chem.AddHs(Chem.RenumberAtoms(mol, ast.literal_eval(mol.GetProp('_smilesAtomOutputOrder'))))
        assert tuple(a.GetSymbol() for a in mol.GetAtoms()) == endpoint.elements
        matrix = np.zeros_like(endpoint.wbo)
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            matrix[i, j] = matrix[j, i] = bond.GetBondTypeAsDouble()
        assert np.array_equal(matrix, endpoint.wbo), 'Display atom order must equal archived search order'
        for atom in mol.GetAtoms():
            atom.SetIntProp('original_index', atom.GetIdx())
        rdDepictor.Compute2DCoords(mol)
        result.append(mol)
    return result


def svg(mol, colors, hydrogens, width, height):
    mol = Chem.Mol(mol) if hydrogens else Chem.RemoveHs(mol)
    # Labels are the original archive indices, including after H removal.
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.annotationFontScale = .65
    options.padding = .10
    atom_colors = {}
    for atom in mol.GetAtoms():
        old = atom.GetIntProp('original_index')
        atom.SetProp('atomNote', str(old))
        if old in colors:
            color = colors[old].lstrip('#')
            atom_colors[atom.GetIdx()] = tuple(int(color[i:i+2], 16)/255 for i in (0, 2, 4))
    bond_colors = {b.GetIdx(): atom_colors[b.GetBeginAtomIdx()] for b in mol.GetBonds()
                   if b.GetBeginAtomIdx() in atom_colors and
                   atom_colors.get(b.GetEndAtomIdx()) == atom_colors[b.GetBeginAtomIdx()]}
    drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors), highlightBonds=list(bond_colors),
                        highlightAtomColors=atom_colors, highlightBondColors=bond_colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().split('<?xml version=\'1.0\' encoding=\'iso-8859-1\'?>')[-1]


def main(args):
    if args.cut_graph:
        from rxn_core import AAMProblem, AAMSearchConfig
        from rxn_core.domain import MolecularEndpoint
        from rxn_core.search_graph import AAMSearchGraph
        raw = json.loads((args.source/'input.json').read_text())
        problem = AAMProblem(MolecularEndpoint(**raw['reactant']), MolecularEndpoint(**raw['product']), raw['name'])
        config = AAMSearchConfig(**json.loads((args.archive.parent/'manifest.json').read_text())['config'])
        graph = AAMSearchGraph.from_record(json.loads(args.archive.read_bytes()), copy=False)
        evaluation = dict(reference_recovery='recovered in this completed cut; full sweep stopped early',
                          top_terminal=args.witness)
    else:
        result = read_aam_checkpoint(args.archive)
        problem, config, graph = result.problem, result.config, result.graph
        evaluation = json.loads((args.archive.parent.parent/'evaluation.json').read_text())
    reference = json.loads((args.source/'reference.json').read_text())
    expected = {int(r): int(p) for r, p in dict(reference['mapping']).items()}
    raw = next(json.loads(line) for line in args.audit.read_text().splitlines()
               if json.loads(line)['index'] == args.index)
    reactant, product = molecules(raw['mapped_reaction'], problem)
    components = Chem.GetMolFrags(reactant)
    parts = Chem.GetMolFrags(reactant, asMols=True)
    owner = {atom: i for i, component in enumerate(components) for atom in component}
    features = reference['features']
    orbits = [pynauty.autgrp(colored_graph([f]))[3][:len(f['heavy'])] for f in features]
    def pairs(mapping):
        return Counter((orbits[0][r], orbits[1][p]) for r, p in project(mapping, features).items())
    goal = pairs(expected)
    def distance(terminal):
        actual = pairs(graph.states[terminal].mapping)
        return sum((actual-goal).values()) + sum((goal-actual).values())
    closest = min(graph.terminals, key=distance)
    if args.cut_graph:
        assert pynauty.certificate(colored_graph(features, project(graph.states[args.witness].mapping, features))) == pynauty.certificate(colored_graph(features, project(expected, features)))
    records = []
    for label, terminal in [('Reference', None), ('Verified recovered witness' if args.cut_graph else 'Top ranked', evaluation['top_terminal']),
                            ('Closest saved representative', closest)]:
        mapping = expected if terminal is None else dict(graph.states[terminal].mapping)
        path = None if terminal is None else next(graph.paths(terminal))
        steps = [] if path is None else [dict(transition=i, seed=graph.transitions[i].seed,
                    fragment=graph.transitions[i].match['fragment'])
                    for i in path.transitions if graph.transitions[i].match is not None]
        counts = Counter(owner[r] for r in mapping if reactant.GetAtomWithIdx(r).GetAtomicNum() != 1)
        records.append(dict(label=label, terminal=terminal, mapping=mapping,
            orbit_distance=None if terminal is None else distance(terminal),
            heavy_donors=dict(counts), context=None if path is None else asdict(path.context), steps=steps))
    args.output.mkdir(parents=True, exist_ok=True)
    payload = dict(index=args.index, archive=str(args.archive), records=records,
        components=[list(c) for c in components], evaluation=evaluation,
        cap_stops=sum(s.reason == 'capped' for s in graph.stops))
    (args.output/'mapping.json').write_text(json.dumps(payload, indent=2)+'\n')
    views = []
    for record in records:
        mapping = record['mapping']
        colors = {r: PALETTE[owner[r] % len(PALETTE)] for r in mapping}
        target_colors = {p: colors[r] for r, p in mapping.items()}
        panels = []
        for hydrogens in (False, True):
            cards = ''.join(f'<article><h3 style="color:{PALETTE[i % len(PALETTE)]}">R{i+1}: '
                f'{record["heavy_donors"].get(i, 0)} heavy atoms donated</h3>'
                + svg(part, colors, hydrogens, 430, 310)+'</article>' for i, part in enumerate(parts))
            panels.append(f'<div class="layout hydrogen-{int(hydrogens)}"><div class="sources">{cards}</div>'
                          f'<article class="target"><h2>P target — assembled atom origins</h2>{svg(product, target_colors, hydrogens, 800, 710)}</article></div>')
        details = '' if record['context'] is None else '<p>Actual sweep cut: '+html.escape(str(record['context']['cuts']))+'</p><ol>'+''.join(
            f'<li>Seed R atom {s["seed"]}: fragment [{", ".join(map(str, s["fragment"]))}]</li>' for s in record['steps'])+'</ol>'
        rows = ''.join(f'<tr><td>{r}</td><td>{p}</td><td>R{owner[r]+1}</td></tr>' for r, p in sorted(mapping.items(), key=lambda item:item[1]))
        views.append('<section>'+''.join(panels)+'<details><summary>Saved mapping and fragment decisions</summary>'+details+
                     '<table><tr><th>Source atom</th><th>Target atom</th><th>Component</th></tr>'+rows+'</table></details></section>')
    legend = ' '.join(f'<span style="color:{PALETTE[i % len(PALETTE)]}">● R{i+1}</span>' for i in range(len(parts)))
    options = ''.join(f'<option>{r["label"]}</option>' for r in records)
    page = '''<!doctype html><meta charset="utf-8"><title>Golden mapping comparison</title>
<style>body{font:16px system-ui;margin:24px;background:#f4f6f8;color:#172432}header{max-width:1200px}.notice{background:#fff0d5;padding:16px;border-left:5px solid #d55e00}select,button{font:inherit;padding:10px;margin:12px 10px 12px 0}.layout{display:grid;grid-template-columns:1fr 1fr;gap:18px}.sources{display:grid;grid-template-columns:1fr 1fr;gap:12px}article{background:white;border:1px solid #ccd5dd;border-radius:10px;padding:10px}h3{margin:5px;font-size:15px}svg{width:100%;height:auto}.target{position:sticky;top:10px;align-self:start}details{background:white;padding:18px;margin:12px 0}td,th{padding:5px 20px}section{display:none}.hydrogen-1{display:none}@media(max-width:1000px){.layout{grid-template-columns:1fr}.target{position:static;grid-row:1}}</style>
'''+f'<header><h1>Case {args.index}: reference versus saved AAM mapping</h1><p>R: {reactant.GetNumAtoms()} atoms; P: {product.GetNumAtoms()} atoms. Search includes explicit H. Tolerance {config.iso_tolerance}; branch cap {config.branch_limit}.</p>'+f'<div class="notice">Reference recovery: {evaluation["reference_recovery"]}. {payload["cap_stops"]:,} cap stops were recorded. This search is not exhaustive. The diagrams show actual saved representatives, not every compressed alternative. Colors identify source components, not AAM fragment boundaries. Exact fragment decisions are listed below.</div><p>Reference specifies heavy-atom mapping; unannotated H have no invented reference assignment. Uncolored source atoms are not donated in the selected mapping. Atom labels are original zero-based archive indices. Layout is generated in 2D for display only.</p><p>“Closest” minimizes the L1 difference in source/target chemical-orbit pair counts over saved representatives; it is diagnostic, not a new ranking or proof of the closest possible symmetry realization.</p></header><select id="mode">{options}</select><button id="hydrogens">Show explicit H</button><p>{legend}</p>'+''.join(views)+'''
<script>let allH=false;function render(){document.querySelectorAll('section').forEach((s,i)=>s.style.display=i===document.querySelector('#mode').selectedIndex?'block':'none');document.querySelectorAll('.hydrogen-0').forEach(x=>x.style.display=allH?'none':'grid');document.querySelectorAll('.hydrogen-1').forEach(x=>x.style.display=allH?'grid':'none');document.querySelector('#hydrogens').textContent=allH?'Hide explicit H':'Show explicit H'}document.querySelector('#mode').onchange=render;document.querySelector('#hydrogens').onclick=()=>{allH=!allH;render()};render();</script>'''
    page = page.replace('Search includes explicit H.',
                        f'Search includes explicit H. {config.seed_count} seed orders per sweep-cut configuration.')
    (args.output/'viewer.html').write_text(page)
    print(json.dumps([{k:v for k,v in r.items() if k in ('label','terminal','orbit_distance','heavy_donors')} for r in records], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cut-graph', action='store_true', help='Read a completed raw cut instead of a full AAM checkpoint')
    parser.add_argument('--witness', type=int, help='Verified terminal to display for --cut-graph')
    main(parser.parse_args())
