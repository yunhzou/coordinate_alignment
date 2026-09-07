"""Offline 2D reference/actual AAM comparison from an existing checkpoint."""
import argparse
import ast
import colorsys
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


def svg(mol, colors, hydrogens, width, height, *, mapping=None, owner=None, side='R', bonds=None, reference=False):
    mol = Chem.Mol(mol) if hydrogens else Chem.RemoveHs(mol)
    mapping = {} if mapping is None else mapping
    owner = {} if owner is None else owner
    inverse = {p:r for r,p in mapping.items()}
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.annotationFontScale = .65
    options.padding = .10
    atom_colors = {}
    for atom in mol.GetAtoms():
        old = atom.GetIntProp('original_index')
        # Both drawings use the SAME source identifier. P's own index remains
        # available in the linked inspector and exact correspondence table.
        source = old if side == 'R' else inverse.get(old)
        atom.SetProp('atomNote', f'r{source}' if source is not None else f'p{old}?')
        if old in colors:
            color = colors[old].lstrip('#')
            atom_colors[atom.GetIdx()] = tuple(int(color[i:i+2], 16)/255 for i in (0, 2, 4))
    bond_colors = {b.GetIdx(): atom_colors[b.GetBeginAtomIdx()] for b in mol.GetBonds()
                   if b.GetBeginAtomIdx() in atom_colors and
                   atom_colors.get(b.GetEndAtomIdx()) == atom_colors[b.GetBeginAtomIdx()] and
                   (bonds is None or tuple(sorted((mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetIntProp('original_index'),
                                                   mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetIntProp('original_index')))) in bonds)}
    drawer.DrawMolecule(mol, highlightAtoms=list(atom_colors), highlightBonds=list(bond_colors),
                        highlightAtomColors=atom_colors, highlightBondColors=bond_colors)
    drawer.FinishDrawing()
    drawing = drawer.GetDrawingText().split('<?xml version=\'1.0\' encoding=\'iso-8859-1\'?>')[-1]
    hits = []
    for atom in mol.GetAtoms():
        old = atom.GetIntProp('original_index')
        source = old if side == 'R' else inverse.get(old)
        target = mapping.get(old) if side == 'R' else old
        position = drawer.GetDrawCoords(atom.GetIdx())
        if source is None:
            description = f'P:p{target} — no source assignment in this mapping'
        elif target is None:
            meaning = 'H identity not reference-annotated' if reference and atom.GetAtomicNum() == 1 else 'not mapped to P'
            description = f'R{owner.get(source, 0)+1}:r{source} — {meaning}'
        else:
            description = f'R{owner.get(source, 0)+1}:r{source} → P:p{target}'
        hits.append(f'<circle class="atom-hit" cx="{position.x:.3f}" cy="{position.y:.3f}" r="10" '
                    f'data-r="{source if source is not None else ""}" data-p="{target if target is not None else ""}" '
                    f'data-side="{side}" data-description="{html.escape(description, quote=True)}" '
                    f'tabindex="0" role="button" aria-label="{html.escape(description, quote=True)}">'
                    f'<title>{html.escape(description)}</title></circle>')
    return drawing.replace('</svg>', '<g class="atom-hits">'+''.join(hits)+'</g></svg>')


def mapping_fragments(record, problem):
    """AAM transition fragments, or explicitly labelled reference-derived regions."""
    mapping = record['mapping']
    if record['context'] is not None:
        groups = [set(step['fragment']) & mapping.keys() for step in record['steps']]
        groups = [group for group in groups if group]
        assert sum(map(len, groups)) == len(set().union(*groups)) == len(mapping)
        return [dict(label=f'F{i+1}', source=sorted(group), kind='AAM fragment',
                     target=sorted(mapping[r] for r in group)) for i,group in enumerate(groups)]
    # The reference contains no AAM search history. Do not invent one.
    remaining = set(mapping)
    groups = []
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        group = {seed}
        queue = [seed]
        while queue:
            r = queue.pop()
            neighbors = {s for s in remaining if problem.reactant.wbo[r,s] > .2 and
                         abs(problem.reactant.wbo[r,s]-problem.product.wbo[mapping[r],mapping[s]]) < 1e-6}
            remaining -= neighbors
            group |= neighbors
            queue.extend(neighbors)
        groups.append(group)
    return [dict(label=f'G{i+1}', source=sorted(group), kind='Reference conserved region',
                 target=sorted(mapping[r] for r in group)) for i,group in enumerate(groups)]


def render_viewer(payload, reactant, product, problem, config):
    records = payload['records']
    parts = Chem.GetMolFrags(reactant, asMols=True)
    owner = {atom:i for i, atoms in enumerate(payload['components']) for atom in atoms}
    expected_inverse = {p:r for r,p in records[0]['mapping'].items()}
    color_by_region = {}
    for record in records:
        record['fragments'] = mapping_fragments(record, problem)
        for fragment in record['fragments']:
            heavy = tuple(p for p in fragment['target'] if problem.product.elements[p] != 'H')
            signature = ('heavy', heavy) if heavy else ('H', tuple(fragment['target']))
            if signature not in color_by_region:
                i = len(color_by_region)
                rgb = colorsys.hsv_to_rgb((i*.61803398875) % 1, .68, .77)
                color_by_region[signature] = PALETTE[i] if i < len(PALETTE) else '#'+''.join(f'{round(v*255):02x}' for v in rgb)
            fragment['color'] = color_by_region[signature]
    views = []
    for record in records:
        mapping = record['mapping']
        inverse = {p:r for r,p in mapping.items()}
        colors = {r:f['color'] for f in record['fragments'] for r in f['source']}
        target_colors = {p:colors[r] for r,p in mapping.items()}
        source_bonds = {tuple(sorted((a,b))) for a in mapping for b in mapping if a < b and
                        problem.reactant.wbo[a,b] > .2 and problem.product.wbo[mapping[a],mapping[b]] > .2}
        target_bonds = {tuple(sorted((mapping[a],mapping[b]))) for a,b in source_bonds}
        legend = ''.join(f'<span class="chip" style="--fragment:{f["color"]}">{f["label"]}: '
                         f'{len([p for p in f["target"] if problem.product.elements[p] != "H"])} heavy + '
                         f'{len([p for p in f["target"] if problem.product.elements[p] == "H"])} H</span>' for f in record['fragments'])
        compare = ''
        if record['label'] == 'Closest saved representative':
            compare = '<p>Selected by minimum source/target orbit-pair count distance among saved representatives, not by globally optimizing symmetry realizations.</p>'
        if record['context'] is not None:
            identical = sum(inverse.get(p) == r for p,r in expected_inverse.items())
            compare += f'<p>{identical}/{len(expected_inverse)} literal reference pairs agree. '
            compare += 'Different raw labels may be symmetry-equivalent; they are not silently relabelled here. See the comparison table.</p>'
        panels = []
        for hydrogens in (False, True):
            cards = []
            contributing = {owner[r] for r in mapping}
            for i in sorted(range(len(parts)), key=lambda i:(i not in contributing, i)):
                tags = ', '.join(f['label'] for f in record['fragments'] if any(owner[r] == i for r in f['source']))
                count = sum(owner[r] == i for r in mapping)
                cards.append(f'<article class="source-card {"unused" if i not in contributing else ""}" data-component="{i}">'
                             f'<h3>Source molecule R{i+1} <span>{tags or "no matched region"}</span></h3>'
                             f'<p>{record["heavy_donors"].get(i,0)} heavy atoms / {count} total mapped</p>'+
                             svg(parts[i], colors, hydrogens, 520, 365, mapping=mapping, owner=owner, bonds=source_bonds, reference=record['context'] is None)+
                             '</article>')
            panels.append(f'<div class="layout hydrogen-{int(hydrogens)}"><div class="sources">'+''.join(cards)+
                          '</div><article class="target"><h2>Product P — labels identify the mapped R atom</h2>'+
                          svg(product, target_colors, hydrogens, 980, 860, mapping=mapping, owner=owner, side='P', bonds=target_bonds)+
                          '<p>Click an atom to locate its source. The same <b>r-number</b> on both sides is one mapped atom.</p></article></div>')
        rows = []
        for p in range(problem.target_atom_count):
            r, ref = inverse.get(p), expected_inverse.get(p)
            identity = lambda a: 'unassigned' if a is None else f'R{owner[a]+1}:r{a}'
            comparison = ('H not reference-annotated' if problem.product.elements[p]=='H'
                          else 'reference-unmatched') if ref is None else 'same pair' if r == ref else 'different raw pair'
            rows.append(f'<tr data-r="{r if r is not None else ""}" data-p="{p}" data-side="P" class="{"raw-difference" if ref is not None and r != ref else ""}">'
                        f'<td>p{p} ({problem.product.elements[p]})</td><td>{identity(r)}</td><td>{identity(ref)}</td><td>{comparison}</td></tr>')
        details = '' if record['context'] is None else '<p>Actual sweep cut on '+record['context'].get('seed_side','r').upper()+': '+html.escape(str(record['context']['cuts']))+'</p><ol>'+''.join(
            f'<li>Seed {record["context"].get("seed_side","r")}{s["seed"]}: mapped R atoms [{", ".join("r"+str(r) for r in s["fragment"])}]</li>' for s in record['steps'])+'</ol>'
        views.append(f'<section><h2>{record["label"]}</h2><p>{"Colors show actual AAM transition fragments (F)." if record["context"] is not None else "Colors show reference-derived conserved regions (G), not AAM fragments."}</p>'+
                     legend+compare+''.join(panels)+'<details><summary>Exact P → R comparison table (click a row to highlight)</summary>'+
                     '<table><tr><th>Product atom</th><th>Selected mapping</th><th>Reference mapping</th><th>Raw-index comparison</th></tr>'+''.join(rows)+
                     '</table></details><details><summary>Saved fragment-growth decisions</summary>'+details+'</details></section>')
    options = ''.join(f'<option>{r["label"]}</option>' for r in records)
    return ('''<!doctype html><meta charset="utf-8"><title>Inspect exact AAM correspondence</title>
<style>
body{font:16px system-ui;margin:24px;background:#f4f6f8;color:#172432}.notice{background:#fff0d5;padding:12px;border-left:5px solid #d55e00}
.controls{position:sticky;top:0;z-index:5;background:#edf2f7;padding:10px;box-shadow:0 2px 4px #bbc}select,button{font:inherit;padding:8px;margin-right:12px}label{margin-right:12px}
#inspector{font-weight:600;padding:8px 0;color:#164a79}.layout{display:grid;grid-template-columns:40% 1fr;gap:18px}.sources{display:grid;gap:12px}article{background:white;border:1px solid #ccd5dd;border-radius:10px;padding:10px}article p{margin:6px 0}h3{margin:5px;font-size:18px}h3 span{float:right}svg{width:100%;height:auto}.target{position:sticky;top:112px;align-self:start}details{background:white;padding:18px;margin:12px 0}td,th{padding:7px 15px;text-align:left}tr[data-p]{cursor:pointer}tr.raw-difference{background:#fff3dc}tr.selected{outline:2px solid #173bdf}section{display:none}.hydrogen-1{display:none}
.chip{display:inline-block;margin:4px 10px 4px 0;padding:5px 10px;border-left:9px solid var(--fragment);background:white}.atom-hit{fill:transparent;stroke:transparent;pointer-events:all;cursor:pointer}.atom-hit.selected,.atom-hit:focus{stroke:#173bdf;stroke-width:3;fill:#b8c8ff;fill-opacity:.3;outline:none}body.hide-unused .unused{display:none}
@media(max-width:1000px){.layout{grid-template-columns:1fr}.target{position:static;grid-row:1}.controls{position:static}}
</style>'''+f'<header><h1>Case {payload["index"]} — inspect the atom correspondence</h1>'+
        '<p><b>r57 means source atom 57 on BOTH drawings.</b> R5 means source molecule 5. Click an atom to see its own product index and source molecule. All indices are zero-based.</p>'+
        f'<p>{config.seed_count} seed orders · cap {config.branch_limit} · tolerance {config.iso_tolerance} · explicit-H search. {payload["cap_stops"]:,} cap stops in this archive. Reference recovery: {html.escape(payload["evaluation"]["reference_recovery"])}</p>'+
        '<p class="notice">Reference and AAM are separate views. No symmetry-equivalent copies are relabelled to hide differences. Same colors across views denote the same target heavy-atom region; compare r-labels to verify the actual source assignment. Reference H identities are unannotated (p…?).</p></header>'+
        f'<div class="controls"><select id="mode" aria-label="Mapping to inspect">{options}</select><button id="hydrogens">Show explicit H</button><button id="clear">Clear selection</button>'+
        '<label><input type="checkbox" id="unused"> Show all input molecules</label><div id="inspector" role="status">Hover or click an atom: its mapped partner will be outlined in blue.</div></div>'+''.join(views)+'''
<script>
let allH=false,pinned=false;const mode=document.querySelector('#mode');mode.selectedIndex=1;
function render(){pinned=false;document.querySelectorAll('section').forEach((s,i)=>s.style.display=i===mode.selectedIndex?'block':'none');document.querySelectorAll('.hydrogen-0').forEach(x=>x.style.display=allH?'none':'grid');document.querySelectorAll('.hydrogen-1').forEach(x=>x.style.display=allH?'grid':'none');document.querySelector('#hydrogens').textContent=allH?'Hide explicit H':'Show explicit H';document.body.classList.toggle('hide-unused',!document.querySelector('#unused').checked);document.querySelectorAll('.selected').forEach(x=>x.classList.remove('selected'));document.querySelector('#inspector').textContent='Hover or click an atom: its mapped partner will be outlined in blue.';}
function highlight(node, locate=false){if(locate)pinned=true;const section=node.closest('section');const r=node.dataset.r,p=node.dataset.p;document.querySelectorAll('.selected').forEach(x=>x.classList.remove('selected'));section.querySelectorAll('[data-r]').forEach(x=>{if((r!==''&&x.dataset.r===r)||(r===''&&x.dataset.r===''&&x.dataset.p===p))x.classList.add('selected')});const partner=Array.from(section.querySelectorAll('.atom-hit')).find(x=>x.dataset.r===r&&x.dataset.p===p);document.querySelector('#inspector').textContent=node.dataset.description||partner?.dataset.description||`P:p${p}: no source assignment`;if(locate&&node.dataset.side==='P'&&r!==''){const source=Array.from(section.querySelectorAll('.atom-hit[data-side="R"]')).find(x=>x.dataset.r===r&&x.getBoundingClientRect().width>0);source?.closest('article').scrollIntoView({behavior:'smooth',block:'center'});}}
document.querySelectorAll('.atom-hit,tr[data-p]').forEach(x=>{x.addEventListener('mouseenter',()=>{if(!pinned)highlight(x)});x.addEventListener('focus',()=>highlight(x));x.addEventListener('click',()=>highlight(x,true));x.addEventListener('keydown',e=>{if(e.key==='Enter')highlight(x,true)});});mode.onchange=render;document.querySelector('#hydrogens').onclick=()=>{allH=!allH;render()};document.querySelector('#unused').onchange=render;document.querySelector('#clear').onclick=render;render();
</script>''')


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
    page = render_viewer(payload, reactant, product, problem, config)
    (args.output/'mapping.json').write_text(json.dumps(payload, indent=2)+'\n')
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
