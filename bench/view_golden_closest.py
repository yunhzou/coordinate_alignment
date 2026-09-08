"""Post-hoc reference-nearest saved representatives; never reruns AAM.

Distance is literal heavy-atom pair disagreement, not family or symmetry distance.
An optional reactant subset must occupy its reference target set collectively.
"""
import argparse
import ast
import colorsys
import json
from pathlib import Path

from rdkit import Chem
from golden_policy_campaign import load_case, save
from view_golden_remaining import graph_for, ranker, draw_record
from view_golden_mapping import molecules, mapping_fragments, PALETTE


def labels(reaction):
    result = []
    for side in reaction.split('>>'):
        mol = Chem.MolFromSmiles(side)
        original = [a.GetAtomMapNum() for a in mol.GetAtoms()]
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.MolToSmiles(mol, canonical=True)
        result.append([original[i] for i in ast.literal_eval(mol.GetProp('_smilesAtomOutputOrder'))])
    return result


def main(args):
    directory, plan = load_case(args.source/'random', args.index)
    problem = plan.input_problem
    reference = dict(json.loads((directory/'reference.json').read_text())['mapping'])
    reference = {int(r):p for r,p in reference.items()}
    audit = next(json.loads(line) for line in args.audit.read_text().splitlines()
                 if json.loads(line)['index'] == args.index)
    rlabels, plabels = labels(audit['mapped_reaction'])
    subset = {r for r,n in enumerate(rlabels) if n in args.subset}
    targets = {reference[r] for r in subset}
    best = None
    archives = []
    for policy in ('random', 'distance'):
        directory, oriented = load_case(args.source/policy, args.index)
        graph, archive, scope = graph_for(directory)
        count = 0
        for terminal in graph.terminals:
            mapping = oriented.to_input_mapping(graph.states[terminal].mapping)
            if {mapping.get(r) for r in subset} != targets:
                continue
            count += 1
            differences = sum(mapping.get(r) != p for r,p in reference.items())
            key = (differences, tuple(sorted(mapping.items())))
            if best is None or key < best[0]:
                best = (key, policy, terminal, mapping, str(archive), scope)
        archives.append(dict(policy=policy, path=str(archive), scope=scope,
            terminals=len(graph.terminals), capped=graph.capped,
            evaluation='unknown', eligible_representatives=count))
        print(policy, 'eligible', count, 'best differences', best[0][0] if best else None, flush=True)
        # Keep only the selected path, with the actual saved fragment decomposition.
        if best is not None and best[1] == policy:
            path = next(graph.paths(best[2]))
            steps = []
            for eid in path.transitions:
                edge = graph.transitions[eid]
                if edge.match is None:
                    continue
                group = edge.match['fragment']
                if oriented.reversed:
                    group = [path.mapping[r] for r in group]
                steps.append(dict(fragment=list(group)))
            selected = dict(label='Reference-nearest saved representative (post-hoc)',
                terminal=best[2], mapping=best[3], archive=best[4], scope=best[5],
                steps=steps, context=dict(cuts=path.context.cuts))
        del graph
    records = [dict(label='Ground truth', mapping=reference, steps=[], context=None), selected]
    reactant, product = molecules(audit['mapped_reaction'], problem)
    owner = {a:i for i,atoms in enumerate(Chem.GetMolFrags(reactant)) for a in atoms}
    colors = {}
    for record in records:
        record['regions'] = mapping_fragments(record, problem)
        for region in record['regions']:
            heavy = tuple(r for r in region['source'] if problem.reactant.elements[r] != 'H')
            signature = ('heavy', heavy) if heavy else ('H', tuple(region['source']))
            if signature not in colors:
                n = len(colors)
                rgb = colorsys.hsv_to_rgb((n*.61803398875)%1, .65, .8)
                colors[signature] = PALETTE[n] if n < len(PALETTE) else '#'+''.join(f'{round(x*255):02x}' for x in rgb)
            region['color'] = colors[signature]
        record['heavy_events'] = ranker(problem, True)(record['mapping'])[1]
        record['mapped_heavy'] = sum(problem.product.elements[p] != 'H' for p in record['mapping'].values())
        record['mapped_total'] = len(record['mapping'])
        draw_record(record, reactant, product, problem, owner)
    differences = [dict(r=r, reference_p=p, detected_p=selected['mapping'].get(r),
        original_r=rlabels[r], original_reference_p=plabels[p],
        original_detected_p=plabels[selected['mapping'][r]])
        for r,p in reference.items() if selected['mapping'].get(r) != p]
    payload = dict(index=args.index, records=records, archives=archives,
        p_elements=problem.product.elements, r_elements=problem.reactant.elements,
        target_heavy=sum(e != 'H' for e in problem.product.elements),
        target_total=problem.target_atom_count, reaction=audit['mapped_reaction'])
    args.output.mkdir(parents=True, exist_ok=True)
    save(args.output/'mapping.json', {**payload, 'records':[{k:v for k,v in r.items() if k != 'drawings'} for r in records],
        'differences':differences, 'required_source_map_labels':args.subset})
    template = Path(__file__).with_name('golden_remaining_template.html').read_text()
    template = template.replace('21 unresolved cases • saved search results only • no reruns or reference-guided results inserted',
        'Post-hoc reference-nearest saved representative • no AAM rerun • not a top-ranked prediction')
    template = template.replace('/21', '/${cases.length}')
    template = template.replace('${a.distinct_raw_heavy.toLocaleString()} raw-heavy representatives',
                                '${a.eligible_representatives.toLocaleString()} eligible representatives')
    start = template.index('The detected choices are up to three')
    end = template.index('</p>', start)
    template = template[:start] + ('Selected across both saved archives by fewest literal heavy-atom reference disagreements, '
        f'with source map-label subset {args.subset} required to occupy its reference target set. '
        'This uses the reference for display selection only. It does not search every realization inside compressed families, '
        'and literal differences can include equivalent symmetry choices.') + template[end:]
    template = template.replace('__DATA__', json.dumps([payload], separators=(',', ':')).replace('</', r'<\/'))
    (args.output/'viewer.html').write_text(template)
    print(json.dumps(dict(differences=differences, events=[r['heavy_events'] for r in records],
                         html=str((args.output/'viewer.html').resolve())), indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index', type=int, required=True)
    p.add_argument('--subset', nargs='*', type=int, default=[])
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--audit', type=Path, default=Path('data/aam_benchmarks/golden_original_20260906/audit.jsonl'))
    main(p.parse_args())
