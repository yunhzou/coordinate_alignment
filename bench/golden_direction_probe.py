"""Compare both search directions without changing the core AAM algorithm."""
import argparse
from dataclasses import asdict, replace
from collections import Counter
import json
from pathlib import Path
import shutil
import subprocess

import pynauty

from golden_evaluation import colored_graph, project, rank_key
from investigate_golden_mapping import save
from rxn_core import AAMProblem, AAMSearchConfig
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.domain import MolecularEndpoint


def original_mapping(pairs, direction):
    """Re-express one saved injection, without choosing new atom assignments."""
    return dict(pairs) if direction == 'R_to_P' else {r:p for p,r in pairs}


def initialize(run, source):
    run.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    raw = json.loads((source/'input.json').read_text())
    reference = json.loads((source/'reference.json').read_text())
    config = replace(AAMSearchConfig(**json.loads((source.parent/'manifest.json').read_text())['config']),
                     seed_count=3, branch_limit=2000, iso_tolerance=1.0)
    for folder in ('src', 'bench'):
        shutil.copytree(root/folder, run/'engine'/folder, ignore=shutil.ignore_patterns('__pycache__'))
    save(run/'design.json', dict(source=str(source), seed_count=3, cap=2000, tolerance=1.0,
        workers=8, directions=['R_to_P','P_to_R'],
        note='Swap endpoints and reference only. Each direction sweeps edges of its own source. Compare inverse P-to-R mappings in original R-to-P orientation.',
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()))
    for direction in ('R_to_P','P_to_R'):
        directory = run/direction
        inputs = directory/'inputs'/'1665'
        inputs.mkdir(parents=True)
        item = dict(raw)
        truth = dict(reference)
        if direction == 'P_to_R':
            item.update(reactant=raw['product'], product=raw['reactant'], name=raw['name']+'_P_to_R')
            truth.update(features=list(reversed(reference['features'])),
                         mapping=[[p,r] for r,p in dict(reference['mapping']).items()])
        save(inputs/'input.json',item)
        save(inputs/'reference.json',truth)
        save(inputs.parent/'manifest.json',dict(config=vars(config)))
        save(directory/'manifest.json',dict(watchdog_seconds=570, jobs=[dict(
            index=1665,seeds=3,cap=2000,tolerance=1.0,workers=8,source=str(inputs))]))
        (directory/'engine').symlink_to(run/'engine',target_is_directory=True)


def compare(run, output, reverse_cut=None):
    from view_golden_mapping import molecules, render_viewer
    source = Path(json.loads((run/'design.json').read_text())['source'])
    raw = json.loads((source/'input.json').read_text())
    problem = AAMProblem(MolecularEndpoint(**raw['reactant']),MolecularEndpoint(**raw['product']),raw['name'])
    ref = json.loads((source/'reference.json').read_text())
    features = ref['features']
    certificate = pynauty.certificate(colored_graph(features,project(ref['mapping'],features)))
    rows = []
    audit = Path(__file__).resolve().parents[1]/'data/aam_benchmarks/golden_original_20260906/audit.jsonl'
    reaction = next(json.loads(line)['mapped_reaction'] for line in audit.read_text().splitlines()
                    if json.loads(line)['index'] == 1665)
    reactant, product = molecules(reaction, problem)
    from rdkit import Chem
    components = Chem.GetMolFrags(reactant)
    owner = {r:i for i,atoms in enumerate(components) for r in atoms}
    for direction in ('R_to_P','P_to_R'):
        directory = run/direction/'case1665_seeds3_cap2000'
        partial = direction == 'P_to_R' and reverse_cut is not None
        if partial:
            from rxn_core.search_graph import AAMSearchGraph
            graph = AAMSearchGraph.from_record(json.loads(reverse_cut.read_bytes()),copy=False)
            config = AAMSearchConfig(**json.loads((reverse_cut.parent/'manifest.json').read_text())['config'])
            archive = reverse_cut
        else:
            result = read_aam_checkpoint(directory/'cuts/aam.pkl.gz')
            graph, config = result.graph, result.config
            archive = directory/'cuts/aam.pkl.gz'
        def mapping(terminal):
            return original_mapping(graph.states[terminal].mapping,direction)
        ranked = sorted(graph.terminals, key=lambda t:rank_key(mapping(t),problem)[0])
        cache = {}
        hit = None
        for terminal in ranked:
            original = mapping(terminal)
            heavy = project(original,features)
            key = tuple(sorted(heavy.items()))
            if key not in cache:
                cache[key] = pynauty.certificate(colored_graph(features,heavy)) == certificate
            if cache[key]:
                hit = terminal
                break
        witness = None
        if hit is not None:
            path = next(graph.paths(hit))
            witness = dict(terminal=hit,original_R_to_P_mapping=sorted(mapping(hit).items()),
                           transitions=list(path.transitions),cuts=path.context.cuts,
                           seed_order=path.context.seed_order)
        search = None if partial else json.loads((directory/'search.json').read_text())
        row = dict(direction=direction,search_seconds=None if partial else search['seconds'],
            cuts=1 if partial else search['metrics']['cut_count'],terminals=len(ranked),capped=graph.capped,
            cap_stops=sum(s.reason=='capped' for s in graph.stops),
            representative_reference_recovered=hit is not None,
            scope='one completed reverse cut; full sweep stopped after verified recovery' if partial else 'full sweep',
            scoring='Representatives checked with identical original-orientation chemical certificate; absence alone is not a compressed-family miss.',
            original_orientation_top_terminal=ranked[0] if ranked else None,witness=witness,
            archive=str(archive))
        save(directory/'original_orientation_comparison.json',row)
        rows.append(row)
        evaluation = dict(reference_recovery='inconclusive partial representatives') if partial else json.loads((directory/'evaluation.json').read_text())
        selected = hit if hit is not None else ranked[0]
        records = []
        for label, terminal in [('Reference',None),('Verified recovered witness' if hit is not None else 'Top ranked',selected)]:
            current = dict(ref['mapping']) if terminal is None else mapping(terminal)
            path = None if terminal is None else next(graph.paths(terminal))
            steps = []
            if path is not None:
                for i in path.transitions:
                    edge = graph.transitions[i]
                    if edge.match is None:
                        continue
                    atoms = edge.match['fragment']
                    fragment = atoms if direction == 'R_to_P' else [path.mapping[p] for p in atoms if p in path.mapping]
                    steps.append(dict(transition=i,seed=edge.seed,fragment=fragment,search_source_fragment=atoms))
            context = None if path is None else dict(asdict(path.context),seed_side='r' if direction == 'R_to_P' else 'p')
            records.append(dict(label=label,terminal=terminal,mapping=current,context=context,steps=steps,
                                heavy_donors=dict(Counter(owner[r] for r in current if problem.reactant.elements[r]!='H'))))
        payload = dict(index=1665,direction=direction,archive=row['archive'],components=components,records=records,
                       cap_stops=row['cap_stops'],evaluation=dict(reference_recovery=(
                           'recovered; search direction '+direction+(' (completed cut only; remaining sweep stopped)' if partial else '') if hit is not None else evaluation['reference_recovery'])))
        viewer_dir = output/direction
        viewer_dir.mkdir(parents=True,exist_ok=True)
        (viewer_dir/'viewer.html').write_text(render_viewer(payload,reactant,product,problem,config))
        save(viewer_dir/'mapping.json',payload)
    save(run/'comparison.json',dict(results=rows))
    save(output/'comparison.json',dict(results=rows))
    print(json.dumps(rows,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['init','compare'])
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--reverse-cut',type=Path,help='Explicitly compare one saved reverse cut; do not claim a full reverse sweep')
    args = parser.parse_args()
    initialize(args.run,args.source) if args.mode == 'init' else compare(args.run,args.output,args.reverse_cut)
