"""Checkpointed independent-fragment AAM experiment; no sequential-core edits."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time

from golden_policy_campaign import load_case, save
from rxn_core.independent_aam import detect_independent, independent_paths, assemble_independent
from rxn_core.artifacts import write_graph_checkpoint, read_graph_checkpoint
from rxn_core.alignment.sweep import cut_sweep_items


def task(payload):
    source, index, out, cut_index, cut, cap, all_atoms = payload
    _, plan = load_case(Path(source), index)
    p = plan.problem
    seeds = tuple(i for i,e in enumerate(p.reactant.elements) if all_atoms or e != 'H')
    start = time.monotonic()
    result = detect_independent(p, seeds=seeds, cuts=cut, branch_limit=cap,
        iso_tolerance=plan.config.iso_tolerance, graph_floor=plan.config.graph_floor)
    out = Path(out)
    write_graph_checkpoint(result.graph, out/f'cut_{cut_index:04d}.pkl.gz')
    record = dict(cut_index=cut_index, cut=cut, seconds=time.monotonic()-start,
        attempts=result.attempts, families=len(result.graph.terminals), capped=result.graph.capped)
    save(out/f'cut_{cut_index:04d}.json', record)
    return record


def detect(args):
    _, plan = load_case(args.source, args.index)
    cuts = cut_sweep_items(plan.problem.reactant.wbo, plan.config.cut_floor) if args.sweep else [()]
    args.output.mkdir(parents=True, exist_ok=False)
    save(args.output/'manifest.json', dict(source=str(args.source.resolve()), index=args.index,
        direction=plan.direction, cuts=cuts, cap=args.cap, all_atom_seeds=args.all_atoms,
        tolerance=plan.config.iso_tolerance, workers=args.workers,
        note='Independent empty-context matching. Explicit H retained. No synthetic completion fragments.'))
    work = [(str(args.source), args.index, str(args.output), i, cut, args.cap, args.all_atoms) for i,cut in enumerate(cuts)]
    start = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for r in pool.map(task, work):
            print(json.dumps({k:r[k] for k in ('cut_index','seconds','families','capped')}), flush=True)
    save(args.output/'search.json', dict(wall_seconds=time.monotonic()-start, cuts=len(cuts)))


def assemble(args):
    manifest = json.loads((args.output/'manifest.json').read_text())
    source, plan = load_case(Path(manifest['source']), manifest['index'])
    paths = tuple(independent_paths(read_graph_checkpoint(p) for p in sorted(args.output.glob('cut_*.pkl.gz'))))
    reference = json.loads((source/'reference.json').read_text())
    required = plan.to_search_mapping(dict(reference['mapping'])) if args.reference else None
    result = assemble_independent(plan.problem, paths, seconds=args.seconds, required_mapping=required)
    result.update(families=len(paths), reference_query='literal original reference' if args.reference else None,
        completed_cuts=len(list(args.output.glob('cut_*.pkl.gz'))), expected_cuts=len(manifest['cuts']))
    if result['status'] == 'covered':
        result['input_mapping'] = sorted(plan.to_input_mapping(dict(result['mapping'])).items())
        result['fragments'] = [dict(mapping=sorted(paths[i].mapping.items()), terminal=paths[i].terminal,
            cuts=paths[i].context.cuts) for i in result['selected']]
    save(args.output/('reference_assembly.json' if args.reference else 'assembly.json'), result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('detect','assemble'))
    p.add_argument('--source', type=Path)
    p.add_argument('--index', type=int, default=590)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sweep', action='store_true')
    p.add_argument('--all-atoms', action='store_true')
    p.add_argument('--reference', action='store_true')
    p.add_argument('--cap', type=int, default=100)
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--seconds', type=float, default=120)
    args = p.parse_args()
    globals()[args.mode](args)
