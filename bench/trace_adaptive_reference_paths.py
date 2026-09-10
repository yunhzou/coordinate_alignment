"""Read-only diagnosis of where saved original successful paths leave a candidate DAG.

Reference-selected paths are diagnostic evidence only, never search inputs.
No mapper is run by this command.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path

from adaptive_full_benchmark import read, save
from rxn_core.artifacts import read_aam_checkpoint


def key(state):
    return (state.mapping, state.islands, state.deferred_edges)


def trace(run, index, direction):
    folder = run/f'results/golden/{index}/{direction}'
    original = read(folder/'original/full_sweep_evaluation.json')
    if not original.get('representative_recovery'):
        return dict(index=index, direction=direction, symbolic_reference=True)
    old = read_aam_checkpoint(folder/'original/cuts/aam.pkl.gz').graph
    row = read(folder/'adaptive/search.json')['rows'][-1]
    new = read_aam_checkpoint(folder/'adaptive'/row['archive']).graph
    terminal = original['witness_terminal']
    path = next(old.paths(terminal))
    cut = path.context.cuts
    context_ids = {i for i,c in enumerate(new.contexts) if c.cuts == cut}
    known = {key(s):s.id for s in new.states if s.context in context_ids}
    outgoing = defaultdict(list)
    for edge in new.transitions:
        if edge.match is not None:
            outgoing[edge.source].append(dict(seed=edge.seed, fragment=edge.match['fragment']))
    steps = []
    for eid in path.transitions:
        edge = old.transitions[eid]
        if edge.match is None:
            continue
        before, after = old.states[edge.source], old.states[edge.target]
        present = known.get(key(before))
        steps.append(dict(seed=edge.seed, fragment=edge.match['fragment'],
            mapping_before=len(before.mapping), mapping_after=len(after.mapping),
            parent_present=present is not None, child_present=key(after) in known,
            other_decisions=outgoing[present] if present is not None else []))
    return dict(index=index, direction=direction, original_cut=cut,
        original_seed_order=path.context.seed_order, steps=steps,
        candidate_stop=row['stop'], candidate_work=row['work'],
        candidate_visited_cuts=row.get('visited_cuts'), candidate_total_cuts=row.get('total_cuts'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--index', type=int, required=True)
    p.add_argument('--direction', choices=['R_to_P','P_to_R'], required=True)
    args = p.parse_args()
    result = trace(args.run, args.index, args.direction)
    save(args.run/f'diagnosis/{args.index}_{args.direction}.json', result)
    print(json.dumps(result), flush=True)
