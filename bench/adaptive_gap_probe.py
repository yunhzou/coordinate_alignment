"""Diagnostic only: separate saved successful cut effects from seed-order effects."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from adaptive_fragment_pilot import save
from compare_elementary_outputs import event_counts
from rxn_core.artifacts import read_graph_checkpoint, write_graph_checkpoint
from rxn_core.frag import build_graph
from rxn_core.matcher import _nauty_orbits
from rxn_core.native_search import find_islands_native
from rxn_core.search_symmetry import finalize_graph_symmetry
from rxn_core.conditioned_symmetry import ConditionedSymmetryWorkspace


def best(graph,r,p):
    terminals=[t for t in graph.terminals if len(graph.states[t].mapping)==len(r)]
    vectors=[list(dict(graph.states[t].mapping).values()) for t in terminals]
    scores=event_counts(r,p,vectors).sum(axis=1) if vectors else []
    return min(zip(map(int,scores),terminals),default=(None,None))


def main(args):
    bank=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
    source=bank/'fragment_native_final_20260909_JGOrLU/results/77/P_to_R/independent_native_dependency'
    raw=json.loads((bank/'conditioned_reuse_full_20260910_0o09px/inputs/77.json').read_text())
    left,right=raw['product'],raw['reactant']
    r,p=[np.asarray(s['wbo']) for s in (left,right)]
    source_graph=build_graph(left['elements'],r,bond_cut=.2)
    target=build_graph(right['elements'],p,bond_cut=.2)
    po=_nauty_orbits(target,wbo_tol=1.)
    workspace=ConditionedSymmetryWorkspace(target,1.)
    rows=[]
    for seed in range(10):
        folder=source/f'seed_{seed:02d}'
        summary=json.loads((folder/'summary.json').read_text())
        ordinal=next(i for i,c in enumerate(summary['cuts']) if c['best']==8)
        saved=read_graph_checkpoint(folder/f'cut_{ordinal:04d}.pkl.gz')
        saved_score,terminal=best(saved,r,p)
        path=next(saved.paths(terminal))
        fragments=[dict(seed=saved.transitions[e].seed,fragment=saved.transitions[e].match['fragment'],
            mapping=saved.states[saved.transitions[e].target].mapping,deferred=saved.transitions[e].match['deferred_edges'])
            for e in path.transitions if saved.transitions[e].match is not None]
        start=time.process_time()
        graph=find_islands_native(source_graph,target,saved.contexts[0].seed_order,
            graph_floor=.2,iso_tol=1.,max_branches=100,p_orbits=po)
        graph,_=finalize_graph_symmetry(graph,target,iso_tolerance=1.,workspace=workspace)
        compute=time.process_time()-start
        write_graph_checkpoint(graph,args.run/f'uncut_seed_{seed:02d}.pkl.gz')
        score,new_terminal=best(graph,r,p)
        row=dict(seed=seed,source=str(folder/f'cut_{ordinal:04d}.pkl.gz'),cut=saved.contexts[0].cuts,
            order=saved.contexts[0].seed_order,saved_best=saved_score,saved_terminal=terminal,
            saved_fragments=fragments,uncut_best=score,uncut_terminal=new_terminal,
            uncut_cpu=compute,uncut_capped=graph.capped)
        rows.append(row)
        save(args.run/'order_vs_cut.json',dict(rows=rows,
            scope='Diagnostic orders selected from saved successful cut runs; NOT a blind recommendation or performance benchmark. Original WBOs used for event scores.'))
        print(json.dumps({k:v for k,v in row.items() if k not in ('order','saved_fragments')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    main(parser.parse_args())
