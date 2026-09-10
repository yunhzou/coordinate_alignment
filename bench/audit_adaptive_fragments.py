"""Validate saved adaptive graphs and inspect correlated families, no search."""
import argparse
import json
from pathlib import Path
import time

import networkx as nx
import numpy as np

from adaptive_fragment_pilot import save
from compare_elementary_outputs import event_counts
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.family_scoring import minimize_events


def audit(args):
    folder=args.run/f'results/{args.slot}'
    summary=json.loads((folder/'summary.json').read_text())
    results=[]
    for row in summary['rows']:
        start=time.perf_counter();result=read_aam_checkpoint(folder/f"{row['label']}.pkl.gz")
        g=result.graph;p=result.problem;edges=0
        assert nx.is_directed_acyclic_graph(nx.DiGraph((e.source,e.target) for e in g.transitions))
        for state in g.states:
            m=dict(state.mapping)
            assert len(set(m.values()))==len(m)
            assert all(p.reactant.elements[r]==p.product.elements[t] for r,t in m.items())
        for edge in g.transitions:
            old=dict(g.states[edge.source].mapping);m=dict(g.states[edge.target].mapping)
            assert old.items()<=m.items()
            if edge.match is None:continue
            for a,b in edge.preserved_bonds:
                assert p.product.wbo[m[a],m[b]]>=result.config.graph_floor
                assert abs(p.reactant.wbo[a,b]-p.product.wbo[m[a],m[b]])<=result.config.iso_tolerance+1e-9
                edges+=1
        record=dict(label=row['label'],states=len(g.states),edges_checked=edges,valid=True,
                    validation_seconds=time.perf_counter()-start)
        results.append(record);save(folder/'validation.json',results)
    # Optional bounded analysis of correlated alternatives on the final snapshot.
    # Order representatives only for evaluation; never prune search families.
    start=time.perf_counter();units=[]
    terminals=[t for t in g.terminals if len(g.states[t].mapping)==p.atom_count]
    vectors=[tuple(dict(g.states[t].mapping)[i] for i in range(p.atom_count)) for t in terminals]
    scores=event_counts(p.reactant.wbo,p.product.wbo,vectors) if vectors else []
    seen=set()
    for _,terminal in sorted(zip((int(sum(v)) for v in scores),terminals)):
        if time.perf_counter()-start>=args.seconds:break
        path=next(g.paths(terminal))
        key=tuple(sorted(path.mapping.items()))
        if key in seen:continue
        seen.add(key)
        score=minimize_events(path,p,seconds=min(.5,args.seconds-(time.perf_counter()-start)))
        units.append(dict(terminal=terminal,transitions=path.transitions,score=score))
        save(folder/'family_diagnostic.json',dict(units=units,seconds=time.perf_counter()-start,
            best=min(u['score']['upper_bound'] for u in units),
            scope='Budgeted diagnostic: one conditional path per selected terminal representative; not exhaustive family traversal or a proof of global optimality. All correlations retained inside each tested path.'))
    print(json.dumps(dict(slot=args.slot,validated=len(results),family_paths=len(units),
        best=min((u['score']['upper_bound'] for u in units),default=None))),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--slot',type=int,required=True)
    parser.add_argument('--seconds',type=float,default=15.)
    audit(parser.parse_args())
