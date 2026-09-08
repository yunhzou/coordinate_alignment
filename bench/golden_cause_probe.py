"""Reference-directed feasibility/fragment diagnostics, never blind accuracy."""
import argparse
import json
from pathlib import Path
import time

from golden_policy_campaign import load_case,save
from rxn_core.alignment.branch import find_islands,_generate_seed_orders
from rxn_core.artifacts import write_graph_checkpoint
from rxn_core.frag import build_graph
from rxn_core.matcher.policy import AttributeNodeMatchPolicy


def main(args):
    args.output.mkdir(parents=True,exist_ok=False)
    directory,plan=load_case(args.source,args.index)
    raw=json.loads((directory/'reference.json').read_text())
    reference=plan.to_search_mapping(dict(raw['mapping']))
    graphs=[build_graph(e.elements,e.wbo,plan.config.graph_floor)
            for e in (plan.problem.reactant,plan.problem.product)]
    for side,g in enumerate(graphs):
        labels=reference if side==0 else {p:r for r,p in reference.items()}
        for atom in g:
            # Equal labels identify only the annotated heavy-atom pairs.
            label=(reference[atom] if side==0 else atom) if atom in labels else f'{side}:{atom}'
            g.nodes[atom]['reference_label']='H' if g.nodes[atom]['element']=='H' else str(label)
    order=_generate_seed_orders(graphs[0],1)[0]
    rows=[]
    for policy in ('ordinary','reference_constrained'):
        events=[];start=time.perf_counter()
        graph=find_islands(*graphs,order,iso_tol=plan.config.iso_tolerance,
            graph_floor=plan.config.graph_floor,max_branches=args.cap,
            node_policy=AttributeNodeMatchPolicy(('element','reference_label')) if policy!='ordinary' else None,
            events=events)
        write_graph_checkpoint(graph,args.output/f'{policy}.pkl.gz')
        save(args.output/f'{policy}.events.json',events)
        hits=[t for t in graph.terminals if all(dict(graph.states[t].mapping).get(r)==p for r,p in reference.items())]
        row=dict(policy=policy,seconds=time.perf_counter()-start,terminals=len(graph.terminals),capped=graph.capped,
            exact_reference_witnesses=hits,
            best_exact_pairs=max((sum(dict(graph.states[t].mapping).get(r)==p for r,p in reference.items()) for t in graph.terminals),default=0))
        rows.append(row);save(args.output/'results.json',dict(index=args.index,direction=plan.direction,
            reference_directed=True,reference=sorted(reference.items()),seed_order=order,cap=args.cap,rows=rows))
        print(json.dumps(row),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--index',type=int,required=True);p.add_argument('--cap',type=int,default=2000)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
