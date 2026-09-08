"""Explicit cut-set cause probe; reference-directed cuts are not blind scores."""
import argparse
from dataclasses import replace,asdict
import json
from pathlib import Path
import time
from golden_policy_campaign import load_case,save
from golden_evaluation import evaluate_planned
from rxn_core.aam import _initialize_search,_search_cut
from rxn_core.domain import AAMResult,AAMSearchMetrics
from rxn_core.artifacts import write_aam_checkpoint
from rxn_core.frag import build_graph
from rxn_core.search_symmetry import finalize_graph_symmetry


def main(args):
    args.output.mkdir(parents=True,exist_ok=False)
    directory,plan=load_case(args.source,1285)
    config=replace(plan.config,seed_count=10,branch_limit=2000)
    plan=replace(plan,config=config)
    reference=json.loads((directory/'reference.json').read_text())
    rows=[]
    cut_sets=[tuple(map(tuple,args.cut))] if args.cut else [(),((7,8),),((9,10),),((7,8),(9,10))]
    for cuts in cut_sets:
        out=args.output/str(len(rows));out.mkdir()
        start=time.perf_counter();_initialize_search(plan.problem,config)
        graph,counts=_search_cut(cuts)
        graph,_=finalize_graph_symmetry(graph,build_graph(plan.problem.product.elements,
            plan.problem.product.wbo,config.graph_floor),iso_tolerance=config.iso_tolerance)
        result=AAMResult(plan.problem,config,graph,AAMSearchMetrics.from_record(counts,time.perf_counter()-start))
        write_aam_checkpoint(result,out/'aam.pkl.gz')
        search_seconds=time.perf_counter()-start
        evaluation=evaluate_planned(result,plan,reference['features'],reference['mapping'],seconds=45,query_timeout_ms=3000)
        save(out/'evaluation.json',evaluation)
        rows.append(dict(cuts=cuts,search_seconds=search_seconds,capped=graph.capped,
            terminals=len(graph.terminals),outcome=evaluation['reference_recovery']))
        save(args.output/'summary.json',dict(reference_directed_cut_selection=True,index=1285,
            direction=plan.direction,config=asdict(config),rows=rows))
        print(json.dumps(rows[-1]),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--cut',type=int,nargs=2,action='append')
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
