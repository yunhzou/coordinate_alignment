"""Controlled diagnostic search; saves a separate full archive, never replaces scores."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import time

from golden_evaluation import evaluate
from investigate_golden_mapping import save
from rxn_core import AAMProblem,AAMSearchConfig,search_aam
from rxn_core.domain import MolecularEndpoint


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seeds',type=int,required=True)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--cap',type=int)
    p.add_argument('--tolerance',type=float)
    p.add_argument('--cut',nargs=2,type=int,action='append',
                   help='Diagnostic explicit cut set; replaces the ordinary sweep, not a blind result')
    p.add_argument('--seed-prefix',nargs='+',type=int,
                   help='Reference-directed diagnostic order; requires --cut')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    raw=json.loads((args.source/'input.json').read_text())
    problem=AAMProblem(MolecularEndpoint(**raw['reactant']),MolecularEndpoint(**raw['product']),raw['name'])
    config=AAMSearchConfig(**json.loads((args.source.parent/'manifest.json').read_text())['config'])
    config=replace(config,seed_count=args.seeds,**({'branch_limit':args.cap} if args.cap else {}))
    if args.tolerance is not None:config=replace(config,iso_tolerance=args.tolerance)
    started=time.perf_counter()
    if args.cut:
        from rxn_core.aam import _initialize_search,_search_cut
        from rxn_core.artifacts import write_aam_checkpoint
        from rxn_core.domain import AAMResult,AAMSearchMetrics
        from rxn_core.frag import build_graph
        from rxn_core.search_symmetry import finalize_graph_symmetry
        cuts=tuple(tuple(sorted(pair)) for pair in args.cut)
        save(args.output/'diagnostic_design.json',dict(reference_directed=True,cuts=cuts,
            seed_prefix=args.seed_prefix,
            note='Explicit cut-set diagnostic; not included in blind benchmark accuracy'))
        _initialize_search(problem,config)
        if args.seed_prefix:
            from rxn_core.alignment.branch import find_islands
            source=build_graph(problem.reactant.elements,problem.reactant.wbo,bond_cut=config.graph_floor)
            source.remove_edges_from(cuts)
            order=list(dict.fromkeys(args.seed_prefix+list(range(problem.source_atom_count))))
            graph=find_islands(source,build_graph(problem.product.elements,problem.product.wbo,
                bond_cut=config.graph_floor),order,graph_floor=config.graph_floor,
                iso_tol=config.iso_tolerance,max_branches=config.branch_limit,cuts=cuts)
            counts={}
        else:graph,counts=_search_cut(cuts)
        graph,groups=finalize_graph_symmetry(graph,build_graph(problem.product.elements,
            problem.product.wbo,bond_cut=config.graph_floor),iso_tolerance=config.iso_tolerance)
        counts.update(groups,cuts=1,raw_result_count=len(graph.terminals),
            retained_branch_count=len(graph.terminals))
        result=AAMResult(problem,config,graph,AAMSearchMetrics.from_record(counts,time.perf_counter()-started))
        (args.output/'cuts').mkdir()
        write_aam_checkpoint(result,args.output/'cuts/aam.pkl.gz')
    else:
        result=search_aam(problem,config,workers=args.workers,
            intermediate_dir=args.output/'cuts',archive_format='checkpoint')
    save(args.output/'search.json',dict(seconds=time.perf_counter()-started,metrics=vars(result.metrics)))
    reference=json.loads((args.source/'reference.json').read_text())
    value=evaluate(result,reference['features'],reference['mapping'],seconds=120)
    save(args.output/'evaluation.json',value)
    print(json.dumps(value),flush=True)


if __name__=='__main__':main()
