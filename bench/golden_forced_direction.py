"""Saved, forced-direction sequential AAM experiment; core search is unchanged."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import time

from golden_policy_campaign import load_case, save
from golden_evaluation import evaluate_planned
from rxn_core import AAMProblem, search_aam
from rxn_core.search_orientation import AAMSearchPlan
from rxn_core.artifacts import read_aam_checkpoint


def main(args):
    directory, original = load_case(args.source, args.index)
    config = replace(original.config, seed_count=args.seeds, branch_limit=args.cap)
    reverse = args.direction == 'P_to_R'
    problem = original.input_problem
    search_problem = AAMProblem(problem.product, problem.reactant, problem.name) if reverse else problem
    if reverse:
        config = replace(config, anchors=tuple((p,r) for r,p in config.anchors))
    plan = AAMSearchPlan(problem, search_problem, config, reverse)
    if args.mode == 'search':
        args.output.mkdir(parents=True, exist_ok=False)
        save(args.output/'design.json', dict(index=args.index, source=str(directory),
            direction=plan.direction, config=asdict(config), workers=args.workers,
            note='Forced direction only; ordinary sequential AAM and single-cut sweep.'))
        start = time.perf_counter()
        result = search_aam(plan.problem, config, workers=args.workers,
            intermediate_dir=args.output/'cuts', archive_format='checkpoint')
        save(args.output/'search.json', dict(seconds=time.perf_counter()-start,
            metrics=asdict(result.metrics), capped=result.graph.capped,
            terminals=len(result.graph.terminals)))
    else:
        design = json.loads((args.output/'design.json').read_text())
        assert design['direction'] == plan.direction and design['config'] == json.loads(json.dumps(asdict(config)))
        result = read_aam_checkpoint(args.output/'cuts/aam.pkl.gz')
        reference = json.loads((directory/'reference.json').read_text())
        report = evaluate_planned(result, plan, reference['features'], reference['mapping'],
                                  seconds=220, query_timeout_ms=5000)
        report.update(direction=plan.direction, search_incomplete=False,
                      archive=str((args.output/'cuts/aam.pkl.gz').resolve()))
        save(args.output/'evaluation.json', report)
        print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('search','score'))
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--index', type=int, default=590)
    p.add_argument('--direction', choices=('R_to_P','P_to_R'), default='P_to_R')
    p.add_argument('--seeds', type=int, default=10)
    p.add_argument('--cap', type=int, default=100)
    p.add_argument('--workers', type=int, default=16)
    main(p.parse_args())
