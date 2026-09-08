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
from rxn_core.artifacts import read_aam_checkpoint, raw_cut_paths, read_raw_cut, read_graph_checkpoint


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
        reference = json.loads((directory/'reference.json').read_text())
        archive = args.output/'cuts/aam.pkl.gz'
        if archive.exists():
            result = read_aam_checkpoint(archive)
            report = evaluate_planned(result, plan, reference['features'], reference['mapping'],
                                      seconds=220, query_timeout_ms=5000)
            report.update(search_incomplete=False, archive=str(archive.resolve()))
        else:
            # A watchdog may interrupt a later cut or parent merge. Saved cuts
            # can certify a positive; incomplete evidence cannot certify absence.
            from rxn_core.domain import AAMResult, AAMSearchMetrics
            from rxn_core.search_symmetry import finalize_graph_symmetry
            from rxn_core.frag import build_graph
            target = build_graph(plan.problem.product.elements, plan.problem.product.wbo, config.graph_floor)
            start = time.monotonic()
            checks = []
            report = dict(reference_recovery='unknown', search_incomplete=True,
                          reason='Full archive not completed')
            for raw in raw_cut_paths(args.output/'cuts'):
                left = 220-(time.monotonic()-start)
                if left <= 0:
                    break
                finalized = raw.with_name(raw.name.replace('.raw.', '.finalized.'))
                if finalized.exists():
                    graph = read_graph_checkpoint(finalized)
                else:
                    graph, _ = finalize_graph_symmetry(read_raw_cut(raw), target,
                                                        iso_tolerance=config.iso_tolerance)
                result = AAMResult(plan.problem, config, graph, AAMSearchMetrics.from_record({},0))
                check = evaluate_planned(result, plan, reference['features'], reference['mapping'],
                    seconds=max(.01,min(25,220-(time.monotonic()-start))), query_timeout_ms=3000)
                check.update(archive=str(raw.resolve()), search_incomplete=True)
                checks.append(check)
                save(args.output/'partial_checks.json', checks)
                if check['reference_recovery'] == 'recovered':
                    report = check
                    break
            report['checked_cuts'] = len(checks)
        report.update(direction=plan.direction)
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
