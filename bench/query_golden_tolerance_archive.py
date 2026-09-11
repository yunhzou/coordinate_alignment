"""Extend verification of one saved Golden output without rerunning mapping."""
import argparse
from pathlib import Path
import time
from types import SimpleNamespace

from golden_tolerance_benchmark import read, save, sha, result_folder


def query(args):
    from adaptive_full_benchmark import problem_plan
    from golden_evaluation import evaluate_planned
    from rxn_core.artifacts import read_aam_checkpoint
    variant = args.run / 'runs' / args.variant
    task = dict(dataset='golden', index=args.index, direction=args.direction)
    _, plan = problem_plan(SimpleNamespace(run=variant, method='original'), task)
    folder = result_folder(variant, task)
    archive = folder / 'cuts/aam.pkl.gz'
    reference_path = variant / f'inputs/golden/{args.index}/reference.json'
    reference = read(reference_path)
    saved = read_aam_checkpoint(archive)
    started = time.monotonic()
    result = evaluate_planned(saved, plan, reference['features'], reference['mapping'],
                              seconds=args.seconds, query_timeout_ms=1500)
    save(args.output, dict(index=args.index, direction=args.direction, variant=args.variant,
         search_rerun=False, archive=str(archive), archive_sha256=sha(archive),
         reference_sha256=sha(reference_path), driver_sha256=sha(Path(__file__)),
         verification_seconds=args.seconds, per_query_timeout_ms=1500,
         elapsed_verification=time.monotonic()-started, evaluation=result))
    print(result['reference_recovery'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--variant', choices=('tol_1p0', 'tol_1p5'), default='tol_1p5')
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--direction', choices=('R_to_P', 'P_to_R'), required=True)
    parser.add_argument('--seconds', type=float, default=240)
    parser.add_argument('--output', type=Path, required=True)
    query(parser.parse_args())
