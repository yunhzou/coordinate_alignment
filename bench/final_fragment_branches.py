"""Export unordered final branches and optionally decode their flat families."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import time


def _read(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        return json.load(stream)


def _save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def run(args):
    from rxn_core import AAMProblem, MolecularEndpoint
    from rxn_core.artifacts import read_aam_checkpoint
    from rxn_core.final_branches import FinalBranchCatalogue
    from rxn_core.event_patterns import SignedEventIndex, extract_path_events
    start = time.perf_counter()
    args.output.mkdir(parents=True, exist_ok=True)
    inputs = {}
    if args.catalogue:
        raw = _read(args.problem)
        problem = AAMProblem(*(MolecularEndpoint(**raw[side]) for side in ('reactant', 'product')))
        catalogue = FinalBranchCatalogue.from_record(problem, _read(args.catalogue))
        inputs[str(args.catalogue)] = hashlib.sha256(args.catalogue.read_bytes()).hexdigest()
    else:
        catalogue = None
        for archive in args.archives:
            if time.perf_counter() - start > args.seconds - 5:
                raise TimeoutError('Final-branch construction exceeded its time budget')
            aam = read_aam_checkpoint(archive)
            if catalogue is None:
                catalogue = FinalBranchCatalogue(aam.problem)
            catalogue.add_aam(aam, str(archive.resolve()))
            inputs[str(archive)] = hashlib.sha256(archive.read_bytes()).hexdigest()
        problem = catalogue.problem
        catalogue.rebuild_all_symmetries()
    destination = args.output / 'final_branches.json.gz'
    with gzip.open(destination, 'wt') as stream:
        json.dump(catalogue.to_record(), stream, separators=(',', ':'))
    result = dict(branches=catalogue.branch_count, flat_families=len(catalogue.families),
                  input_paths=catalogue.path_count, incomplete_paths=catalogue.incomplete_paths,
                  intrinsic_fragment_groups=len(catalogue.symmetries),
                  problem_sha256=catalogue.problem_sha256, input_sha256=inputs,
                  decode_requested=args.decode, max_events=args.max_events,
                  complete_window=False, decoded_families=0, unfinished=[], patterns={})
    if args.decode:
        index = SignedEventIndex(problem)
        last_save = time.perf_counter()
        for fid, family in enumerate(catalogue.families):
            remaining = args.seconds - 5 - (time.perf_counter() - start)
            if remaining <= 0:
                break
            def record(pattern):
                result['patterns'].setdefault(pattern['id'], dict(pattern, family=fid))
            decoded = extract_path_events(family.as_path(problem), problem, index,
                max_events=args.max_events, seconds=min(args.family_seconds, remaining),
                max_patterns=None, on_pattern=record)
            result['decoded_families'] += 1
            if not decoded['complete']:
                result['unfinished'].append(dict(family=fid, reason=decoded['reason']))
            if time.perf_counter() - last_save >= 2:
                _save(args.output / 'summary.json', result)
                last_save = time.perf_counter()
        result['complete_window'] = (result['decoded_families'] == len(catalogue.families)
                                     and not result['unfinished'])
    result['wall_seconds'] = time.perf_counter() - start
    _save(args.output / 'summary.json', result)
    print(json.dumps({k: v for k, v in result.items() if k not in ('patterns', 'input_sha256')}))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--archives', nargs='+', type=Path)
    source.add_argument('--catalogue', type=Path)
    parser.add_argument('--problem', type=Path, help='Endpoint input JSON when reloading a flat catalogue')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--decode', action='store_true')
    parser.add_argument('--max-events', type=int)
    parser.add_argument('--seconds', type=float, default=300)
    parser.add_argument('--family-seconds', type=float, default=10)
    parser.add_argument('--child', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.catalogue and not args.problem:
        parser.error('--catalogue requires --problem')
    if args.decode and (args.max_events is None or args.max_events < 0):
        parser.error('--decode requires a nonnegative --max-events window')
    if args.seconds <= 5 or args.family_seconds <= 0:
        parser.error('time budgets must be positive, with --seconds greater than five')
    if args.child:
        run(args)
        return
    from run_event_campaign import MemoryGuard, bounded_process
    args.output.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    with (args.output / 'run.log').open('w') as log:
        status, peak = bounded_process([sys.executable, str(Path(__file__).resolve()),
            *sys.argv[1:], '--child'], env, log, args.seconds, MemoryGuard(3072, 6144, 6144))
    _save(args.output / 'execution.json', dict(status=status, peak_mib=peak,
                                              watchdog_seconds=args.seconds))
    if status != 'passed':
        raise SystemExit(f'Final-branch job stopped: {status}; inspect {args.output / "run.log"}')
    print(f'Wrote {args.output / "summary.json"}')


if __name__ == '__main__':
    main()
