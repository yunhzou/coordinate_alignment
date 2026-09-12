"""Fully decode a saved AAM archive, with exact certificates and resumable progress.

No AAM search runs here. Each invocation has a process watchdog. Omitting
--max-events requests all event counts; otherwise completeness is windowed.
"""
import argparse
from dataclasses import asdict
import math
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'bench')]
from run_event_campaign import read, save, sha, MemoryGuard, bounded_process


def verified_resume(folder, digest, policy, window):
    record = read(folder / 'comparison.json')
    previous_digest = record.get('archive_sha256') or record.get('resume', {}).get('archive_sha256')
    if previous_digest is None and (folder / 'search.json').exists():
        previous_digest = read(folder / 'search.json')['archive_sha256']
    if previous_digest != digest or record['event_policy'] != policy or record['max_events'] != window:
        raise ValueError('Resume archive, event policy, or event window differs')
    return record


def same_config(recorded, config):
    # JSON preserves sequence contents, but not Python tuple/list distinctions.
    return recorded == json.loads(json.dumps(asdict(config)))


def certificate_rows(path):
    with path.open() as stream:
        for line in stream:
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                if line.endswith('\n'):
                    raise
                # A watchdog can interrupt the final line; only complete
                # certificates before that line are reusable.
                break


def child(args):
    import json
    from rxn_core.artifacts import read_aam_checkpoint
    from rxn_core.event_patterns import SignedEventIndex, extract_path_events
    from validate_event_examples import compare_saved_slap
    start, cpu = time.perf_counter(), time.process_time()
    deadline = start + max(0., args.seconds - 2.)
    digest = sha(args.archive)
    aam = read_aam_checkpoint(args.archive)
    index = SignedEventIndex(aam.problem)
    completed, patterns, record = set(), {}, {}
    for folder in args.resume:
        previous = verified_resume(folder, digest, index.policy, args.max_events)
        if not same_config(previous.get('config'), aam.config):
            raise ValueError('Resume search configuration differs')
        record.update({k: v for k, v in previous.items() if k not in ('path_results', 'resume', 'patterns')})
        for pattern in previous['patterns'].values():
            checked = index.describe(pattern['mapping'])
            if checked['id'] != pattern['id'] or checked['total'] != pattern['total']:
                raise ValueError('Resume event witness differs from its recorded identity or count')
            if args.max_events is not None and checked['total'] > args.max_events:
                raise ValueError('Resume event witness lies outside the requested window')
            patterns.setdefault(pattern['id'], dict(pattern, **checked))
        for row in previous.get('path_results', []):
            if row['complete']:
                completed.add((row['terminal'], tuple(row['transitions'])))
        if (folder / 'families.jsonl').exists():
            for row in certificate_rows(folder / 'families.jsonl'):
                if row['complete']:
                    completed.add((row['terminal'], tuple(row['transitions'])))
    record.update(archive_sha256=digest, config=asdict(aam.config), event_policy=index.policy,
                  max_events=args.max_events, search_capped=aam.graph.capped)
    initial_ids = set(patterns)
    visited = skipped = finished = incomplete = queries = 0
    active, last_save, reasons = None, 0., {}

    def checkpoint(force=False, complete=False, reason='running'):
        nonlocal last_save
        if not force and time.perf_counter() - last_save < 3:
            return
        record.update(patterns=patterns, saved_graph_window_complete=complete, reason=reason,
            best_found=min((p['total'] for p in patterns.values()), default=None),
            resume=dict(archive_sha256=digest, previous_directories=[str(p) for p in args.resume],
                paths_visited=visited, previously_complete_skipped=skipped, new_complete_families=finished,
                incomplete_families=incomplete, active_path=active, solver_queries=queries, reasons=reasons,
                additional_event_classes=len(patterns.keys() - initial_ids),
                wall_seconds=time.perf_counter() - start, cpu_seconds=time.process_time() - cpu,
                hard_watchdog_seconds=args.seconds, max_patterns=None, search_rerun=False))
        if 'slap' in record:
            record['comparisons'] = compare_saved_slap(record['slap'], patterns,
                float('inf') if args.max_events is None else args.max_events, complete)
        save(args.output / 'comparison.json', record)
        last_save = time.perf_counter()

    checkpoint(True)
    exhausted = True
    with (args.output / 'families.jsonl').open('w', buffering=1) as journal:
        for path in aam.graph.paths():
            if len(path.mapping) != index.n:
                continue
            visited += 1
            if time.perf_counter() >= deadline:
                exhausted = False
                break
            if (path.terminal, tuple(path.transitions)) in completed:
                skipped += 1
                checkpoint()
                continue
            active = dict(terminal=path.terminal, transitions=list(path.transitions))
            checkpoint()

            def on_pattern(pattern):
                if pattern['id'] not in patterns:
                    patterns[pattern['id']] = dict(pattern, terminal=path.terminal,
                        transitions=list(path.transitions), provenance='compressed_family_full_decode')
                    checkpoint(True)

            result = extract_path_events(path, aam.problem, index, max_events=args.max_events,
                seconds=max(0., deadline - time.perf_counter()), max_patterns=None, on_pattern=on_pattern)
            queries += result['solver_queries']
            reasons[result['reason']] = reasons.get(result['reason'], 0) + 1
            if result['complete']:
                finished += 1
            else:
                incomplete += 1
            journal.write(json.dumps({k: v for k, v in result.items() if k != 'patterns'}) + '\n')
            active = None
            checkpoint()
    complete = exhausted and incomplete == 0
    checkpoint(True, complete, 'all_retained_full_paths_checked' if complete
               else 'watchdog_budget' if not exhausted else 'solver_unknown')


def main(args):
    args.output.mkdir(parents=True, exist_ok=False)
    sources = [ROOT / 'src/rxn_core' / name for name in
               ('event_patterns.py', 'event_certificates.py', 'family_query.py')]
    save(args.output / 'manifest.json', dict(archive=str(args.archive), archive_sha256=sha(args.archive),
        max_events=args.max_events, hard_watchdog_seconds=args.seconds, workers=1,
        resume=[str(p) for p in args.resume], sources={str(p): sha(p) for p in sources + [Path(__file__)]},
        memory_limits_mib=dict(worker=3072, total=4096, reserve=6144), search_rerun=False))
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0')
    command = [sys.executable, str(Path(__file__).resolve()), '--child', '--archive', str(args.archive),
               '--output', str(args.output), '--seconds', str(args.seconds)]
    if args.max_events is not None:
        command += ['--max-events', str(args.max_events)]
    for folder in args.resume:
        command += ['--resume', str(folder)]
    start = time.perf_counter()
    with (args.output / 'decode.log').open('w') as log:
        status, peak = bounded_process(command, env, log, args.seconds, MemoryGuard(3072, 4096, 6144))
    result = dict(status=status, wall_seconds=time.perf_counter() - start, peak_sampled_rss_mib=peak, complete=False)
    if (args.output / 'comparison.json').exists():
        record = read(args.output / 'comparison.json')
        if status != 'passed':
            record['saved_graph_window_complete'] = False
            record['reason'] = status
            save(args.output / 'comparison.json', record)
        result.update(complete=record['saved_graph_window_complete'], resume=record['resume'])
    save(args.output / 'execution.json', result)
    print(__import__('json').dumps(result), flush=True)
    if not result['complete']:
        raise SystemExit(1)


def parse():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', type=Path, action='append', default=[])
    parser.add_argument('--max-events', type=int)
    parser.add_argument('--seconds', type=float, default=300)
    parser.add_argument('--child', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not math.isfinite(args.seconds) or args.seconds <= 0 or (args.max_events is not None and args.max_events < 0):
        parser.error('seconds must be positive and event window nonnegative')
    args.archive, args.output = args.archive.resolve(), args.output.resolve()
    args.resume = [p.resolve() for p in args.resume]
    return args


if __name__ == '__main__':
    args = parse()
    (child if args.child else main)(args)
