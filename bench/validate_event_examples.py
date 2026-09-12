"""Small real AAM/SLAP comparison with hard per-process limits.

Fresh forward AAM searches, saved SLAP witnesses rescored under one raw-WBO
policy, and bounded compressed-family decoding. No permutation expansion.
This is a selected-case regression experiment, not a chemical accuracy study.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'bench')]


def read(path):
    with gzip.open(path, 'rt') if path.suffix == '.gz' else path.open() as stream:
        return json.load(stream)


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2) + '\n')
    tmp.replace(path)


def compare_saved_slap(slap, patterns, window, complete):
    comparisons = {}
    for method, record in slap.items():
        wanted = set(record['minimum_ids'])
        missing = wanted - patterns.keys()
        comparisons[method] = dict(slap_saved_minimum=record['minimum'],
            minimum_pattern_count=len(wanted), recovered=sorted(wanted & patterns.keys()),
            unrecovered=sorted(missing), all_recovered=not missing,
            unrecovered_status={key: ('outside_event_window' if record['patterns'][key]['total'] > window
                else 'absent_from_this_saved_forward_graph' if complete else 'unresolved') for key in sorted(missing)})
    return comparisons


def finish_saved(args, aam, canonical, folder):
    """Retry only unfinished families and unvisited paths from the same archive."""
    from rxn_core.event_patterns import extract_path_events
    record = read(folder / 'comparison.json')
    if (folder / 'comparison_initial.json').exists():
        raise ValueError('finish is a single bounded follow-up, not an unlimited retry loop')
    save(folder / 'comparison_initial.json', record)
    start, cpu = time.perf_counter(), time.process_time()
    deadline = start + 5.0
    rows, patterns = record['path_results'], record['patterns']
    visited = 0
    reached_end = True
    for path in aam.graph.paths():
        if len(path.mapping) != canonical.n:
            continue
        ordinal = visited
        visited += 1
        if ordinal < len(rows) and rows[ordinal]['complete']:
            continue
        if time.perf_counter() >= deadline or visited > 5000:
            reached_end = False
            break
        result = extract_path_events(path, aam.problem, canonical, max_events=record['max_events'],
                                     seconds=min(2., deadline - time.perf_counter()), max_patterns=32)
        record['solver_queries'] += result['solver_queries']
        for pattern in result['patterns']:
            patterns.setdefault(pattern['id'], dict(pattern, provenance='compressed_path',
                terminal=path.terminal, transitions=list(path.transitions)))
        row = {k: v for k, v in result.items() if k != 'patterns'}
        if ordinal < len(rows):
            rows[ordinal] = row
        else:
            rows.append(row)
    record['path_count'] = len(rows)
    record['saved_graph_window_complete'] = reached_end and all(row['complete'] for row in rows)
    record['best_found'] = min((p['total'] for p in patterns.values()), default=None)
    record['additional_classes_from_families'] = len(patterns) - record['representative_classes_in_window']
    record['comparisons'] = compare_saved_slap(record['slap'], patterns, record['max_events'], record['saved_graph_window_complete'])
    record['finish'] = dict(cpu_seconds=time.process_time() - cpu, wall_seconds=time.perf_counter() - start,
                           mapping_search_rerun=False)
    record['final_source_sha256'] = hashlib.sha256((ROOT / 'src/rxn_core/event_patterns.py').read_bytes()).hexdigest()
    save(folder / 'comparison.json', record)


def child(args):
    import random
    import numpy as np
    from rxn_core import AAMProblem, AAMSearchConfig, MolecularEndpoint, search_aam
    from rxn_core.artifacts import read_aam_checkpoint, write_aam_checkpoint
    from rxn_core.event_patterns import SignedEventIndex, extract_path_events
    raw = read(args.repo / f'manuscript/evidence/case{args.case}_input.json')
    problem = AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant', 'product')))
    folder = args.output / f'case{args.case}'
    folder.mkdir(parents=True, exist_ok=True)
    archive = folder / 'aam.pkl.gz'
    if args.phase == 'search':
        random.seed(42)
        np.random.seed(42)
        config = AAMSearchConfig(**read(args.repo / 'reports/holdout_cap1000_seed1_20260910/manifest.json')['original_config'])
        start, cpu = time.perf_counter(), time.process_time()
        aam = search_aam(problem, config, execution='reused_native', workers=1)
        wall, elapsed_cpu = time.perf_counter() - start, time.process_time() - cpu
        write_aam_checkpoint(aam, archive)
        save(folder / 'search.json', dict(case=args.case, direction='R_to_P', config=asdict(config),
            cpu_seconds=elapsed_cpu, wall_seconds=wall, capped=aam.graph.capped,
            terminals=len(aam.graph.terminals), metrics=asdict(aam.metrics),
            archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest()))
        return
    aam = read_aam_checkpoint(archive)
    canonical = SignedEventIndex(problem)
    if args.phase == 'finish':
        finish_saved(args, aam, canonical, folder)
        return
    start, cpu = time.perf_counter(), time.process_time()
    sources = {
        'native_slap': ('reports/holdout_cap1000_seed1_20260910/slap_scoring.json.gz', str(args.case)),
        'slap_sweep': ('reports/holdout_slap_xyz_sweep_20260910/scored_candidates.json.gz', f'{args.case}/scored')}
    slap = {}
    for method, (filename, key) in sources.items():
        rows = read(args.repo / filename)[key]['slap']
        if len(rows) > 512:
            raise ValueError('bounded pilot allows at most 512 saved SLAP witnesses per method')
        patterns = {}
        for row in rows:
            mapping = dict(row['mapping'])
            p = canonical.describe([mapping[a] for a in range(canonical.n)])
            patterns.setdefault(p['id'], p)
        minimum = min(p['total'] for p in patterns.values())
        slap[method] = dict(saved_witnesses=len(rows), minimum=minimum,
                            patterns=patterns, minimum_ids=sorted(k for k, p in patterns.items() if p['total'] == minimum))
    window = max(v['minimum'] for v in slap.values()) + 1
    # Score every retained representative; no path or atom-permutation expansion.
    vectors = sorted({tuple(dict(aam.graph.states[t].mapping)[a] for a in range(canonical.n))
                      for t in aam.graph.terminals if len(aam.graph.states[t].mapping) == canonical.n})
    patterns = {}
    for offset in range(0, len(vectors), 64):
        batch = vectors[offset:offset + 64]
        for vector, counts in zip(batch, canonical.counts(batch)):
            if int(counts.sum()) <= window:
                p = canonical.describe(vector)
                patterns.setdefault(p['id'], dict(p, provenance='saved_terminal'))
    representative_ids = set(patterns)
    paths, complete, queries = [], True, 0
    deadline = time.perf_counter() + 5.0
    for path in aam.graph.paths():
        if len(path.mapping) != canonical.n:
            continue
        if time.perf_counter() >= deadline or len(paths) >= 1500:
            complete = False
            break
        result = extract_path_events(path, aam.problem, canonical, max_events=window,
                                     seconds=min(.2, max(0., deadline - time.perf_counter())), max_patterns=32)
        complete &= result['complete']
        queries += result['solver_queries']
        for p in result['patterns']:
            patterns.setdefault(p['id'], dict(p, provenance='compressed_path',
                                            terminal=path.terminal, transitions=list(path.transitions)))
        paths.append({k: v for k, v in result.items() if k != 'patterns'})
    comparisons = compare_saved_slap(slap, patterns, window, complete)
    historical = read(args.repo / 'reports/holdout_minimum_event_patterns_20260910/aam_enumeration.json.gz')[str(args.case)]
    save(folder / 'comparison.json', dict(case=args.case, atoms=canonical.n,
        direction='R_to_P', event_policy=canonical.policy, max_events=window,
        full_representative_count=len(vectors), representative_classes_in_window=len(representative_ids),
        patterns=patterns, additional_classes_from_families=len(patterns.keys() - representative_ids),
        best_found=min((p['total'] for p in patterns.values()), default=None),
        saved_graph_window_complete=complete, search_capped=aam.graph.capped,
        path_count=len(paths), solver_queries=queries, path_results=paths,
        comparisons=comparisons, slap=slap,
        cpu_seconds=time.process_time() - cpu, wall_seconds=time.perf_counter() - start,
        historical=dict(complete=historical['complete'], seconds=historical['seconds'],
                        event_policy='historical floor-based four-kind events; NOT directly comparable timing'),
        scope='Exact signed events on fresh retained forward AAM families versus saved SLAP mapping witnesses. '
              'SLAP label families are not fully decoded. A capped search or a bounded decode cannot prove global absence.'))


def main(args):
    if (os.cpu_count() or 1) < 2:
        raise RuntimeError('two logical CPUs required to leave one available')
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONHASHSEED='0', PYTHONDONTWRITEBYTECODE='1')
    def run(case):
        phases = {}
        for phase in ('search', 'decode'):
            command = [sys.executable, str(Path(__file__).resolve()), '--repo', str(args.repo),
                       '--output', str(args.output), '--case', str(case), '--phase', phase]
            start = time.perf_counter()
            try:
                result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=12)
                phases[phase] = dict(status='passed' if result.returncode == 0 else 'failed',
                                     wall_seconds=time.perf_counter() - start, stderr=result.stderr[-3000:])
            except subprocess.TimeoutExpired:
                phases[phase] = dict(status='timeout', wall_seconds=time.perf_counter() - start)
            if phases[phase]['status'] != 'passed':
                break
        return dict(case=case, phases=phases)
    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=min(2, (os.cpu_count() or 1) - 1)) as pool:
        rows = list(pool.map(run, (1, 64, 135)))
    files = [ROOT / 'src/rxn_core/event_patterns.py', Path(__file__)]
    summary = dict(cases=rows, elapsed_seconds=time.perf_counter() - start, worker_limit=2,
                   numerical_threads_per_worker=1, process_timeout_seconds=12,
                   sources={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    save(args.output / 'summary.json', summary)
    print(json.dumps(summary, indent=2))
    if any(p['status'] != 'passed' for row in rows for p in row['phases'].values()):
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', type=int)
    parser.add_argument('--phase', choices=('search', 'decode', 'finish'))
    args = parser.parse_args()
    (child if args.phase else main)(args)
