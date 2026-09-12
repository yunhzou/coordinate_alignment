"""Bounded local replay of saved event witnesses; never reruns AAM search.

Two single-threaded workers by default, always below the logical CPU count.
Each worker has a process timeout. Existing inputs and reports are read-only.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'bench'))


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + '\n')


def child(args):
    import pynauty
    from rxn_core import AAMProblem, MolecularEndpoint
    from rxn_core.event_patterns import SignedEventIndex
    from reference_dense_events import DeltaPatterns
    raw = json.loads((args.repo / f'manuscript/evidence/case{args.case}_input.json').read_text())
    catalogue = json.loads((args.repo / 'reports/holdout_minimum_event_patterns_20260910/per_case.json').read_text())
    vectors = [record['mapping'] for record in catalogue[args.case]['patterns'].values()]
    problem = AAMProblem(*(MolecularEndpoint(**raw[side]) for side in ('reactant', 'product')))
    if len(vectors) > 16:
        raise ValueError('pilot deliberately limits each case to 16 saved witnesses')
    samples = []
    reference_certificate = pynauty.certificate
    for repeat in range(2):
        for backend in (('legacy', 'sparse') if repeat == 0 else ('sparse', 'legacy')):
            canonical_cpu, calls = 0.0, 0
            def timed_certificate(graph):
                nonlocal canonical_cpu, calls
                start = time.process_time()
                try:
                    return reference_certificate(graph)
                finally:
                    canonical_cpu += time.process_time() - start
                    calls += 1
            pynauty.certificate = timed_certificate
            start, wall = time.process_time(), time.perf_counter()
            engine = DeltaPatterns(raw) if backend == 'legacy' else SignedEventIndex(problem)
            counts = engine.counts(vectors).tolist()
            rows = [engine.describe(vector) for vector in vectors]
            cpu_seconds, wall_seconds = time.process_time() - start, time.perf_counter() - wall
            # Group validation is separate from the replay timing.
            generators = engine.generators if backend == 'legacy' else engine.source_generators
            samples.append(dict(backend=backend, repeat=repeat, cpu_seconds=cpu_seconds,
                wall_seconds=wall_seconds, canonical_cpu_seconds=canonical_cpu, certificate_calls=calls,
                graph_vertices=engine.base.number_of_vertices, generators=generators,
                counts=counts, patterns=[dict(id=r['id'], total=r['total'], events=r['events']) for r in rows]))
            pynauty.certificate = reference_certificate
    record = dict(index=args.case, atoms=problem.source_atom_count, witnesses=len(vectors),
                  samples=samples, mapping_search_rerun=False)
    # Include the exact group checks inside the same process watchdog.
    record['validation'] = validate_case(record)
    save(args.output, record)


def validate_case(record):
    from sympy.combinatorics import Permutation, PermutationGroup
    samples = record['samples']
    reference = samples[0]
    n = record['atoms']
    def group(sample):
        generators = sample['generators'] or [list(range(n))]
        return PermutationGroup([Permutation(g, size=n) for g in generators])
    reference_group = group(reference)
    def partition(rows):
        return [[a['id'] == b['id'] for b in rows] for a in rows]
    for sample in samples:
        assert sample['counts'] == reference['counts']
        assert [(p['total'], p['events']) for p in sample['patterns']] == [(p['total'], p['events']) for p in reference['patterns']]
        assert partition(sample['patterns']) == partition(reference['patterns'])
        candidate_group = group(sample)
        assert all(candidate_group.contains(g) for g in reference_group.generators)
        assert all(reference_group.contains(g) for g in candidate_group.generators)
    return dict(event_counts_equal=True, literal_events_equal=True,
                equivalence_partition_equal=True, exact_source_groups_equal=True,
                group_verification='mutual generator containment; no group-element enumeration')


def main(args):
    logical = os.cpu_count() or 1
    if logical < 2:
        raise RuntimeError('cannot reserve one logical core on a one-core host')
    workers = min(max(1, args.workers), logical - 1)
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    selected = (1, 64, 135)
    def run(index):
        output = args.output / f'case{index}.json'
        command = [sys.executable, str(Path(__file__).resolve()), '--child', '--repo', str(args.repo),
                   '--case', str(index), '--output', str(output)]
        try:
            result = subprocess.run(command, env=env, text=True, capture_output=True, timeout=args.timeout)
        except subprocess.TimeoutExpired:
            return dict(index=index, status='timeout', timeout_seconds=args.timeout)
        if result.returncode:
            return dict(index=index, status='failed', stderr=result.stderr[-4000:])
        record = json.loads(output.read_text())
        times = {backend: sum(s['cpu_seconds'] for s in record['samples'] if s['backend'] == backend) / 2
                 for backend in ('legacy', 'sparse')}
        return dict(index=index, status='passed', atoms=record['atoms'], witnesses=record['witnesses'],
                    mean_cpu_seconds=times, sparse_over_legacy_cpu=times['sparse'] / times['legacy'],
                    graph_vertices={s['backend']: s['graph_vertices'] for s in record['samples']},
                    validation=record['validation'])
    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        cases = list(pool.map(run, selected))
    source_files = [ROOT / 'src/rxn_core/event_patterns.py', ROOT / 'bench/reference_dense_events.py', Path(__file__)]
    inputs = [args.repo / f'manuscript/evidence/case{i}_input.json' for i in selected]
    inputs.append(args.repo / 'reports/holdout_minimum_event_patterns_20260910/per_case.json')
    summary = dict(logical_cpus=logical, worker_limit=workers, logical_cpus_not_assigned_to_workers=logical - workers,
        numerical_threads_per_worker=1, process_timeout_seconds=args.timeout,
        elapsed_seconds=time.perf_counter() - start, mapping_searches=0,
        mapping_bijections_enumerated=0, scope='Three selected saved-input cases, two paired repeats; postprocessing only.',
        sources={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
        inputs={str(p.relative_to(args.repo)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        versions={package: version(package) for package in ('numpy', 'pynauty', 'sympy', 'z3-solver', 'rdkit')},
        python=sys.version, cases=cases)
    save(args.output / 'summary.json', summary)
    print(json.dumps(summary, indent=2))
    if any(case['status'] != 'passed' for case in cases):
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True, help='Checkout containing saved manuscript/report inputs')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=20.0)
    parser.add_argument('--child', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--case', type=int, help=argparse.SUPPRESS)
    options = parser.parse_args()
    if not math.isfinite(options.timeout) or options.timeout <= 0:
        parser.error('timeout must be positive')
    (child if options.child else main)(options)
