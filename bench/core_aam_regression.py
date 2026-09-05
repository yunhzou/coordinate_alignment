#!/usr/bin/env python3
"""Replay identical saved full-AAM inputs, keeping search and validation separate."""
import argparse
import cProfile
import gc
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import pickle
import platform
import pstats
import resource
import subprocess
import sys
import time


def normalized(states):
    result = []
    for mapping, islands, cuts in states:
        groups = {}
        for atom, label in islands:
            groups.setdefault(label, []).append(atom)
        result.append((tuple(map(tuple, mapping)),
                       tuple(sorted(tuple(sorted(atoms)) for atoms in groups.values())),
                       tuple(map(tuple, cuts))))
    return sorted(result)


def replay(args):
    sys.path.insert(0, str(args.package_root / 'src'))
    from rxn_core.alignment.branch import find_islands
    data = (args.inputs / f'{args.family}.input.pkl').read_bytes()
    positional, keywords = pickle.loads(data)
    profiler = cProfile.Profile() if args.profile else None
    collections = []
    collection_start = 0.
    def gc_time(phase, info):
        nonlocal collection_start
        if phase == 'start':
            collection_start = time.perf_counter()
        else:
            collections.append((info['generation'], time.perf_counter() - collection_start))
    gc.callbacks.append(gc_time)
    if profiler:
        profiler.enable()
    start = time.perf_counter()
    graph = find_islands(*positional, **keywords)
    elapsed = time.perf_counter() - start
    gc.callbacks.remove(gc_time)
    search_peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    if profiler:
        profiler.disable()
        profiler.dump_stats(str(args.output.with_suffix('.prof')))
        with args.output.with_suffix('.profile.txt').open('w') as stream:
            pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats('cumulative').print_stats(70)
    # Neither loading inputs nor checking/persisting results is timed as AAM.
    if args.legacy:
        states = [(tuple(sorted(b.mapping.items())), tuple(sorted(b.islands_R.items())),
                   tuple(sorted(b.deferred_edges))) for b in graph]
        exact_graph = None
        details = dict(branches=len(graph))
    else:
        states = [(graph.states[t].mapping, graph.states[t].islands,
                   graph.states[t].deferred_edges) for t in graph.terminals]
        with (args.inputs / f'{args.family}.graph.pkl').open('rb') as stream:
            expected = pickle.load(stream)
        exact_graph = graph == expected
        details = dict(branches=len(graph.terminals), states=len(graph.states),
                       transitions=len(graph.transitions), capped=graph.capped)
        if args.save_graph or not exact_graph:
            with args.output.with_suffix('.graph.pkl').open('wb') as stream:
                pickle.dump(graph, stream, protocol=5)
    old = json.loads((args.inputs / f'{args.family}.old.states.json').read_text())
    record = dict(family=args.family, seconds=elapsed,
                  hostname=platform.node(),
                  input_sha256=hashlib.sha256(data).hexdigest(),
                  terminal_relations_equal=normalized(states) == normalized(old),
                  full_graph_equal=exact_graph,
                  search_peak_mb=search_peak_mb,
                  gc_seconds=sum(seconds for _, seconds in collections),
                  gc_collections=len(collections),
                  validation_peak_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
                  **details)
    args.output.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record), flush=True)
    if not record['terminal_relations_equal'] or exact_graph is False:
        raise SystemExit(1)


def suite(args):
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.inputs / 'inputs.json').read_text())
    families = args.families or [f['family'] for f in manifest['families']]
    jobs = [(family, repeat) for repeat in range(args.repeats) for family in families]

    def pair(job):
        family, repeat = job
        rows = []
        # Alternate order, in the same worker allocation, to avoid warm-up bias.
        versions = ['old', 'current', 'before'] if args.before else ['old', 'current']
        offset = (family + repeat) % len(versions)
        versions = versions[offset:] + versions[:offset]
        for version in versions:
            output = args.output / f'{family}.{repeat}.{version}.json'
            root = dict(old=args.baseline, current=Path.cwd(), before=args.before)[version]
            command = [sys.executable, str(Path(__file__).resolve()), 'replay',
                       '--inputs', str(args.inputs), '--family', str(family),
                       '--package-root', str(root), '--output', str(output)]
            if version == 'old':
                command.append('--legacy')
            if args.profile:
                command.append('--profile')
            with output.with_suffix('.log').open('w') as stream:
                completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                           timeout=300, env=os.environ.copy())
            if completed.returncode:
                raise RuntimeError(f'Failed replay: {output.with_suffix(".log")}')
            row = dict(json.loads(output.read_text()), version=version, repeat=repeat)
            rows.append(row)
            print(json.dumps(row), flush=True)
        return rows

    with ThreadPoolExecutor(args.workers) as pool:
        records = [row for rows in pool.map(pair, jobs) for row in rows]
    (args.output / 'results.json').write_text(json.dumps(records, indent=2) + '\n')


def workflows(args):
    """Exercise full multi-seed/cut workflows, including optional grouping."""
    args.output.mkdir(parents=True, exist_ok=True)
    def pair(case):
        rows = {}
        for version, root in [('before', args.baseline), ('current', Path.cwd())]:
            output = args.output / f'{case}.{version}.json'
            command = [sys.executable, str(Path(__file__).with_name('search_graph_regression.py')),
                       '--case', case, '--output', str(output), '--workers', '1',
                       '--package-root', str(root), '--repeats', str(args.repeats)]
            with output.with_suffix('.log').open('w') as stream:
                subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                               timeout=300, check=True)
            rows[version] = json.loads(output.read_text())
        graphs = [json.loads((args.output / f'{case}.{v}.aam.json').read_text())['graph']
                  for v in ('before', 'current')]
        row = dict(case=case, full_graph_equal=graphs[0] == graphs[1],
                   mechanisms_equal=rows['before']['mechanisms'] == rows['current']['mechanisms'],
                   before_seconds=rows['before']['all_seconds'],
                   current_seconds=rows['current']['all_seconds'])
        print(json.dumps(row), flush=True)
        if not row['full_graph_equal'] or not row['mechanisms_equal']:
            raise RuntimeError(f'Workflow regression: {case}')
        return row
    with ThreadPoolExecutor(args.workers) as pool:
        records = list(pool.map(pair, args.cases))
    (args.output / 'results.json').write_text(json.dumps(records, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    for name in ('replay', 'suite', 'workflows'):
        command = sub.add_parser(name)
        if name != 'workflows':
            command.add_argument('--inputs', type=Path, required=True)
        command.add_argument('--output', type=Path, required=True)
        command.add_argument('--profile', action='store_true')
        if name == 'replay':
            command.add_argument('--family', type=int, required=True)
            command.add_argument('--package-root', type=Path, required=True)
            command.add_argument('--legacy', action='store_true')
            command.add_argument('--save-graph', action='store_true')
        else:
            command.add_argument('--baseline', type=Path, required=True)
            command.add_argument('--before', type=Path)
            command.add_argument('--families', type=int, nargs='+')
            command.add_argument('--workers', type=int, default=4)
            command.add_argument('--repeats', type=int, default=3)
            if name == 'workflows':
                command.add_argument('--cases', nargs='+', default=['tempo', 'tetraphenyl', 'tetratbu'])
    args = parser.parse_args()
    dict(replay=replay, suite=suite, workflows=workflows)[args.mode](args)


if __name__ == '__main__':
    main()
