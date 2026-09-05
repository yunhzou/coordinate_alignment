#!/usr/bin/env python3
"""Stage-isolated AAM memory diagnosis; no old retro baseline or search changes."""
import argparse
import csv
import gc
import gzip
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import subprocess
import sys
import time


def memory():
    with open('/proc/self/status') as stream:
        values = {line.split(':')[0]: line.split(':')[1].strip() for line in stream if ':' in line}
    return dict(rss_mb=int(values['VmRSS'].split()[0]) / 1024,
                peak_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)


def emit(stage, **data):
    print(json.dumps(dict(stage=stage, **memory(), **data)), flush=True)


def prepare(args):
    from rxn_core.smiles import smiles_to_weighted_graph
    from rxn_core.fragment_matching import FragmentDetectionConfig
    from rxn_core.subgraph import _coerce_graph
    from rxn_core.fragment_matching.detection import _initial_fragment_placements, _augment_initial_family
    import rxn_core.fragment_matching.augmentation as augmentation
    with gzip.open('data/inventory/processed/inventory_structure_bank.csv.gz', 'rt') as stream:
        row = next(r for r in csv.DictReader(stream) if r['Inventory ID'] == args.source_id)
    summary = json.loads(Path('data/retro_runs/native_exact_20260905/full_bank/parts/part_0.jsonl.gz.summary.json').read_text())
    config = FragmentDetectionConfig(**summary['config'])
    source = _coerce_graph(smiles_to_weighted_graph(row['SMILES'], expand_hydrogens=True), config.graph_floor)
    target = _coerce_graph(smiles_to_weighted_graph(summary['target_smiles'], expand_hydrogens=True), config.graph_floor)
    initial = _initial_fragment_placements(source, target, config)
    # Native graph caches are execution objects, not part of the molecular input.
    # Both versions rebuild this same cache from the saved graph/WBO data.
    source.graph.pop('_native_source', None)
    args.directory.mkdir(parents=True, exist_ok=True)
    class Captured(Exception): pass
    rows = []
    for i, placement in enumerate(initial[0]):
        # Capture the *actual* arguments supplied by the current consumer, then
        # stop before it runs AAM. Both core versions receive these exact bytes.
        def capture(*positional, **keywords):
            data = pickle.dumps((positional, keywords), protocol=5)
            (args.directory / f'{i}.input.pkl').write_bytes(data)
            rows.append(dict(family=i, retained=len(placement.retained_atoms),
                r_atoms=len(positional[0]), p_aug_atoms=len(positional[1]),
                input_sha256=hashlib.sha256(data).hexdigest()))
            raise Captured
        augmentation.find_islands = capture
        with (args.directory / f'{i}.family.pkl').open('wb') as stream:
            pickle.dump((source, target, placement, config), stream, protocol=5)
        try:
            _augment_initial_family(source, target, placement, config, None)
        except Captured:
            pass
    manifest = dict(source_id=args.source_id, smiles=row['SMILES'], families=rows)
    (args.directory / 'inputs.json').write_text(json.dumps(manifest, indent=2) + '\n')
    emit('prepared', **manifest)


def core(args):
    from rxn_core.alignment.branch import find_islands
    data = (args.directory / f'{args.family}.input.pkl').read_bytes()
    positional, keywords = pickle.loads(data)
    emit('core_start', family=args.family, version=args.version,
         input_sha256=hashlib.sha256(data).hexdigest())
    started = time.perf_counter()
    result = find_islands(*positional, **keywords)
    elapsed = time.perf_counter() - started
    if args.version == 'old':
        states = [(tuple(sorted(b.mapping.items())), tuple(sorted(b.islands_R.items())),
                   tuple(sorted(b.deferred_edges))) for b in result]
        extra = dict(branches=len(result), histories=sum(len(b.symmetry_paths) for b in result))
    else:
        states = [(result.states[t].mapping, result.states[t].islands,
                   result.states[t].deferred_edges) for t in result.terminals]
        extra = dict(branches=len(result.terminals), states=len(result.states),
                     transitions=len(result.transitions), capped=result.capped)
    metrics = dict(stage='core_end', version=args.version, family=args.family, seconds=elapsed,
        **memory(), **extra,
        state_sha256=hashlib.sha256(repr(sorted(states)).encode()).hexdigest())
    (args.directory / f'{args.family}.{args.version}.states.json').write_text(
        json.dumps(sorted(states), default=int) + '\n')
    print(json.dumps(metrics), flush=True)
    (args.directory / f'{args.family}.{args.version}.metrics.json').write_text(json.dumps(metrics) + '\n')
    if args.version == 'current':
        with (args.directory / f'{args.family}.graph.pkl').open('wb') as stream:
            pickle.dump(result, stream, protocol=5)
        emit('core_saved', family=args.family)


def post(args):
    import rxn_core.fragment_matching.augmentation as augmentation
    import rxn_core.fragment_matching.symmetry as symmetry
    from rxn_core.fragment_matching.detection import _augment_initial_family
    with (args.directory / f'{args.family}.family.pkl').open('rb') as stream:
        source, target, placement, config = pickle.load(stream)
    with (args.directory / f'{args.family}.graph.pkl').open('rb') as stream:
        graph = pickle.load(stream)
    # Reuse the already saved core result: this measurement performs no AAM.
    augmentation.find_islands = lambda *a, **k: graph
    finalize = augmentation.finalize_graph_symmetry
    orbit = symmetry.materialize_target_coverage_orbit
    serial = 0
    def measured_finalize(*a, **k):
        emit('symmetry_start')
        started = time.perf_counter()
        result = finalize(*a, **k)
        emit('symmetry_end', seconds=time.perf_counter()-started, **result[1])
        return result
    def measured_orbit(*a, **k):
        nonlocal serial
        serial += 1
        emit('occupation_start', path=serial)
        started = time.perf_counter()
        result = orbit(*a, **k)
        emit('occupation_end', path=serial, occupations=len(result), seconds=time.perf_counter()-started)
        return result
    augmentation.finalize_graph_symmetry = measured_finalize
    symmetry.materialize_target_coverage_orbit = measured_orbit
    emit('post_start', family=args.family)
    started = time.perf_counter()
    result = _augment_initial_family(source, target, placement, config, None)
    emit('post_end', family=args.family, candidates=len(result[0]), paths=serial,
         seconds=time.perf_counter()-started)
    # Save this evidence independently; do not collect other families in RAM.
    destination = args.output_directory or args.directory
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / f'{args.family}.post.pkl').open('wb') as stream:
        pickle.dump(result, stream, protocol=5)
    emit('post_saved', family=args.family)


def suite(args):
    from concurrent.futures import ThreadPoolExecutor
    manifest = json.loads((args.directory / 'inputs.json').read_text())
    versions = (args.version,) if args.version else ('old', 'current')
    jobs = [(family['family'], version) for family in manifest['families'] for version in versions]
    def run(job):
        family, version = job
        env = dict(os.environ, PYTHONPATH=str(args.baseline / 'src') if version == 'old' else str(Path('src').resolve()))
        command = [sys.executable, __file__, 'core', '--directory', str(args.directory),
                   '--family', str(family), '--version', version]
        with (args.directory / f'{family}.{version}.log').open('w') as stream:
            try:
                result = subprocess.run(command, env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=300)
                status = result.returncode
            except subprocess.TimeoutExpired:
                status = 'diagnostic_timeout_300s'
        print(json.dumps(dict(family=family, version=version, exit=status)), flush=True)
    with ThreadPoolExecutor(args.workers) as pool:
        list(pool.map(run, jobs))


def audit(args):
    from collections import Counter, defaultdict
    with (args.directory / f'{args.family}.post.pkl').open('rb') as stream:
        candidates, capped, maximum, graphs = pickle.load(stream)
    emit('audit_loaded', candidates=len(candidates))
    graph = graphs[0]
    incoming = defaultdict(list)
    for edge in graph.transitions: incoming[edge.target].append(edge)
    seen, pending, reachable_edges = set(), list(graph.terminals), set()
    while pending:
        node = pending.pop()
        if node in seen: continue
        seen.add(node)
        for edge in incoming[node]:
            reachable_edges.add(edge.id)
            pending.append(edge.source)
    generators = Counter()
    references = 0
    image_slots = 0
    for edge in graph.transitions:
        if edge.match is None: continue
        for generator in edge.match['symmetry'].get('automorph_generators', ()):
            generators[tuple(generator)] += 1
            references += 1
            image_slots += len(generator)
    stats = dict(family=args.family, branches=len(graph.terminals),
        transitions=len(graph.transitions), transitions_reaching_terminals=len(reachable_edges),
        stored_generator_occurrences=references, distinct_generator_values=len(generators),
        stored_image_slots=image_slots, unique_image_slots=sum(map(len,generators)),
        candidates=len(candidates), **memory())
    (args.directory / f'{args.family}.storage_audit.json').write_text(json.dumps(stats, indent=2) + '\n')
    print(json.dumps(stats), flush=True)


def compare(args):
    def normalized(states):
        output = []
        for mapping, islands, cuts in states:
            parts = {}
            for atom, label in islands: parts.setdefault(label, []).append(atom)
            output.append((mapping, sorted(sorted(p) for p in parts.values()), cuts))
        return json.dumps(sorted(output), default=int)
    rows = []
    for family in json.loads((args.directory / 'inputs.json').read_text())['families']:
        i = family['family']
        old = json.loads((args.directory / f'{i}.old.states.json').read_text())
        with (args.directory / f'{i}.graph.pkl').open('rb') as stream:
            graph = pickle.load(stream)
        current = [(s.mapping, s.islands, s.deferred_edges)
                   for s in (graph.states[t] for t in graph.terminals)]
        rows.append(dict(family=i, identical_terminal_relations=normalized(old)==normalized(current)))
        del graph, current
        gc.collect()
    (args.directory / 'core_comparison.json').write_text(json.dumps(rows, indent=2) + '\n')
    print(json.dumps(dict(compared=len(rows), differing=[r for r in rows if not r['identical_terminal_relations']])))


def compare_post(args):
    """Compare every occupation and its correlated evidence, not just counts.

    Histories not leading to returned terminals retain their raw matching
    records, but need not have their unused exact groups eagerly calculated.
    """
    from dataclasses import fields
    from rxn_core.search_graph import frozen_value
    def fingerprint(directory):
        with (directory / f'{args.family}.post.pkl').open('rb') as stream:
            candidates, capped, maximum, graphs = pickle.load(stream)
        bases, paths = {}, {}
        def digest(value):
            return hashlib.sha256(repr(frozen_value(value)).encode()).hexdigest()
        def hierarchy(value):
            parts = []
            for base, action in value.segments:
                if id(base) not in bases:
                    bases[id(base)] = digest(base.to_record())
                parts.append((bases[id(base)], action))
            return tuple(parts)
        def path(value):
            if id(value) not in paths:
                paths[id(value)] = digest((value.terminal, value.transitions,
                    value.mapping, value.deferred_edges, hierarchy(value.hierarchy)))
            return paths[id(value)]
        occupations = []
        for candidate in candidates:
            raw = {field.name: getattr(candidate, field.name) for field in fields(candidate)
                   if field.name not in {'aam_hierarchy', 'derivations'}}
            raw['hierarchy'] = hierarchy(candidate.aam_hierarchy)
            raw['derivations'] = tuple((tuple(map(path, d.initial_paths)),
                tuple(map(path, d.residual_paths)), d.target_action, d.occupation_projected)
                for d in candidate.derivations)
            occupations.append(digest(raw))
        structure = []
        exact_groups = []
        for graph in graphs:
            retained = graph.ancestor_transitions(graph.terminals)
            matches = []
            for edge in graph.transitions:
                match = edge.match
                if match is not None:
                    state = match['symmetry']
                    if edge.id in retained:
                        exact_groups.append(digest(state))
                    match = {**match, 'symmetry': {k: v for k, v in state.items()
                        if k not in {'automorph_generators', 'automorph_group_source'}}}
                matches.append((edge.id, edge.source, edge.target, edge.seed,
                                edge.step, match, edge.preserved_bonds))
            structure.append(digest((graph.contexts, graph.roots, graph.states,
                                     graph.stops, matches)))
        result = dict(capped=capped, maximum=maximum, occupations=occupations,
                      structure=structure, exact_groups=exact_groups)
        del candidates, graphs
        gc.collect()
        return result
    old = fingerprint(args.directory)
    new = fingerprint(args.output_directory)
    result = dict(family=args.family, candidates=len(new['occupations']),
                  identical={key: old[key] == new[key] for key in old})
    (args.output_directory / f'{args.family}.comparison.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)
    if not all(result['identical'].values()):
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['prepare', 'core', 'post', 'suite', 'audit', 'compare', 'compare_post'])
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--output-directory', type=Path)
    parser.add_argument('--source-id', default='INVENTORY-000400')
    parser.add_argument('--family', type=int)
    parser.add_argument('--version', choices=['old', 'current'])
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    globals()[args.mode](args)
