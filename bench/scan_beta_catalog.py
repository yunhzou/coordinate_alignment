#!/usr/bin/env python3
"""Sharded connected-query execution for the opt-in beta workflow."""
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import csv
import gzip
import json
from pathlib import Path
import pickle
import time

import numpy as np

from rxn_core.frag import WeightedGraph
from rxn_core.fragment_matching import FragmentDetectionConfig, prepare_fragment_target
from rxn_core.fragment_matching.connected import find_connected_fragments
from rxn_core.fragment_matching.symmetry import materialize_target_coverage_orbit
from rxn_core.retrosynthesis.beta import BetaPlacement
from rxn_core.smiles import smiles_to_weighted_graph


def initialize(manifest, query):
    global CONTEXT, CONFIG, REGION, QUERY
    CONFIG = FragmentDetectionConfig(**manifest['config'])
    REGION = tuple(manifest['region'])
    QUERY = Path(query)
    target = smiles_to_weighted_graph(manifest['target_smiles'], expand_hydrogens=True)
    local = WeightedGraph([target.nodes[a] for a in REGION],
        np.asarray(target.weights)[np.ix_(REGION, REGION)])
    CONTEXT = prepare_fragment_target(local, config=CONFIG)


def scan(row):
    index, source_id, smiles = row
    prefix = QUERY / 'sources' / str(index)
    summary = prefix.with_suffix('.json')
    if summary.exists():
        return dict(json.loads(summary.read_text()), reused=True)
    started = time.perf_counter()
    evidence = prefix.with_suffix('.connected.pkl.gz')
    reused_detection = evidence.exists()
    if reused_detection:
        with gzip.open(evidence, 'rb') as stream:
            result, detection_seconds = pickle.load(stream)
    else:
        source = smiles_to_weighted_graph(smiles, expand_hydrogens=True)
        result = find_connected_fragments(source, CONTEXT, source_id=source_id, config=CONFIG)
        detection_seconds = time.perf_counter() - started
        temporary = evidence.with_suffix('.tmp')
        with gzip.open(temporary, 'wb', compresslevel=1) as stream:
            pickle.dump((result, detection_seconds), stream, protocol=pickle.HIGHEST_PROTOCOL)
        temporary.replace(evidence)
    blocks = {}
    for candidate in result.candidates:
        for variant in materialize_target_coverage_orbit(candidate, CONTEXT.graph,
                iso_tolerance=CONFIG.iso_tolerance,
                generators=CONTEXT.automorphism_generators):
            block = BetaPlacement(variant, REGION)
            blocks.setdefault(block.key, block)
    archive = prefix.with_suffix('.blocks.pkl.gz')
    temporary = archive.with_suffix('.tmp')
    with gzip.open(temporary, 'wb', compresslevel=1) as stream:
        pickle.dump(tuple(sorted(blocks.values(), key=lambda p: p.key)), stream,
                    protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(archive)
    record = dict(row_index=index, source_id=source_id, block_file=str(archive),
        connected_file=str(evidence), candidates=len(result.candidates), blocks=len(blocks),
        capped=bool(result.capped_seed_count), complete=result.complete,
        detection_seconds=detection_seconds, pair_seconds=time.perf_counter()-started,
        reused_detection=reused_detection)
    temporary = summary.with_suffix('.tmp')
    temporary.write_text(json.dumps(record)+'\n')
    temporary.replace(summary)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--query', type=Path, required=True)
    parser.add_argument('--shard', type=int, required=True)
    parser.add_argument('--workers', type=int, required=True)
    parser.add_argument('--work-seconds', type=float)
    args = parser.parse_args()
    manifest = json.loads((args.query/'manifest.json').read_text())
    with gzip.open(manifest['catalog'], 'rt') as stream:
        rows = [(i, r[manifest['id_column']], r['SMILES'])
                for i, r in enumerate(csv.DictReader(stream))
                if i % manifest['shards'] == args.shard]
    started = time.perf_counter()
    deadline = started + args.work_seconds if args.work_seconds is not None else float('inf')
    records = []
    remaining = iter(rows)
    with ProcessPoolExecutor(args.workers, initializer=initialize,
            initargs=(manifest, str(args.query))) as pool:
        pending = set()
        for _ in range(args.workers):
            row = next(remaining, None)
            if row is not None:
                pending.add(pool.submit(scan, row))
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                record = future.result()
                records.append(record)
                print(json.dumps(record), flush=True)
                row = next(remaining, None) if time.perf_counter() < deadline else None
                if row is not None:
                    pending.add(pool.submit(scan, row))
    if len(records) < len(rows):
        print(json.dumps(dict(stage='yield_checkpoint',shard=args.shard,
            completed=len(records),expected=len(rows))),flush=True)
        raise SystemExit(75)
    report = dict(shard=args.shard, rows=len(records),
        seconds=time.perf_counter()-started,
        blocks=sum(r['blocks'] for r in records),
        candidates=sum(r['candidates'] for r in records),
        capped=sum(r['capped'] for r in records))
    (args.query/'parts'/f'{args.shard}.json').write_text(json.dumps(report)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
