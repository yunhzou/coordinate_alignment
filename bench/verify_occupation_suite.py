#!/usr/bin/env python3
"""Validate the observation quotient on captured calls without repeating AAM."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import gzip
import json
from pathlib import Path
import pickle
import sys
import time

from rxn_core._group_ops import occupation_orbit
from verify_observed_occupations import relation


def verify(task):
    index, path, output = task
    from _occupation_probe import probe
    with gzip.open(path, 'rb') as stream:
        captured = pickle.load(stream)
    parameters = captured['parameters']
    baseline = probe(*parameters[:6], False)
    observed = set(range(captured['context']['target_count']))
    started = time.perf_counter()
    compact = occupation_orbit(*parameters, sorted(observed))
    seconds = time.perf_counter() - started
    _, _, _, attachments, fragments, bonds, _ = parameters
    key = lambda images: relation(images, observed, attachments, fragments, bonds)
    full = baseline['states']
    assert {key(x) for x, _ in compact} == {key(x) for x, _ in full}, str(path)
    pairs = lambda images: tuple((i, p) for i, p in enumerate(images) if p in observed)
    row = dict(index=index, input=str(path), full_states=len(full), compact_states=len(compact),
               baseline_seconds=sum(s['seconds'] for s in baseline['stages']),
               compact_seconds=seconds, observed_relations_identical=True,
               exact_source_target_pairs_identical={pairs(x) for x, _ in compact} ==
                                                  {pairs(x) for x, _ in full})
    with gzip.open(output / f'{index}.states.pkl.gz', 'xb', compresslevel=1) as stream:
        pickle.dump(dict(full=full, compact=compact), stream, protocol=5)
    (output / f'{index}.json').write_text(json.dumps(row) + '\n')
    return row


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--captures', type=Path, nargs='+', required=True)
    parser.add_argument('--reference-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.reference_directory.resolve()))
    args.output.mkdir(parents=True, exist_ok=False)
    paths = [p for directory in args.captures for p in sorted(directory.glob('*.pkl.gz'))]
    started = time.perf_counter()
    rows = []
    with ProcessPoolExecutor(args.workers) as pool:
        futures = [pool.submit(verify, (i, path, args.output)) for i, path in enumerate(paths)]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(json.dumps(dict(completed=len(rows), expected=len(paths), **row)), flush=True)
    report = dict(cases=len(rows), elapsed_seconds=time.perf_counter()-started,
                  full_states=sum(r['full_states'] for r in rows),
                  compact_states=sum(r['compact_states'] for r in rows),
                  observed_relations_identical=all(r['observed_relations_identical'] for r in rows),
                  exact_source_target_pairs_identical=all(r['exact_source_target_pairs_identical'] for r in rows))
    (args.output / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    main()
