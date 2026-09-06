#!/usr/bin/env python3
"""Compare the production observation quotient with saved pre-fix native logic."""
import argparse
import gzip
import json
from pathlib import Path
import pickle
import sys
import time

from rxn_core._group_ops import occupation_orbit


def relation(images, observed, attachments, fragments, bonds):
    values = tuple(p if p in observed else -1 for p in images)
    return (tuple(sorted(values)), tuple(sorted(values[i] for i in attachments)),
            tuple(sorted((label, tuple(sorted(values[i] for i in positions)))
                         for label, positions in fragments)),
            tuple(sorted(tuple(sorted((values[a], values[b]))) for a, b in bonds)))


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--input', type=Path, nargs='+', required=True)
    parser.add_argument('--reference-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.reference_directory.resolve()))
    from _occupation_probe import probe
    args.output.mkdir(parents=True, exist_ok=False)
    reports = []
    for index, path in enumerate(args.input):
        with gzip.open(path, 'rb') as stream:
            captured = pickle.load(stream)
        parameters = captured['parameters']
        baseline = probe(*parameters[:6], False)
        observed = set(range(captured['context']['target_count']))
        started = time.perf_counter()
        full = occupation_orbit(*parameters)
        full_seconds = time.perf_counter() - started
        assert full == baseline['states'], 'unobserved output/witness/order changed'
        started = time.perf_counter()
        compact = occupation_orbit(*parameters, sorted(observed))
        compact_seconds = time.perf_counter() - started
        _, _, _, attachments, fragments, bonds, _ = parameters
        key = lambda images: relation(images, observed, attachments, fragments, bonds)
        assert {key(x) for x, _ in full} == {key(x) for x, _ in compact}
        exact_pairs = lambda images: tuple((i, p) for i, p in enumerate(images) if p in observed)
        pairs_equal = {exact_pairs(x) for x, _ in full} == {exact_pairs(x) for x, _ in compact}
        with gzip.open(args.output / f'{index}.states.pkl.gz', 'xb', compresslevel=1) as stream:
            pickle.dump(dict(full=full, compact=compact), stream, protocol=5)
        row = dict(input=str(path), baseline_seconds=sum(s['seconds'] for s in baseline['stages']),
                   full_seconds=full_seconds, compact_seconds=compact_seconds,
                   full_states=len(full), compact_states=len(compact),
                   full_output_identical=True, observed_relations_identical=True,
                   exact_source_target_pairs_identical=pairs_equal)
        reports.append(row)
        (args.output / 'results.json').write_text(json.dumps(reports, indent=2) + '\n')
        print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
