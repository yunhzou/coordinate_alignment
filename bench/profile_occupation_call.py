#!/usr/bin/env python3
"""Replay an isolated captured call; compare profiling/caching replicas exactly."""
import argparse
import gzip
import json
from pathlib import Path
import pickle
import sys
import time

from rxn_core._group_ops import occupation_orbit


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--probe-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.probe_directory.resolve()))
    from _occupation_probe import probe
    with gzip.open(args.input, 'rb') as stream:
        captured = pickle.load(stream)
    parameters = captured['parameters']
    started = time.perf_counter()
    actual = occupation_orbit(*parameters)
    native_seconds = time.perf_counter() - started
    with gzip.open(args.output.with_suffix('.states.pkl.gz'), 'xb', compresslevel=1) as stream:
        pickle.dump(actual, stream, protocol=5)
    print(json.dumps(dict(phase='native_done', seconds=native_seconds, states=len(actual))), flush=True)
    measurements = {}
    for label, cached in [('baseline', False), ('edge_cache', True)]:
        result = probe(*parameters[:6], cached)
        assert result['states'] == actual, 'changed occupation, witness, action, or order'
        measurements[label] = result['stages']
        print(json.dumps(dict(phase=label, seconds=sum(s['seconds'] for s in result['stages']),
                              states=len(actual), identical=True)), flush=True)
    _, _, stages, attachments, fragments, bonds, _ = parameters
    target_count = captured['context']['target_count']
    def projected_key(images):
        return (tuple(sorted(p for p in images if p < target_count)),
                tuple(sorted(images[i] for i in attachments if images[i] < target_count)),
                tuple(sorted((label, tuple(sorted(images[i] for i in positions
                            if images[i] < target_count))) for label, positions in fragments)),
                tuple(sorted(tuple(sorted((images[a], images[b]))) for a, b in bonds
                             if images[a] < target_count and images[b] < target_count)))
    report = dict(input=str(args.input), native_seconds=native_seconds, states=len(actual),
                  projected_target_relations=len({projected_key(images) for images, _ in actual}),
                  target_atom_sets=len({tuple(sorted(p for p in images if p < target_count))
                                        for images, _ in actual}),
                  generator_entries=sum(map(len, stages)),
                  distinct_generators=len({g for stage in stages for g in stage}),
                  exact_states_witnesses_actions_order_equal=True, measurements=measurements)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'measurements'}), flush=True)


if __name__ == '__main__':
    main()
