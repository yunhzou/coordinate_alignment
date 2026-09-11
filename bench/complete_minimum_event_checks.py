"""Reproduce targeted SLAP enumeration retries and enriched AAM membership checks.

Write to a new destination; preserve the original result and query timeouts.
"""
import argparse
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
import time

import pynauty
from golden_evaluation import colored_graph, exact_action
from holdout_minimum_events import AAM, SWEEP, EventPatterns, read, save, sha
from holdout_minimum_event_families import aam_membership
from enumerate_slap_minimum_events import enumerate_family


def retry_slap(args):
    source = args.run / f'enumeration/{args.index}.json'
    result = read(source)
    raw = read(AAM / f'inputs/{args.index}/input.json')
    canonical = EventPatterns(raw)

    @lru_cache(None)
    def twin(side, a, b):
        images = list(range(canonical.n))
        images[a], images[b] = b, a
        return exact_action(images, canonical.features[side])

    deadline = time.perf_counter() + args.seconds
    for method, root in [('native_slap', AAM), ('slap_sweep', SWEEP)]:
        output = result['methods'][method]
        native = read(root / f'slap_xyz/{args.index}.json')['candidates']
        scores = read(root / f'slap_scoring/{args.index}.refined.json')['slap']
        for family in output['families']:
            if family['complete']:
                continue
            ordinal = family['candidate']
            retry = enumerate_family(raw, canonical, native[ordinal], scores[ordinal],
                                     twin, deadline, check_ms=args.query_ms)
            for key, pattern in retry['patterns'].items():
                output['patterns'].setdefault(key, dict(pattern, candidate=ordinal))
            family.update({key: value for key, value in retry.items() if key != 'patterns'})
            family['pattern_ids'] = sorted(retry['patterns'])
        output['complete'] = all(family['complete'] for family in output['families'])
    result['extended_analysis'] = dict(source=str(source), source_sha256=sha(source),
                                       query_timeout_ms=args.query_ms)
    return result


def enrich_aam(args):
    source = args.run / f'families/{args.index}.json'
    initial = read(source)
    enum_path = args.run / f'enumeration_extended/{args.index}.json'
    if not enum_path.exists():
        enum_path = args.run / f'enumeration/{args.index}.json'
    enumeration = read(enum_path)
    snapshot = deepcopy(read(args.run / f'snapshots/{args.index}.json'))
    for key, value in initial['methods']['aam'].items():
        if value['status'] == 'represented':
            snapshot['methods']['aam']['patterns'].setdefault(key, value['witness'])
    patterns = dict(initial['patterns'])
    for output in enumeration['methods'].values():
        patterns.update(output['patterns'])
    raw = read(AAM / f'inputs/{args.index}/input.json')
    canonical = EventPatterns(raw)
    generators = tuple(tuple(g[:canonical.n]) for g in
                       pynauty.autgrp(colored_graph([canonical.features[0]]))[0])
    started = time.perf_counter()
    found, metrics = aam_membership(raw, canonical, patterns, snapshot, generators,
                                    started + args.seconds)
    return dict(index=args.index, name=raw['name'], patterns=patterns,
                methods={'aam': found}, query_metrics={'aam': metrics},
                seconds=time.perf_counter() - started,
                scope='Additional AAM family checks for the complete SLAP minimum-event '
                      'catalogue; initial positive witnesses retained.',
                sources_sha256={str(p): sha(p) for p in (source, enum_path)})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('retry_slap', 'enrich_aam'))
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--seconds', type=int, default=240)
    parser.add_argument('--query-ms', type=int, default=30000)
    args = parser.parse_args()
    assert not args.destination.exists(), 'Preserve existing analysis artifacts.'
    save(args.destination, globals()[args.command](args))
