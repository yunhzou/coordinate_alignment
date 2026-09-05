#!/usr/bin/env python3
"""Replay saved exact occupations through archive construction and persistence."""
import argparse
import cProfile
import gzip
import hashlib
import json
from pathlib import Path
import pickle
import pstats
import resource
import time
from dataclasses import replace
from types import SimpleNamespace

from rxn_core.fragment_matching.models import FragmentDetectionResult
from rxn_core.fragment_matching.serialization import fragment_detection_to_record


def memory():
    values = dict(line.split(':', 1) for line in Path('/proc/self/status').read_text().splitlines()
                  if ':' in line)
    return dict(rss_mb=int(values['VmRSS'].split()[0]) / 1024,
                peak_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--profile', action='store_true')
    parser.add_argument('--codec', choices=['stdlib', 'rapidjson', 'msgspec'], default='stdlib')
    parser.add_argument('--compression-level', type=int, default=9)
    parser.add_argument('--no-circular-check', action='store_true')
    parser.add_argument('--merge-only', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.input.open('rb') as stream:
        candidates, capped, maximum, graphs = pickle.load(stream)
    if args.merge_only:
        from rxn_core.fragment_matching.detection import _candidate_identity, _detect_fragments_from_initial
        from rxn_core.fragment_matching.models import FragmentDetectionConfig
        started = time.perf_counter()
        old, seen = [], {}
        for raw in candidates:
            candidate = replace(raw, source_id='saved-family')
            identity = _candidate_identity(candidate)
            if identity in seen:
                index = seen[identity]
                old[index] = replace(old[index], derivations=old[index].derivations + candidate.derivations)
            else:
                seen[identity] = len(old)
                old.append(candidate)
        old.sort(key=lambda c: (c.covered_target_atoms, c.attachment_atoms_target, c.mapping))
        before = time.perf_counter()-started
        placement = SimpleNamespace(retained_atoms=(), representative_mapping=(), encounter_count=1)
        initial = ((placement,), 0, 0, False, False, 1, 0, False, ())
        started = time.perf_counter()
        new = _detect_fragments_from_initial(None, SimpleNamespace(graph=None), initial,
            source_id='saved-family', config=FragmentDetectionConfig(), region=None,
            augmentation_runner=lambda _: ((candidates, capped, maximum, graphs),))
        after = time.perf_counter()-started
        stats = dict(input_candidates=len(candidates), retained_candidates=len(old),
                     before_seconds=before, after_seconds=after,
                     identical_candidates=tuple(old) == new.candidates, **memory())
        (args.output / 'merge.json').write_text(json.dumps(stats, indent=2) + '\n')
        print(json.dumps(stats), flush=True)
        if not stats['identical_candidates']:
            raise SystemExit(1)
        return
    result = FragmentDetectionResult(source_id='saved-family', candidates=candidates,
        status='capped' if capped else 'matched', complete=not capped,
        branch_limit=100, maximum_branch_count=maximum, capped_seed_count=int(capped),
        best_fragment_size=max(c.retained_size for c in candidates),
        initial_placement_encounters=1, initial_family_count=1, best_initial_family_count=1,
        seed_attempt_count=1, seed_pruned_count=0, rough_stop_hit=False, search_graphs=graphs)
    stats = dict(candidates=len(candidates), input=str(args.input), profiled=args.profile,
                 codec=args.codec, compression_level=args.compression_level,
                 check_circular=not args.no_circular_check,
                 loaded=memory())
    print(json.dumps(stats), flush=True)
    profiler = cProfile.Profile() if args.profile else None
    if profiler is not None:
        profiler.enable()
    started = time.perf_counter()
    record = fragment_detection_to_record(result, row_index=0, representation='saved-family')
    stats['archive'] = dict(seconds=time.perf_counter()-started, **memory())
    print(json.dumps({'archive': stats['archive']}), flush=True)
    started = time.perf_counter()
    if args.codec == 'msgspec':
        import msgspec
        encoded = msgspec.json.encode(record) + b'\n'
        payload_size = len(encoded)
    elif args.codec == 'rapidjson':
        import rapidjson
        payload = rapidjson.dumps(record, mapping_mode=rapidjson.MM_COERCE_KEYS_TO_STRINGS) + '\n'
    else:
        payload = json.dumps(record, separators=(',', ':'),
                             check_circular=not args.no_circular_check) + '\n'
    if args.codec != 'msgspec':
        payload_size = len(payload)
    stats['json'] = dict(seconds=time.perf_counter()-started, bytes=payload_size, **memory())
    print(json.dumps({'json': stats['json']}), flush=True)
    started = time.perf_counter()
    if args.codec != 'msgspec':
        encoded = payload.encode('utf-8')
    compressed = gzip.compress(encoded, compresslevel=args.compression_level, mtime=0)
    stats['gzip'] = dict(seconds=time.perf_counter()-started, bytes=len(compressed), **memory())
    stats['sha256'] = hashlib.sha256(encoded).hexdigest()
    if profiler is not None:
        profiler.disable()
        profiler.dump_stats(str(args.output / 'profile.pstats'))
        with (args.output / 'profile.txt').open('w') as stream:
            pstats.Stats(profiler, stream=stream).sort_stats('cumulative').print_stats(35)
    (args.output / 'result.json.gz').write_bytes(compressed)
    (args.output / 'metrics.json').write_text(json.dumps(stats, indent=2) + '\n')
    print(json.dumps(stats), flush=True)


if __name__ == '__main__':
    main()
