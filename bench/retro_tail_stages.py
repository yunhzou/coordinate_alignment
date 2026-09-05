#!/usr/bin/env python3
"""Checkpoint and time detection, archive construction, and encoding separately."""
import argparse
import faulthandler
import json
from pathlib import Path
import pickle
import resource
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from search_mcule_retro import _worker_init, _detect_one, _record_detection, _encode_records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior-run', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--input-dir', type=Path)
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--workers', type=int, default=48)
    parser.add_argument('--stage', choices=['detect', 'archive', 'encode', 'all'], default='all')
    args = parser.parse_args()
    directory = args.output_dir / str(args.index)
    inputs = (args.input_dir or args.output_dir) / str(args.index)
    directory.mkdir(parents=True, exist_ok=True)
    # Trace the process from inside: external ptrace is disabled on these nodes.
    trace = (directory / f'{args.stage}.stacks.txt').open('w')
    faulthandler.dump_traceback_later(60, repeat=True, file=trace)
    config = json.loads(next((args.prior_run / 'parts').glob('*.summary.json')).read_text())
    source = json.loads((args.prior_run / 'progress_audit.json').read_text())['unfinished_sources'][args.index]
    row = source['row_index'], source['smiles'], source['source_id']
    _worker_init(config['target_smiles'], config['config'],
                 config['minimum_target_coverage_fraction'], True, None)
    metrics = dict(source_id=source['source_id'],
                   workers=args.workers if args.stage in {'detect', 'all'} else 1,
                   stage=args.stage, config=config['config'])
    def checkpoint(stage, started):
        metrics[stage] = dict(seconds=time.perf_counter()-started,
            peak_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
        (directory / f'{args.stage}.metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
        print(json.dumps({stage: metrics[stage], 'source_id': source['source_id']}), flush=True)
    if args.stage in {'detect', 'all'}:
        started = time.perf_counter()
        detection = _detect_one(row, seed_workers=args.workers-1)
        checkpoint('detection', started)
        started = time.perf_counter()
        with (directory / 'detection.pkl').open('xb') as stream:
            pickle.dump(detection, stream, protocol=5)
        checkpoint('save_detection', started)
    if args.stage in {'archive', 'all'}:
        if args.stage == 'archive':
            started = time.perf_counter()
            with (inputs / 'detection.pkl').open('rb') as stream:
                detection = pickle.load(stream)
            checkpoint('load_detection', started)
        started = time.perf_counter()
        counts, record = _record_detection(row, *detection)
        checkpoint('archive', started)
        started = time.perf_counter()
        with (directory / 'archive.pkl').open('xb') as stream:
            pickle.dump((counts, record), stream, protocol=5)
        checkpoint('save_archive', started)
    if args.stage in {'encode', 'all'}:
        if args.stage == 'encode':
            started = time.perf_counter()
            with (inputs / 'archive.pkl').open('rb') as stream:
                counts, record = pickle.load(stream)
            checkpoint('load_archive', started)
        print(json.dumps({'candidates': 0 if record is None else len(record['candidates']),
                          'counts': dict(counts)}), flush=True)
        started = time.perf_counter()
        encoded = _encode_records([] if record is None else [record])
        checkpoint('encode', started)
        with (directory / 'result.jsonl.gz').open('xb') as stream:
            stream.write(encoded)
        metrics['counts'] = dict(counts)
        metrics['compressed_bytes'] = len(encoded)
        checkpoint('saved', started)
    faulthandler.cancel_dump_traceback_later()
    trace.close()


if __name__ == '__main__':
    main()
