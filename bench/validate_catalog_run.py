#!/usr/bin/env python3
"""Validate saved whole-bank records and index their locations without matching."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import hashlib
import json
from pathlib import Path

from msgspec import Raw, Struct, json as fast_json


class _Record(Struct):
    source_id: str
    row_index: int
    representation: str
    status: str
    complete: bool
    candidates: list[Raw]
    best_fragment_size: int
    capped_seed_count: int
    branch_limit: int
    timing: dict[str, float]


def scan_file(path):
    rows = []
    # Validate JSON syntax without rebuilding millions of nested archive
    # objects. This is an integrity/index check, not a chemical result replay.
    decoder = fast_json.Decoder(_Record)
    with gzip.open(path, 'rb') as stream:
        for line_number, line in enumerate(stream, 1):
            record = decoder.decode(line)
            rows.append(dict(
                source_id=record.source_id, row_index=record.row_index,
                representation=record.representation, status=record.status,
                complete=record.complete, candidates=len(record.candidates),
                best_fragment_size=record.best_fragment_size,
                capped_seed_count=record.capped_seed_count,
                branch_limit=record.branch_limit,
                detection_seconds=record.timing['detection_seconds'],
                file=str(path), line=line_number,
                record_sha256=hashlib.sha256(line).hexdigest()))
    return rows


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--catalog', type=Path, required=True)
    parser.add_argument('--id-column', default='Inventory ID')
    parser.add_argument('--parts', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    with gzip.open(args.catalog, 'rt') as stream:
        bank = list(csv.DictReader(stream))
    expected = {i: row for i, row in enumerate(bank)}
    files = [path for directory in args.parts for path in sorted(directory.glob('*.jsonl.gz'))]
    found = {}
    with ProcessPoolExecutor(args.workers) as pool:
        for rows in pool.map(scan_file, files):
            for row in rows:
                index = row['row_index']
                if index in found:
                    raise ValueError(f'duplicate bank row: {index}')
                if index not in expected:
                    raise ValueError(f'unknown bank row: {index}')
                source = expected[index]
                if (row['source_id'], row['representation']) != (source[args.id_column], source['SMILES']):
                    raise ValueError(f'source identity mismatch: {index}')
                found[index] = row
    rows = [found[index] for index in sorted(found)]
    report = dict(expected_rows=len(bank), saved_rows=len(rows),
                  complete=set(found) == set(expected),
                  missing_rows=sorted(set(expected) - set(found)),
                  statuses=dict(Counter(row['status'] for row in rows)),
                  capped_sources=sum(row['capped_seed_count'] > 0 for row in rows),
                  exact_complete_sources=sum(row['complete'] for row in rows),
                  matched_sources=sum(row['candidates'] > 0 for row in rows),
                  candidates=sum(row['candidates'] for row in rows),
                  branch_limits=sorted({row['branch_limit'] for row in rows}),
                  files=len(files), records=rows)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'records'}), flush=True)
    if not report['complete']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
