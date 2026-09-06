#!/usr/bin/env python3
"""Build a disk-sorted proposal stream from saved connected placements."""
import argparse
import gzip
import json
from pathlib import Path
import pickle
import sqlite3
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--query', type=Path, required=True)
    p.add_argument('--shard', type=int, required=True)
    a = p.parse_args()
    manifest = json.loads((a.query/'manifest.json').read_text())
    output = a.query/'indexes'/f'{a.shard}.sqlite'
    if output.exists():
        print(json.dumps({'shard':a.shard, 'reused': True}), flush=True)
        return
    started = time.perf_counter()
    temp = output.with_suffix('.tmp')
    with sqlite3.connect(temp) as db:
        db.execute('PRAGMA journal_mode=OFF')
        db.execute('PRAGMA synchronous=OFF')
        db.execute('DROP TABLE IF EXISTS blocks')
        db.execute('CREATE TABLE blocks (coverage INTEGER, fragments INTEGER, atoms INTEGER, source TEXT, row INTEGER, ordinal INTEGER, mapping TEXT, partition TEXT)')
        sources = blocks = 0
        for row in range(a.shard, manifest['bank_rows'], manifest['shards']):
            path = a.query/'sources'/f'{row}.blocks.pkl.gz'
            with gzip.open(path,'rb') as stream:
                values = pickle.load(stream)
            db.executemany('INSERT INTO blocks VALUES (?,?,?,?,?,?,?,?)',
                ((len(v.covered_atoms),v.fragment_count,v.input_atom_count,v.source_id,
                  row,i,json.dumps(v.mapping),json.dumps(v.candidate.retained_fragments))
                 for i,v in enumerate(values)))
            sources += 1
            blocks += len(values)
        db.execute('CREATE INDEX ordering ON blocks (coverage DESC, fragments, atoms, source, row, ordinal)')
        db.commit()
    temp.replace(output)
    report=dict(shard=a.shard,sources=sources,blocks=blocks,seconds=time.perf_counter()-started)
    (a.query/'indexes'/f'{a.shard}.json').write_text(json.dumps(report)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
