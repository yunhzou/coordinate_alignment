"""Validate and compress a completed raw archive after a wrapper timeout.

No AAM search is invoked. Original timeout logs and checkpoints are preserved.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import time

from golden_campaign import save
from rxn_core.artifacts import aam_from_record


def recover(run,index):
    directory=run/str(index)
    if (directory/'aam.json.gz').exists():
        raise ValueError('Complete compressed archive already exists')
    source=directory/'cuts/aam.json'
    started=time.perf_counter()
    # A killed raw write is not a complete result. Require a valid typed
    # archive before publishing it; a parse/validation error remains explicit.
    with source.open() as stream:result=aam_from_record(json.load(stream))
    with source.open('rb') as stream:
        digest=hashlib.file_digest(stream,'sha256').hexdigest()
    temporary=directory/'aam.recovered.json.gz.tmp'
    with source.open('rb') as raw,gzip.open(temporary,'wb',compresslevel=1) as compressed:
        shutil.copyfileobj(raw,compressed,length=1024*1024)
    temporary.replace(directory/'aam.json.gz')
    save(directory/'archive_recovery.json',dict(source=str(source),sha256=digest,
        elapsed_seconds=time.perf_counter()-started,search_rerun=False,
        search_complete=True,metrics=vars(result.metrics),
        previous_status=json.loads((directory/'status.json').read_text())))
    # Preserve partial evidence separately before full-archive rescoring.
    if (directory/'evaluation.json').exists():
        shutil.copy2(directory/'evaluation.json',directory/'partial_evaluation.json')
    print(json.dumps(dict(index=index,complete_archive_recovered=True)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int,required=True)
    args=parser.parse_args();recover(args.run,args.index)
