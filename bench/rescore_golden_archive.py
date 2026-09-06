"""Rescore saved evidence, preserving the previous evaluation and task status."""
import argparse
import os
from pathlib import Path
import shutil
import time

from golden_campaign import score,save
from recover_golden_archive import recover


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--index',type=int,required=True)
    args=parser.parse_args();directory=args.run/str(args.index)
    attempt=directory/'rescoring'/str(time.time_ns());attempt.mkdir(parents=True)
    for name in ('status.json','evaluation.json'):
        if (directory/name).exists():shutil.copy2(directory/name,attempt/name)
    save(attempt/'attempt.json',dict(index=args.index,slurm_job=os.environ.get('SLURM_JOB_ID'),
        started=time.time(),search_rerun=False))
    if not (directory/'aam.json.gz').exists():recover(args.run,args.index)
    score(args)
    shutil.copy2(directory/'evaluation.json',attempt/'new_evaluation.json')
    save(directory/'status.json',dict(index=args.index,stage='complete',rescored=True,
        evidence=str(attempt),updated=time.time()))


if __name__=='__main__':main()
