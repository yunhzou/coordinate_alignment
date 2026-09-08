"""Dispatch a recorded list of comparison workers in one Slurm allocation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--tasks',type=Path,required=True)
    p.add_argument('--workers',type=int,default=8)
    a=p.parse_args()
    indices=json.loads(a.tasks.read_text())
    def work(index):
        status=a.run/'status'/f'{index}_union_combine.json'
        if status.exists():
            raise RuntimeError(f'Comparison already started: {index}')
        started=time.time()
        command=[sys.executable,str(a.run/'engine/bench/golden_publication.py'),
                 'worker','--run',str(a.run),'--phase','combine','--offset','0','--slot',str(index)]
        result=subprocess.run(command,check=False)
        print(json.dumps(dict(index=index,worker_exit=result.returncode,
                              started=started,finished=time.time())),flush=True)
        return result.returncode
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        exits=list(pool.map(work,indices))
    if any(exits):
        raise SystemExit(1)


if __name__=='__main__':
    main()
