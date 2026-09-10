"""Finish confirmed unstarted full-benchmark tasks in one bounded allocation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main(args):
    from adaptive_full_benchmark import save
    tasks = json.loads(args.tasks.read_text())
    cpus = sorted(os.sched_getaffinity(0))
    needed = sum(8 if t['method'] == 'original' else 1 for t in tasks)
    if len(cpus) < needed:
        raise ValueError(f'Need {needed} allocated CPUs, received {len(cpus)}')
    assigned = []
    for task in tasks:
        n = 8 if task['method'] == 'original' else 1
        selected, cpus = cpus[:n], cpus[n:]
        status = args.run/f"status/{task['method']}_{task['slot']}.json"
        if status.exists():
            raise ValueError(f'Worker already started: {status}')
        assigned.append(dict(**task, affinity=selected))
    save(args.run/'completion_batch_started.json', dict(tasks=assigned, time=time.time(),
        job=os.environ.get('SLURM_JOB_ID'), host=os.uname().nodename))
    def work(task):
        command = ['taskset', '-c', ','.join(map(str, task['affinity'])),
            'timeout', '--kill-after=5s', '550', sys.executable,
            str(args.run/'engine/bench/adaptive_full_benchmark.py'), 'worker',
            '--run', str(args.run), '--method', task['method'], '--slot', str(task['slot'])]
        result = subprocess.run(command)
        return dict(**task, exit=result.returncode)
    with ThreadPoolExecutor(max_workers=len(tasks)) as pool:
        results = list(pool.map(work, assigned))
    save(args.run/'completion_batch_finished.json', dict(tasks=results, time=time.time()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--tasks', type=Path, required=True)
    main(parser.parse_args())
