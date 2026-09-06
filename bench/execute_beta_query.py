"""Reuse a Slurm allocation for scanning and indexing; record actual CPU work."""
import argparse
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


def usage():
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return (own.ru_utime + children.ru_utime, own.ru_stime + children.ru_stime)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--query', type=Path, required=True)
    parser.add_argument('--shard', type=int, required=True)
    parser.add_argument('--workers', type=int, required=True)
    parser.add_argument('--wall-seconds', type=float, default=580)
    args = parser.parse_args()
    started = time.monotonic()
    stages = []
    code = 0
    for label, script, extra in (
            ('scan', 'scan_beta_catalog.py', ['--workers', str(args.workers),
                '--work-seconds', str(args.wall_seconds * .75)]),
            ('index', 'index_beta_query.py', [])):
        before = usage()
        stage_start = time.monotonic()
        result = subprocess.run([sys.executable, str(Path(__file__).with_name(script)),
            '--query', str(args.query), '--shard', str(args.shard), *extra])
        after = usage()
        stages.append(dict(stage=label, wall_seconds=time.monotonic()-stage_start,
            user_cpu_seconds=after[0]-before[0], system_cpu_seconds=after[1]-before[1],
            exit_code=result.returncode))
        if result.returncode:
            code = result.returncode
            break
    elapsed = time.monotonic()-started
    cpu_seconds = sum(s['user_cpu_seconds']+s['system_cpu_seconds'] for s in stages)
    report = dict(shard=args.shard, job_id=os.environ.get('SLURM_JOB_ID'),
        allocated_cpus=args.workers, wall_seconds=elapsed, actual_cpu_seconds=cpu_seconds,
        allocated_cpu_seconds=elapsed*args.workers,
        utilization=cpu_seconds/(elapsed*args.workers), stages=stages, exit_code=code)
    directory = args.query/'resources'
    directory.mkdir(exist_ok=True)
    path = directory/f"{args.shard}.{os.environ.get('SLURM_JOB_ID', 'local')}.json"
    path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(stage='resource_usage', **report)), flush=True)
    raise SystemExit(code)


if __name__ == '__main__':
    main()
