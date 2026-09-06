"""Operator recovery of this suite's Slurm startup stalls; retain all chemistry checkpoints."""
import argparse
import json
from pathlib import Path
import subprocess
import time


def seconds(value):
    days=0
    if '-' in value:
        days,value=value.split('-');days=int(days)
    return days*86400+sum(int(v)*60**i for i,v in enumerate(reversed(value.split(':'))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    items=json.loads((args.root/'submissions.json').read_text())
    workers={}
    active=[]
    for item in items:
        path=Path(item['manifest'])
        status=json.loads(path.with_name('status.json').read_text())
        if status['stage'] in ('complete','failed'):
            continue
        active.append(item)
        for line in path.with_name('blind.out').read_text().splitlines():
            if '"submitted"' in line:
                workers[json.loads(line)['job']]=item
    bad=set('bosque38,bosque39,bosque48,bosque56,bosque61,bosque75'.split(','))
    stalls=[]
    queue=subprocess.check_output(['squeue','-h','-u','yunhengzou','-o','%i|%T|%M|%N'],text=True)
    for line in queue.splitlines():
        job,state,elapsed,node=line.split('|')
        if job.split('_')[0] in workers and state=='CONFIGURING' and seconds(elapsed)>120:
            bad.add(node);stalls.append(dict(job=job,node=node,elapsed=elapsed))
    if not stalls:
        raise RuntimeError('No measured startup stalls; refusing an unnecessary restart')
    record=dict(time=time.time(),stalls=stalls,excluded=sorted(bad),
        retired_coordinators=[i['job'] for i in active],retired_workers=list(workers))
    (args.root/f'startup_recovery_{int(time.time())}.json').write_text(json.dumps(record,indent=2)+'\n')
    subprocess.run(['scancel',*[i['job'] for i in active],*workers],check=True)
    # Give Slurm time to stop tasks before any new worker can touch their checkpoints.
    time.sleep(10)
    retired=set(workers)
    states=subprocess.check_output(['squeue','-h','-u','yunhengzou','-o','%i|%T|%N'],text=True)
    if any(j.split('_')[0] in retired and state=='RUNNING'
           for j,state,node in (line.split('|') for line in states.splitlines())):
        raise RuntimeError('Old workers still running; do not launch overlapping source writers')
    excluded=','.join(sorted(bad))
    for item in active:
        path=Path(item['manifest'])
        case=json.loads(path.read_text());case['exclude']=excluded
        path.write_text(json.dumps(case,indent=2)+'\n')
        (Path(case['run'])/'scheduler_exclusions.json').write_text(json.dumps(sorted(bad))+'\n')
        job=subprocess.check_output(['sbatch','--parsable',f'--exclude={excluded}',
            f'--output={path.parent}/coordinator_retry.out',f'--error={path.parent}/coordinator_retry.err',
            'hpc/beta_case.sbatch',str(path),str(item['cpu_budget'])],text=True).strip().split(';')[0]
        item.setdefault('previous_jobs',[]).append(item['job']);item['job']=job
        print(json.dumps(item),flush=True)
    (args.root/'submissions.json').write_text(json.dumps(items,indent=2)+'\n')


if __name__=='__main__':
    main()
