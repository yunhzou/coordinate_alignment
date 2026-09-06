"""Monitor only owned Slurm tasks; enforce a ten-minute outer watchdog."""
import argparse
import json
from pathlib import Path
import subprocess
import time

from golden_campaign import save,status


def retry_unstarted(run, manifest, submission, queue):
    """Retry only proven launch failures, never a search already attempted."""
    active={line.split('|')[0].split('_')[0] for line in queue.splitlines()}
    if any(job.get('infrastructure_retry') and job['job'] in active for job in submission['jobs']):
        return
    slots=[]
    for slot,record in enumerate(manifest['records']):
        directory=run/str(record['index'])
        state=json.loads((directory/'status.json').read_text())
        unstarted=(state.get('slurm_stage')=='CONFIGURING' or
                   (state['stage'].startswith('scheduler_') and state.get('previous_stage')=='pending'))
        if not unstarted or (directory/'infrastructure_retry.json').exists():continue
        if (directory/'search.log').exists() or (directory/'aam.json.gz').exists():continue
        slots.append(slot)
    if not slots:return
    config=subprocess.check_output(['scontrol','show','config'],text=True)
    limit=int(next(line.split('=')[1] for line in config.splitlines() if line.startswith('MaxArraySize')))
    chunks=sorted({slot//limit for slot in slots})
    for chunk in chunks:
        offset=chunk*limit;selected=[slot for slot in slots if slot//limit==chunk]
        for slot in selected:
            directory=run/str(manifest['records'][slot]['index'])
            save(directory/'infrastructure_retry.json',dict(previous_status=json.loads((directory/'status.json').read_text()),
                 reason='Slurm failed to start search; retry once on cpunodes partition',created=time.time()))
            save(directory/'status.json',dict(index=manifest['records'][slot]['index'],stage='pending'))
        array=','.join(str(slot-offset) for slot in selected)+'%'+str(32//len(chunks))
        job=subprocess.check_output(['sbatch','--parsable','--partition=cpunodes','--cpus-per-task=2',
            f'--array={array}',f'--output={run}/slurm_%A_%a.out',f'--error={run}/slurm_%A_%a.err',
            'hpc/golden_benchmark.sbatch',str(run),str(offset)],text=True).strip().split(';')[0]
        submission['jobs'].append(dict(job=job,offset=offset,infrastructure_retry=True,slots=selected))
        save(run/'submission.json',submission)


def elapsed_seconds(value):
    days=0
    if '-' in value:
        prefix,value=value.split('-');days=int(prefix)
    return days*86400+sum(int(v)*60**i for i,v in enumerate(reversed(value.split(':'))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--retry-infrastructure',action='store_true',
                        help='Reserve 64 of the 256 CPU budget for one batch of launch-failure retries')
    args=parser.parse_args()
    manifest=json.loads((args.run/'manifest.json').read_text())
    submission=json.loads((args.run/'submission.json').read_text())
    offsets={item['job']:item['offset'] for item in submission['jobs']}
    completed_seen={}
    while True:
        offsets={item['job']:item['offset'] for item in submission['jobs']}
        queue=subprocess.check_output(['squeue','-r','-h','-j',','.join(offsets),
            '-o','%i|%T|%M|%N'],text=True)
        for line in queue.splitlines():
            job,state,elapsed,node=line.split('|')
            if '_' not in job:continue
            parent,slot=job.split('_')
            if parent not in offsets or not slot.isdigit():continue
            if state not in ('RUNNING','CONFIGURING') or elapsed_seconds(elapsed)<=600:continue
            index=manifest['records'][offsets[parent]+int(slot)]['index']
            subprocess.run(['scancel',job],check=True)
            save(args.run/str(index)/'status.json',dict(index=index,stage='supervisor_timeout',
                slurm_job=job,node=node,slurm_stage=state,seconds=elapsed_seconds(elapsed)))
        accounting=subprocess.check_output(['sacct','-n','-P','-j',','.join(offsets),
            '--format=JobID,State,ExitCode,ElapsedRaw,TotalCPU,MaxRSS,AllocCPUS,NodeList'],text=True)
        (args.run/'slurm_accounting.psv').write_text(accounting)
        for line in accounting.splitlines():
            fields=line.split('|');job,state=fields[:2]
            if '_' not in job or '.' in job:continue
            parent,slot=job.split('_')
            if parent not in offsets or not slot.isdigit():continue
            if state in ('RUNNING','PENDING','CONFIGURING','COMPLETING'):continue
            first=completed_seen.setdefault(job,time.time())
            if time.time()-first<60:continue  # Allow shared-filesystem status visibility.
            index=manifest['records'][offsets[parent]+int(slot)]['index']
            path=args.run/str(index)/'status.json'
            recorded=json.loads(path.read_text())
            retry=path.parent/'infrastructure_retry.json'
            if retry.exists() and not next(item for item in submission['jobs'] if item['job']==parent).get('infrastructure_retry'):
                continue  # The original failed launch must not overwrite its retry.
            if recorded['stage'] in ('pending','search','score'):
                save(path,dict(index=index,stage='scheduler_'+state.lower(),slurm_job=job,
                               exit_code=fields[2],previous_stage=recorded['stage']))
        if args.retry_infrastructure:
            retry_unstarted(args.run,manifest,submission,queue)
        status(args)
        progress=json.loads((args.run/'progress.json').read_text())
        if not any(progress['stages'].get(stage) for stage in ('pending','search','score')):
            from report_golden_campaign import export
            export(args.run,args.run/'report')
            break
        time.sleep(30)


if __name__=='__main__':main()
