"""Monitor only owned Slurm tasks; enforce a ten-minute outer watchdog."""
import argparse
import json
from pathlib import Path
import subprocess
import time

from golden_campaign import save,status


def elapsed_seconds(value):
    days=0
    if '-' in value:
        prefix,value=value.split('-');days=int(prefix)
    return days*86400+sum(int(v)*60**i for i,v in enumerate(reversed(value.split(':'))))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args()
    manifest=json.loads((args.run/'manifest.json').read_text())
    submission=json.loads((args.run/'submission.json').read_text())
    offsets={item['job']:item['offset'] for item in submission['jobs']}
    completed_seen={}
    while True:
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
            if recorded['stage'] in ('pending','search','score'):
                save(path,dict(index=index,stage='scheduler_'+state.lower(),slurm_job=job,
                               exit_code=fields[2],previous_stage=recorded['stage']))
        status(args)
        progress=json.loads((args.run/'progress.json').read_text())
        if not any(progress['stages'].get(stage) for stage in ('pending','search','score')):
            break
        time.sleep(30)


if __name__=='__main__':main()
