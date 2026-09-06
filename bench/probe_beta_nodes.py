"""Probe previously stalled nodes and restore only successful nodes to suite scheduling."""
import argparse
import json
from pathlib import Path
import subprocess


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--apply',action='store_true')
    args=p.parse_args()
    recovery=json.loads(sorted(args.root.glob('startup_recovery_*.json'))[-1].read_text())
    original=set('bosque38,bosque39,bosque48,bosque56,bosque61,bosque75'.split(','))
    probe_nodes=set(recovery['excluded'])-original
    if not args.apply:
        probes=[]
        for node in sorted(probe_nodes):
            job=subprocess.check_output(['sbatch','--parsable','--job-name=beta_startup_probe',
                '--partition=cpunodes_nia',f'--nodelist={node}','--cpus-per-task=2',
                '--mem=1G','--time=00:01:00',f'--output={args.root}/probe_{node}.out',
                '--wrap=hostname'],text=True).strip().split(';')[0]
            probes.append(dict(node=node,job=job))
        (args.root/'node_probes.json').write_text(json.dumps(probes,indent=2)+'\n')
        print(json.dumps(probes),flush=True)
        return
    probes=json.loads((args.root/'node_probes.json').read_text())
    completed=set()
    for row in probes:
        result=subprocess.check_output(['sacct','-X','-n','-P','-j',row['job'],
            '--format=State,ExitCode'],text=True).strip()
        if result=='COMPLETED|0:0':
            completed.add(row['node'])
    excluded=sorted(set(recovery['excluded'])-completed)
    updates=[]
    for item in json.loads((args.root/'submissions.json').read_text()):
        path=Path(item['manifest']);case=json.loads(path.read_text())
        status=json.loads(path.with_name('status.json').read_text())
        if status['stage']=='complete':continue
        (Path(case['run'])/'scheduler_exclusions.json').write_text(json.dumps(excluded)+'\n')
        for line in path.with_name('blind.out').read_text().splitlines():
            if '"submitted"' in line:
                job=json.loads(line)['job']
                update=subprocess.run(['scontrol','update',f'JobId={job}',
                    'ExcNodeList='+','.join(excluded)],capture_output=True,text=True)
                updates.append(dict(job=job,returncode=update.returncode,output=update.stdout+update.stderr))
    (args.root/'restored_nodes.json').write_text(json.dumps(dict(
        verified=sorted(completed),excluded=excluded,updates=updates),indent=2)+'\n')
    print(json.dumps(dict(verified=sorted(completed),excluded=excluded)),flush=True)


if __name__=='__main__':
    main()
