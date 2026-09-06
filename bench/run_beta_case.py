"""Execute one checkpointed blind case, then independent reference validation."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,required=True)
    parser.add_argument('--cpu-budget',type=int,required=True)
    parser.add_argument('--worker-cpus',type=int,default=8)
    parser.add_argument('--shards',type=int,default=256)
    args=parser.parse_args()
    case=json.loads(args.case.read_text())
    run=Path(case['run'])
    status_path=args.case.with_name('status.json')
    started=time.time()
    if status_path.exists():
        previous=json.loads(status_path.read_text())
        history=args.case.parent/'attempts'/str(int(time.time()))
        history.mkdir(parents=True)
        for path in (*args.case.parent.glob('*.out'), *args.case.parent.glob('*.err'), status_path):
            # Slurm may already have opened the coordinator logs for this attempt.
            if not path.name.startswith('coordinator'):
                path.rename(history/path.name)
        started=previous.get('started',started)
    def status(stage,**extra):
        status_path.write_text(json.dumps(dict(stage=stage,started=started,
            elapsed=time.time()-started,case=case['id'],**extra),indent=2)+'\n')
    def execute(command,label,timeout=None):
        status(label)
        with args.case.with_name(label+'.out').open('w') as out, args.case.with_name(label+'.err').open('w') as err:
            return subprocess.run([sys.executable,*command],stdout=out,stderr=err,timeout=timeout)
    status('starting')
    search=execute(['bench/run_beta_distributed.py','--run',str(run),
        '--catalog',case['catalog'],'--target-smiles',case['target_smiles'],
        '--shards',str(args.shards),'--cpu-budget',str(args.cpu_budget),
        '--worker-cpus',str(args.worker_cpus),'--recommendations','20','--patterns','4',
        '--exclude',case.get('exclude','bosque38,bosque39,bosque48,bosque56,bosque61,bosque75'),
        '--result-stem','assemblies_eight'],'blind')
    if search.returncode:
        status('failed',phase='blind',exit_code=search.returncode)
        return
    scan_summary=args.case.with_name('scan_summary.json')
    if case.get('scan_summary'):
        scan_summary=Path(case['scan_summary'])
    else:
        parts=[json.loads(p.read_text()) for p in (run/'query_full/parts').glob('*.json')]
        scan_summary.write_text(json.dumps(dict(scan_counts=dict(
            matched_precursors=sum(p['matched_sources'] for p in parts))))+'\n')
    view=['bench/view_beta_result.py','--run',str(run),'--result-stem','assemblies_eight',
        '--case',str(args.case),'--scan-summary',str(scan_summary),
        '--html-output',str(args.case.with_name('viewer.html')),'--title',case['id']]
    rendered=execute(view,'viewer',580)
    if rendered.returncode:
        status('failed',phase='viewer',exit_code=rendered.returncode)
        return
    # A preview exists even if the independent reference-set check exceeds its watchdog.
    try:
        validation=execute(['bench/assemble_beta_known_sources.py','--run',str(run),
            '--case',str(args.case),'--result-stem','reference_eight'],'validation',580)
    except subprocess.TimeoutExpired:
        status('complete',html=str(args.case.with_name('viewer.html')),validation='timeout; blind viewer available')
        return
    if validation.returncode == 0:
        rendered=execute(view+['--validation-stem','reference_eight'],'viewer_validated',580)
        if rendered.returncode:
            status('failed',phase='viewer_validated',exit_code=rendered.returncode)
            return
    status('complete',html=str(args.case.with_name('viewer.html')),validation_exit=validation.returncode)


if __name__=='__main__':
    main()
