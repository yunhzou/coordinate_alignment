"""Keep fixed-archive coverage distinct from a union with search ablations."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import subprocess
from golden_policy_campaign import save


def main(args):
    manifest=json.loads((args.run/'manifest.json').read_text());rows=[];witnesses=[]
    for index in manifest['indices']:
        directory=args.run/str(index);e=json.loads((directory/'evaluation.json').read_text())
        s=json.loads((directory/'status.json').read_text())
        rows.append(dict(index=index,reference_recovery=e['reference_recovery'],
            seconds=e.get('total_seconds'),queries=e.get('symbolic_queries'),
            rejected_heavy_projections=e.get('rejected_heavy_projections'),
            search_incomplete=e.get('search_incomplete'),exit_code=s['exit_code'],
            evaluation=str(directory/'evaluation.json')))
        if e['reference_recovery']=='recovered':witnesses.append(dict(index=index,**e))
    positives=[r['index'] for r in rows if r['reference_recovery']=='recovered']
    old=manifest['cached_complete_reference_positives'];denominator=manifest['complete_reference_total']
    opposite=json.loads(args.diagnostics.read_text())['opposite_recovered']
    combined=old+len(set(positives)|set(opposite))
    summary=dict(cached_verified_positives=old,rescored=len(rows),new_recovered=positives,
        fixed_archive_recovered=old+len(positives),denominator=denominator,
        fixed_archive_percent=100*(old+len(positives))/denominator,
        outcomes=dict(Counter(r['reference_recovery'] for r in rows)),
        not_recovered=[r['index'] for r in rows if r['reference_recovery']=='not_recovered'],
        unknown=[r['index'] for r in rows if r['reference_recovery']=='unknown'],
        combined_with_opposite_recovered=combined,combined_with_opposite_percent=100*combined/denominator,
        warning='Combined coverage is an explicitly separate union of different search directions, not fixed-policy benchmark accuracy.')
    args.output.mkdir(parents=True,exist_ok=True)
    save(args.output/'summary.json',summary);save(args.output/'manifest.json',manifest)
    save(args.output/'witnesses.json',witnesses)
    with (args.output/'cases.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    if args.jobs:
        accounting=subprocess.check_output(['sacct','-j',args.jobs,'-P',
            '--format=JobID,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,Start,End'],text=True)
        (args.output/'slurm_accounting.psv').write_text(accounting)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--diagnostics',type=Path,required=True)
    p.add_argument('--jobs');main(p.parse_args())
