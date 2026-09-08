"""Separate observed candidate recovery from exhaustive extraction completion."""
import argparse
from collections import Counter
import json
from pathlib import Path
import shutil
import subprocess
from golden_policy_campaign import save


def report(run):
    manifest=json.loads((run/'manifest.json').read_text());indices=manifest['indices']
    states=Counter();directions={};rows=[]
    for i in indices:
        entries={}
        for d in ('R_to_P','P_to_R'):
            folder=run/'results'/str(i)/d;path=folder/'summary.json'
            s=json.loads(path.read_text()) if path.exists() else dict(stage='pending')
            status_path=run/'logs'/f'{i}_{d}.status.json'
            if status_path.exists():
                status=json.loads(status_path.read_text())['status']
                if status!='ok':s=dict(s,stage=status)
            states[s['stage']]+=1
            if s['stage']=='finished':entries[d]=(s,json.loads((folder/'candidates.json').read_text()))
            directions[i,d]=s
        default=next((s['default_direction'] for s,_ in entries.values()),None)
        for mode,wanted in [('single',[default]),('bidirectional',['R_to_P','P_to_R'])]:
            available=[entries[d] for d in wanted if d in entries]
            candidates=[r for _,rs in available for r in rs]
            best=min(((-r['heavy'],-r['total'],r['events']['total']) for r in candidates),default=None)
            eligible=[r for r in candidates if best and (-r['heavy'],-r['total'])==best[:2]]
            rows.append(dict(index=i,mode=mode,directions_finished=len(available),
                all_directions_finished=len(available)==len(wanted),
                exhaustive=len(available)==len(wanted) and all(s['complete'] for s,_ in available),
                baseline_detected=any(s['baseline_reference'] for s,_ in available),
                detected=any(r['reference_equivalent'] for r in candidates),
                event_windows={str(n):any(r['reference_equivalent'] and r['events']['total']<=best[2]+n
                                         for r in eligible) for n in range(21)}))
    result=dict(denominator=len(indices),directional_stages=dict(states),modes={},
                attempted_directions=sum(n for stage,n in states.items()
                                         if stage not in ('pending','loading','extracting')))
    for mode in ('single','bidirectional'):
        selected=[r for r in rows if r['mode']==mode]
        result['modes'][mode]=dict(finished=sum(r['all_directions_finished'] for r in selected),
            exhaustive=sum(r['exhaustive'] for r in selected),baseline_detected=sum(r['baseline_detected'] for r in selected),
            detected=sum(r['detected'] for r in selected),
            event_windows={str(n):sum(r['event_windows'][str(n)] for r in selected) for n in range(21)})
    done=[s for s in directions.values() if s['stage']=='finished']
    result['totals']={k:sum(s.get(k,0) for s in done) for k in ('visited_families','unresolved_families',
        'duplicate_histories','new_candidates','archive_loading_seconds','extraction_wall_seconds',
        'persistence_wall_seconds','reference_verification_wall_seconds')}
    result['interpretation']='Detection is a lower bound until extraction completes; unfinished directions remain in the 1851 denominator. No new AAM search.'
    save(run/'report.json',result);save(run/'per_case.json',rows);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path);p.add_argument('--jobs',help='Comma-separated Slurm IDs for accounting export')
    args=p.parse_args();result=report(args.run)
    if args.output:
        args.output.mkdir(parents=True,exist_ok=True)
        save(args.output/'benchmark_report.json',result)
        for name in ('per_case.json','manifest.json','repair_manifest.json','parallel_repair_manifest.json','exclusion_comparison.json'):
            shutil.copy2(args.run/name,args.output/name)
        if args.jobs:
            accounting=subprocess.check_output(['sacct','-j',args.jobs,'--parsable2',
                '--format=JobID,JobName,State,Start,End,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS'],text=True)
            (args.output/'slurm_accounting.tsv').write_text(accounting)
    print(json.dumps(result,indent=2))
