"""Report all Golden records, retaining missing/partial outcomes in denominators."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import statistics
import subprocess


def main(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    old={int(r['index']):r for r in csv.DictReader(args.baseline.open())}
    rows=[]
    for item in manifest['records']:
        directory=args.run/str(item['index'])
        status=json.loads((directory/'status.json').read_text())
        evaluation=json.loads((directory/'evaluation.json').read_text()) if (directory/'evaluation.json').exists() else {}
        search=json.loads((directory/'search.json').read_text()) if (directory/'search.json').exists() else {}
        previous=old.get(item['index'],{})
        recovered=evaluation.get('reference_recovery','pending' if status['stage'] in ('pending','search','score','partial') else 'unresolved')
        baseline=previous.get('verified_reference_recovery','')
        change='gained' if recovered=='recovered' and baseline=='not_recovered' else 'lost' if recovered=='not_recovered' and baseline=='recovered' else ''
        rows.append(dict(index=item['index'],reference_complete=item['reference_annotation_complete'],
            direction=item['direction'],stage=status['stage'],reference_recovery=recovered,
            top1_correct=evaluation.get('top1_correct'),search_incomplete=evaluation.get('search_incomplete'),
            capped=evaluation.get('capped'),search_seconds=search.get('wall_seconds'),
            evaluation_seconds=evaluation.get('evaluation_seconds'),cpu_seconds=search.get('cpu_seconds'),
            peak_rss_mb=search.get('peak_rss_mb'),old_reference_recovery=baseline,
            old_top1_correct=previous.get('original_top1_correct',''),recovery_change=change,
            archive=str(directory/'cuts/aam.pkl.gz') if (directory/'cuts/aam.pkl.gz').exists() else str(directory/'cuts'),
            evaluation=str(directory/'evaluation.json')))
    complete=[r for r in rows if r['reference_complete']]
    times=[r['search_seconds'] for r in rows if r['search_seconds'] is not None]
    summary=dict(total=len(rows),states=dict(Counter(r['stage'] for r in rows)),
        settled=all(r['stage']=='complete' for r in rows),
        complete_reference_total=len(complete),complete_reference_outcomes=dict(Counter(r['reference_recovery'] for r in complete)),
        complete_reference_recovered=sum(r['reference_recovery']=='recovered' for r in complete),
        complete_reference_top1=sum(r['top1_correct'] is True for r in complete),
        recovered_percent=100*sum(r['reference_recovery']=='recovered' for r in complete)/len(complete),
        top1_percent=100*sum(r['top1_correct'] is True for r in complete)/len(complete),
        gains=[r['index'] for r in complete if r['recovery_change']=='gained'],
        losses=[r['index'] for r in complete if r['recovery_change']=='lost'],
        not_recovered=[r['index'] for r in complete if r['reference_recovery']=='not_recovered'],
        unknown=[r['index'] for r in complete if r['reference_recovery']=='unknown'],
        partial_search_indices=[r['index'] for r in rows if r['search_incomplete'] is True],
        partial_reference_records=len(rows)-len(complete),
        incomplete_searches=sum(r['search_incomplete'] is True for r in rows),
        completed_searches=len(times),search_median_seconds=statistics.median(times) if times else None,
        search_max_seconds=max(times,default=None),
        completed_search_cpu_hours=sum(r['cpu_seconds'] or 0 for r in rows)/3600,
        completed_search_seconds=sum(times),
        evaluation_seconds=sum(r['evaluation_seconds'] or 0 for r in rows),
        peak_rss_mb=max((r['peak_rss_mb'] or 0 for r in rows),default=0),
        note='Accuracy uses all complete-reference records, including unresolved cases. Recovered partial-search witnesses count as positive evidence; their full-sweep top-1 is unknown. An unsettled report is interim. Completed-search CPU totals omit interrupted searches; Slurm allocation totals are separate.')
    args.output.mkdir(parents=True,exist_ok=True)
    if args.accounting:
        submission=json.loads((args.run/'submission.json').read_text())
        recovery=args.run/'resubmission_unstarted.json'
        if recovery.exists():
            resubmission=json.loads(recovery.read_text())
            submission['jobs']+=resubmission['jobs']
            (args.output/recovery.name).write_text(json.dumps(resubmission,indent=2)+'\n')
        jobs=','.join(j['job'] for j in submission['jobs'])
        accounting=subprocess.check_output(['sacct','-j',jobs,'-P',
            '--format=JobID,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,Start,End'],text=True)
        (args.output/'slurm_accounting.psv').write_text(accounting)
        allocations=[r for r in csv.DictReader(accounting.splitlines(),delimiter='|') if '.' not in r['JobID']]
        summary['allocation_states']=dict(Counter(r['State'] for r in allocations))
        summary['allocated_cpu_hours']=sum(int(r['ElapsedRaw'])*int(r['AllocCPUS']) for r in allocations)/3600
        summary['allocation_note']='Main campaign allocations include interrupted searches and process/startup overhead; separate saved-cut reevaluation jobs are excluded.'
    for name,data in [('summary.json',summary),('manifest.json',manifest)]:
        (args.output/name).write_text(json.dumps(data,indent=2)+'\n')
    with (args.output/'cases.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,default=Path('reports/golden_mapping_diagnosis_20260907/unchanged_search_reference_status.csv'))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--accounting',action='store_true')
    main(parser.parse_args())
