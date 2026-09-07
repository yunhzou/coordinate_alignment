"""Report all Golden records, retaining missing/partial outcomes in denominators."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import statistics


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
        complete_reference_total=len(complete),complete_reference_outcomes=dict(Counter(r['reference_recovery'] for r in complete)),
        complete_reference_recovered=sum(r['reference_recovery']=='recovered' for r in complete),
        complete_reference_top1=sum(r['top1_correct'] is True for r in complete),
        recovered_percent=100*sum(r['reference_recovery']=='recovered' for r in complete)/len(complete),
        top1_percent=100*sum(r['top1_correct'] is True for r in complete)/len(complete),
        gains=[r['index'] for r in complete if r['recovery_change']=='gained'],
        losses=[r['index'] for r in complete if r['recovery_change']=='lost'],
        partial_reference_records=len(rows)-len(complete),
        incomplete_searches=sum(r['search_incomplete'] is True for r in rows),
        completed_searches=len(times),search_median_seconds=statistics.median(times) if times else None,
        search_max_seconds=max(times,default=None),
        completed_search_cpu_hours=sum(r['cpu_seconds'] or 0 for r in rows)/3600,
        peak_rss_mb=max((r['peak_rss_mb'] or 0 for r in rows),default=0),
        note='Accuracy uses all complete-reference records, including pending/errors/unknown. Until the run settles these percentages are lower bounds, not final scores. CPU totals omit interrupted searches; use Slurm accounting for allocation totals.')
    args.output.mkdir(parents=True,exist_ok=True)
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
    main(parser.parse_args())
