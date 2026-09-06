"""Export every Golden case, including failures, without loading/rerunning AAM."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import statistics


def read(path):
    return json.loads(path.read_text()) if path.exists() else {}


def duration(value):
    days=0
    if '-' in value:
        days,value=value.split('-');days=int(days)
    return days*86400+sum(float(v)*60**i for i,v in enumerate(reversed(value.split(':'))))


def stats(values):
    values=sorted(v for v in values if v is not None)
    return dict(n=len(values),total=sum(values),median=statistics.median(values),
                p95=values[min(len(values)-1,int(.95*len(values)))],maximum=values[-1]) if values else dict(n=0)


def summarize(rows):
    return dict(records=len(rows),evaluated=sum(r['evaluated'] for r in rows),
        family_recovery=dict(Counter(r.get('reference_recovery','unscored') for r in rows)),
        top1_eligible=sum(r.get('top1_correct') is not None for r in rows),
        top1_correct=sum(r.get('top1_correct') is True for r in rows),
        representative_recovered=sum(r.get('representative_recovery') is True for r in rows),
        full_P_heavy=sum(r.get('best_target_heavy_coverage')==1 for r in rows),
        full_P_explicit_H_included=sum(r.get('best_target_all_atom_coverage')==1 for r in rows),
        capped=sum(r.get('capped') is True for r in rows),
        incomplete_searches=sum(r['search_incomplete'] for r in rows),
        stages=dict(Counter(r['stage'] for r in rows)))


def export(run, output):
    output.mkdir(parents=True,exist_ok=True)
    manifest=read(run/'manifest.json');rows=[]
    for item in manifest['records']:
        directory=run/str(item['index'])
        state=read(directory/'status.json');search=read(directory/'search.json')
        evaluation=read(directory/'evaluation.json');reference=read(directory/'reference.json')
        heavy=set(reference['features'][1]['heavy'])
        annotated={p for r,p in reference['mapping'] if p in heavy}
        row=dict(index=item['index'],balanced=item['balanced'],source_atoms=item['source_atoms'],
            target_atoms=item['target_atoms'],target_heavy_atoms=len(heavy),
            reference_complete=len(annotated)==len(heavy),reference_target_coverage=len(annotated)/len(heavy),
            stage=state['stage'],evaluated=bool(evaluation),
            search_incomplete=bool(evaluation.get('search_incomplete') or (directory/'search_timeout.json').exists()),
            reused_search=(directory/'reuse.json').exists(),
            archive=str(directory/'aam.json.gz') if (directory/'aam.json.gz').exists() else '',
            partial_archive=str(directory/'partial_archive.json') if (directory/'partial_archive.json').exists() else '',
            evaluation=str(directory/'evaluation.json') if evaluation else '')
        row.update({k:v for k,v in evaluation.items() if not isinstance(v,(dict,list))})
        row.update({'search_'+k:v for k,v in search.items() if k!='metrics'})
        row.update({'metric_'+k:v for k,v in search.get('metrics',{}).items()})
        row.update({'top_'+k:v for k,v in evaluation.get('top_events',{}).items()})
        rows.append(row)
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with (output/'cases.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    summary=dict(run=str(run),search_revision=manifest['revision'],configuration=manifest['config'],
        all_records=summarize(rows),
        complete_reference=summarize([r for r in rows if r['reference_complete']]),
        partial_reference=summarize([r for r in rows if not r['reference_complete']]),
        balanced=summarize([r for r in rows if r['balanced']]),
        unbalanced=summarize([r for r in rows if not r['balanced']]),
        completed_searches=summarize([r for r in rows if r['archive']]),
        interrupted_searches=summarize([r for r in rows if r['search_incomplete']]),
        reused_searches=sum(r['reused_search'] for r in rows),
        search_wall_seconds=stats([r.get('search_wall_seconds') for r in rows]),
        archive_seconds=stats([r.get('search_archive_seconds') for r in rows]),
        evaluation_seconds=stats([r.get('evaluation_seconds') for r in rows]),
        completed_search_cpu_seconds=stats([r.get('search_cpu_seconds') for r in rows]),
        completed_search_peak_rss_mb=stats([r.get('search_peak_rss_mb') for r in rows]),
        nodes=sorted({r['search_hostname'] for r in rows if 'search_hostname' in r}),
        cpu_models=sorted({r['search_cpu_model'] for r in rows if 'search_cpu_model' in r}))
    for key in ('worker_search_seconds','checkpoint_seconds','parent_merge_seconds','symmetry_finalization_seconds'):
        summary[key]=stats([r.get('metric_'+key) for r in rows])
    accounting=run/'slurm_accounting.psv'
    if accounting.exists():
        tasks=[line.split('|') for line in accounting.read_text().splitlines()
               if '_' in line.split('|')[0] and '.' not in line.split('|')[0]]
        summary['slurm']=dict(tasks=len(tasks),states=dict(Counter(t[1] for t in tasks)),
            measured_cpu_seconds=sum(duration(t[4]) for t in tasks if t[4]),
            allocated_cpu_seconds=sum(float(t[3])*int(t[6]) for t in tasks if t[3] and t[6]),
            note='Parent array tasks only: includes search, scoring, serialization and timed-out jobs; excludes cached pilot searches and supervisor. Accounting can lag task completion.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    all_=summary['all_records'];complete=summary['complete_reference'];partial=summary['partial_reference']
    lines=['# Golden AAM benchmark — measured campaign', '',
        f"Run: `{run}`. Search revision: `{manifest['revision']}`. Frozen source hashes and inputs: `manifest.json`.", '',
        '## Recovery and coverage', '',
        '| Measure | Count |','|---|---:|',
        f"| Records / evaluated | {len(rows)} / {all_['evaluated']} |",
        f"| Complete heavy-atom reference annotations | {complete['records']} |",
        f"| Complete-reference family recovery | {complete['family_recovery'].get('recovered',0)} / {complete['records']} |",
        f"| Selected-representative correct / scored eligible | {complete['top1_correct']} / {complete['top1_eligible']} |",
        f"| Partial-reference consistency (not full mapping accuracy) | {partial['family_recovery'].get('recovered',0)} / {partial['records']} |",
        f"| One mapping covers all P heavy atoms | {all_['full_P_heavy']} / {len(rows)} |",
        f"| One mapping covers all P atoms, including explicit H | {all_['full_P_explicit_H_included']} / {len(rows)} |",
        f"| Incomplete searches | {all_['incomplete_searches']} |", '',
        'Unresolved scoring and unscored records remain in `summary.json` and `cases.csv`; they are not counted as successes. '
        'A positive recovery from a completed cut is sound even when later search timed out. Its coverage is a lower bound. '
        'Cap-limited and seed-limited searches do not prove global exhaustiveness. '
        'Coverage is not a union of incompatible branches and does not establish chemical correctness.', '',
        '## Measured performance', '',
        '| Phase | Completed measurements | Median (s) | p95 (s) | Maximum (s) |',
        '|---|---:|---:|---:|---:|']
    for label,key in [('Search including checkpoints','search_wall_seconds'),('Additional gzip archive','archive_seconds'),('Evaluation','evaluation_seconds')]:
        s=summary[key]
        if s['n']:lines.append(f"| {label} | {s['n']} | {s['median']:.3f} | {s['p95']:.3f} | {s['maximum']:.3f} |")
    lines+=['',f"Search watchdog: {manifest['search_timeout']} s; outer task watchdog: 600 s. CPU allocation ceiling: {manifest['cpu_budget']}; one AAM thread per reaction. {summary['reused_searches']} searches reused from the saved pilot.", '',
        'Completed-search latency excludes timed-out searches and additional gzip serialization; report those separately, not as a cold full-run speed claim. '
        'Phase metrics are diagnostic, can overlap, and should not be blindly added. Slurm measured CPU and allocated CPU-seconds are different quantities. '
        'Peak RSS in the case table is the search process high-water mark; raw scheduler accounting also retains scoring/task memory.', '',
        '## Evaluation scope', '',
        'Full explicit-H graphs are mapped, with original reference map labels removed before canonical input ordering. '
        'Heavy-atom reference agreement is assessed modulo joint, chemically colored endpoint automorphisms. '
        'Family queries preserve correlated group actions and do not enumerate all bijections. '
        'Top-1 refers to the reference-blind selected terminal representative, not the best possible member of every compressed family. '
        'Reference annotations do not establish identities for automatically added H atoms; all-atom coverage is not H-mapping accuracy. '
        'This is graph AAM with discrete bond orders, not a 3D stereo-validation or mechanism-grouping benchmark.', '',
        'No external baseline has been rerun here. Published Golden subsets/reference corrections differ; these results are not yet a like-for-like paper comparison.', '',
        '## Saved evidence', '',
        '- Every input, separate reference, status, search log and available complete compressed AAM archive.',
        '- Atomic per-cut checkpoints, including for interrupted searches; partial manifests identify their exact files.',
        '- Per-case evaluations and symbolic witness actions when required; source/configuration hashes and Slurm accounting.',
        '- `cases.csv`: every dataset row, including failures. `summary.json`: stratified denominators and timings.', '',
        'Recreate this report without remapping: `python bench/report_golden_campaign.py --run RUN --output OUTPUT`.',
        'Rescore a completed archive without remapping: `PYTHONPATH=RUN/engine/src python RUN/engine/bench/golden_campaign.py score --run RUN --index INDEX`. '
        'Preserve the previous evaluation before intentionally rescoring; use the frozen engine for exact reproducibility.', '']
    (output/'README.md').write_text('\n'.join(lines))
    print(json.dumps(summary['all_records']))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();export(args.run.resolve(),args.output)
