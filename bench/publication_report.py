"""Publication tables from frozen benchmark outputs; never reruns AAM."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import statistics
import subprocess

from golden_policy_campaign import save
from golden_publication import report as refresh


def distribution(values):
    values=sorted(v for v in values if v is not None)
    if not values:return dict(count=0)
    def percentile(q):
        index=(len(values)-1)*q
        lo=int(index);hi=min(lo+1,len(values)-1)
        return values[lo]+(values[hi]-values[lo])*(index-lo)
    return dict(count=len(values),sum=sum(values),mean=statistics.mean(values),median=statistics.median(values),
                p90=percentile(.9),p95=percentile(.95),maximum=max(values))


def main(args):
    refresh(args)
    manifest=json.loads((args.run/'manifest.json').read_text())
    cases=json.loads((args.run/'cases.json').read_text())
    # A finished mode with no ranking is unresolved, not still running.
    summary=json.loads((args.run/'summary.json').read_text())
    for name in ('single','bidirectional'):
        records=[c.get('modes',{}).get(name,{}) for c in cases]
        summary['modes'][name]['top5_family_outcomes']=dict(Counter(
            m.get('top5_family',{}).get('outcome','unknown' if m else 'pending') for m in records))
    save(args.run/'summary.json',summary)
    out=args.run/'publication';out.mkdir(exist_ok=True)
    table=[]
    for case in cases:
        row=dict(index=case['index'],default_direction=case.get('default_direction'))
        for name in ('single','bidirectional'):
            m=case.get('modes',{}).get(name,{})
            row[name+'_recovery']=m.get('reference_recovery','pending')
            row[name+'_ranking_complete']=m.get('ranking_complete',False)
            row[name+'_search_complete']=m.get('search_complete',False)
            row[name+'_top5_family']=m.get('top5_family',{}).get('outcome','pending' if not m else 'unknown')
            row[name+'_compute_cpu_seconds']=m.get('compute_cpu_excluding_io_seconds')
            row[name+'_directional_elapsed_sum_seconds']=m.get('directional_elapsed_including_io_sum_seconds')
            representative=m.get('representative',{})
            row[name+'_class_count']=representative.get('class_count')
            for k in (1,3,5,10):row[f'{name}_top{k}']=representative.get('topk',{}).get(str(k))
            for tolerance in (0,1,2,3,5):
                item=representative.get('event_windows',{}).get(str(tolerance),{})
                row[f'{name}_event_window_{tolerance}_recovered']=item.get('recovered')
                row[f'{name}_event_window_{tolerance}_candidates']=item.get('count')
        table.append(row)
    with (out/'per_case.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    denominator=len(cases)
    analysis=dict(denominator=denominator,modes={},paired={})
    for name in ('single','bidirectional'):
        counts=Counter(r[name+'_recovery'] for r in table)
        unranked=sum(not r[name+'_ranking_complete'] for r in table)
        modes=dict(recovery_outcomes=dict(counts),certified_recovery_percent=100*counts['recovered']/denominator,
            upper_bound_percent=100*(counts['recovered']+counts['unknown']+counts['pending'])/denominator,
            incomplete_rankings=unranked,
            top5_family_outcomes=dict(Counter(r[name+'_top5_family'] for r in table)),
            compute_cpu_completed_cases=distribution([r[name+'_compute_cpu_seconds'] for r in table]),
            class_counts=distribution([r[name+'_class_count'] for r in table]),topk={},event_windows={})
        for k in (1,3,5,10):
            count=sum(r[f'{name}_top{k}'] is True for r in table)
            modes['topk'][str(k)]=dict(correct=count,total=denominator,percent=100*count/denominator,
                upper_bound_percent=100*(count+unranked)/denominator)
        for tolerance in (0,1,2,3,5):
            count=sum(r[f'{name}_event_window_{tolerance}_recovered'] is True for r in table)
            modes['event_windows'][str(tolerance)]=dict(correct=count,total=denominator,percent=100*count/denominator,
                candidates=distribution([r[f'{name}_event_window_{tolerance}_candidates'] for r in table]))
        analysis['modes'][name]=modes
    analysis['paired']=dict(
        additional_certified_cases=[r['index'] for r in table if r['single_recovery']!='recovered' and r['bidirectional_recovery']=='recovered'],
        impossible_union_regressions=[r['index'] for r in table if r['single_recovery']=='recovered' and r['bidirectional_recovery'] not in ('recovered','pending')],
        top1_gains=[r['index'] for r in table if r['single_top1'] is False and r['bidirectional_top1'] is True],
        top1_losses=[r['index'] for r in table if r['single_top1'] is True and r['bidirectional_top1'] is False])
    directions=[]
    for task in manifest['tasks']:
        base=args.run/'directions'/str(task['index'])/task['direction']
        def read(name):
            path=base/name
            return json.loads(path.read_text()) if path.exists() else {}
        search,ranking,evaluation=read('search.json'),read('ranking.json'),read('evaluation.json')
        env=read('environment.json')
        status_path=args.run/'status'/f'{task["index"]}_{task["direction"]}_search.json'
        status=json.loads(status_path.read_text()) if status_path.exists() else {}
        phases=search.get('phases',{})
        directions.append(dict(**task,host=env.get('host'),cpu_models=';'.join(env.get('cpu_models',[])),
            search_process_exit=status.get('exit'),search_complete=bool(search),capped=search.get('capped'),
            terminals=search.get('terminals'),cut_count=search.get('metrics',{}).get('cut_count'),
            compute_cpu_seconds=search.get('compute_cpu_excluding_persistence_and_loading_seconds'),
            elapsed_wall_including_io=search.get('elapsed_wall_including_io_seconds'),
            persistence_cpu_seconds=sum(v['cpu_seconds'] for k,v in phases.items() if k.startswith('persistence/')) if search else None,
            persistence_summed_wall_seconds=sum(v['summed_wall_seconds'] for k,v in phases.items() if k.startswith('persistence/')) if search else None,
            loading_cpu_seconds=sum(v['cpu_seconds'] for k,v in phases.items() if k.startswith('loading/')) if search else None,
            ranking_cpu_seconds=ranking.get('cpu_seconds'),ranking_wall_seconds=ranking.get('wall_seconds'),
            ranking_saving_wall_seconds=ranking.get('saving_wall_seconds'),
            verification_cpu_seconds=evaluation.get('verification_cpu_seconds'),verification_wall_seconds=evaluation.get('verification_wall_seconds'),
            reference_recovery=evaluation.get('reference_recovery','pending')))
    with (out/'per_direction.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(directions[0]));writer.writeheader();writer.writerows(directions)
    analysis['directional_searches']=dict(total=len(directions),completed=sum(r['search_complete'] for r in directions),
        process_exits=dict(Counter(str(r['search_process_exit']) for r in directions)),
        cap_hits=sum(r['capped'] is True for r in directions),
        compute_cpu_completed=distribution([r['compute_cpu_seconds'] for r in directions]),
        elapsed_wall_completed=distribution([r['elapsed_wall_including_io'] for r in directions]),
        persistence_cpu=distribution([r['persistence_cpu_seconds'] for r in directions]),
        hardware=dict(Counter(r['cpu_models'] for r in directions)))
    analysis['interpretation']='Fixed full benchmark, not union with prior experiments. Every record remains in denominators. '
    analysis['interpretation']+='Timing distributions describe completed cases and list their count; timeouts are censored, not zero. '
    analysis['interpretation']+='Net CPU removes measured archive persistence/loading and retains IPC/metadata/profiling overhead. '
    analysis['interpretation']+='Directional elapsed sums are NOT measured concurrent bidirectional latency. '
    analysis['interpretation']+='Representative and family-recovery metrics are distinct. Event-window thresholds were fixed before the run.'
    jobs=json.loads((args.run/'jobs.json').read_text())
    recovery_path=args.run/'reporting/recovery.json'
    recovery=json.loads(recovery_path.read_text()) if recovery_path.exists() else {}
    analysis['scheduler_recovery']=recovery
    accounting_jobs=[str(j['job']) for j in jobs]+[str(j['job']) for j in recovery.get('replacement_jobs',[])]
    accounting=subprocess.run(['sacct','-j',','.join(accounting_jobs),'-P',
        '--format=JobID,State,ElapsedRaw,TotalCPU,CPUTimeRAW,AllocCPUS,MaxRSS,NodeList,ExitCode,Submit,Start,End'],
        text=True,capture_output=True,timeout=60)
    (out/'slurm_accounting.psv').write_text(accounting.stdout)
    analysis['accounting_export']=dict(exit=accounting.returncode,stderr=accounting.stderr)
    save(out/'analytics.json',analysis)
    text=['# Frozen Golden publication benchmark','',analysis['interpretation'],'',
          '| Mode | Certified reference recovery | Unknown | Pending | Representative top-1 |',
          '|---|---:|---:|---:|---:|']
    for name,m in analysis['modes'].items():
        text.append(f"| {name} | {m['recovery_outcomes'].get('recovered',0)}/{denominator} ({m['certified_recovery_percent']:.3f}%) | "
                    f"{m['recovery_outcomes'].get('unknown',0)} | {m['recovery_outcomes'].get('pending',0)} | {m['topk']['1']['percent']:.3f}% |")
    text += ['',f"Frozen engine commit: `{manifest['git_commit']}`.",'',
             'See `analytics.json`, `per_case.csv`, and `per_direction.csv` for full counts, timing, caps, event windows, and paired comparisons.',
             '', 'The parent directory contains the full input/source/dependency manifests, checkpoint archives, evaluation witnesses, and per-process timing logs.']
    (out/'README.md').write_text('\n'.join(text)+'\n')
    print(str(out.resolve()),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    main(p.parse_args())
