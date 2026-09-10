"""Cheap live recovery/compute summary from saved evaluations, without class loading."""
import argparse
from collections import Counter
import json
from pathlib import Path
import time

from adaptive_full_benchmark import read, save, DIRECTIONS
from publication_analysis import union_outcome


def summarize(run):
    tasks = read(run/'tasks.json')
    specs = {(t['dataset'],t['index']):t for t in tasks}
    rows = []
    for (dataset, index), spec in specs.items():
        methods = {}
        for method in ('original','adaptive'):
            directed = {}
            for direction in DIRECTIONS:
                folder = run/f'results/{dataset}/{index}/{direction}/{method}'
                path = folder/'search.json'
                search = read(path) if path.exists() else {}
                last = search.get('rows', [{}])[-1]
                path = folder/f"{last.get('label','missing')}_evaluation.json"
                evaluation = read(path) if path.exists() else {}
                directed[direction] = dict(search=last, complete=search.get('complete',False), evaluation=evaluation)
            modes = {}
            for mode, directions in (('single',[spec['smaller_first']]), ('bidirectional',DIRECTIONS)):
                values = [directed[d] for d in directions]
                modes[mode] = dict(complete=all(v['complete'] and v['evaluation'] for v in values),
                    cpu=sum(v['search'].get('compute_cpu_excluding_persistence_and_loading_seconds',0)
                            for v in values if v['complete']))
                if dataset=='golden':
                    modes[mode]['outcome'] = union_outcome([v['evaluation'].get('reference_recovery','unknown')
                                                          for v in values])
                else:
                    scores = [v['evaluation'].get('best_events') for v in values]
                    modes[mode]['best'] = min((s for s in scores if s is not None), default=None)
            methods[method] = modes
        rows.append(dict(dataset=dataset,index=index,methods=methods))
    totals = {}
    for dataset in ('golden','holdout'):
        cases = [r for r in rows if r['dataset']==dataset]
        for mode in ('single','bidirectional'):
            summary = {}
            for method in ('original','adaptive'):
                values = [r['methods'][method][mode] for r in cases]
                summary[method] = dict(complete=sum(bool(v['complete']) for v in values),
                    cpu=sum(v['cpu'] for v in values))
                if dataset=='golden':
                    summary[method]['outcomes'] = dict(Counter(v['outcome'] for v in values))
            if dataset=='golden':
                for outcome in ('not_recovered','unknown'):
                    summary['previously_recovered_now_'+outcome] = [r['index'] for r in cases
                        if r['methods']['original'][mode]['outcome']=='recovered'
                        and r['methods']['adaptive'][mode]['outcome']==outcome]
            else:
                summary['worse'] = [r['index'] for r in cases if None not in (
                    a:=r['methods']['original'][mode]['best'], b:=r['methods']['adaptive'][mode]['best']) and b>a]
                summary['equal'] = sum(r['methods']['original'][mode]['best']==r['methods']['adaptive'][mode]['best']
                                      and r['methods']['original'][mode]['best'] is not None for r in cases)
            totals[dataset+'_'+mode] = summary
    result = dict(updated=time.time(), totals=totals,
        scope='Recovery and saved representative best events only. Full class-equivalence report is separate.')
    save(run/'recovery_progress.json', result)
    save(run/'recovery_per_case.json', rows)
    return result


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    args = parser.parse_args()
    print(json.dumps(summarize(args.run),indent=2))
