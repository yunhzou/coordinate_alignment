"""Summarize the saved, paired family-scoring experiments."""
import argparse
import csv
import json
from pathlib import Path
import statistics


def stats(rows):
    times=sorted(r['seconds'] for r in rows)
    return dict(count=len(rows),proven=sum(r['optimal'] for r in rows),
        median_seconds=statistics.median(times),maximum_seconds=max(times),
        summed_path_wall_seconds=sum(times),
        invariant_certificates=sum(r['method']=='invariance_certificate' for r in rows),
        zero_solver_queries=sum(r['solver_queries']==0 for r in rows))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    roots={name:a.base/f'family_event_scoring_{name+"_" if name!="initial" else ""}20260908'
           for name in ('initial','blind','pruned','profile')}
    experiments={name:[json.loads(p.read_text()) for p in sorted(root.glob('*_to_*.json'))]
                 for name,root in roots.items()}
    summary={'roots':{k:str(v) for k,v in roots.items()},'experiments':{},'comparison':{},
        'scope':'Five reactions, six directional archives; not a full Golden reranking. '
        'Explicit H; matching tolerance 1.0, scoring tolerance 0.5. '
        'Scoring excludes archive loading, path indexing and first library import. '
        'No full bijection enumeration in production; tiny exhaustive oracles only in tests.'}
    table=[]
    for name,cases in experiments.items():
        rows=[r for case in cases for r in case['results']]
        summary['experiments'][name]=stats(rows)
        if name in ('blind','pruned'):
            summary['experiments'][name]['top8']=stats([r for c in cases for r in c['results'][:8]])
            summary['experiments'][name]['random24']=stats([r for c in cases for r in c['results'][8:]])
        for case in cases:
            for row in case['results']:
                table.append(dict(experiment=name,index=case['index'],direction=case['direction'],
                    terminal=row.get('terminal'),repeat=row.get('repeat'),method=row['method'],
                    representative=row['representative_score'],lower=row['lower_bound'],upper=row['upper_bound'],
                    optimal=row['optimal'],seconds=row['seconds'],solver_queries=row['solver_queries'],
                    archive_loading_seconds=case['archive_loading_seconds']))
    old={(c['index'],c['direction'],r['terminal']):r for c in experiments['blind'] for r in c['results']}
    new={(c['index'],c['direction'],r['terminal']):r for c in experiments['pruned'] for r in c['results']}
    assert old.keys()==new.keys()
    for key,previous in old.items():
        current=new[key]
        if previous['optimal']:assert current['lower_bound']<=previous['upper_bound']<=current['upper_bound']
        if current['optimal']:assert previous['lower_bound']<=current['upper_bound']<=previous['upper_bound']
    summary['comparison']=dict(paired_families=len(old),contradictory_bounds=0,
        pruning_retained=False,reason='No demonstrated speed or proof-completion benefit.')
    from golden_publication import plans
    from publication_analysis import certificate_id
    source=Path(json.loads((roots['initial']/'manifest.json').read_text())['source'])
    checks=[]
    for case in experiments['initial']:
        pair,_=plans(source,case['index']);plan=pair[case['direction']]
        reference=json.loads((source/'inputs'/str(case['index'])/'reference.json').read_text())
        expected=certificate_id(reference['features'],reference['mapping'])
        checks.append(dict(index=case['index'],direction=case['direction'],
            representative_score=case['results'][0]['representative_score'],
            optimal_score=case['results'][0]['upper_bound'],
            reference_equivalent=[certificate_id(reference['features'],plan.to_input_mapping(dict(r['mapping'])))==expected
                                  for r in case['results']]))
    summary['targeted_witness_checks']=checks
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (a.output/'per_family.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    print(json.dumps(summary,indent=2))
