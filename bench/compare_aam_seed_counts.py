"""Compare saved seed-count ablations; never launches a mapper or evaluator."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def aggregate(variants,slap):
    labels=sorted(variants,key=lambda label:int(label.removeprefix('seeds')))
    keys=set(variants[labels[0]])
    assert all(set(variants[label])==keys for label in labels)
    golden=sorted(index for dataset,index in keys if dataset=='golden')
    holdout=sorted(index for dataset,index in keys if dataset=='holdout')
    assert set(golden)==set(slap)
    common=[index for index in golden if all(variants[label]['golden',index]['searches_complete'] for label in labels)
            and not slap[index]['incomplete_variants'] and not slap[index]['mapping_errors']
            and not slap[index]['invalid_predictions']]
    methods={}
    baseline=labels[-1]
    for label in labels:
        values=[variants[label]['golden',i] for i in golden]
        outcomes=Counter(v['recovery'] for v in values)
        comparable=[i for i in holdout if variants[label]['holdout',i]['best_events'] is not None
                    and variants[baseline]['holdout',i]['best_events'] is not None]
        cpu=sum(variants[label]['golden',i]['compute_cpu'] for i in common)
        methods[label]=dict(golden_cases=len(golden),golden_outcomes=dict(outcomes),
            golden_recovery_percent=100*outcomes['recovered']/len(golden),
            common_complete_cases=len(common),common_cpu_seconds=cpu,common_mean_cpu_seconds=cpu/len(common),
            full_golden_completed_cases=sum(v['searches_complete'] for v in values),
            full_golden_completed_cpu_seconds=sum(v['compute_cpu'] for v in values if v['searches_complete']),
            holdout_cases=len(holdout),holdout_full_mapping_cases=sum(variants[label]['holdout',i]['best_events'] is not None for i in holdout),
            holdout_best_event_comparison=dict(Counter(
                'equal' if variants[label]['holdout',i]['best_events']==variants[baseline]['holdout',i]['best_events'] else
                'worse' if variants[label]['holdout',i]['best_events']>variants[baseline]['holdout',i]['best_events'] else 'better'
                for i in comparable)),
            missing_relative_to_ten=[i for i in golden if variants[baseline]['golden',i]['recovery']=='recovered'
                                     and variants[label]['golden',i]['recovery']!='recovered'])
    cpu=sum(slap[i]['workflow_cpu_excluding_io'] for i in common)
    methods['slap_sweep']=dict(golden_cases=len(golden),golden_recovered=sum(slap[i]['recovered'] for i in golden),
        golden_recovery_percent=100*sum(slap[i]['recovered'] for i in golden)/len(golden),
        common_complete_cases=len(common),common_cpu_seconds=cpu,common_mean_cpu_seconds=cpu/len(common))
    for label in labels:
        methods[label]['cpu_relative_to_slap']=methods[label]['common_cpu_seconds']/cpu
        methods[label]['cpu_reduction_factor_vs_ten']=methods[baseline]['common_cpu_seconds']/methods[label]['common_cpu_seconds']
    return dict(methods=methods,common_case_indices=common,
        scope='Recovery uses all Golden cases, including unknowns. CPU comparisons use the same mutually completed cases, '
              'excluding measured saving/loading. AAM includes process/IPC setup; SLAP workflow timing excludes worker startup. '
              'Not a full-cost or equal-resource cluster latency comparison. Holdout has no annotated ground truth.')


def main(args):
    from aam_seed_ablation import read,save,sha
    variants={};sources=[];manifests=[]
    for count in args.seeds:
        directory=ROOT/f'reports/aam_seed{count}_20260910'
        rows=read(directory/'case_metrics.json')
        manifests.append(read(directory/'manifest.json'))
        for row in rows:
            key=(row['dataset'],row['index'])
            for label,value in row['variants'].items():
                previous=variants.setdefault(label,{})
                assert key not in previous or previous[key]==value, 'Saved baseline changed between ablations'
                previous[key]=value
        sources.append(dict(path=str(directory/'case_metrics.json'),sha256=sha(directory/'case_metrics.json')))
    reference=manifests[0]
    for manifest in manifests[1:]:
        assert manifest['original_commit']==reference['original_commit']
        assert manifest['original_sha256']==reference['original_sha256']
        assert manifest['adapter_sha256']==reference['adapter_sha256']
        assert {k:v for k,v in manifest['original_config'].items() if k!='seed_count'}=={
            k:v for k,v in reference['original_config'].items() if k!='seed_count'}
    slap_path=ROOT/'reports/slap_sweep_cut_20260910/case_metrics.json'
    slap={row['index']:row for row in read(slap_path)}
    result=aggregate(variants,slap)
    result['sources']=sources+[dict(path=str(slap_path),sha256=sha(slap_path))]
    result['aam_engine_commit']=reference['original_commit']
    result['only_search_config_variation']='seed_count'
    save(args.output/'comparison.json',result)
    with (args.output/'comparison.csv').open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['method','golden_cases','golden_recovery_percent','common_complete_cases','common_mean_cpu_seconds'])
        for label,value in result['methods'].items():
            writer.writerow([label,*(value[k] for k in ('golden_cases','golden_recovery_percent','common_complete_cases','common_mean_cpu_seconds'))])
    print(json.dumps(result['methods'],indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seeds',type=int,nargs='+',default=[1,2,3])
    parser.add_argument('--output',type=Path,default=ROOT/'reports/aam_seed1_2_20260910')
    main(parser.parse_args())
