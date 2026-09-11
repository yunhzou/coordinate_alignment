"""Collect or publish the binary-metal experiment, keeping unknowns explicit."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import statistics
import time

from metal_binary_benchmark import read,save,sha,folder,VARIANTS,GOLDEN
from metal_binary_events import binary_metal_input


def collect(run):
    manifest=read(run/'manifest.json');tasks=read(run/'tasks.json');rows=[]
    for slot,task in enumerate(tasks):
        row=dict(slot=slot,**task,variants={});evaluations={}
        for variant in VARIANTS:
            f=folder(run,task,variant);search=read(f/'search.json') if (f/'search.json').exists() else {}
            evaluation=read(f/'raw_evaluation.json') if (f/'raw_evaluation.json').exists() else {}
            evaluations[variant]=evaluation
            timing=search.get('rows',[{}])[-1]
            row['variants'][variant]=dict(search_complete=search.get('complete',False),analysis_complete=bool(evaluation),
                mapping_cpu=timing.get('compute_cpu_excluding_persistence_and_loading_seconds'),
                mapping_wall=timing.get('elapsed_wall_including_io_seconds'),capped=timing.get('capped'),
                states=timing.get('states'),terminals=timing.get('terminals'),
                full_mapping=evaluation.get('full_mapping'),full_representatives=evaluation.get('full_representatives'),
                best_representative_events=evaluation.get('best_representative_events'),
                minimum_patterns=len(evaluation.get('minimum_patterns',{})),
                minimum_pattern_ids=sorted(evaluation.get('minimum_patterns',{})),
                analysis_cpu=evaluation.get('analysis_cpu'),analysis_wall=evaluation.get('analysis_wall'),
                recovery=evaluation.get('reference_recovery','unknown'),audit=evaluation.get('audit'),
                target_score=evaluation.get('target',{}).get('total') if evaluation.get('target') else None,
                target_id=evaluation.get('target',{}).get('id') if evaluation.get('target') else None)
        if task['dataset']=='holdout' and all(evaluations.values()):
            a,b=[set(evaluations[v]['minimum_patterns']) for v in VARIANTS]
            row['representative_pattern_comparison']=dict(shared=len(a&b),original_only=len(a-b),binary_only=len(b-a))
        q=run/f"family_queries/{task['index']}.json"
        if task['dataset']=='holdout' and q.exists():
            query=read(q);row['family_check_complete']=query.get('complete',False)
            for variant in VARIANTS:
                data=row['variants'][variant];found=query['variants'].get(variant,{}).get('results',{})
                scores=[query['patterns'][k]['total'] for k,v in found.items() if v['status']=='represented']
                baseline=data['best_representative_events']
                data['best_verified_events']=min(([baseline] if baseline is not None else [])+scores,default=None)
                data['target_status']=found.get(data['target_id'],{}).get('status')
                data['queried_pattern_statuses']=dict(Counter(v['status'] for v in found.values()))
            comparison=Counter()
            for key in query['patterns']:
                a,b=[query['variants'].get(v,{}).get('results',{}).get(key,{}).get('status','unresolved') for v in VARIANTS]
                comparison['shared' if a==b=='represented' else
                           'binary_gain' if a=='excluded_from_saved_families' and b=='represented' else
                           'binary_loss' if a=='represented' and b=='excluded_from_saved_families' else
                           'neither' if a==b=='excluded_from_saved_families' else 'unresolved']+=1
            row['queried_pattern_comparison']=dict(comparison)
        rows.append(row)
    holdout=[r for r in rows if r['dataset']=='holdout']
    paired=[r for r in holdout if all(v['search_complete'] for v in r['variants'].values())]
    metrics={}
    for variant in VARIANTS:
        values=[r['variants'][variant] for r in holdout]
        metrics[variant]=dict(searches_complete=sum(v['search_complete'] for v in values),
            analyses_complete=sum(v['analysis_complete'] for v in values),
            full_mapping_cases=sum(bool(v['full_mapping']) for v in values),
            paired_mapping_cpu=sum(r['variants'][variant]['mapping_cpu'] for r in paired),
            paired_median_cpu=statistics.median(r['variants'][variant]['mapping_cpu'] for r in paired) if paired else None,
            capped_cases=sum(bool(v['capped']) for v in values),
            saved_full_representatives=sum(v['full_representatives'] or 0 for v in values),
            analysis_cpu=sum(v['analysis_cpu'] or 0 for v in values),
            audit_violations=sum(v['audit']['violations'] for v in values if v['audit']))
    changes=[];counts=Counter();verified_counts=Counter()
    for row in holdout:
        a,b=[row['variants'][v]['best_representative_events'] for v in VARIANTS]
        if a is None or b is None:counts['unavailable']+=1;continue
        counts['lower' if b<a else 'higher' if b>a else 'equal']+=1
        av,bv=[row['variants'][v].get('best_verified_events',row['variants'][v]['best_representative_events']) for v in VARIANTS]
        verified_counts['lower' if bv<av else 'higher' if bv>av else 'equal']+=1
        if a!=b:changes.append(dict(index=row['index'],original=a,metal_binary=b,
                                   original_best_verified=av,metal_binary_best_verified=bv))
    selected=[r for r in rows if r['dataset']=='golden'];by={(r['index'],r['direction']):r for r in selected}
    archived=read(GOLDEN/'comparison/cases.json');archived={r['index']:r['variants']['seeds1']['recovery'] for r in archived if r['dataset']=='golden'}
    golden=[]
    for index in manifest['golden']['cases']:
        result=dict(index=index,archived=archived[index])
        for variant in VARIANTS:
            outcomes=[by[index,d]['variants'][variant]['recovery'] for d in ('R_to_P','P_to_R')]
            result[variant]='recovered' if 'recovered' in outcomes else 'not_recovered' if outcomes==['not_recovered','not_recovered'] else 'unknown'
        golden.append(result)
    summary=dict(holdout=dict(cases=140,paired_completed_cases=len(paired),metrics=metrics,
        cpu_ratio_binary_over_original=metrics['metal_binary']['paired_mapping_cpu']/metrics['original']['paired_mapping_cpu'] if paired else None,
        representative_event_comparison=dict(counts),best_verified_event_comparison=dict(verified_counts),score_changes=changes,
        representative_pattern_changes=[dict(index=r['index'],**r['representative_pattern_comparison']) for r in holdout
            if r.get('representative_pattern_comparison',{}).get('original_only') or r.get('representative_pattern_comparison',{}).get('binary_only')],
        targeted_family_comparison=[dict(index=r['index'],**r['queried_pattern_comparison']) for r in holdout if 'queried_pattern_comparison' in r],
        score_scope='Best saved representative / best verified supplied pattern; these are achievable scores, not global optimality certificates.'),
        golden=dict(changed_cases=golden,coverage_with_unchanged_archives={variant:dict(Counter(
            {**archived,**{r['index']:r[variant] for r in golden}}.values())) for variant in VARIANTS},
            scope='Only 11 input-changing Golden cases rerun; other 1840 exact unchanged search problems reuse archived reference checks. No full-Golden fresh timing claim.'),
        scoring=manifest['scoring'],historical_score_caveat=manifest['historical_score_caveat'])
    save(run/'comparison/per_case.json',rows);save(run/'comparison/summary.json',summary)
    return rows,summary


def publish(args):
    args.run=args.run.resolve()
    rows,summary=collect(args.run);manifest=read(args.run/'manifest.json')
    statuses={}
    for phase in ('search','analyze'):
        data=[read(args.run/f'status/{phase}/{slot}.json') for slot in range(manifest['tasks'])]
        assert all('finished' in r for r in data)
        statuses[phase]=dict(outcomes=dict(Counter(str(v.get('exit',v.get('status'))) for r in data for v in r['variants'].values())),
                            exceptions=[dict(slot=r['slot'],variant=v,**x) for r in data for v,x in r['variants'].items() if x.get('exit')!=0])
    assert all(v['search_complete'] and v['analysis_complete'] and v['audit']['violations']==0 for r in rows for v in r['variants'].values())
    required={11,64,101}|{c['index'] for c in summary['holdout']['representative_pattern_changes']}
    assert all(r.get('family_check_complete') for r in rows if r['dataset']=='holdout' and r['index'] in required)
    for relative,expected in manifest['frozen_sha256'].items():assert sha(args.run/relative)==expected,relative
    for name,expected in manifest['drivers_sha256'].items():assert sha(args.run/name)==expected,name
    repair=read(args.run/'analysis_repair.json')
    assert sha(Path(repair['corrected_driver']))==repair['corrected_driver_sha256']
    assert repair['initial_driver_sha256']==manifest['drivers_sha256']['metal_binary_benchmark.py']
    summary['analysis_repair']=dict(reason=repair['reason'],retried_variants=sum(len(s['slots']) for s in repair['submissions']),
        scope='Postprocessing only; initial frozen drivers and all saved search results remain unchanged. Failed attempt logs are retained.')
    summary['scheduler_accounting']=read(args.run/'followup_accounting.json')
    # Measure the input transform separately; raw data loading is outside this timer.
    conversion_cpu=conversion_wall=0.
    for row in rows:
        if row['dataset']!='holdout':continue
        path=args.run/f"raw/holdout/{row['index']}/input.json";raw=read(path)
        expected=read(args.run/f"runs/holdout/metal_binary/inputs/holdout/{row['index']}/input.json")
        cpu,wall=time.process_time(),time.perf_counter();actual=binary_metal_input(raw)
        conversion_cpu+=time.process_time()-cpu;conversion_wall+=time.perf_counter()-wall
        assert actual==expected
    hashes=read(args.run/'input_audit.json')['hashes']
    for dataset,entries in hashes.items():
        source=Path(manifest['frozen_source'] if dataset=='holdout' else manifest['golden_baseline'])
        for index,entry in entries.items():assert sha(source/f'inputs/{dataset}/{index}/input.json')==entry['original_input_sha256']
    unchanged=[]
    for row in rows:
        dataset,index=row['dataset'],row['index']
        assert sha(args.run/f'raw/{dataset}/{index}/input.json')==hashes[dataset][str(index)]['original_input_sha256']
        if dataset=='holdout' and not hashes[dataset][str(index)]['search_input_changes']:
            a,b=[row['variants'][v] for v in VARIANTS]
            for field in ('states','terminals','capped','full_representatives','best_representative_events','minimum_pattern_ids'):
                assert a[field]==b[field],(index,field)
            unchanged.append(index)
    summary['unchanged_holdout_controls']=dict(cases=len(unchanged),indices=unchanged,
        checks='Identical state/terminal counts, cap flags, full representative counts, lowest representative scores and minimum pattern IDs.')
    summary['phase_statuses']=statuses
    summary['input_conversion']=dict(cpu_seconds=conversion_cpu,wall_seconds=conversion_wall,cases=140,
        scope='One post-run measurement of the same conversion, excluding raw input loading. Not included in mapping CPU.')
    summary['source_and_input_hashes_verified']=True
    args.destination.mkdir(parents=True,exist_ok=False)
    for name in ('manifest.json','input_audit.json','submissions.json','historical_targets.json','analysis_repair.json','followup_accounting.json'):
        shutil.copy2(args.run/name,args.destination/name)
    save(args.destination/'summary.json',summary);save(args.destination/'per_case.json',rows)
    if (args.run/'family_queries').exists():shutil.copytree(args.run/'family_queries',args.destination/'family_queries')
    for path in args.run.glob('*query*submission*.json'):shutil.copy2(path,args.destination/path.name)
    fields=['dataset','index','direction','variant','mapping_cpu','mapping_wall','capped','best_representative_events',
            'best_verified_events','full_representatives','target_score','target_status','recovery']
    with (args.destination/'per_case.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore',lineterminator='\n');writer.writeheader()
        for r in rows:
            for variant,value in r['variants'].items():writer.writerow(dict(dataset=r['dataset'],index=r['index'],direction=r['direction'],variant=variant,**value))
    a,b=[summary['holdout']['metrics'][v] for v in VARIANTS]
    lines=['# Binary metal connections for search; raw-WBO events downstream','',
        'Binary metal search recovers the previously missing PR7 four-event alternative with one seed. '
        'The lowest saved event count is unchanged in all 140 holdout cases, and Golden reference recovery is unchanged. '
        'The observed mapping CPU reduction is small (2.6%). This supports a targeted benefit for metal-sensitive matching, '
        'without establishing general chemical accuracy or a consistent speed advantage.','',
        'The frozen engine `98b01b1` searches original weights or weights in which every metal-containing pair is 1 '
        'when raw WBO >= 0.2 and 0 otherwise. Nonmetal weights remain unchanged. One seed, tolerance 1.0, explicit H, '
        'and the original uncut/single-edge sweep are retained. The holdout uses R→P and cap 1000; Golden uses both directions and cap 100.','',
        '**Both arms use the original raw WBOs for event scoring:** decrease >= 0.5 means breaking/weakening and increase >= 0.5 means forming/strengthening; '
        'pairs involving a metal use 0.3. This is the existing `classify_bonds` rule, with inclusive thresholds and no extra 0.2 event-presence gate. '
        'The older holdout benchmark helper used a different uniform-0.5/floor-0.2 rule, so its event totals cannot be directly compared to these.','',
        '| 140-case forward holdout | Original WBO search | Binary metal search |','|---|---:|---:|',
        f"| Full saved mappings | {a['full_mapping_cases']}/140 | {b['full_mapping_cases']}/140 |",
        f"| Mapping CPU, seconds | {a['paired_mapping_cpu']:.3f} | {b['paired_mapping_cpu']:.3f} |",
        f"| Median mapping CPU, seconds/case | {a['paired_median_cpu']:.3f} | {b['paired_median_cpu']:.3f} |",
        f"| Cases with cap flags | {a['capped_cases']} | {b['capped_cases']} |",
        f"| Unique full representatives, summed | {a['saved_full_representatives']} | {b['saved_full_representatives']} |",'',
        f"CPU ratio binary/original: {summary['holdout']['cpu_ratio_binary_over_original']:.4f}. Fresh pairs use the same CPU allocation and alternate order. "
        'CPU includes parent plus workers and excludes measured archive persistence/loading. Conversion and analysis are separate. '
        'Single paired measurements do not establish a statistical speed guarantee.','',
        f"Best representative event comparison (binary versus original): {summary['holdout']['representative_event_comparison']}. "
        f"After targeted saved-family checks, best verified scores: {summary['holdout']['best_verified_event_comparison']}. "
        'These are achievable scores, not proofs of the global minimum. The holdout has no annotated chemical ground truth.']
    if summary['holdout']['score_changes']:
        lines+=['','## Changed event scores','', '| Index | Original representative | Binary representative | Original best verified | Binary best verified |', '|---|---:|---:|---:|---:|']
    for c in summary['holdout']['score_changes']:
        lines.append(f"| {c['index']} | {c['original']} | {c['metal_binary']} | {c['original_best_verified']} | {c['metal_binary_best_verified']} |")
    lines+=['','## Previously missing alternatives','', '| Index | Target events at 0.5/0.3 | Original membership | Binary membership |', '|---|---:|---|---|']
    for row in rows:
        if row['dataset']=='holdout' and row['index'] in (11,64,101):
            a,b=[row['variants'][v] for v in VARIANTS]
            lines.append(f"| {row['index']} | {a['target_score']} | {a.get('target_status','not_checked')} | {b.get('target_status','not_checked')} |")
    pr7=read(args.run/'family_queries/101.json')
    target=next(r['variants']['original']['target_id'] for r in rows if r['dataset']=='holdout' and r['index']==101)
    witness=pr7['variants']['metal_binary']['results'][target]['witness']
    raw=read(args.run/'raw/holdout/101/input.json');mapping=witness['mapping']
    lines+=['','PR7 gained witness: two breaking/weakening and two forming/strengthening events. '
        'The table uses zero-based input atom indices; all weights are the original WBOs.','',
        '| Event | R pair | P pair | Raw R WBO | Raw P WBO | ΔWBO |','|---|---|---|---:|---:|---:|']
    for kind,edges in witness['events'].items():
        for x,y in edges:
            rw=raw['reactant']['wbo'][x][y];pw=raw['product']['wbo'][mapping[x]][mapping[y]]
            lines.append(f'| {kind} | {x}–{y} | {mapping[x]}–{mapping[y]} | {rw:.6f} | {pw:.6f} | {pw-rw:+.6f} |')
    lines+=['','Targets are the historical missing witness mappings, reclassified under the requested delta thresholds. '
        'Membership includes compressed families, with search constraints checked on search weights and events evaluated on raw WBO. '
        'Unresolved checks are not counted as exclusions. Pattern identity uses signed event pairs modulo exact raw-WBO score-response symmetry. '
        'This can distinguish patterns merged by an unweighted connectivity convention.','',
        f"Cases with differing minimum representative pattern sets: {summary['holdout']['representative_pattern_changes']}", '',
        f"Targeted family comparison: {summary['holdout']['targeted_family_comparison']}", '',
        'The query set includes both arms’ minimum representative patterns wherever their sets differ, plus the three historical targets. '
        'Identical representative sets do not prove equality of the entire compressed-family event spectrum; unseen patterns are not enumerated.','', '## Golden','',
        'Normalization changes only 11 Golden inputs. Those are rerun in both directions; the other 1840 search inputs are exactly unchanged '
        'and retain their archived reference checks. This is not a fresh timing run of all 1851 cases.','',
        f"Reference recovery including unchanged archives: {summary['golden']['coverage_with_unchanged_archives']}", '',
        '| Index | Archived original | Fresh original | Binary |','|---|---|---|---|']
    for c in summary['golden']['changed_cases']:lines.append(f"| {c['index']} | {c['archived']} | {c['original']} | {c['metal_binary']} |")
    lines+=['','## Validation and reproduction','',
        'Small exhaustive tests cover pair-specific inclusive thresholds, input immutability, symmetry-equivalent events, '
        'and a compressed family whose raw-WBO event scores differ despite identical binary metal weights. '
        'Every emitted event witness is independently rescored with `classify_bonds`. All saved states are checked for element preservation '
        'and injectivity; preserved transition bonds and stored-generator target adjacency are checked. '
        'Frozen source and original input hashes are verified. Engine defaults and earlier reports remain unchanged.','',
        f"All {len(unchanged)} holdout inputs unaffected by normalization reproduce the same checked search counts and minimum representative patterns. "
        f"The initial auxiliary audit failed on {summary['analysis_repair']['retried_variants']} analyses because discarded branches legitimately lack finalized symmetry groups. "
        'Only those analyses were retried after correcting the ancestry check; every returned path still requires finalized groups. '
        'No mapping searches were rerun. Original failure statuses, repair hashes and scheduler accounting are retained.','',
        f"Input conversion: {conversion_cpu:.3f} CPU seconds for 140 cases (separate post-run measurement). "
        f"Successful holdout event analysis CPU: original {summary['holdout']['metrics']['original']['analysis_cpu']:.3f} s; "
        f"binary {summary['holdout']['metrics']['metal_binary']['analysis_cpu']:.3f} s. "
        'These analysis totals exclude archive loading, auxiliary structural audits, failed attempts and optional family queries.','',
        f'Raw outputs and logs: `{args.run}`. Drivers: `bench/metal_binary_benchmark.py`, `bench/metal_binary_events.py`, '
        '`bench/publish_metal_binary.py`. The run manifest records source hashes, configurations and scoring semantics.','']
    (args.destination/'README.md').write_text('\n'.join(lines))
    save(args.destination/'artifacts.json',dict(raw_run=str(args.run),publisher_sha256=sha(Path(__file__))))
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('collect','publish'))
    p.add_argument('--run',type=Path,required=True);p.add_argument('--destination',type=Path)
    args=p.parse_args()
    if args.command=='collect':print(json.dumps(collect(args.run)[1],indent=2))
    else:publish(args)
