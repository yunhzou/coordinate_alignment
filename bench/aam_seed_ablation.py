"""Seed-count-only AAM ablation against a saved frozen native-reuse baseline.

The mapper and native binary are copied unchanged. Search sees input structures
only; evaluation is a separate, watchdog-limited phase over saved archives.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
BASE=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks')
BASELINE=BASE/'adaptive_full_20260910'
PYTHON='/h/399/yunhengzou/coordinate_alignment/.venv/bin/python'


def read(path):return json.loads(path.read_text())


def save(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n')
    temporary.replace(path)


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def result_folder(run,task):
    return run/f"results/{task['dataset']}/{task['index']}/{task['direction']}/original"


def prepare(args):
    source=args.baseline.resolve()
    old=read(source/'manifest.json')
    assert old['original_config']['seed_count']==10
    assert old['original_commit'].startswith('98b01b1')
    args.run.mkdir(parents=True,exist_ok=False)
    shutil.copytree(source/'original',args.run/'original',ignore=shutil.ignore_patterns('__pycache__'))
    # Keep the baseline's exact search, profiling and evaluation adapters too.
    shutil.copytree(source/'engine/bench',args.run/'engine/bench',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(__file__,args.run/'aam_seed_ablation.py')
    (args.run/'inputs').symlink_to(source/'inputs',target_is_directory=True)
    for name in ('tasks.json','input_hashes.json'):shutil.copy2(source/name,args.run/name)
    tasks=read(args.run/'tasks.json')
    weights=[]
    for slot,task in enumerate(tasks):
        raw=read(args.run/f"inputs/{task['dataset']}/{task['index']}/input.json")
        n=max(len(raw[side]['elements']) for side in ('reactant','product'))
        weights.append((n,slot))
    save(args.run/'work_order.json',[slot for _,slot in sorted(weights,reverse=True)])
    for phase in ('search','analyze'):
        (args.run/f'status/{phase}').mkdir(parents=True)
        save(args.run/f'{phase}_queue.json',dict(next=0))
    config=dict(old['original_config'],seed_count=args.seeds)
    manifest=dict(schema='aam_seed_ablation/v1',baseline=str(source),counts=old['counts'],
        original_commit=old['original_commit'],original_config=config,original_workers=old['original_workers'],
        original_execution=old['original_execution'],tasks=len(tasks),root_seed=old['root_seed'],
        python_hash_seed=old['python_hash_seed'],explicit_H=old['explicit_H'],modes=old['modes'],
        baseline_config=old['original_config'],changed_config_fields=['seed_count'],
        seed_policy='First N of the same deterministic per-cut seed orders; no new seed selection',
        search_watchdog=300,analysis_watchdog=240,cpu_budget=args.cpu_budget,node_cpus=args.node_cpus,
        analysis_node_cpus=args.analysis_cpus,
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        original_sha256={str(p.relative_to(args.run/'original')):sha(p) for p in (args.run/'original').rglob('*') if p.is_file()},
        adapter_sha256={str(p.relative_to(args.run/'engine/bench')):sha(p) for p in (args.run/'engine/bench').rglob('*.py')},
        driver_sha256=sha(Path(__file__)),
        timing='Unchanged SearchProfiler: parent plus child CPU minus measured checkpoint writing/loading. '
               'Search, evaluation, watchdogs, scheduling and persistence are separate.',
        scheduling='Global locked task queue; allocations claim fresh cases dynamically, longest structures first.')
    save(args.run/'manifest.json',manifest)
    print(json.dumps(dict(run=str(args.run),tasks=len(tasks),seeds=args.seeds,config_changes=['seed_count'])),flush=True)


def claim(run,phase):
    # Stable lock inode; only the small queue counter is atomically replaced.
    with (run/f'{phase}_queue.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        counter=read(run/f'{phase}_queue.json')
        order=read(run/'work_order.json')
        if counter['next']==len(order):return None
        slot=order[counter['next']]
        save(run/f'{phase}_queue.json',dict(next=counter['next']+1))
        return slot


def batch(args):
    manifest=read(args.run/'manifest.json');tasks=read(args.run/'tasks.json')
    cores=sorted(os.sched_getaffinity(0))
    width=manifest['original_workers'] if args.phase=='search' else 1
    allocation=manifest['node_cpus'] if args.phase=='search' else manifest['analysis_node_cpus']
    assert len(cores)>=allocation and allocation%width==0
    env=dict(os.environ,PYTHONPATH=f'{args.run}/original/src:{args.run}/engine/bench',RXN_CORE_NATIVE='1',
        PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def work(affinity):
        while (slot:=claim(args.run,args.phase)) is not None:
            task=tasks[slot];folder=result_folder(args.run,task)
            status=args.run/f'status/{args.phase}/{slot}.json'
            assert not status.exists(), 'A claimed task must not be rerun'
            if args.phase=='analyze' and not (folder/'search.json').exists():
                save(status,dict(slot=slot,**task,status='search_incomplete',finished=time.time()))
                continue
            row=dict(slot=slot,**task,phase=args.phase,started=time.time(),host=socket.gethostname(),
                     affinity=affinity,slurm_job=os.environ.get('SLURM_JOB_ID'))
            save(status,row)
            command=['taskset','-c',','.join(map(str,affinity)),'timeout','--kill-after=5s',
                str(manifest['search_watchdog'] if args.phase=='search' else manifest['analysis_watchdog']),
                PYTHON,str(args.run/'engine/bench/adaptive_full_benchmark.py'),args.phase,
                '--run',str(args.run),'--slot',str(slot),'--method','original']
            start=time.perf_counter()
            with status.with_suffix('.log').open('w') as log:
                code=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
            save(status,dict(**row,exit=code,finished=time.time(),elapsed_including_startup_io=time.perf_counter()-start))
    groups=[cores[i:i+width] for i in range(0,allocation,width)]
    with ThreadPoolExecutor(max_workers=len(groups)) as pool:list(pool.map(work,groups))


def submit(args):
    manifest=read(args.run/'manifest.json')
    assert not (args.run/'submissions.json').exists(), 'Do not resubmit a started campaign'
    manifest['analysis_node_cpus']=args.analysis_cpus
    # Finalize the orchestration snapshot before any worker starts. The frozen
    # AAM engine, configuration and baseline adapters are never modified.
    shutil.copy2(__file__,args.run/'aam_seed_ablation.py')
    manifest['driver_sha256']=sha(Path(__file__))
    save(args.run/'manifest.json',manifest)
    allocations=manifest['cpu_budget']//manifest['node_cpus']
    jobs=[]
    for phase in ('search','analyze'):
        cpus=manifest['node_cpus'] if phase=='search' else manifest['analysis_node_cpus']
        command=['sbatch','--parsable','--partition=cpunodes_nia','--exclude=bosque49,bosque56',
            '--nodes=1',f'--cpus-per-task={cpus}','--mem=64G','--time=00:20:00','--no-requeue',
            f'--array=0-{allocations-1}%{allocations}',f'--job-name=aam_seed{manifest["original_config"]["seed_count"]}_{phase}',
            f'--output={args.run}/status/{phase}_%A_%a.out']
        if jobs:command.append('--dependency=afterany:'+jobs[0]['job'].split(';')[0])
        command+=['--wrap',shlex.join(['env','PYTHONDONTWRITEBYTECODE=1',PYTHON,str(args.run/'aam_seed_ablation.py'),
                   'batch','--run',str(args.run),'--phase',phase])]
        job=subprocess.check_output(command,text=True).strip()
        jobs.append(dict(phase=phase,job=job,command=command))
        save(args.run/'submissions.json',jobs)
        print(phase,job,flush=True)


def collect(args):
    """Summarize saved results only; unknown evaluations are never negatives."""
    manifest=read(args.run/'manifest.json')
    trial=f"seeds{manifest['original_config']['seed_count']}"
    baseline=f"seeds{manifest['baseline_config']['seed_count']}"
    rows=[]
    for dataset,count in manifest['counts'].items():
        for index in range(count):
            variants={}
            for name,root in ((trial,args.run),(baseline,Path(manifest['baseline']))):
                directions={}
                for direction in ('R_to_P','P_to_R'):
                    folder=result_folder(root,dict(dataset=dataset,index=index,direction=direction))
                    path=folder/'search.json'
                    search=read(path) if path.exists() else {}
                    path=folder/'full_sweep_evaluation.json'
                    evaluation=read(path) if path.exists() else {}
                    directions[direction]=dict(search_complete=search.get('complete',False),
                        compute_cpu=search['rows'][-1]['compute_cpu_excluding_persistence_and_loading_seconds']
                                    if search.get('complete') else None,
                        reference_recovery=evaluation.get('reference_recovery','unknown'),
                        best_events=evaluation.get('best_events'),valid_full_representatives=evaluation.get('valid_full_representatives'),
                        evaluation_complete=bool(evaluation))
                outcomes=[d['reference_recovery'] for d in directions.values()]
                variants[name]=dict(directions=directions,
                    recovery='recovered' if 'recovered' in outcomes else
                             'not_recovered' if outcomes==['not_recovered','not_recovered'] else 'unknown',
                    searches_complete=all(d['search_complete'] for d in directions.values()),
                    compute_cpu=sum(d['compute_cpu'] for d in directions.values())
                                if all(d['compute_cpu'] is not None for d in directions.values()) else None,
                    best_events=min((d['best_events'] for d in directions.values() if d['best_events'] is not None),default=None))
            rows.append(dict(dataset=dataset,index=index,variants=variants))
    summaries={}
    for dataset,count in manifest['counts'].items():
        selected=[r for r in rows if r['dataset']==dataset]
        paired=[r for r in selected if all(v['searches_complete'] for v in r['variants'].values())]
        totals={name:sum(r['variants'][name]['compute_cpu'] for r in paired) for name in (trial,baseline)}
        summary=dict(cases=count,paired_complete_cases=len(paired),paired_compute_cpu=totals,
                     paired_cpu_speedup=totals[baseline]/totals[trial] if totals[trial] else None)
        summary['evaluation_complete_cases']={name:sum(all(d['evaluation_complete'] for d in r['variants'][name]['directions'].values())
            for r in selected) for name in totals}
        if dataset=='golden':
            summary['recovery']={name:dict(Counter(r['variants'][name]['recovery'] for r in selected)) for name in totals}
            summary['previously_recovered_now_not_verified']=[r['index'] for r in selected
                if r['variants'][baseline]['recovery']=='recovered' and r['variants'][trial]['recovery']!='recovered']
            summary['newly_verified_from_previous_unknown_or_miss']=[r['index'] for r in selected
                if r['variants'][baseline]['recovery']!='recovered' and r['variants'][trial]['recovery']=='recovered']
        else:
            summary['scope']='No annotated holdout ground truth. Event minima describe saved representatives only.'
            summary['full_mapping_cases']={name:sum(r['variants'][name]['best_events'] is not None for r in selected) for name in totals}
            summary['same_best_events']=sum(r['variants'][trial]['best_events']==r['variants'][baseline]['best_events']
                for r in selected if all(v['best_events'] is not None for v in r['variants'].values()))
        summaries[dataset]=summary
    save(args.run/'comparison/cases.json',rows)
    save(args.run/'comparison/summary.json',summaries)
    print(json.dumps(summaries,indent=2),flush=True)


def publish(args):
    """Publish compact comparisons without repeating any mapping or evaluation."""
    import csv
    manifest=read(args.run/'manifest.json')
    trial=f"seeds{manifest['original_config']['seed_count']}"
    baseline=f"seeds{manifest['baseline_config']['seed_count']}"
    statuses={phase:[read(p) for p in (args.run/f'status/{phase}').glob('*.json')]
              for phase in ('search','analyze')}
    assert all(len(rows)==manifest['tasks'] and all('finished' in r for r in rows) for rows in statuses.values())
    rows=read(args.run/'comparison/cases.json')
    destination=ROOT/f"reports/aam_seed{manifest['original_config']['seed_count']}_20260910"
    destination.mkdir(exist_ok=True)
    for name in ('manifest.json','submissions.json'):
        shutil.copy2(args.run/name,destination/name)
    shutil.copy2(args.run/'comparison/summary.json',destination/'summary.json')
    shutil.copy2(args.run/'comparison/cases.json',destination/'case_metrics.json')
    with (destination/'per_case.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['dataset','index','variant','recovery','searches_complete','compute_cpu_seconds','best_events'])
        writer.writeheader()
        for row in rows:
            for variant,value in row['variants'].items():
                writer.writerow(dict(dataset=row['dataset'],index=row['index'],variant=variant,
                    recovery=value['recovery'] if row['dataset']=='golden' else 'no_ground_truth',
                    searches_complete=value['searches_complete'],compute_cpu_seconds=value['compute_cpu'],best_events=value['best_events']))
    timing=dict(phase_statuses={phase:dict(Counter(str(r.get('exit',r.get('status'))) for r in values))
                               for phase,values in statuses.items()},
        phase_execution_spans_seconds={phase:max(r['finished'] for r in values)-min(r['started'] for r in values if 'started' in r)
                                       if any('started' in r for r in values) else None
                                       for phase,values in statuses.items()},
        full_search_cpu_by_dataset={dataset:sum(r['variants'][trial]['compute_cpu'] for r in rows if r['dataset']==dataset)
                                    for dataset in manifest['counts']},
        scope='CPU excludes measured persistence/loading. Execution spans include persistence and dispatch, not initial Slurm queue. '
              'Bidirectional calls were independently scheduled; no equal-resource concurrent bidirectional latency claim.')
    save(destination/'timing.json',timing)
    slap_path=ROOT/'reports/slap_sweep_cut_20260910/case_metrics.json'
    slap={r['index']:r for r in read(slap_path)}
    common=[r for r in rows if r['dataset']=='golden' and all(v['searches_complete'] for v in r['variants'].values())
            and not slap[r['index']]['incomplete_variants'] and not slap[r['index']]['mapping_errors']
            and not slap[r['index']]['invalid_predictions']]
    sums={name:sum(r['variants'][name]['compute_cpu'] for r in common) for name in (trial,baseline)}
    sums['slap_sweep']=sum(slap[r['index']]['workflow_cpu_excluding_io'] for r in common)
    save(destination/'slap_cpu_comparison.json',dict(paired_complete_cases=len(common),cpu_seconds=sums,
        mean_cpu_seconds={k:v/len(common) for k,v in sums.items()},
        trial_cpu_over_slap=sums[trial]/sums['slap_sweep'],
        scope='Matched completed cases only. AAM includes IPC/process setup; SLAP workflow timer excludes worker startup. '
              'Not an exact end-to-end latency or full-cost comparison; accuracy keeps all 1851 cases.'))
    jobs=[r['job'].split(';')[0] for r in read(args.run/'submissions.json')]
    fields=['JobID','State','ElapsedRaw','TotalCPU','AllocCPUS','Submit','Start','End','MaxRSS']
    command=['sacct','-j',','.join(jobs),'--array','--noheader','--parsable2','--format='+','.join(fields)]
    accounting=subprocess.check_output(command,text=True)
    save(destination/'slurm_accounting.json',dict(command=command,
        rows=[dict(zip(fields,line.split('|'))) for line in accounting.splitlines()]))
    save(destination/'artifacts.json',dict(full_run=str(args.run),
        native_engine=str(args.run/'original'),full_archives=str(args.run/'results'),
        frozen_driver=str(args.run/'aam_seed_ablation.py'),baseline=manifest['baseline']))
    print(json.dumps(dict(report=str(destination),timing=timing,slap_comparison=sums),indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','batch','submit','collect','publish'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--baseline',type=Path,default=BASELINE)
    parser.add_argument('--seeds',type=int,default=3)
    parser.add_argument('--cpu-budget',type=int,default=1024)
    parser.add_argument('--node-cpus',type=int,default=32)
    parser.add_argument('--analysis-cpus',type=int,default=8)
    parser.add_argument('--phase',choices=('search','analyze'))
    args=parser.parse_args();globals()[args.command](args)
