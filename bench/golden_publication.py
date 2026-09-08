"""Frozen full-Golden benchmark: smaller-first and bidirectional modes.

All subprocesses and timeouts are explicit. Computation, persistence, ranking,
and reference verification have separate artifacts and timing records.
"""
import argparse
from collections import Counter
from dataclasses import asdict, replace
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import random
import resource
import shlex
import shutil
import subprocess
import sys
import time

import numpy as np
from golden_policy_campaign import load_case, save, guarded
from golden_evaluation import evaluate_planned
from publication_timing import SearchProfiler
from publication_analysis import rank_archive, merge_classes, representative_metrics, certificate_id, union_outcome, top_five_family
from rxn_core import AAMProblem, AAMSearchConfig, search_aam
from rxn_core.search_orientation import AAMSearchPlan
from rxn_core.artifacts import read_aam_checkpoint

DIRECTIONS = ('R_to_P','P_to_R')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def metadata():
    return dict(host=platform.node(),platform=platform.platform(),python=sys.version,
        cpu_models=sorted({line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                           if line.startswith('model name')}),
        affinity=sorted(os.sched_getaffinity(0)),
        slurm={k:os.environ.get(k) for k in ('SLURM_JOB_ID','SLURM_ARRAY_JOB_ID','SLURM_ARRAY_TASK_ID','SLURM_CPUS_PER_TASK','SLURM_JOB_PARTITION')})


def plans(run,index):
    _, base = load_case(run/'inputs',index)
    p, c = base.input_problem, base.config
    return {d:AAMSearchPlan(p,AAMProblem(p.product,p.reactant,p.name) if d=='P_to_R' else p,
                c,d=='P_to_R') for d in DIRECTIONS}, base.direction


def initialize(args):
    start=time.perf_counter()
    args.run.mkdir(parents=True,exist_ok=False)
    for name in ('inputs','directions','results','status'):
        (args.run/name).mkdir()
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    original=json.loads((args.source/'manifest.json').read_text())
    indices=sorted(r['index'] for r in original['records'])
    assert indices==list(range(1851)), 'Publication denominator must contain every Golden record'
    config=AAMSearchConfig(seed_count=10,branch_limit=100,iso_tolerance=1.0)
    inputs=[]
    for index in indices:
        out=args.run/'inputs'/str(index);out.mkdir()
        for name in ('input.json','reference.json'):
            shutil.copy2(args.source/str(index)/name,out/name)
        inputs.append(dict(index=index,input_sha256=sha256(out/'input.json'),
                           reference_sha256=sha256(out/'reference.json')))
    save(args.run/'inputs/manifest.json',dict(config=asdict(config),records=inputs))
    tasks=[dict(index=i,direction=d) for i in indices for d in DIRECTIONS]
    dataset=Path('data/aam_benchmarks/golden_original_20260906')
    provenance=json.loads((dataset/'manifest.json').read_text())
    assert sha256(dataset/'golden_dataset.rdf')==provenance['rdf_sha256']
    manifest=dict(schema='golden_publication/v1',indices=indices,tasks=tasks,config=asdict(config),
        workers=16,search_watchdog=300,rank_watchdog=240,verification_watchdog=260,
        union_watchdog=300,root_seed=42,python_hash_seed=0,z3_default_random_seed=0,
        seed_policy='Existing cut_seed: no-cut 42; canonical cut edges plus root 42 hashed with BLAKE2b AAM-cut-seeds-v1. '
                    'Ten existing random seed orders per cut; actual orders retained in graph contexts.',
        modes=dict(single='smaller explicit endpoint first; R_to_P on ties',
                   bidirectional='unconditional union of both directed compressed archives; canonical representative classes merged in original R/P space'),
        source_inputs=str(args.source.resolve()),dataset=provenance,
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        dependencies={d.metadata['Name']:d.version for d in importlib.metadata.distributions()},
        engine_sha256={str(p.relative_to(args.run/'engine')):sha256(p)
                       for p in (args.run/'engine').rglob('*') if p.is_file()},
        preparation_seconds=time.perf_counter()-start,
        scoring='All 1851 in denominator. Reference heavy-atom relation including unmatched atoms, modulo endpoint chemical symmetry. '
                'Explicit H in search and event ranking, not reference H-identity accuracy. No energy or geometric chirality post-filter. '
                'Top-k and event-window representative metrics separate from compressed-family recovery.',
        timing='Primary additive compute CPU excludes checkpoint encoding/compression/writes and reads/decoding. '
               'Raw elapsed includes IO. Ranking and reference evaluation separately timed. No summed-IO subtraction from parallel elapsed. '
               'Offline input preparation, archive hashing and scheduling excluded from matching compute. '
               'Mode compute CPU is single directional CPU plus ranking, or both directional CPU/ranking plus union indexing. '
               'No claim of directly measured equal-resource bidirectional latency from independently scheduled tasks.')
    save(args.run/'manifest.json',manifest)
    save(args.run/'environment.json',metadata())
    print(json.dumps(dict(run=str(args.run),records=len(indices),directional_searches=len(tasks))),flush=True)


def search(args):
    pair,_=plans(args.run,args.index);plan=pair[args.direction]
    out=args.run/'directions'/str(args.index)/args.direction
    out.mkdir(parents=True,exist_ok=False)
    save(out/'environment.json',metadata())
    with SearchProfiler(out/'timing_events') as profiler:
        result=search_aam(plan.problem,plan.config,workers=16,
            intermediate_dir=out/'cuts',archive_format='checkpoint')
    report=profiler.summary()
    report.update(index=args.index,direction=args.direction,metrics=asdict(result.metrics),
        capped=result.graph.capped,terminals=len(result.graph.terminals),
        peak_parent_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        seed_orders_sha256=hashlib.sha256(json.dumps([(c.cuts,c.seed_order) for c in result.graph.contexts],
                                        separators=(',',':')).encode()).hexdigest())
    save(out/'search.json',report)
    start=time.perf_counter()
    archive=out/'cuts/aam.pkl.gz'
    save(out/'artifact.json',dict(path=str(archive.resolve()),sha256=sha256(archive),
        bytes=archive.stat().st_size,integrity_wall_seconds=time.perf_counter()-start))


def rank(args):
    pair,_=plans(args.run,args.index);plan=pair[args.direction]
    out=args.run/'directions'/str(args.index)/args.direction
    start=time.perf_counter();aam=read_aam_checkpoint(out/'cuts/aam.pkl.gz')
    loading=time.perf_counter()-start
    wall,cpu=time.perf_counter(),time.process_time()
    classes=rank_archive(aam,plan)
    timing=dict(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
                loading_wall_seconds=loading,classes=len(classes))
    start=time.perf_counter();save(out/'classes.json',classes)
    timing['saving_wall_seconds']=time.perf_counter()-start
    save(out/'ranking.json',timing)


def verify(args):
    pair,_=plans(args.run,args.index);plan=pair[args.direction]
    out=args.run/'directions'/str(args.index)/args.direction
    reference=json.loads((args.run/'inputs'/str(args.index)/'reference.json').read_text())
    start=time.perf_counter();aam=read_aam_checkpoint(out/'cuts/aam.pkl.gz')
    loading=time.perf_counter()-start
    wall,cpu=time.perf_counter(),time.process_time()
    result=evaluate_planned(aam,plan,reference['features'],reference['mapping'],seconds=220,query_timeout_ms=5000)
    result.update(verification_wall_seconds=time.perf_counter()-wall,
        verification_cpu_seconds=time.process_time()-cpu,loading_wall_seconds=loading,
        archive=str((out/'cuts/aam.pkl.gz').resolve()),search_incomplete=False)
    save(out/'evaluation.json',result)


def combine(args):
    pair,default=plans(args.run,args.index)
    out=args.run/'results'/str(args.index);out.mkdir(exist_ok=True)
    reference=json.loads((args.run/'inputs'/str(args.index)/'reference.json').read_text())
    expected=certificate_id(reference['features'],reference['mapping'])
    records={}
    for direction in DIRECTIONS:
        base=args.run/'directions'/str(args.index)/direction
        def read(name):
            path=base/name
            return json.loads(path.read_text()) if path.exists() else {}
        records[direction]=dict(classes=read('classes.json'),search=read('search.json'),
            ranking=read('ranking.json'),evaluation=read('evaluation.json'))
    modes={}
    for name,directions in (('single',[default]),('bidirectional',list(DIRECTIONS))):
        full=all(isinstance(records[d]['classes'],list) for d in directions)
        outcomes=[records[d]['evaluation'].get('reference_recovery','unknown') for d in directions]
        mode=dict(directions=directions,reference_recovery=union_outcome(outcomes),ranking_complete=full,
                  search_complete=all(bool(records[d]['search']) for d in directions))
        if full:
            wall,cpu=time.perf_counter(),time.process_time()
            classes=merge_classes([records[d]['classes'] for d in directions])
            mode['class_merge_wall_seconds']=time.perf_counter()-wall
            mode['class_merge_cpu_seconds']=time.process_time()-cpu
            mode['representative']=representative_metrics(classes,expected)
            mode['compute_cpu_excluding_io_seconds']=sum(
                records[d]['search']['compute_cpu_excluding_persistence_and_loading_seconds']+
                records[d]['ranking']['cpu_seconds'] for d in directions)+mode['class_merge_cpu_seconds']
            mode['directional_elapsed_including_io_sum_seconds']=sum(
                records[d]['search']['elapsed_wall_including_io_seconds']+
                records[d]['ranking']['wall_seconds'] for d in directions)+mode['class_merge_wall_seconds']
            save(out/f'{name}_classes.json',classes)
            mode['top5_family']=dict(outcome='unknown',reason='verification pending')
        modes[name]=mode
    result=dict(index=args.index,default_direction=default,modes=modes)
    save(out/'result.json',result)
    # Save cheap, complete ranking metrics before optional bounded family checks.
    start=time.monotonic()
    for name,mode in modes.items():
        if not mode['ranking_complete']:
            continue
        if mode['reference_recovery']=='not_recovered':
            mode['top5_family']=dict(outcome='not_recovered',reason='whole returned family union verified absent')
        elif time.monotonic()-start < 220:
            classes=json.loads((out/f'{name}_classes.json').read_text())
            mode['top5_family']=top_five_family(args.run,args.index,classes,pair,reference,
                                               seconds=220-(time.monotonic()-start))
        save(out/'result.json',result)


def worker(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    if args.phase=='combine':
        index=args.offset+args.slot;direction=None
    else:
        task=manifest['tasks'][args.offset+args.slot]
        index,direction=task['index'],task['direction']
    status=args.run/'status'/f'{index}_{direction or "union"}_{args.phase}.json'
    prefix=[sys.executable,str(args.run/'engine/bench/golden_publication.py')]
    common=['--run',str(args.run),'--index',str(index)]
    if direction:common+=['--direction',direction]
    record=dict(index=index,direction=direction,phase=args.phase,started=time.time(),environment=metadata())
    save(status,record)
    if args.phase=='search':
        record['exit']=guarded([*prefix,'search',*common],status.with_suffix('.log'),300)
    elif args.phase=='analyze':
        base=args.run/'directions'/str(index)/direction
        if (base/'cuts/aam.pkl.gz').exists():
            record['rank_exit']=guarded([*prefix,'rank',*common],status.with_suffix('.rank.log'),240)
            record['verify_exit']=guarded([*prefix,'verify',*common],status.with_suffix('.verify.log'),260)
        else:
            record['verify_exit']='incomplete_search'
        if not (base/'evaluation.json').exists():
            base.mkdir(parents=True,exist_ok=True)
            save(base/'evaluation.json',dict(reference_recovery='unknown',
                reason=record.get('verify_exit'),search_incomplete=not (base/'cuts/aam.pkl.gz').exists()))
    else:
        record['exit']=guarded([*prefix,'combine',*common],status.with_suffix('.log'),300)
    record.update(finished=time.time(),complete=True)
    save(status,record)


def submit(args):
    root=args.run.resolve();engine=root/'engine'
    manifest=json.loads((root/'manifest.json').read_text())
    jobs=[]
    def launch(phase,offset,count,dependency=None):
        cpus=16 if phase=='search' else 1
        memory='32G' if phase=='search' else '24G'
        concurrency=40 if phase=='search' else 24
        command=shlex.join(['env',f'PYTHONPATH={engine}/src:{engine}/bench','PYTHONHASHSEED=0',
            'RXN_CORE_NATIVE=1','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1','MKL_NUM_THREADS=1',
            sys.executable,str(engine/'bench/golden_publication.py'),'worker','--run',str(root),
            '--phase',phase,'--offset',str(offset),'--slot'])+' "$SLURM_ARRAY_TASK_ID"'
        options=['sbatch','--parsable',f'--partition={args.partition}',f'--cpus-per-task={cpus}',
            f'--mem={memory}','--time=00:10:00',f'--job-name=gold_paper_{phase}',
            f'--array=0-{count-1}%{concurrency}',f'--output={root}/status/slurm_%A_%a.log']
        if dependency:options.append(f'--dependency=afterany:{dependency}')
        job=subprocess.check_output([*options,'--wrap',command],text=True).strip()
        jobs.append(dict(job=job,phase=phase,offset=offset,count=count,dependency=dependency,
                         partition=args.partition,cpus=cpus,memory=memory))
        save(root/'jobs.json',jobs)
        return job
    analysis=[]
    for offset in range(0,len(manifest['tasks']),1000):
        count=min(1000,len(manifest['tasks'])-offset)
        search_job=launch('search',offset,count)
        analysis.append(launch('analyze',offset,count,search_job))
    for offset in range(0,len(manifest['indices']),1000):
        launch('combine',offset,min(1000,len(manifest['indices'])-offset),':'.join(analysis))
    print(json.dumps(jobs,indent=2),flush=True)


def report(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    rows=[]
    for index in manifest['indices']:
        path=args.run/'results'/str(index)/'result.json'
        rows.append(json.loads(path.read_text()) if path.exists() else dict(index=index,modes={}))
    summary=dict(total=len(rows),modes={})
    for mode in ('single','bidirectional'):
        records=[r['modes'].get(mode,{}) for r in rows]
        counts=Counter(r.get('reference_recovery','pending') for r in records)
        summary['modes'][mode]=dict(outcomes=dict(counts),
            certified_recovery_percent=100*counts['recovered']/len(rows),
            upper_bound_percent=100*(counts['recovered']+counts['unknown']+counts['pending'])/len(rows),
            representative_topk={str(k):sum(r.get('representative',{}).get('topk',{}).get(str(k),False) for r in records)
                                 for k in (1,3,5,10)},
            top5_family_outcomes=dict(Counter(r.get('top5_family',{}).get('outcome','unknown' if r else 'pending') for r in records)),
            incomplete_rankings=sum(not r.get('ranking_complete',False) for r in records))
    summary['completed_mode_records']=sum(bool(r['modes']) for r in rows)
    save(args.run/'summary.json',summary);save(args.run/'cases.json',rows)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    random.seed(42);np.random.seed(42)
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('initialize','search','rank','verify','combine','worker','submit','report'))
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--source',type=Path)
    p.add_argument('--index',type=int)
    p.add_argument('--direction',choices=DIRECTIONS)
    p.add_argument('--phase',choices=('search','analyze','combine'))
    p.add_argument('--offset',type=int,default=0)
    p.add_argument('--slot',type=int,default=0)
    p.add_argument('--partition',default='cpunodes_nia')
    a=p.parse_args();globals()[a.mode](a)
