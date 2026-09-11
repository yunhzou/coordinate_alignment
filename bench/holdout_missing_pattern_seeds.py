"""Test seed budgets on the three known forward AAM alternative omissions only.

The frozen mapper reads inputs/configuration only. Exact target patterns are used
after search for membership checks, never to guide fragment growth or cuts.
"""
import argparse
from collections import Counter
import csv
import gzip
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import time

import numpy as np
import pynauty
from holdout_minimum_events import AAM,DATA,PYTHON,EventPatterns,read,save,sha,event_counts,read_aam_checkpoint,colored_graph
from holdout_minimum_event_families import aam_membership
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events

CASES=(11,64,101)
SEEDS=(1,3,10,30)
CACHED=DATA/'adaptive_full_20260910'
FORWARD=Path(__file__).resolve().parents[1]/'reports/holdout_forward_cap1000_seed1_20260910'


def prepare(args):
    args.run.mkdir(parents=True,exist_ok=False)
    (args.run/'status').mkdir()
    baseline=read(AAM/'manifest.json')
    rows=read(FORWARD/'per_case.json')
    targets={}
    for index in CASES:
        row=rows[index];keys=row['comparisons']['slap_sweep']['slap_only_proven']
        assert len(keys)==1 and row['catalogue_complete']['aam']
        targets[str(index)]=dict(index=index,name=row['name'],pattern=row['patterns'][keys[0]],
            baseline_minimum=row['minima']['aam'],source=str(FORWARD/'per_case.json'))
    save(args.run/'targets.json',targets)
    tasks=[dict(dataset='holdout',index=i,direction='R_to_P') for i in CASES]
    variants={}
    for seeds in SEEDS:
        name=f'seed_{seeds}';folder=args.run/'runs'/name;folder.mkdir(parents=True)
        (folder/'inputs').mkdir();(folder/'inputs/holdout').symlink_to(AAM/'inputs/holdout',target_is_directory=True)
        config=dict(baseline['original_config'],seed_count=seeds)
        assert config['branch_limit']==1000
        save(folder/'manifest.json',dict(original_config=config,original_execution=baseline['original_execution'],
            original_workers=8,root_seed=baseline['root_seed'],original_commit=baseline['original_commit']))
        save(folder/'tasks.json',tasks);variants[name]=dict(run=str(folder),seeds=seeds,cap=1000,fresh=True)
    variants['cached_seed10_cap100']=dict(run=str(CACHED),seeds=10,cap=100,fresh=False)
    save(args.run/'variants.json',variants)
    names=['holdout_missing_pattern_seeds.py','holdout_minimum_events.py','holdout_minimum_event_families.py']
    for name in names:shutil.copy2(Path(__file__).with_name(name),args.run/name)
    save(args.run/'manifest.json',dict(cases=CASES,seeds=SEEDS,direction='R_to_P',branch_cap=1000,
        original_engine_commit=baseline['original_commit'],frozen_source=str(AAM),
        config=baseline['original_config'],scope='Targeted diagnostic on three known omissions; not a rerun of the 140-case benchmark.',
        target_policy='Targets withheld from the frozen search; only post-search exact event-family membership uses them.',
        seed_policy='Original deterministic random seed-order prefixes, independent streams per cut, root seed 42.',
        workers=8,search_watchdog=300,analysis_watchdog=360,
        frozen_sha256={name:sha(args.run/name) for name in names},
        sources_sha256={str(p):sha(p) for p in [AAM/'manifest.json',FORWARD/'per_case.json',CACHED/'manifest.json']},
        input_sha256={str(i):sha(AAM/f'inputs/{i}/input.json') for i in CASES}))
    print(args.run,flush=True)


def analyze(args):
    assert args.index in CASES
    variants=read(args.run/'variants.json');variant=variants[args.variant];source=Path(variant['run'])
    out=args.output if args.output is not None else args.run/f'analysis/{args.variant}/{args.index}.json'
    assert not out.exists()
    target=read(args.run/'targets.json')[str(args.index)]['pattern'];target_score=target['total']
    raw=read(source/f'inputs/holdout/{args.index}/input.json');canonical=EventPatterns(raw)
    assert sha(source/f'inputs/holdout/{args.index}/input.json')==sha(AAM/f'inputs/{args.index}/input.json')
    folder=source/f'results/holdout/{args.index}/R_to_P/original'
    search=read(folder/'search.json');assert search['complete']
    started,cpu=time.perf_counter(),time.process_time()
    graph=read_aam_checkpoint(folder/'cuts/aam.pkl.gz').graph
    full=[t for t in graph.terminals if len(graph.states[t].mapping)==canonical.n]
    representatives={};minimum=None;min_patterns={}
    for offset in range(0,len(full),256):
        terminals=full[offset:offset+256]
        vectors=[[dict(graph.states[t].mapping)[a] for a in range(canonical.n)] for t in terminals]
        scores=event_counts(canonical.r,canonical.p,vectors).sum(axis=1)
        for terminal,vector,score in zip(terminals,vectors,scores,strict=True):
            if minimum is None or score<minimum:minimum=int(score);min_patterns={}
            if score==minimum or score==target_score:
                pattern=canonical.describe(vector);value=dict(pattern,terminal=terminal,direction='R_to_P')
                if score==minimum:min_patterns.setdefault(pattern['id'],value)
                if score==target_score:representatives.setdefault(pattern['id'],value)
    assert minimum is not None
    del graph
    snapshot=dict(index=args.index,methods={'aam':dict(minimum=target_score,patterns=representatives)})
    generators=tuple(tuple(g[:canonical.n]) for g in pynauty.autgrp(colored_graph([canonical.features[0]]))[0])
    found,metrics=aam_membership(raw,canonical,{target['id']:target},snapshot,generators,
        time.perf_counter()+args.query_seconds,directions=('R_to_P',),archive_root=source)
    result=found[target['id']]
    if result['status']=='represented':
        actual=canonical.describe(result['witness']['mapping'])
        assert actual['id']==target['id'] and actual['total']==target_score
        problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
        scalar=bond_events(problem,dict(enumerate(actual['mapping'])))
        assert scalar['total']==target_score
    timing=search['rows'][-1]
    output=dict(index=args.index,name=raw['name'],variant=args.variant,seeds=variant['seeds'],cap=variant['cap'],
        fresh=variant['fresh'],direction='R_to_P',minimum_events=minimum,target_score=target_score,
        target_pattern=target['id'],target_result=result,target_in_terminal_representatives=target['id'] in representatives,
        full_terminals_scanned=len(full),minimum_representative_patterns=len(min_patterns),query_metrics=metrics,
        search_cpu=timing['compute_cpu_excluding_persistence_and_loading_seconds'],
        search_wall=timing['elapsed_wall_including_io_seconds'],capped=timing['capped'],search_metrics=timing['metrics'],
        analysis_cpu=time.process_time()-cpu,analysis_wall=time.perf_counter()-started,
        query_seconds=args.query_seconds,
        minimum_witness=next(iter(min_patterns.values())),archive=str(folder/'cuts/aam.pkl.gz'),
        archive_sha256=sha(folder/'cuts/aam.pkl.gz'),search_sha256=sha(folder/'search.json'))
    save(out,output)
    print(dict(index=args.index,seeds=variant['seeds'],cap=variant['cap'],target=result['status'],
               minimum=minimum,cpu=output['search_cpu']),flush=True)


def worker(args):
    index=CASES[args.slot]
    assert index in CASES
    path=args.run/f'status/case_{index}.json';assert not path.exists()
    status=dict(index=index,host=socket.gethostname(),started=time.time(),phases={})
    save(path,status)
    env=dict(os.environ,PYTHONPATH=f'{AAM}/original/src:{AAM}/engine/bench',RXN_CORE_NATIVE='1',
        PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
    variants=read(args.run/'variants.json')
    for name,variant in variants.items():
        commands=[]
        if variant['fresh']:
            commands.append(('search',300,[PYTHON,str(AAM/'engine/bench/adaptive_full_benchmark.py'),'search',
                '--run',variant['run'],'--slot',str(args.slot),'--method','original']))
        commands.append(('analyze',360,[PYTHON,str(args.run/'holdout_missing_pattern_seeds.py'),'analyze',
            '--run',str(args.run),'--variant',name,'--index',str(index)]))
        for phase,seconds,command in commands:
            tick=time.perf_counter()
            with (args.run/f'status/{index}_{name}_{phase}.log').open('w') as stream:
                code=subprocess.run(['timeout','--kill-after=5s',str(seconds),*command],env=env,stdout=stream,stderr=subprocess.STDOUT).returncode
            status['phases'][f'{name}/{phase}']=dict(exit=code,elapsed=time.perf_counter()-tick)
            save(path,status)
            if code!=0:break
    save(path,dict(status,finished=time.time()))


def submit(args):
    assert not (args.run/'submission.json').exists()
    command=['sbatch','--parsable','--partition=cpunodes_nia','--nodelist=bosque83','--nodes=1','--cpus-per-task=8',
        '--mem=32G','--time=00:45:00','--no-requeue','--array=0-2%3','--job-name=three_case_seeds',
        f'--output={args.run}/status/slurm_%A_%a.out','--wrap',shlex.join(['env',
            f'PYTHONPATH={AAM}/original/src:{AAM}/engine/bench','PYTHONDONTWRITEBYTECODE=1','OPENBLAS_NUM_THREADS=1',
            PYTHON,str(args.run/'holdout_missing_pattern_seeds.py'),'worker','--run',str(args.run),'--slot'])+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(command,text=True).strip();save(args.run/'submission.json',dict(job=job,command=command));print(job,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','worker','analyze','submit'))
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--slot',type=int,choices=range(len(CASES)))
    parser.add_argument('--index',type=int,choices=CASES)
    parser.add_argument('--variant',choices=tuple(f'seed_{s}' for s in SEEDS)+('cached_seed10_cap100',))
    parser.add_argument('--query-seconds',type=int,default=240)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();globals()[args.command](args)
