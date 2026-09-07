"""Complete interrupted Golden searches from immutable cut checkpoints."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import faulthandler
import gzip
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time

from golden_campaign import initialize,save
from golden_evaluation import evaluate
from rxn_core import AAMProblem,AAMSearchConfig,search_aam
from rxn_core.aam import checkpoint_manifest
from rxn_core.artifacts import read_aam,read_aam_checkpoint,write_aam_checkpoint
from rxn_core.domain import MolecularEndpoint


def problem_at(directory):
    raw=json.loads((directory/'input.json').read_text())
    return AAMProblem(MolecularEndpoint(**raw['reactant']),MolecularEndpoint(**raw['product']),raw['name'])


def init(args):
    with (args.source/'report/cases.csv').open() as stream:
        args.indices=[int(row['index']) for row in csv.DictReader(stream) if row['search_incomplete']=='True']
    original=json.loads((args.source/'manifest.json').read_text())
    args.audit=Path(original['audit']);args.reuse=None;args.cpu_budget=96
    initialize(args)
    shutil.copy2(__file__,args.run/'engine/bench'/Path(__file__).name)
    for name in ('report_golden_campaign.py',):
        shutil.copy2(Path(__file__).with_name(name),args.run/'engine/bench'/name)
    manifest=json.loads((args.run/'manifest.json').read_text())
    for path in (args.run/'engine/bench').glob('*.py'):
        manifest['engine_sha256'][str(path.relative_to(args.run/'engine'))]=hashlib.sha256(path.read_bytes()).hexdigest()
    manifest['source_campaign']=str(args.source)
    save(args.run/'manifest.json',manifest)
    config=AAMSearchConfig(**original['config']);plans=[]
    for index in args.indices:
        old=args.source/str(index);new=args.run/str(index)
        if (old/'input.json').read_bytes()!=(new/'input.json').read_bytes():
            raise ValueError(f'Input changed for {index}')
        if json.loads((args.run/'manifest.json').read_text())['config']!=original['config']:
            raise ValueError('Search configuration changed')
        cuts=new/'cuts';cuts.mkdir()
        raw=old/'cuts/aam.json';complete_raw=False
        if raw.exists() and raw.stat().st_size:
            with raw.open('rb') as stream:
                stream.seek(-1,2);complete_raw=stream.read()==b'\n'
        if (old/'partial_archive.json').exists():
            chunks=[Path(p) for p in json.loads((old/'partial_archive.json').read_text())['cuts']]
        else:chunks=list((old/'cuts').glob('cut_*.json'))
        for chunk in chunks:(cuts/chunk.name).symlink_to(chunk.resolve())
        save(cuts/'manifest.json',checkpoint_manifest(problem_at(new),config))
        plan=dict(index=index,original=str(old),saved_cuts=len(chunks),
                  complete_raw=str(raw) if complete_raw else None)
        save(new/'completion_plan.json',plan);plans.append(plan)
    save(args.run/'completion_plan.json',plans)


def finish(args):
    faulthandler.enable()
    directory=args.run/str(args.index);started=time.perf_counter()
    plan=json.loads((directory/'completion_plan.json').read_text())
    config=AAMSearchConfig(**json.loads((args.run/'manifest.json').read_text())['config'])
    raw=directory/'cuts/aam.json'
    snapshot=directory/'cuts/aam.pkl.gz'
    if snapshot.exists():
        result=read_aam_checkpoint(snapshot);mode='complete_snapshot_reused'
    elif (directory/'aam.json.gz').exists():
        with gzip.open(directory/'aam.json.gz','rt') as stream:result=read_aam(stream)
        mode='complete_archive_reused'
    elif raw.exists() or plan['complete_raw']:
        if not raw.exists():raw=Path(plan['complete_raw'])
        with raw.open() as stream:result=read_aam(stream)
        mode='complete_archive_reused'
    else:
        result=search_aam(problem_at(directory),config,workers=args.cut_workers,
                          intermediate_dir=directory/'cuts',resume=True,archive_format='checkpoint')
        mode='unfinished_cuts_resumed'
    search_seconds=time.perf_counter()-started;archive_started=time.perf_counter()
    if snapshot.exists():(directory/'aam.pkl.gz').hardlink_to(snapshot)
    else:write_aam_checkpoint(result,directory/'aam.pkl.gz')
    save(directory/'search.json',dict(wall_seconds=search_seconds,
        archive_seconds=time.perf_counter()-archive_started,metrics=vars(result.metrics),
        mode=mode,saved_cuts_reused=plan['saved_cuts'],cold_search=False,
        cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,
        cpu_seconds_scope='parent_process_only',
        child_cpu_seconds=resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,
        peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        hostname=os.uname().nodename,slurm_job=os.environ.get('SLURM_JOB_ID'),cpu_model=''))


def score(args):
    directory=args.run/str(args.index)
    result=read_aam_checkpoint(directory/'aam.pkl.gz')
    reference=json.loads((directory/'reference.json').read_text())
    evaluated=evaluate(result,reference['features'],reference['mapping'],seconds=120)
    evaluated.update(search_incomplete=False,
        evaluator_sha256=hashlib.sha256(Path(__file__).with_name('golden_evaluation.py').read_bytes()).hexdigest())
    save(directory/'evaluation.json',evaluated)


def worker(args):
    directory=args.run/str(args.index)
    for phase,artifact,seconds in [('finish','aam.pkl.gz',300),('score','evaluation.json',270)]:
        if (directory/artifact).exists():continue
        save(directory/'status.json',dict(index=args.index,stage=phase,updated=time.time()))
        with (directory/f'{phase}.log').open('a') as log:
            child=subprocess.Popen([sys.executable,__file__,phase,'--run',str(args.run),
                '--index',str(args.index),'--cut-workers',str(args.cut_workers)],
                stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                code=child.wait(timeout=seconds)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGKILL);child.wait()
                save(directory/'status.json',dict(index=args.index,stage=phase+'_timeout'));return
        if code:
            save(directory/'status.json',dict(index=args.index,stage=phase+'_error',code=code));return
    save(directory/'status.json',dict(index=args.index,stage='complete'))


def pool(args):
    records=json.loads((args.run/'manifest.json').read_text())['records']
    indices=[row['index'] for row in records][args.shard::args.shards]
    def run(index):
        directory=args.run/str(index)
        if json.loads((directory/'status.json').read_text())['stage']=='complete':return
        attempt=directory/'attempts'/str(time.time_ns());attempt.mkdir(parents=True)
        shutil.copy2(directory/'status.json',attempt/'previous_status.json')
        for phase in ('finish','score','worker'):
            path=directory/f'{phase}.log'
            if path.exists():shutil.copy2(path,attempt/path.name)
        with (directory/'worker.log').open('a') as log:
            child=subprocess.Popen([sys.executable,__file__,'worker','--run',str(args.run),
                '--index',str(index),'--cut-workers',str(args.cut_workers)],
                stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:child.wait(timeout=600)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGKILL);child.wait()
                save(directory/'status.json',dict(index=index,stage='outer_timeout'))
        print(index,(directory/'status.json').read_text(),flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:list(executor.map(run,indices))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['init','finish','score','worker','pool'])
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--source',type=Path)
    parser.add_argument('--index',type=int)
    parser.add_argument('--shard',type=int,default=0)
    parser.add_argument('--shards',type=int,default=1)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--cut-workers',type=int,default=1)
    args=parser.parse_args();dict(init=init,finish=finish,score=score,worker=worker,pool=pool)[args.mode](args)
