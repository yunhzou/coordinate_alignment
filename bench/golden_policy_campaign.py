"""Frozen Golden rerun: smaller-first planning and independent cut RNG streams."""
import argparse
from collections import Counter
from dataclasses import asdict
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

from golden_evaluation import prepare, evaluate_planned, project
from rxn_core import AAMProblem,AAMSearchConfig,plan_aam_search,search_aam
from rxn_core.artifacts import read_aam_checkpoint,raw_cut_paths,read_raw_cut
from rxn_core.domain import MolecularEndpoint


def save(path,value):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2)+'\n')
    temp.replace(path)


def load_case(run,index):
    directory=run/str(index)
    raw=json.loads((directory/'input.json').read_text())
    problem=AAMProblem(MolecularEndpoint(**raw['reactant']),MolecularEndpoint(**raw['product']),raw['name'])
    config=AAMSearchConfig(**json.loads((run/'manifest.json').read_text())['config'])
    return directory,plan_aam_search(problem,config)


def initialize(args):
    args.run.mkdir(parents=True,exist_ok=False)
    engine=args.run/'engine'
    shutil.copytree('src',engine/'src',ignore=shutil.ignore_patterns('__pycache__'))
    (engine/'bench').mkdir()
    for name in ('golden_policy_campaign.py','golden_evaluation.py'):
        shutil.copy2(Path(__file__).with_name(name),engine/'bench'/name)
    config=AAMSearchConfig()
    rows=[]
    for line in args.audit.read_text().splitlines():
        row=json.loads(line)
        if args.indices is not None and row['index'] not in args.indices:
            continue
        problem,features,reference=prepare(row['mapped_reaction'])
        plan=plan_aam_search(problem,config)
        directory=args.run/str(row['index']);directory.mkdir()
        def endpoint(e):
            return dict(elements=e.elements,coordinates=e.coordinates.tolist(),wbo=e.wbo.tolist(),
                        label=e.label,metadata=dict(e.metadata))
        save(directory/'input.json',dict(reactant=endpoint(problem.reactant),product=endpoint(problem.product),name=f'golden_{row["index"]}'))
        save(directory/'reference.json',dict(features=features,mapping=sorted(reference.items())))
        save(directory/'orientation.json',dict(direction=plan.direction,reversed=plan.reversed,
            input_R_atoms=problem.source_atom_count,input_P_atoms=problem.target_atom_count,
            graph_space='search source/target; do not interpret reversed graphs as input R/P',
            rule='fewer explicit atoms first; input R first on ties'))
        save(directory/'status.json',dict(stage='pending',index=row['index']))
        rows.append(dict(index=row['index'],direction=plan.direction,
                        reference_annotation_complete=len(project(reference,features))==len(features[1]['heavy']),
                        input_sha256=hashlib.sha256((directory/'input.json').read_bytes()).hexdigest()))
    save(args.run/'manifest.json',dict(records=rows,config=asdict(config),
        search_timeout=300,score_timeout=180,partial_timeout=90,cpu_budget=args.cpu_budget,
        workers_per_reaction=1,memory_per_job='12G',
        seed_policy='independent_per_cut_blake2b_v1',orientation_policy='smaller_explicit_endpoint_first',
        created=time.time(),audit=str(args.audit.resolve()),
        audit_sha256=hashlib.sha256(args.audit.read_bytes()).hexdigest(),
        revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        engine_sha256={str(p.relative_to(engine)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in engine.rglob('*') if p.is_file()}))
    print(json.dumps(dict(records=len(rows),directions=dict(Counter(r['direction'] for r in rows)))))


def search(args):
    directory,plan=load_case(args.run,args.index)
    started=time.perf_counter()
    result=search_aam(plan.problem,plan.config,workers=1,intermediate_dir=directory/'cuts',archive_format='checkpoint')
    usage=resource.getrusage(resource.RUSAGE_SELF)
    save(directory/'search.json',dict(wall_seconds=time.perf_counter()-started,metrics=asdict(result.metrics),
        cpu_seconds=usage.ru_utime+usage.ru_stime,peak_rss_mb=usage.ru_maxrss/1024,
        direction=plan.direction,hostname=os.uname().nodename,slurm_job=os.environ.get('SLURM_JOB_ID')))


def score(args):
    directory,plan=load_case(args.run,args.index)
    result=read_aam_checkpoint(directory/'cuts/aam.pkl.gz')
    ref=json.loads((directory/'reference.json').read_text())
    report=evaluate_planned(result,plan,ref['features'],ref['mapping'],seconds=120)
    report['search_incomplete']=False
    save(directory/'evaluation.json',report)


def partial(args):
    """Explicitly labelled positive-only evidence from completed cut records."""
    import pynauty
    from golden_evaluation import colored_graph
    directory,plan=load_case(args.run,args.index)
    ref=json.loads((directory/'reference.json').read_text())
    chunks=raw_cut_paths(directory/'cuts')
    save(directory/'partial_archive.json',dict(cuts=[str(p) for p in chunks],search_incomplete=True))
    def finish(report):
        save(directory/'evaluation.json',report)
        save(directory/'status.json',dict(stage='complete',index=args.index,
            search_interrupted=True,evaluation_scope='saved partial cuts',updated=time.time()))
    started=time.perf_counter()
    features=list(reversed(ref['features'])) if plan.reversed else ref['features']
    expected=pynauty.certificate(colored_graph(features,project(plan.to_search_mapping(ref['mapping']),features)))
    reports=[]
    for path in chunks:
        if time.perf_counter()-started>75:
            break
        graph=read_raw_cut(path)
        seen=set();witness=None
        for terminal in graph.terminals:
            if time.perf_counter()-started>75:break
            mapping=graph.states[terminal].mapping
            heavy=project(mapping,features);key=tuple(sorted(heavy.items()))
            if key in seen:continue
            seen.add(key)
            if pynauty.certificate(colored_graph(features,heavy))==expected:
                witness=terminal;break
        report=dict(cut=str(path),checked_heavy_representatives=len(seen),
            reference_recovery='recovered' if witness is not None else 'unknown',
            reference_annotation_complete=len(project(ref['mapping'],ref['features']))==len(ref['features'][1]['heavy']),
            witness_terminal=witness,search_direction=plan.direction,
            evaluation_seconds=time.perf_counter()-started)
        reports.append(report)
        save(directory/'partial_checks.json',reports)
        if report['reference_recovery']=='recovered':
            report['input_orientation_witness']=sorted(plan.to_input_mapping(graph.states[witness].mapping).items())
            report.update(search_incomplete=True,top1_correct=None,evaluation_scope='positive-only completed cut, not full sweep')
            finish(report)
            return
    finish(dict(reference_recovery='unknown',top1_correct=None,
        reference_annotation_complete=len(project(ref['mapping'],ref['features']))==len(ref['features'][1]['heavy']),
        search_incomplete=True,search_direction=plan.direction,evaluated_cuts=len(reports),
        completed_cuts=len(chunks),evaluation_scope='no positive witness found in checked partial representatives; not proof of absence'))


def guarded(command,log,timeout):
    cgroup=Path('/sys/fs/cgroup')/Path('/proc/self/cgroup').read_text().strip().split('::')[-1].lstrip('/')
    events=next((p/'memory.events' for p in (cgroup,*cgroup.parents) if p.name.startswith('job_')),None)
    def oom():
        return int(dict(line.split() for line in events.read_text().splitlines())['oom_kill']) if events else 0
    initial=oom();started=time.monotonic()
    with log.open('a') as stream:
        child=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        while True:
            try:
                return child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                reason='oom' if oom()>initial else 'timeout' if time.monotonic()-started>=timeout else None
                if reason:
                    os.killpg(child.pid,signal.SIGTERM)
                    try:child.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        os.killpg(child.pid,signal.SIGKILL);child.wait()
                    return reason


def worker(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    index=manifest['records'][args.slot]['index'];directory=args.run/str(index)
    if (directory/'evaluation.json').exists():return
    def phase(name,timeout):
        save(directory/'status.json',dict(stage=name,index=index,updated=time.time(),job=os.environ.get('SLURM_JOB_ID')))
        code=guarded([sys.executable,__file__,name,'--run',str(args.run),'--index',str(index)],directory/f'{name}.log',timeout)
        save(directory/f'{name}_status.json',dict(exit_code=code,ended=time.time()))
        return code
    if not (directory/'cuts/aam.pkl.gz').exists():
        code=phase('search',manifest['search_timeout'])
        if code!=0:
            if (directory/'cuts/aam.pkl.gz').exists():
                code=phase('score',manifest['score_timeout'])
            else:
                code=phase('partial',manifest['partial_timeout'])
            save(directory/'status.json',dict(stage='complete' if code==0 else 'partial_error',index=index,
                search_interrupted=True,updated=time.time()))
            return
    code=phase('score',manifest['score_timeout'])
    save(directory/'status.json',dict(stage='complete' if code==0 else 'score_error',index=index,updated=time.time()))


def submit(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    settings=subprocess.check_output(['scontrol','show','config'],text=True)
    limit=int(next(s.split('=')[1] for s in settings.splitlines() if s.startswith('MaxArraySize')))
    slots=[i for i,r in enumerate(manifest['records']) if not (args.run/str(r['index'])/'evaluation.json').exists()]
    chunks=sorted({s//limit for s in slots});jobs=[]
    concurrency=max(1,manifest['cpu_budget']//len(chunks))
    for chunk in chunks:
        offset=chunk*limit
        array=','.join(str(s-offset) for s in slots if s//limit==chunk)+f'%{concurrency}'
        job=subprocess.check_output(['sbatch','--parsable',f'--array={array}',
            f'--output={args.run}/slurm_%A_%a.log','hpc/golden_policy.sbatch',str(args.run.resolve()),str(offset)],text=True).strip().split(';')[0]
        jobs.append(dict(job=job,offset=offset));print(job,flush=True)
    save(args.run/'submission.json',dict(jobs=jobs,submitted=time.time(),slots=slots))


def status(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    states=Counter();evaluations=[];unfinished=[]
    for row in manifest['records']:
        directory=args.run/str(row['index']);state=json.loads((directory/'status.json').read_text());states[state['stage']]+=1
        if (directory/'evaluation.json').exists():
            evaluations.append(dict(index=row['index'],**json.loads((directory/'evaluation.json').read_text())))
        else:unfinished.append(row['index'])
    complete=[r for r in evaluations if r.get('reference_annotation_complete')]
    report=dict(total=len(manifest['records']),states=dict(states),evaluated=len(evaluations),unfinished=unfinished,
        complete_reference_total=sum(r['reference_annotation_complete'] for r in manifest['records']),
        reference_recovered=sum(r['reference_recovery']=='recovered' for r in evaluations),
        complete_reference_evaluated=len(complete),complete_reference_recovered=sum(r['reference_recovery']=='recovered' for r in complete),
        top1_correct=sum(r.get('top1_correct') is True for r in complete),
        incomplete_searches=sum(r.get('search_incomplete',False) for r in evaluations),
        unknown=sum(r['reference_recovery']=='unknown' for r in evaluations))
    save(args.run/'progress.json',report);save(args.run/'evaluations.json',evaluations)
    print(json.dumps({k:v for k,v in report.items() if k!='unfinished'},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['init','search','score','partial','worker','submit','status'])
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--audit',type=Path)
    parser.add_argument('--indices',nargs='+',type=int)
    parser.add_argument('--index',type=int)
    parser.add_argument('--slot',type=int)
    parser.add_argument('--cpu-budget',type=int,default=128)
    args=parser.parse_args()
    dict(init=initialize,search=search,score=score,partial=partial,worker=worker,submit=submit,status=status)[args.mode](args)
