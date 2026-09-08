"""Paired, repeated core and pipeline resource regression experiments."""
import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import pickle
import resource
import signal
import shutil
import statistics
import subprocess
import sys
import time


def save(path,value):
    path=Path(path);temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


def tree_rss(pid):
    pending=[pid];total=0
    while pending:
        current=pending.pop();root=Path('/proc')/str(current)
        try:
            status=(root/'status').read_text()
            pending.extend(map(int,(root/'task'/str(current)/'children').read_text().split()))
        except FileNotFoundError:continue  # process exited between samples
        total+=next((int(line.split()[1]) for line in status.splitlines() if line.startswith('VmRSS:')),0)
    return total/1024


def init(args):
    from cases import CASES
    from rxn_core import AAMProblem,AAMSearchConfig,MolecularEndpoint
    from rxn_core.smiles import smiles_to_weighted_graph
    import numpy as np
    args.run.mkdir(parents=True,exist_ok=False)
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'after'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    problems={name:AAMProblem(*factory(),name) for name,factory in CASES.items()}
    def endpoint(smiles):
        graph=smiles_to_weighted_graph(smiles,expand_hydrogens=True)
        from rxn_core.subgraph import _coerce_graph
        graph=_coerce_graph(graph,.2)
        atoms=list(graph);weights=np.zeros((len(atoms),len(atoms)))
        for a,b,d in graph.edges(data=True):weights[a,b]=weights[b,a]=d['wbo']
        return MolecularEndpoint(tuple(graph.nodes[i]['element'] for i in atoms),np.zeros((len(atoms),3)),weights)
    problems['unequal']=AAMProblem(endpoint('CCBr.O'),endpoint('CCO'),'unequal')
    raw=json.loads(Path('data/retro_runs/aam_memory_boundary_20260905/inputs/inputs.json').read_text())
    tasks=[dict(kind='core',family=f['family'],name=f"augmented_{f['family']}") for f in raw['families']]
    core_count=len(tasks)
    inputs=args.run/'inputs';inputs.mkdir()
    for name,problem in problems.items():
        def record(e):return dict(elements=e.elements,coordinates=e.coordinates.tolist(),wbo=e.wbo.tolist(),label=e.label)
        save(inputs/f'{name}.json',dict(reactant=record(problem.reactant),product=record(problem.product),name=name))
        for workers in (1,4):
            for fmt in ('memory','json','checkpoint'):
                tasks.append(dict(kind='pipeline',name=name,workers=workers,format=fmt))
    save(args.run/'manifest.json',dict(before=str(args.before.resolve()),after=str((args.run/'after').resolve()),
        inputs=str(Path('data/retro_runs/aam_memory_boundary_20260905/inputs').resolve()),
        config=asdict(AAMSearchConfig()),tasks=tasks,core_count=core_count,repeats=args.repeats,
        after_base_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        note='Configured repetitions, alternating version order on the same allocation; core inputs are 59 augmented inputs from one prior experiment, not 59 independent reactions. Cold means no saved search intermediates, not a flushed OS page cache. Frozen source hashes in the report identify the exact tested code, including uncommitted edits.'))
    if args.reuse_core:
        old=json.loads((args.reuse_core/'manifest.json').read_text())
        for p in (args.run/'after/src').rglob('*.py'):
            assert p.read_bytes()==(Path(old['after'])/'src'/p.relative_to(args.run/'after/src')).read_bytes()
        for slot in range(core_count):(args.run/str(slot)).symlink_to((args.reuse_core/str(slot)).resolve())
    print(json.dumps(dict(core=core_count,pipeline=len(tasks)-core_count,total=len(tasks))))


def replay(args):
    m=json.loads((args.run/'manifest.json').read_text());task=m['tasks'][args.slot]
    sys.path.insert(0,str(Path(m[args.version])/'src'))
    from rxn_core import AAMProblem,AAMSearchConfig,MolecularEndpoint,search_aam
    from rxn_core.artifacts import write_aam_checkpoint,read_aam_checkpoint,aam_json,read_aam,write_graph_checkpoint
    out=args.run/str(args.slot)/f'{args.repeat}.{args.version}'
    out.mkdir(parents=True,exist_ok=False);times={};details={}
    def measure(name,fn):
        save(out/'stage.json',dict(stage=name));start=time.perf_counter();result=fn()
        times[name]=time.perf_counter()-start
        save(out/'stage.json',dict(stage='validation'))
        save(out/'partial_timings.json',times)
        return result
    def canonical(graph):return json.dumps(graph.to_record(copy=False),sort_keys=True)
    if task['kind']=='core':
        from rxn_core.alignment.branch import find_islands
        data=(Path(m['inputs'])/f"{task['family']}.input.pkl").read_bytes()
        positional,keywords=pickle.loads(data)
        graph=measure('core',lambda:find_islands(*positional,**keywords))
        details['input_sha256']=hashlib.sha256(data).hexdigest()
    else:
        raw=json.loads((args.run/'inputs'/f"{task['name']}.json").read_text())
        problem=AAMProblem(MolecularEndpoint(**raw['reactant']),MolecularEndpoint(**raw['product']),raw['name'])
        config=AAMSearchConfig(**m['config']);options=dict(workers=task['workers'])
        if task['format']!='memory':options.update(intermediate_dir=out/'cuts',archive_format=task['format'])
        result=measure('cold_search',lambda:search_aam(problem,config,**options));graph=result.graph
        details['cold_metrics']=asdict(result.metrics)
        if task['format']!='memory':
            resumed=measure('warm_resume',lambda:search_aam(problem,config,resume=True,**options))
            assert canonical(resumed.graph)==canonical(graph)
            assert resumed.metrics.worker_search_seconds==0
            details['resume_equal']=True
        measure('json_write',lambda:(out/'result.json').write_text(aam_json(result)))
        def read_json():
            with (out/'result.json').open() as stream:return read_aam(stream)
        restored=measure('json_read',read_json);assert canonical(restored.graph)==canonical(graph);del restored
        measure('binary_write',lambda:write_aam_checkpoint(result,out/'result.pkl.gz'))
        restored=measure('binary_read',lambda:read_aam_checkpoint(out/'result.pkl.gz'))
        assert canonical(restored.graph)==canonical(graph);del restored
        details['json_bytes']=(out/'result.json').stat().st_size
        details['binary_bytes']=(out/'result.pkl.gz').stat().st_size
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    save(out/'stage.json',dict(stage='validation'))
    digest=hashlib.sha256(json.dumps(graph.to_record(copy=False),sort_keys=True).encode()).hexdigest()
    if task['kind']=='core':write_graph_checkpoint(graph,out/'graph.pkl.gz')
    save(out/'result_metrics.json',dict(task=task,version=args.version,repeat=args.repeat,times=times,
        peak_process_mib=peak,graph_sha256=digest,states=len(graph.states),transitions=len(graph.transitions),
        terminals=len(graph.terminals),capped=graph.capped,**details))


def worker(args):
    m=json.loads((args.run/'manifest.json').read_text());out=args.run/str(args.slot);out.mkdir()
    statuses=[];start=time.monotonic()
    for repeat in range(m['repeats']):
        versions=['before','after'] if (args.slot+repeat)%2==0 else ['after','before']
        for version in versions:
            if time.monotonic()-start>540:
                statuses.append(dict(repeat=repeat,version=version,status='job_watchdog'));continue
            command=[sys.executable,__file__,'replay','--run',str(args.run),'--slot',str(args.slot),
                     '--repeat',str(repeat),'--version',version]
            with (out/f'{repeat}.{version}.log').open('w') as stream:
                proc=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
                deadline=time.monotonic()+min(180,570-(time.monotonic()-start));peaks={}
                while True:
                    stage_file=out/f'{repeat}.{version}'/'stage.json'
                    stage=json.loads(stage_file.read_text())['stage'] if stage_file.exists() else 'startup'
                    peaks[stage]=max(peaks.get(stage,0),tree_rss(proc.pid))
                    try:code=proc.wait(timeout=.1);break
                    except subprocess.TimeoutExpired:
                        if time.monotonic()>deadline:
                            os.killpg(proc.pid,signal.SIGKILL);proc.wait();code='timeout';break
                save(out/f'{repeat}.{version}'/'sampled_memory.json',dict(
                    peak_tree_rss_mib=peaks,interval_seconds=.1,
                    note='Sum of process RSS, sampled; shared pages can be counted in multiple processes.'))
            statuses.append(dict(repeat=repeat,version=version,status=code));save(out/'statuses.json',statuses)


def report(args):
    m=json.loads((args.run/'manifest.json').read_text());rows=[]
    provenance={version:{str(p.relative_to(Path(m[version]))):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((Path(m[version])/'src').rglob('*'))
        if p.is_file() and p.suffix in ('.py','.so','.cpp','.h')}
        for version in ('before','after')}
    for slot,task in enumerate(m['tasks']):
        records={v:[json.loads(p.read_text()) for p in sorted((args.run/str(slot)).glob(f'*.{v}/result_metrics.json'))]
                 for v in ('before','after')}
        hashes={r['graph_sha256'] for rs in records.values() for r in rs}
        stages=sorted({stage for rs in records.values() for r in rs for stage in r['times']})
        times={v:{stage:statistics.median(r['times'][stage] for r in rs) for stage in stages} for v,rs in records.items() if rs}
        rows.append(dict(slot=slot,task=task,counts={v:len(rs) for v,rs in records.items()},
            equal=len(hashes)==1 and all(records.values()),median_seconds=times,
            ratios={stage:times['after'][stage]/times['before'][stage] for stage in stages} if len(times)==2 else {},
            records=records))
    args.output.mkdir(parents=True,exist_ok=True);save(args.output/'results.json',rows);save(args.output/'manifest.json',m)
    save(args.output/'source_sha256.json',provenance)
    print(json.dumps(dict(tasks=len(rows),complete=sum(all(n==m['repeats'] for n in r['counts'].values()) for r in rows),
        mismatches=[r['slot'] for r in rows if all(r['counts'].values()) and not r['equal']],
        slower=[dict(slot=r['slot'],task=r['task'],ratio=r['ratios']) for r in rows if any(v>1.1 for v in r['ratios'].values())]),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['init','replay','worker','report'])
    p.add_argument('--run',type=Path,required=True);p.add_argument('--before',type=Path);p.add_argument('--slot',type=int)
    p.add_argument('--repeat',type=int);p.add_argument('--version',choices=['before','after']);p.add_argument('--output',type=Path)
    p.add_argument('--reuse-core',type=Path)
    p.add_argument('--repeats',type=int,default=3)
    a=p.parse_args();dict(init=init,replay=replay,worker=worker,report=report)[a.mode](a)
