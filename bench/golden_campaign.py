"""Checkpointed Golden campaign. Search and evaluation are independent jobs."""
import argparse
from dataclasses import asdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import sys
import time

from golden_evaluation import prepare, evaluate
from rxn_core import AAMProblem, AAMSearchConfig, search_aam
from rxn_core.artifacts import aam_from_record, aam_record
from rxn_core.domain import MolecularEndpoint


def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n');temporary.replace(path)


def initialize(args):
    args.run.mkdir(parents=True,exist_ok=False)
    engine=args.run/'engine'
    shutil.copytree(Path('src'),engine/'src',ignore=shutil.ignore_patterns('__pycache__'))
    (engine/'bench').mkdir()
    for name in ('golden_campaign.py','golden_evaluation.py'):
        shutil.copy2(Path(__file__).with_name(name),engine/'bench'/name)
    rows=[json.loads(s) for s in args.audit.read_text().splitlines()]
    indices=args.indices if args.indices is not None else [r['index'] for r in rows]
    manifests=[]
    for row in rows:
        if row['index'] not in indices:continue
        directory=args.run/str(row['index']);directory.mkdir()
        problem,features,reference=prepare(row['mapped_reaction'])
        def endpoint(e):
            return dict(elements=e.elements,coordinates=e.coordinates.tolist(),wbo=e.wbo.tolist(),
                        label=e.label,metadata=dict(e.metadata))
        payload=dict(reactant=endpoint(problem.reactant),product=endpoint(problem.product),
                     name=f'golden_{row["index"]}')
        save(directory/'input.json',payload)
        save(directory/'reference.json',dict(mapping=sorted(reference.items()),features=features))
        save(directory/'status.json',dict(stage='pending',index=row['index']))
        if args.reuse:
            previous=args.reuse/str(row['index'])
            if (previous/'aam.json.gz').exists() or (previous/'partial_archive.json').exists():
                old_manifest=json.loads((args.reuse/'manifest.json').read_text())
                if old_manifest['config'] != json.loads(json.dumps(asdict(AAMSearchConfig()))) or (previous/'input.json').read_bytes()!=(directory/'input.json').read_bytes():
                    raise ValueError('Cannot reuse a different input or AAM configuration')
                for name in ('aam.json.gz','search.json','partial.json','search_timeout.json'):
                    if (previous/name).exists():shutil.copy2(previous/name,directory/name)
                save(directory/'reuse.json',dict(source=str(previous),search_reused=True))
                if (previous/'partial_archive.json').exists():
                    for name in ('partial_archive.json','evaluation.json'):
                        shutil.copy2(previous/name,directory/name)
                    save(directory/'status.json',dict(stage='complete',index=row['index'],reused=True,incomplete_search=True))
        manifests.append(dict(index=row['index'],source_atoms=problem.source_atom_count,
                              target_atoms=problem.target_atom_count,balanced=problem.balanced,
                              input_sha256=hashlib.sha256((directory/'input.json').read_bytes()).hexdigest()))
    manifest=dict(created=time.time(),revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        engine_sha256={str(p.relative_to(engine)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in engine.rglob('*') if p.is_file()},
        worktree_diff_sha256=hashlib.sha256(subprocess.check_output(['git','diff'])).hexdigest(),
        audit=str(args.audit.resolve()),audit_sha256=hashlib.sha256(args.audit.read_bytes()).hexdigest(),
        config=asdict(AAMSearchConfig()),search_timeout=300,evaluation_timeout=120,
        native=True,workers_per_reaction=1,allocated_cpus_per_task=2,cpu_budget=args.cpu_budget,records=manifests,
        evaluation_policy='heavy-atom paired-graph equality with charge/isotope/H-count/CIP/EZ labels; '
          'joint endpoint automorphisms; compressed path group products queried symbolically',
        rank_policy='unmapped P heavy, unmapped P all, bond edits, lexicographic representative mapping')
    save(args.run/'manifest.json',manifest)
    print(json.dumps(dict(run=str(args.run),records=len(manifests))),flush=True)


def search(args):
    directory=args.run/str(args.index)
    raw=json.loads((directory/'input.json').read_text())
    problem=AAMProblem(MolecularEndpoint(**raw['reactant']),MolecularEndpoint(**raw['product']),raw['name'])
    manifest=json.loads((args.run/'manifest.json').read_text())
    started=time.perf_counter()
    result=search_aam(problem,AAMSearchConfig(**manifest['config']),workers=1,
                      intermediate_dir=directory/'cuts')
    search_elapsed=time.perf_counter()-started
    archive_started=time.perf_counter()
    temporary=directory/'aam.json.gz.tmp'
    with gzip.open(temporary,'wt',compresslevel=1) as stream:json.dump(aam_record(result),stream)
    temporary.replace(directory/'aam.json.gz')
    save(directory/'search.json',dict(wall_seconds=search_elapsed,metrics=vars(result.metrics),
        archive_seconds=time.perf_counter()-archive_started,
        peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,
        hostname=platform.node(),slurm_job=os.environ.get('SLURM_JOB_ID'),
        cpu_model=next((s.split(':',1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines()
                        if s.startswith('model name')),'')))


def score(args):
    directory=args.run/str(args.index)
    with gzip.open(directory/'aam.json.gz','rt') as stream:archive=json.load(stream)
    result=aam_from_record(archive)
    reference=json.loads((directory/'reference.json').read_text())
    evaluation=evaluate(result,reference['features'],reference['mapping'],seconds=120)
    evaluation['search_incomplete']=archive.get('benchmark_completion')=='search_timeout'
    save(directory/'evaluation.json',evaluation)


def partial_score(args):
    """Stream completed cuts; a positive recovery remains valid after timeout.

    Do not merge/rebuild a giant partial archive merely to ask whether one
    recorded family contains the reference. Original compressed cuts persist.
    """
    from rxn_core.domain import AAMResult,AAMSearchMetrics
    from rxn_core.search_graph import AAMSearchGraph
    directory=args.run/str(args.index)
    chunks=sorted((directory/'cuts').glob('cut_*.json'))
    if (directory/'partial_archive.json').exists():
        chunks=[Path(p) for p in json.loads((directory/'partial_archive.json').read_text())['cuts']]
    save(directory/'partial_archive.json',dict(cuts=[str(p.resolve()) for p in chunks],
        search_incomplete=True,scope='All completed cut checkpoints; no fabricated missing cut results'))
    raw=json.loads((directory/'input.json').read_text())
    reference=json.loads((directory/'reference.json').read_text())
    problem=AAMProblem(MolecularEndpoint(**raw['reactant']),MolecularEndpoint(**raw['product']),raw['name'])
    config=AAMSearchConfig(**json.loads((args.run/'manifest.json').read_text())['config'])
    started=time.perf_counter();reports=[]
    for chunk in chunks:
        if time.perf_counter()-started>110:break
        graph=AAMSearchGraph.from_record(json.loads(chunk.read_text()))
        result=AAMResult(problem,config,graph,AAMSearchMetrics.from_record({},300))
        report=evaluate(result,reference['features'],reference['mapping'],symbolic=False)
        report['cut']=str(chunk)
        reports.append(report)
        if report['reference_recovery']=='recovered':break
    if not reports:raise ValueError('No completed cut available to score')
    selected=next((r for r in reports if r['reference_recovery']=='recovered'),reports[-1]).copy()
    selected.update(search_incomplete=True,top1_correct=None,reference_recovery=(
        'recovered' if any(r['reference_recovery']=='recovered' for r in reports) else 'unknown'),
        best_target_heavy_coverage=max(r['best_target_heavy_coverage'] for r in reports),
        best_target_all_atom_coverage=max(r['best_target_all_atom_coverage'] for r in reports),
        evaluated_cuts=len(reports),completed_cuts=len(chunks),evaluation_seconds=time.perf_counter()-started,
        evaluation_scope='Representative lower bound from completed cuts of a timed-out search; '
                         'negative means unknown, not absence from the compressed families')
    save(directory/'evaluation.json',selected)


def worker(args):
    directory=args.run/str(args.index)
    started=time.time()
    if (directory/'partial_archive.json').exists() and (directory/'evaluation.json').exists():
        save(directory/'status.json',dict(stage='complete',index=args.index,incomplete_search=True));return
    if not (directory/'aam.json.gz').exists() and json.loads((directory/'status.json').read_text())['stage']=='search_timeout':
        with (directory/'partial.log').open('a') as stream:
            try:
                completed=subprocess.run([sys.executable,__file__,'partial_score','--run',str(args.run),
                    '--index',str(args.index)],stdout=stream,stderr=subprocess.STDOUT,timeout=120)
            except subprocess.TimeoutExpired:return
            if completed.returncode:return
        save(directory/'status.json',dict(stage='complete',index=args.index,incomplete_search=True));return
    for phase,artifact,timeout in [('search','aam.json.gz',300),('score','evaluation.json',135)]:
        if (directory/artifact).exists():continue
        save(directory/'status.json',dict(index=args.index,stage=phase,started=started,updated=time.time()))
        command=[sys.executable,__file__,phase,'--run',str(args.run),'--index',str(args.index)]
        with (directory/f'{phase}.log').open('a') as stream:
            try:
                completed=subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout,
                    env=dict(os.environ,RXN_CORE_NATIVE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'))
                if completed.returncode:
                    save(directory/'status.json',dict(stage=f'{phase}_error',exit_code=completed.returncode,
                                                     index=args.index,elapsed=time.time()-started));return
            except subprocess.TimeoutExpired:
                save(directory/'status.json',dict(stage=f'{phase}_timeout',seconds=timeout,index=args.index,
                                                  elapsed=time.time()-started))
                if phase!='search':return
                save(directory/'search_timeout.json',dict(seconds=300,index=args.index))
                with (directory/'partial.log').open('a') as partial_log:
                    try:
                        completed=subprocess.run([sys.executable,__file__,'partial_score','--run',str(args.run),
                            '--index',str(args.index)],stdout=partial_log,stderr=subprocess.STDOUT,timeout=120)
                    except subprocess.TimeoutExpired:return
                    if completed.returncode:return
                save(directory/'status.json',dict(stage='complete',index=args.index,incomplete_search=True));return
    save(directory/'status.json',dict(stage='complete',index=args.index,elapsed=time.time()-started))


def submit(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    slots=[i for i,r in enumerate(manifest['records'])
           if not (args.run/str(r['index'])/'evaluation.json').exists()]
    config=subprocess.check_output(['scontrol','show','config'],text=True)
    limit=int(next(line.split('=')[1] for line in config.splitlines() if line.startswith('MaxArraySize')))
    chunks=sorted({slot//limit for slot in slots})
    concurrency=max(1,manifest['cpu_budget']//manifest['allocated_cpus_per_task']//len(chunks))
    jobs=[]
    for chunk in chunks:
        offset=chunk*limit
        array=','.join(str(slot-offset) for slot in slots if slot//limit==chunk)+'%'+str(concurrency)
        job=subprocess.check_output(['sbatch','--parsable',f'--array={array}','--cpus-per-task=2',
            f'--output={args.run}/slurm_%A_%a.out',f'--error={args.run}/slurm_%A_%a.err',
            'hpc/golden_benchmark.sbatch',str(args.run.resolve()),str(offset)],text=True).strip().split(';')[0]
        jobs.append(dict(job=job,offset=offset));print(job,flush=True)
    save(args.run/'submission.json',dict(jobs=jobs,submitted=time.time(),slots=slots))


def status(args):
    from collections import Counter
    import statistics
    manifest=json.loads((args.run/'manifest.json').read_text())
    states=Counter();evaluations=[];searches=[];slow=[]
    for row in manifest['records']:
        directory=args.run/str(row['index'])
        state=json.loads((directory/'status.json').read_text());states[state['stage']]+=1
        if state['stage'] in ('search','score') and time.time()-state['updated']>600:slow.append(row['index'])
        if (directory/'evaluation.json').exists():evaluations.append(json.loads((directory/'evaluation.json').read_text()))
        if (directory/'search.json').exists():searches.append(json.loads((directory/'search.json').read_text()))
    times=sorted(s['wall_seconds'] for s in searches)
    report=dict(total=len(manifest['records']),stages=dict(states),evaluated=len(evaluations),
        reference_recovered=sum(e['reference_recovery']=='recovered' for e in evaluations),
        reference_not_recovered=sum(e['reference_recovery']=='not_recovered' for e in evaluations),
        reference_unknown=sum(e['reference_recovery']=='unknown' for e in evaluations),
        top1_correct=sum(e['top1_correct'] is True for e in evaluations),
        incomplete_reference_annotations=sum(not e.get('reference_annotation_complete',True) for e in evaluations),
        full_P_heavy_coverage=sum(e['best_target_heavy_coverage']==1 for e in evaluations),
        full_P_all_atom_coverage=sum(e['best_target_all_atom_coverage']==1 for e in evaluations),
        capped=sum(e['capped'] for e in evaluations),supervisor_overdue=slow,
        incomplete_searches=sum(e.get('search_incomplete',False) for e in evaluations),
        search_median_seconds=statistics.median(times) if times else None,
        search_max_seconds=max(times) if times else None,
        search_p95_seconds=times[min(len(times)-1,int(.95*len(times)))] if times else None,
        cpu_seconds=sum(s['cpu_seconds'] for s in searches),
        peak_rss_mb=max((s['peak_rss_mb'] for s in searches),default=0))
    save(args.run/'progress.json',report);print(json.dumps(report),flush=True)
    from html import escape
    lines=[]
    for row in manifest['records']:
        directory=args.run/str(row['index'])
        state=json.loads((directory/'status.json').read_text())['stage']
        evaluation=json.loads((directory/'evaluation.json').read_text()) if (directory/'evaluation.json').exists() else {}
        coverage=(str(round(100*evaluation['best_target_heavy_coverage'],1))+'%'
                  if 'best_target_heavy_coverage' in evaluation else '—')
        cells=[row['index'],state,evaluation.get('reference_recovery','—'),
               evaluation.get('top1_correct','—'),coverage,evaluation.get('capped','—'),
               evaluation.get('search_incomplete','—')]
        lines.append('<tr>'+''.join('<td>'+escape(str(cell))+'</td>' for cell in cells)+
            '<td><a href="'+str(row['index'])+'/evaluation.json">scores</a> · <a href="'+
            str(row['index'])+'/input.json">input</a></td></tr>')
    page='<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="30"><title>Golden AAM progress</title>'
    page+='<style>body{font:14px system-ui;margin:24px}td,th{padding:7px;border-bottom:1px solid #ddd;text-align:left}pre{background:#f1f5f9;padding:16px}</style>'
    page+='<h1>Golden AAM · '+str(len(manifest['records']))+' reactions</h1><p>Reference recovery ≠ top-1. Coverage is achieved by one returned mapping, not a union of incompatible branches. Search/evaluation failures remain in the denominator.</p>'
    page+='<pre>'+escape(json.dumps(report,indent=2))+'</pre><table><tr><th>Index</th><th>Stage</th><th>Reference recovered</th><th>Top-1 representative</th><th>P heavy coverage</th><th>Cap hit</th><th>Incomplete search</th><th>Evidence</th></tr>'+''.join(lines)+'</table>'
    (args.run/'index.html').write_text(page)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['init','search','score','partial_score','worker','submit','status'])
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--audit',type=Path)
    parser.add_argument('--indices',type=int,nargs='+')
    parser.add_argument('--index',type=int)
    parser.add_argument('--slot',type=int)
    parser.add_argument('--cpu-budget',type=int,default=32)
    parser.add_argument('--reuse',type=Path)
    args=parser.parse_args()
    if args.slot is not None:
        args.index=json.loads((args.run/'manifest.json').read_text())['records'][args.slot]['index']
    dict(init=initialize,search=search,score=score,partial_score=partial_score,worker=worker,submit=submit,status=status)[args.mode](args)


if __name__=='__main__':main()
