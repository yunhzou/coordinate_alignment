"""Full saved-Golden pattern extraction, with immutable inputs and explicit budgets."""
import argparse
from concurrent.futures import ThreadPoolExecutor,ProcessPoolExecutor
from dataclasses import asdict
import gzip
import hashlib
from itertools import chain
import json
import os
import multiprocessing
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time

from golden_policy_campaign import save
from golden_publication import plans,metadata,sha256
from collect_golden_patterns import first_paths
from view_golden_mapping import molecules
from publication_analysis import certificate_id
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.pattern_collection import PatternEquivalence,extract_path_patterns
from rxn_core.search_graph import frozen_value


def _baseline_initialize(equivalence):
    global _baseline_equivalence
    _baseline_equivalence=equivalence


def _baseline_key(mapping):
    return _baseline_equivalence.key(mapping)


def baseline_keys(equivalence,mappings,workers):
    """Independent exact certificates; ordered output keeps ranking reproducible."""
    if workers==1:
        yield from map(equivalence.key,mappings)
    else:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('fork'),
                initializer=_baseline_initialize,initargs=(equivalence,)) as pool:
            yield from pool.map(_baseline_key,mappings,chunksize=64)


def family_key(path):
    """Deduplicate identical query inputs, not merely identical witnesses."""
    steps=[]
    for eid in path.transitions:
        edge=path.graph.transitions[eid]
        if edge.match is not None:
            steps.append((path.graph.states[edge.source].mapping,frozen_value(edge.match)))
    context=path.context
    value=(tuple(sorted(path.mapping.items())),context.cuts,context.graph_floor,context.iso_tolerance,steps)
    return hashlib.sha256(repr(value).encode()).hexdigest()


def worker(args):
    root=args.run;manifest=json.loads((root/'manifest.json').read_text())
    source=Path(manifest['source']);index=args.index;direction=args.direction
    out=root/'results'/str(index)/direction;out.mkdir(parents=True,exist_ok=True)
    pair,default=plans(source,index);plan=pair[direction]
    base=source/'directions'/str(index)/direction;archive=base/'cuts/aam.pkl.gz'
    summary=dict(index=index,direction=direction,default_direction=default,complete=False,
                 environment=metadata(),stage='loading',source_archive=str(archive),baseline_workers=args.baseline_workers)
    save(out/'summary.json',summary)
    if not archive.exists():
        summary.update(stage='missing_archive');save(out/'summary.json',summary);return
    t=time.perf_counter();aam=read_aam_checkpoint(archive)
    summary['archive_loading_seconds']=time.perf_counter()-t
    t=time.perf_counter();summary['archive_sha256']=sha256(archive)
    summary['archive_hashing_seconds']=time.perf_counter()-t
    reaction=json.loads((root/'inputs'/f'{index}.json').read_text())['mapped_reaction']
    mols=molecules(reaction,plan.input_problem)
    atoms=tuple({a.GetIdx():(a.GetFormalCharge(),a.GetIsotope(),
        a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else '') for a in mol.GetAtoms()} for mol in mols)
    bonds=tuple({tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))):str(b.GetStereo())
        for b in mol.GetBonds()} for mol in mols)
    if plan.reversed:atoms=atoms[::-1];bonds=bonds[::-1]
    eq=PatternEquivalence(plan.problem,atom_tags=atoms,bond_tags=bonds)
    first=first_paths(aam.graph)
    ranked=json.loads((base/'classes.json').read_text()) if (base/'classes.json').exists() else []
    # Missing prior ranking is recorded; it never removes terminal families.
    summary['prior_ranking_available']=(base/'classes.json').exists()
    reference=json.loads((source/'inputs'/str(index)/'reference.json').read_text())
    expected=certificate_id(reference['features'],reference['mapping'])
    candidate_rows={};persist_seconds=0.;verification_seconds=0.;visited=0;duplicates=0;unresolved=0
    start=time.perf_counter();cpu=time.process_time();seen=set()
    pattern_file=gzip.open(out/'patterns.jsonl.gz','wt')
    path_file=gzip.open(out/'paths.jsonl.gz','wt')
    def write(stream,record):
        nonlocal persist_seconds
        t=time.perf_counter();stream.write(json.dumps(record,separators=(',',':'))+'\n');stream.flush()
        persist_seconds+=time.perf_counter()-t
    def candidate(record,certified_baseline=None):
        nonlocal verification_seconds
        key=record['key']
        if key in candidate_rows:return
        # Reference labels are computed only after reference-blind extraction.
        t=time.perf_counter();mapping=plan.to_input_mapping(record['mapping'])
        correct=(certificate_id(reference['features'],mapping)==expected if certified_baseline is None
                 else certified_baseline['id']==expected)
        verification_seconds+=time.perf_counter()-t
        heavy=sum(plan.input_problem.reactant.elements[r]!='H' for r in mapping)
        row=dict(key=key,events=record['events'],heavy=heavy,total=len(mapping),reference_equivalent=correct)
        candidate_rows[key]=row;write(pattern_file,dict(record,**{k:v for k,v in row.items() if k!='key'}))
    try:
        # Preserve all prior displayed representatives before extending families.
        mappings=(dict(aam.graph.states[row['terminal']].mapping) for row in ranked)
        for row,key in zip(ranked,baseline_keys(eq,mappings,args.baseline_workers)):
            path=first(row['terminal']);mapping=path.mapping
            events=dict(row['events']);events['total']=sum(events.values())
            candidate(dict(key=key,mapping=sorted(mapping.items()),actions=[],terminal=path.terminal,
                transitions=list(path.transitions),events=events,
                origin='saved_representative'),certified_baseline=row)
        summary['baseline_candidates']=len(candidate_rows)
        summary['baseline_reference']=any(r['reference_equivalent'] for r in candidate_rows.values())
        save(out/'summary.json',dict(summary,stage='extracting'))
        for path in chain((first(row['terminal']) for row in ranked),aam.graph.paths()):
            remaining=manifest['case_seconds']-(time.perf_counter()-start-persist_seconds-verification_seconds)
            if remaining<=0:break
            key=family_key(path)
            if key in seen:duplicates+=1;continue
            seen.add(key);visited+=1
            result=extract_path_patterns(path,plan.problem,eq,seconds=min(remaining,manifest['path_seconds']),
                reverse=plan.reversed,on_pattern=candidate)
            unresolved+=not result['complete']
            write(path_file,dict(family_key=key,terminal=path.terminal,transitions=list(path.transitions),
                **{k:v for k,v in result.items() if k!='patterns'},pattern_keys=[p['key'] for p in result['patterns']]))
        else:summary['complete']=unresolved==0
    finally:
        t=time.perf_counter();pattern_file.close();path_file.close();persist_seconds+=time.perf_counter()-t
    rows=sorted(candidate_rows.values(),key=lambda r:(-r['heavy'],-r['total'],r['events']['total'],r['key']))
    best=rows[0] if rows else None
    eligible=[r for r in rows if best and (r['heavy'],r['total'])==(best['heavy'],best['total'])]
    summary.update(stage='finished',visited_families=visited,duplicate_histories=duplicates,
        unresolved_families=unresolved,candidates=len(rows),new_candidates=len(rows)-summary['baseline_candidates'],
        reference_recovered=any(r['reference_equivalent'] for r in rows),
        event_windows={str(n):any(r['reference_equivalent'] and r['events']['total']<=best['events']['total']+n
                                 for r in eligible) for n in range(21)},
        extraction_wall_seconds=time.perf_counter()-start-persist_seconds-verification_seconds,
        worker_cpu_seconds=time.process_time()-cpu,persistence_wall_seconds=persist_seconds,
        reference_verification_wall_seconds=verification_seconds,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        archive_terminals=len(aam.graph.terminals))
    save(out/'candidates.json',rows);save(out/'summary.json',summary)


def initialize(args):
    args.run.mkdir(parents=True,exist_ok=False)
    for folder in ('inputs','results','logs'): (args.run/folder).mkdir()
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    audit=Path('data/aam_benchmarks/golden_original_20260906/audit.jsonl')
    records=list(map(json.loads,audit.read_text().splitlines()))
    for row in records:save(args.run/'inputs'/f'{row["index"]}.json',row)
    save(args.run/'manifest.json',dict(source=str(args.source.resolve()),indices=sorted(r['index'] for r in records),
        git_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        engine_sha256={str(p.relative_to(args.run/'engine')):sha256(p)
                       for p in (args.run/'engine').rglob('*.py')},audit_sha256=sha256(audit),
        case_seconds=args.case_seconds,path_seconds=args.path_seconds,watchdog_seconds=300,
        scoring_tolerance=.5,explicit_hydrogen=True,random_seed=0,
        scope='All saved terminal histories, exact duplicate query inputs merged; unresolved work is reported.',
        search='No AAM rerun: frozen seed10 cap100 tol1 bidirectional archives.'))


def dispatch(args):
    manifest=json.loads((args.run/'manifest.json').read_text())
    tasks=[(i,d) for i in manifest['indices'] for d in ('R_to_P','P_to_R')]
    tasks=[t for n,t in enumerate(tasks) if n%args.shards==args.shard]
    def run(task):
        i,d=task;directory=args.run/'results'/str(i)/d
        summary=directory/'summary.json'
        if summary.exists() and json.loads(summary.read_text()).get('stage')=='finished':return
        cmd=[sys.executable,str(Path(__file__).resolve()),'worker','--run',str(args.run),'--index',str(i),'--direction',d,
             '--baseline-workers',str(args.baseline_workers)]
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONHASHSEED='0')
        with (args.run/'logs'/f'{i}_{d}.log').open('w') as log:
            process=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:
                code=process.wait(timeout=manifest['watchdog_seconds'])
                status='ok' if code==0 else f'exit_{code}'
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL);process.wait();status='watchdog_timeout'
        save(args.run/'logs'/f'{i}_{d}.status.json',dict(index=i,direction=d,status=status))
        print(i,d,status,flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:list(pool.map(run,tasks))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('init','worker','dispatch'))
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--source',type=Path)
    parser.add_argument('--index',type=int);parser.add_argument('--direction',choices=('R_to_P','P_to_R'))
    parser.add_argument('--case-seconds',type=float,default=60);parser.add_argument('--path-seconds',type=float,default=3)
    parser.add_argument('--workers',type=int,default=16);parser.add_argument('--shard',type=int,default=0);parser.add_argument('--shards',type=int,default=1)
    parser.add_argument('--baseline-workers',type=int,default=1,help='Processes per archive for independent saved baseline certificates')
    args=parser.parse_args();{'init':initialize,'worker':worker,'dispatch':dispatch}[args.action](args)
