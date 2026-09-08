"""Reference-blind anchor policy pilot outside the production AAM pipeline."""
import argparse
from collections import Counter
from dataclasses import replace,asdict
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

from golden_policy_campaign import load_case,save,guarded


def proposals(problem,budget,seed=42):
    """Read endpoint graphs only; never reference mappings or saved outcomes."""
    endpoints=(problem.reactant,problem.product)
    counts=[Counter(e.elements) for e in endpoints]
    neighborhoods=[[Counter((e.elements[j],float(e.wbo[i,j]))
        for j in range(e.atom_count) if e.wbo[i,j]>.2)
        for i in range(e.atom_count)] for e in endpoints]
    pairs=[(r,p) for r,e in enumerate(endpoints[0].elements)
           for p,f in enumerate(endpoints[1].elements) if e==f]
    rng=random.Random(seed);rng.shuffle(pairs)
    def rarity(pair):
        element=endpoints[0].elements[pair[0]]
        return counts[0][element]*counts[1][element]
    def similarity(pair):
        a,b=neighborhoods[0][pair[0]],neighborhoods[1][pair[1]]
        total=sum((a|b).values())
        return sum((a&b).values())/total if total else 1.
    result={'random_anchor':pairs[:budget],
        'stable_anchor':sorted(pairs,key=lambda q:(rarity(q),-similarity(q)))[:budget],
        'change_anchor':sorted(pairs,key=lambda q:(rarity(q),similarity(q)))[:budget]}
    edges=[(a,b) for a in range(endpoints[0].atom_count) for b in range(a+1,endpoints[0].atom_count)
           if endpoints[0].wbo[a,b]>.2]
    rng.shuffle(edges);result['random_single_cut']=edges[:budget]
    return result


def init(args):
    args.run.mkdir(parents=True,exist_ok=False);tasks=[];selections={}
    for index in args.indices:
        _,plan=load_case(args.source,index)
        selected=proposals(plan.problem,args.budget);selections[index]=selected
        shared={}
        for policy,pairs in selected.items():
            for rank,pair in enumerate(pairs):
                kind='cut' if policy=='random_single_cut' else 'anchor'
                key=(kind,tuple(pair))
                if key not in shared:
                    shared[key]=dict(index=index,kind=kind,pair=pair,policies={})
                shared[key]['policies'][policy]=rank+1
        tasks.append(dict(index=index,kind='baseline',pair=None,policies={}))
        tasks.extend(shared.values())
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    save(args.run/'manifest.json',dict(source=str(args.source.resolve()),budget=args.budget,
        indices=args.indices,selections=selections,tasks=tasks,seeds=10,cap=100,
        source_hashes={str(p.relative_to(args.run/'engine')):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (args.run/'engine').rglob('*.py')},
        note='Independent single anchors versus independent single cuts; no joint anchors or reference-guided selection. Each policy receives the same maximum number of tasks. Shared proposals are computed once.'))
    print(len(tasks))


def search(args):
    from rxn_core.aam import _initialize_search,_search_cut
    from rxn_core.artifacts import write_aam_checkpoint
    from rxn_core.domain import AAMResult,AAMSearchMetrics
    from rxn_core.frag import build_graph
    from rxn_core.search_symmetry import finalize_graph_symmetry
    m=json.loads((args.run/'manifest.json').read_text());task=m['tasks'][args.slot]
    _,plan=load_case(Path(m['source']),task['index']);out=args.run/str(args.slot)
    anchors=(tuple(task['pair']),) if task['kind']=='anchor' else ()
    cuts=(tuple(task['pair']),) if task['kind']=='cut' else ()
    config=replace(plan.config,seed_count=m['seeds'],branch_limit=m['cap'],anchors=anchors)
    start=time.perf_counter();_initialize_search(plan.problem,config)
    graph,counts=_search_cut(cuts)
    search_seconds=time.perf_counter()-start
    graph,_=finalize_graph_symmetry(graph,build_graph(plan.problem.product.elements,
        plan.problem.product.wbo,config.graph_floor),iso_tolerance=config.iso_tolerance)
    assert all(all(dict(graph.states[t].mapping).get(r)==p for r,p in anchors) for t in graph.terminals)
    result=AAMResult(plan.problem,config,graph,AAMSearchMetrics.from_record(counts,time.perf_counter()-start))
    write_aam_checkpoint(result,out/'aam.pkl.gz')
    save(out/'search.json',dict(search_seconds=search_seconds,total_seconds=time.perf_counter()-start,
        terminals=len(graph.terminals),capped=graph.capped,anchors_preserved=True))


def score(args):
    from rxn_core.artifacts import read_aam_checkpoint
    from golden_evaluation import evaluate_planned
    m=json.loads((args.run/'manifest.json').read_text());task=m['tasks'][args.slot]
    directory,plan=load_case(Path(m['source']),task['index']);out=args.run/str(args.slot)
    result=read_aam_checkpoint(out/'aam.pkl.gz');plan=replace(plan,config=result.config)
    reference=json.loads((directory/'reference.json').read_text())
    evaluation=evaluate_planned(result,plan,reference['features'],reference['mapping'],seconds=30,query_timeout_ms=2000)
    save(out/'evaluation.json',evaluation)


def worker(args):
    out=args.run/str(args.slot);out.mkdir()
    statuses={}
    for phase,limit in [('search',90),('score',45)]:
        statuses[phase]=guarded([sys.executable,__file__,phase,'--run',str(args.run),'--slot',str(args.slot)],
            out/f'{phase}.log',limit)
        save(out/'status.json',statuses)
        if phase=='search' and statuses[phase]!=0:break
    if not (out/'evaluation.json').exists():
        save(out/'evaluation.json',dict(reference_recovery='unknown',reason=statuses))


def report(args):
    m=json.loads((args.run/'manifest.json').read_text());rows=[]
    for slot,task in enumerate(m['tasks']):
        out=args.run/str(slot)
        def read(name):
            path=out/name
            return json.loads(path.read_text()) if path.exists() else {}
        rows.append(dict(slot=slot,task=task,evaluation=read('evaluation.json'),search=read('search.json'),status=read('status.json')))
    summary={}
    for index in m['indices']:
        summary[index]={}
        for policy in ('baseline','random_anchor','stable_anchor','change_anchor','random_single_cut'):
            selected=[r for r in rows if r['task']['index']==index and
                      (r['task']['kind']=='baseline' if policy=='baseline' else policy in r['task']['policies'])]
            outcomes=Counter(r['evaluation'].get('reference_recovery','pending') for r in selected)
            summary[index][policy]=dict(outcomes=dict(outcomes),
                recovered_slots=[r['slot'] for r in selected if r['evaluation'].get('reference_recovery')=='recovered'],
                completed_searches=sum(bool(r['search']) for r in selected),
                search_seconds=sum(r['search'].get('search_seconds',0) for r in selected))
    args.output.mkdir(parents=True,exist_ok=True);save(args.output/'summary.json',summary)
    save(args.output/'cases.json',rows);save(args.output/'manifest.json',m)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['init','worker','search','score','report'])
    p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path)
    p.add_argument('--indices',type=int,nargs='+',default=[986,1285,1553,1739,300,600])
    p.add_argument('--budget',type=int,default=16);p.add_argument('--slot',type=int);p.add_argument('--output',type=Path)
    args=p.parse_args();globals()[args.mode](args)
