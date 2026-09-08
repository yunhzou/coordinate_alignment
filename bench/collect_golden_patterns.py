"""Resumable, reference-blind pattern extraction from an unchanged archive."""
import hashlib
import json
import time
from itertools import chain

from golden_policy_campaign import save
from rxn_core.pattern_collection import PatternEquivalence,extract_path_patterns
from rxn_core.search_graph import SearchPath


def first_paths(graph):
    """One shared incoming index, rather than rebuilding it for each terminal."""
    incoming={}
    for edge in graph.transitions:incoming.setdefault(edge.target,edge)
    roots=set(graph.roots)
    def path(terminal):
        state=terminal;edges=[]
        while state not in roots:
            edge=incoming[state];edges.append(edge.id);state=edge.source
        return SearchPath(graph,terminal,tuple(reversed(edges)))
    return path


def collect(archive,aam,plan,mols,priority,out,seconds,per_path_seconds):
    atom_tags=tuple({a.GetIdx():(a.GetFormalCharge(),a.GetIsotope(),
        a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else '') for a in mol.GetAtoms()} for mol in mols)
    bond_tags=tuple({tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))):str(b.GetStereo())
        for b in mol.GetBonds()} for mol in mols)
    if plan.reversed:atom_tags=atom_tags[::-1];bond_tags=bond_tags[::-1]
    eq=PatternEquivalence(plan.problem,atom_tags=atom_tags,bond_tags=bond_tags)
    digest=hashlib.sha256()
    with archive.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    identity=dict(archive_sha256=digest.hexdigest(),direction=plan.direction,
        pattern_version=1,event_tolerance=.5,endpoint_tags=repr((atom_tags,bond_tags)))
    state=json.loads(out.read_text()) if out.exists() else dict(identity=identity,paths={},complete=False)
    if state['identity']!=identity:raise ValueError('Pattern checkpoint belongs to different input/settings')
    first=first_paths(aam.graph);start=time.perf_counter();all_complete=True
    visited=set()
    for path in chain((first(t) for t in priority),aam.graph.paths()):
        key=','.join(map(str,path.transitions))
        if key in visited:continue
        visited.add(key)
        previous=state['paths'].get(key,{})
        if previous.get('complete'):continue
        remaining=seconds-(time.perf_counter()-start)
        if remaining<=0:break
        def checkpoint(record):
            item=state['paths'].setdefault(key,dict(patterns=[],complete=False))
            if not any(p['key']==record['key'] for p in item['patterns']):item['patterns'].append(record)
            save(out,state)
        result=extract_path_patterns(path,plan.problem,eq,previous=previous.get('patterns',()),
            seconds=min(remaining,per_path_seconds),reverse=plan.reversed,on_pattern=checkpoint)
        state['paths'][key]=result
        save(out,state)
        if not result['complete']:all_complete=False
    else:state['complete']=all_complete
    if 'last_extraction_seconds' in state:
        state.setdefault('previous_pass_seconds',[]).append(state['last_extraction_seconds'])
    state['last_extraction_seconds']=time.perf_counter()-start
    save(out,state)
    patterns={}
    for result in state['paths'].values():
        for record in result['patterns']:
            item=patterns.setdefault(record['key'],dict(record,origins=[]))
            item['origins'].append(dict(terminal=record['terminal'],transitions=record['transitions'],
                actions=record['actions'],mapping=record['mapping'],orientation=plan.direction))
    return state,list(patterns.values()),eq
