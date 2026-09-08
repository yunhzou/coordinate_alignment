"""Bounded scoring experiment on saved Golden paths; never reruns AAM."""
import argparse
import json
from pathlib import Path
import resource
import random
from array import array
import time

from golden_publication import plans
from golden_policy_campaign import save
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.search_graph import SearchPath
from rxn_core.family_scoring import minimize_events,bond_events

CASES=((0,'P_to_R'),(44,'R_to_P'),(44,'P_to_R'),(865,'P_to_R'),
       (1033,'P_to_R'),(1786,'P_to_R'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--slot',type=int,required=True)
    p.add_argument('--blind',action='store_true')
    a=p.parse_args();index,direction=CASES[a.slot]
    out=a.output/f'{index}_{direction}.json'
    pair,_=plans(a.source,index);plan=pair[direction]
    base=a.source/'directions'/str(index)/direction
    evaluation=json.loads((base/'evaluation.json').read_text())
    start=time.perf_counter();aam=read_aam_checkpoint(base/'cuts/aam.pkl.gz')
    record=dict(index=index,direction=direction,archive_loading_seconds=time.perf_counter()-start,
                explicit_atoms=[aam.problem.source_atom_count,aam.problem.target_atom_count],results=[])
    save(out,record)
    if a.blind:
        # Reference-blind selection: ranked representatives plus reproducible
        # random terminals. Select one recorded path per terminal, not orbits.
        ranked=json.loads((base/'classes.json').read_text())
        terminals=list(dict.fromkeys([c['terminal'] for c in ranked[:8]]+
            random.Random(42).sample(list(aam.graph.terminals),min(24,len(aam.graph.terminals)))))
        started=time.perf_counter();incoming=array('q',[-1])*len(aam.graph.states)
        for edge in aam.graph.transitions:
            if incoming[edge.target]<0:incoming[edge.target]=edge.id
        record['path_index_seconds']=time.perf_counter()-started
        record['selection']='first 8 reference-blind ranked classes plus 24 random terminals, random seed 42'
        record['selected_terminals']=terminals;save(out,record)
        for terminal in terminals:
            transitions=[];state=terminal
            while incoming[state]>=0:
                edge=aam.graph.transitions[incoming[state]];transitions.append(edge.id);state=edge.source
            path=SearchPath(aam.graph,terminal,tuple(reversed(transitions)))
            result=minimize_events(path,aam.problem,reverse=plan.reversed,seconds=3.)
            result['terminal']=terminal;record['results'].append(result)
            record['peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            save(out,record)
            print(json.dumps({k:v for k,v in result.items() if k!='mapping'}),flush=True)
        raise SystemExit(0)
    terminal=evaluation['witness_terminal']
    start=time.perf_counter()
    path=(SearchPath(aam.graph,terminal,tuple(evaluation['witness_path']))
          if evaluation.get('witness_path') else next(aam.graph.paths(terminal)))
    record['path_access_seconds']=time.perf_counter()-start
    if evaluation.get('input_orientation_witness'):
        witness=plan.to_search_mapping(dict(evaluation['input_orientation_witness']))
        record['saved_reference_witness_events']=bond_events(aam.problem,witness,reverse=plan.reversed)
    save(out,record)
    for repeat in range(3):
        result=minimize_events(path,aam.problem,reverse=plan.reversed,seconds=10.)
        result['repeat']=repeat
        record['results'].append(result)
        record['peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        save(out,record)
        print(json.dumps({k:v for k,v in result.items() if k!='mapping'}),flush=True)
