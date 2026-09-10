"""Batched historical-pattern membership in saved compressed families; no search."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import pynauty
import z3

from adaptive_fragment_pilot import save
from compare_elementary_outputs import features, certificate
from golden_evaluation import colored_graph
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.family_query import compile_path, SymbolicActions
from rxn_core.family_scoring import event_objective, bond_events
from rxn_core.search_graph import frozen_value


def main(args):
    folder=args.run/f'results/{args.slot}'
    summary=json.loads((folder/'summary.json').read_text())
    patterns=json.loads((folder/'pattern_witnesses.json').read_text())
    comparison=json.loads((folder/'pattern_comparison.json').read_text())['rows'][0]
    limit=comparison['event_limit']
    expected={k:v for k,v in patterns['baseline'].items() if v['events']<=limit}
    found={k:dict(method='saved_representative',**patterns['adaptive'][k])
           for k in expected if k in patterns['adaptive'] and patterns['adaptive'][k]['events']<=limit}
    result=read_aam_checkpoint(folder/f"{summary['rows'][-1]['label']}.pkl.gz")
    raw=json.loads((args.run/f"inputs/{summary['index']}.json").read_text())
    if summary['direction']=='P_to_R':raw['reactant'],raw['product']=raw['product'],raw['reactant']
    feat=features(raw)
    heavy=[i for i,e in enumerate(raw['reactant']['elements']) if e!='H']
    generators=[tuple(tuple(g[:len(f['colors'])]) for g in pynauty.autgrp(colored_graph([f]))[0]) for f in feat]
    scored=json.loads((folder/f"{summary['rows'][-1]['label']}_witnesses.json").read_text())
    scores={tuple(v):sum(e) for v,e in zip(scored['mappings'],scored['events'],strict=True)}
    terminals=sorted(result.graph.terminals,key=lambda t:scores.get(tuple(dict(result.graph.states[t].mapping).values()),float('inf')))
    started=time.perf_counter();deadline=started+args.seconds
    seen=set();queries=unknown=encoded=0;exhausted=True
    def checkpoint():
        save(folder/'pattern_family_queries.json',dict(event_limit=limit,baseline_patterns=len(expected),
            recovered_patterns=len(found),unresolved=sorted(set(expected)-set(found)),witnesses=found,
            encoded_paths=encoded,queries=queries,unknown_queries=unknown,exhausted=exhausted,
            seconds=time.perf_counter()-started,
            scope='Positive full-H witnesses with required event bound, checked modulo exact score-response endpoint symmetry. A time-limited unresolved class is not proven absent. No AAM rerun or permutation enumeration.'))
    for terminal in terminals:
        if len(found)==len(expected):break
        for path in result.graph.paths(terminal):
            if time.perf_counter()>=deadline:exhausted=False;break
            key=(tuple(sorted(path.mapping.items())),frozen_value(path.fragments))
            if key in seen:continue
            seen.add(key)
            compiled=compile_path(path,result.problem,{},source_atoms=(),complete_reference=False)
            objective,lower,_=event_objective(compiled)
            encoded += 1
            if lower>limit:continue
            solver=compiled.solver;solver.add(objective<=limit)
            encoder=SymbolicActions(solver,fresh=True)
            sources=encoder.act([encoder.constant(r) for r in heavy],generators[0],'source_equivalence')
            normalized=[encoder.lookup(v,dict(enumerate(compiled.values))) for v in sources]
            normalized=encoder.act(normalized,generators[1],'target_equivalence')
            for key,reference in expected.items():
                if key in found:continue
                if time.perf_counter()>=deadline:exhausted=False;break
                vector=reference['mapping']
                if any(vector[r] not in value[1] for r,value in zip(heavy,normalized)):continue
                solver.push();solver.add(*(value[0]==vector[r] for r,value in zip(heavy,normalized)))
                solver.set(timeout=max(1,min(1000,int(1000*(deadline-time.perf_counter())))))
                status=solver.check();queries += 1
                if status==z3.sat:
                    witness=compiled.realize(solver.model())
                    mapping=dict(witness['mapping']);events=bond_events(result.problem,mapping)
                    actual=[mapping[i] for i in range(result.problem.atom_count)]
                    assert events['total']<=limit
                    assert hashlib.sha256(certificate(feat,actual,heavy)).hexdigest()==key
                    found[key]=dict(method='compressed_family_query',events=events['total'],mapping=actual,
                        terminal=terminal,transitions=path.transitions,checked_fragment_edges=witness['checked_fragment_edges'])
                elif status==z3.unknown:unknown += 1
                solver.pop()
            checkpoint()
            if not exhausted or len(found)==len(expected):break
        if not exhausted:break
    checkpoint()
    print(json.dumps(dict(slot=args.slot,recovered=len(found),expected=len(expected),encoded=encoded,
        queries=queries,unknown=unknown,exhausted=exhausted,seconds=time.perf_counter()-started)),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--slot',type=int,required=True)
    parser.add_argument('--seconds',type=float,default=60.)
    main(parser.parse_args())
