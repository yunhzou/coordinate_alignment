"""Locate one unrecovered historical family and isolate its cut/order effects."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from adaptive_fragment_pilot import save
from compare_elementary_outputs import features,certificate,event_counts
from rxn_core.artifacts import read_graph_checkpoint,write_graph_checkpoint
from rxn_core.frag import build_graph
from rxn_core.matcher import _nauty_orbits
from rxn_core.native_search import find_islands_native
from rxn_core.search_symmetry import finalize_graph_symmetry


def main(args):
    folder=args.run/f'results/{args.slot}'
    summary=json.loads((folder/'summary.json').read_text())
    families=json.loads((folder/'pattern_family_queries.json').read_text())
    refs=json.loads((folder/'pattern_witnesses.json').read_text())['baseline']
    key=families['unresolved'][0];reference=refs[key]
    raw=json.loads((args.run/f"inputs/{summary['index']}.json").read_text())
    if summary['direction']=='P_to_R':raw['reactant'],raw['product']=raw['product'],raw['reactant']
    feat=features(raw);heavy=[i for i,e in enumerate(raw['reactant']['elements']) if e!='H']
    baseline=Path(json.loads((folder/'pattern_comparison.json').read_text())['baseline'])
    saved=None
    for seed in sorted(baseline.glob('seed_*')):
        witness=json.loads((seed/'witnesses.json').read_text())
        needle=reference['mapping']
        if summary['direction']=='P_to_R':needle=np.argsort(needle).tolist()
        if needle not in witness['mappings']:continue
        cuts=json.loads((seed/'summary.json').read_text())['cuts']
        for ordinal,cut in enumerate(cuts):
            if cut['best']>families['event_limit']:continue
            path=seed/f'cut_{ordinal:04d}.pkl.gz';g=read_graph_checkpoint(path)
            terminal=next((t for t in g.terminals if dict(g.states[t].mapping)==dict(enumerate(reference['mapping']))),None)
            if terminal is not None:saved=(path,g,terminal);break
        if saved is not None:break
    assert saved is not None
    path,g,terminal=saved
    left,right=raw['reactant'],raw['product'];r,p=[np.asarray(v['wbo']) for v in (left,right)]
    source=build_graph(left['elements'],r,bond_cut=.2);target=build_graph(right['elements'],p,bond_cut=.2)
    po=_nauty_orbits(target,wbo_tol=1.)
    start=time.process_time()
    uncut=find_islands_native(source,target,g.contexts[0].seed_order,iso_tol=1.,max_branches=100,p_orbits=po)
    uncut,_=finalize_graph_symmetry(uncut,target,iso_tolerance=1.)
    compute=time.process_time()-start
    write_graph_checkpoint(uncut,folder/'missing_pattern_uncut.pkl.gz')
    vectors=[list(dict(uncut.states[t].mapping).values()) for t in uncut.terminals if len(uncut.states[t].mapping)==len(r)]
    scores=event_counts(r,p,vectors).sum(axis=1)
    recovered=any(s<=families['event_limit'] and hashlib.sha256(certificate(feat,v,heavy)).hexdigest()==key
                  for v,s in zip(vectors,scores))
    selected=next(g.paths(terminal))
    fragments=[g.transitions[e].match for e in selected.transitions if g.transitions[e].match]
    mapping=dict(g.states[terminal].mapping)
    cut_effect=[dict(atoms=(a,b),source_wbo=float(r[a,b]),target_wbo=float(p[mapping[a],mapping[b]]),
        together_in_fragment=any(a in f['fragment'] and b in f['fragment'] for f in fragments))
        for a,b in g.contexts[0].cuts]
    save(folder/'missing_pattern_trace.json',dict(reference_key=key,source=str(path),terminal=terminal,
        seed_order=g.contexts[0].seed_order,cut_effect=cut_effect,fragments=fragments,
        uncut_best=int(min(scores)),uncut_capped=uncut.capped,uncut_cpu=compute,
        uncut_representative_recovers_same_pattern=recovered,
        scope='Historical-reference-guided diagnosis only, never a blind search policy. Negative representative membership does not exclude compressed alternatives.'))
    print(json.dumps(dict(source=str(path),cut_effect=cut_effect,uncut_best=int(min(scores)),
        uncut_recovers_pattern=recovered,cpu=compute)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--slot',type=int,required=True)
    main(p.parse_args())
