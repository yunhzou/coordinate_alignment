"""Compare saved low-event heavy-atom patterns; no new matching or permutations."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from adaptive_fragment_pilot import save
from compare_elementary_outputs import features, certificate
from score_equivalence import element_pair_features


def main(args):
    started=time.perf_counter()
    folder=args.run/f'results/{args.slot}'
    summary=json.loads((folder/'summary.json').read_text())
    baseline=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/fragment_native_final_20260909_JGOrLU/results')
    baseline=baseline/str(summary['index'])/summary['direction']/'independent_native_dependency'
    raw=json.loads((args.run/f"inputs/{summary['index']}.json").read_text())
    if summary['direction']=='P_to_R':raw['reactant'],raw['product']=raw['product'],raw['reactant']
    feat=element_pair_features(raw) if args.element_pair else features(raw)
    prefix='pair_' if args.element_pair else ''
    if args.checkpoint_label:prefix+=args.checkpoint_label+'_'
    heavy=[i for i,e in enumerate(raw['reactant']['elements']) if e!='H']
    previous=[json.loads(p.read_text()) for p in sorted(baseline.glob('seed_*/witnesses.json'))]
    # Historical witness exports are normalized to R->P even when their search
    # ran P->R. New checkpoint witnesses retain the search direction.
    if summary['direction']=='P_to_R':
        for data in previous:data['mappings']=[np.argsort(v).tolist() for v in data['mappings']]
    best=min(sum(events) for data in previous for events in data['events'])
    label=args.checkpoint_label or summary['rows'][-1]['label']
    current=[json.loads((folder/f'{label}_witnesses.json').read_text())]
    classes={}
    for name,datasets in (('baseline',previous),('adaptive',current)):
        representatives={}
        for data in datasets:
            for vector,events in zip(data['mappings'],data['events'],strict=True):
                score=sum(events)
                if score>best+2:continue
                assert len(set(vector))==len(vector)
                assert all(raw['reactant']['elements'][i]==raw['product']['elements'][j]
                           for i,j in enumerate(vector))
                relation=tuple(vector[a] for a in heavy)
                old=representatives.get(relation)
                if old is None or score<old[0]:representatives[relation]=(score,vector)
        found={}
        for score,vector in representatives.values():
            key=hashlib.sha256(certificate(feat,vector,heavy)).hexdigest()
            if key not in found or score<found[key]['events']:
                found[key]=dict(events=score,mapping=vector)
        classes[name]=found
    rows=[]
    for delta in range(3):
        expected={k for k,v in classes['baseline'].items() if v['events']<=best+delta}
        observed={k for k,v in classes['adaptive'].items() if v['events']<=best+delta}
        rows.append(dict(delta=delta,event_limit=best+delta,baseline_patterns=len(expected),
            recovered_patterns=len(expected & observed),adaptive_patterns=len(observed),
            missing=sorted(expected-observed),new=sorted(observed-expected)))
    save(folder/f'{prefix}pattern_comparison.json',dict(rows=rows,baseline=str(baseline),
        equivalence='element_pair_score_response' if args.element_pair else 'all_bond_weights_score_response',
        seconds=time.perf_counter()-started,
        scope='Saved representative heavy-atom relations modulo exact event-response endpoint graph automorphisms. Full-H event counts. Missing representative class does not prove absence from compressed families. No independent chemical ground truth.'))
    save(folder/f'{prefix}pattern_witnesses.json',classes)
    print(json.dumps([{k:v for k,v in row.items() if k not in ('missing','new')} for row in rows]),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--slot',type=int,required=True)
    p.add_argument('--element-pair',action='store_true')
    p.add_argument('--checkpoint-label')
    main(p.parse_args())
