"""Explain cached ranking failures without changing candidates or scoring."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

from golden_policy_campaign import save


def diagnose(rows):
    counts=Counter();gaps=Counter();cases=[]
    for case in rows:
        top=case['classes'][0]
        reference=next((c for c in case['classes'] if c['reference_equivalent']),None)
        counts['completed_cases']+=1
        counts['top_one_center_component']+=top['features']['center_components']==1
        if top['reference_equivalent']:
            counts['top1_correct']+=1;continue
        counts['top1_wrong']+=1
        row=dict(index=case['index'],top=top,reference=reference)
        cases.append(row)
        if reference is None:
            counts['no_reference_representative']+=1;row['category']='no_reference_representative';continue
        counts['wrong_with_reference_representative']+=1
        if reference['score'][:2]!=top['score'][:2]:
            counts['reference_lower_coverage']+=1;row['category']='reference_lower_coverage';continue
        gap=reference['score'][2]-top['score'][2];gaps[gap]+=1
        row['category']='event_tie' if gap==0 else 'reference_more_events';row['event_gap']=gap
        a,b=top['features'],reference['features']
        row['identical_cached_features']=a==b
        counts['same_all_cached_features']+=a==b
        counts['equal_center_components']+=a['center_components']==b['center_components']
        if a['weighted_edits'] is not None and b['weighted_edits'] is not None:
            relation='lower' if b['weighted_edits']<a['weighted_edits']-1e-9 else 'higher' if b['weighted_edits']>a['weighted_edits']+1e-9 else 'equal'
            counts['reference_energy_'+relation]+=1;row['reference_energy_relation']=relation
    return dict(counts=dict(counts),reference_event_gaps=dict(sorted(gaps.items()))),cases


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--features',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    rows=json.loads(gzip.decompress(args.features.read_bytes()));summary,cases=diagnose(rows)
    args.output.mkdir(parents=True,exist_ok=True);save(args.output/'summary.json',summary);save(args.output/'cases.json',cases)
    print(json.dumps(summary,indent=2))
