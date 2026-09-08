"""Extend representative bond-event windows from saved metrics, without search."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path


def summarize(rows):
    result={'total':len(rows),'modes':{},'semantics':
        'Post-hoc windows on saved representatives in the best heavy/all-atom coverage tier. '
        'Not verification of every compressed symmetry alternative. No search rerun.'}
    for name in ('single','bidirectional'):
        histogram=Counter();excluded=defaultdict(list)
        for row in rows:
            mode=row['modes'][name];rep=mode.get('representative',{})
            gap=rep.get('reference_event_gap')
            if gap is not None:
                histogram[gap]+=1
            else:
                reason=('incomplete_ranking' if not mode['ranking_complete'] else
                        'outside_best_coverage_tier' if rep.get('reference_class_rank') is not None else
                        'no_matching_representative:'+mode['reference_recovery'])
                excluded[reason].append(row['index'])
        windows={str(n):sum(v for gap,v in histogram.items() if gap<=n)
                 for n in (*range(11),15,20)}
        for n in (0,1,2,3,5):
            assert windows[str(n)]==sum(row['modes'][name].get('representative',{}).
                get('event_windows',{}).get(str(n),{}).get('recovered',False) for row in rows)
        windows['unlimited']=sum(histogram.values())
        result['modes'][name]={'windows':{n:{'recovered':v,'percent':100*v/len(rows)}
                                        for n,v in windows.items()},
            'gap_histogram':dict(sorted(histogram.items())), 'excluded_cases':dict(excluded)}
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    a=p.parse_args()
    result=summarize(json.loads((a.run/'cases.json').read_text()))
    out=a.run/'publication/event_windows_extended.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
