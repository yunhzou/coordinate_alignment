"""Compare saved native mappings on shared WBOs, without new mapping searches."""
import argparse
from collections import defaultdict, Counter
import gzip
import hashlib
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
import time

import numpy as np
import pynauty
from golden_competitors import save
from golden_evaluation import colored_graph


def event_counts(r,p,mappings):
    mappings=np.asarray(mappings,dtype=int)
    ra,rb=np.where(np.triu(r,1)>.2);pa,pb=np.where(np.triu(p,1)>.2)
    matched=p[mappings[:,ra],mappings[:,rb]]
    broken=(matched<=.2).sum(axis=1)
    changed=((matched>.2)&(np.abs(matched-r[ra,rb])>.5)).sum(axis=1)
    inverse=np.argsort(mappings,axis=1)
    formed=(r[inverse[:,pa],inverse[:,pb]]<=.2).sum(axis=1)
    return np.stack((broken,formed,changed),axis=1)


def features(raw):
    endpoints=[raw[k] for k in ('reactant','product')];result=[]
    # Exact score-response colors: no rounding of WBOs or arbitrary binary
    # equivalence. Allowed endpoint automorphisms preserve event classification.
    for side,endpoint in enumerate(endpoints):
        w=np.asarray(endpoint['wbo']);other=np.asarray(endpoints[1-side]['wbo'])
        values=np.unique(other[np.triu(other,1)>.2])
        bonds=[]
        for a,b in zip(*np.where(np.triu(w,1)>.2)):
            bonds.append((int(a),int(b),tuple(map(int,np.abs(w[a,b]-values)>.5))))
        result.append(dict(colors=[(e,) for e in endpoint['elements']],bonds=bonds))
    return result


def certificate(feat,mapping,heavy):
    relation={i:int(mapping[i]) for i in heavy}
    return pynauty.certificate(colored_graph(feat,relation))


def events_row(counts):
    a,b,c=map(int,counts)
    return dict(broken=a,formed=b,order_changed=c,total=a+b+c)


def compare(args):
    start=time.perf_counter();source=args.source
    raw=json.loads((source/'inputs'/str(args.index)/'input.json').read_text())
    r,p=[np.asarray(raw[k]['wbo']) for k in ('reactant','product')]
    heavy=[i for i,e in enumerate(raw['reactant']['elements']) if e!='H']
    feat=features(raw)
    native=json.loads((source/'slap_xyz'/f'{args.index}.json').read_text())['candidates']
    slap=[];slap_keys=set();slap_certificates=set()
    for ordinal,candidate in enumerate(native):
        groups=[]
        for graph in candidate['graphs']:
            partition=defaultdict(list)
            for i,label in enumerate(graph['labels']):partition[label].append(i)
            groups.append(partition)
        mapping=dict(pair for label,left in groups[0].items() for pair in zip(left,groups[1][label]))
        assert len(mapping)==len(r) and len(set(mapping.values()))==len(p)
        assert all(raw['reactant']['elements'][a]==raw['product']['elements'][b] for a,b in mapping.items())
        assert all(sum(raw['reactant']['elements'][i]!='H' for i in group)<=1 for group in groups[0].values())
        vector=[mapping[i] for i in range(len(r))]
        key=tuple(vector[i] for i in heavy);slap_keys.add(key)
        cert=certificate(feat,vector,heavy);slap_certificates.add(cert)
        # All within-label alternatives arise through target permutations.
        # Certify event invariance using the score-response colored target graph.
        labels={(a,b):c for a,b,c in feat[1]['bonds']}
        def edge(a,b):return labels.get(tuple(sorted((a,b))))
        invariant=True
        for group in groups[1].values():
            for a,b in zip(group,group[1:]):
                perm=list(range(len(p)));perm[a],perm[b]=b,a
                invariant &= all(edge(i,j)==edge(perm[i],perm[j]) for i in (a,b) for j in range(len(p)) if i!=j)
        slap.append(dict(candidate=ordinal,mapping=sorted(mapping.items()),
            events=events_row(event_counts(r,p,[vector])[0]),
            all_h_label_permutations_score_invariant=bool(invariant),
            heavy_pattern_certificate=hashlib.sha256(cert).hexdigest()))
    best=10**9;minimum={};seen_slap_keys=set();witnesses=0
    run=Path(args.cap200) if args.index==123 else source
    for direction in ('R_to_P','P_to_R'):
        path=run/'directions'/str(args.index)/direction/'terminal_mappings.jsonl.gz'
        with gzip.open(path,'rt') as stream:
            exhausted=False
            while not exhausted:
                batch=[];metadata=[]
                for _ in range(2048):
                    line=stream.readline()
                    if not line:exhausted=True;break
                    row=json.loads(line);pairs=dict(row['mapping'])
                    if not row['full_element_bijection']:continue
                    batch.append([pairs[i] for i in range(len(r))]);metadata.append(row['terminal'])
                if not batch:continue
                counts=event_counts(r,p,batch);totals=counts.sum(axis=1)
                for i,vector in enumerate(batch):
                    witnesses+=1;key=tuple(vector[j] for j in heavy)
                    if key in slap_keys:seen_slap_keys.add(key)
                    score=int(totals[i])
                    if score<best:best=score;minimum.clear()
                    if score==best and key not in minimum:
                        minimum[key]=dict(mapping=list(enumerate(vector)),events=events_row(counts[i]),
                            terminal=metadata[i],direction=direction)
    best_classes=set();matching_best=None;representative=None
    for entry in minimum.values():
        vector=[p for r,p in entry['mapping']];cert=certificate(feat,vector,heavy)
        best_classes.add(cert)
        if representative is None:representative=entry
        if cert in slap_certificates:matching_best=entry
    result=dict(index=args.index,name=raw['name'],atoms=len(r),
        aam_source=str(run),aam_witnesses_scored=witnesses,
        aam_min_saved_events=best,aam_min_distinct_indexed_heavy_maps=len(minimum),
        aam_min_symmetry_classes=len(best_classes),aam_best_example=representative,
        slap=slap,slap_min_representative_events=min(x['events']['total'] for x in slap),
        any_slap_heavy_map_exactly_seen_in_aam_witnesses=bool(seen_slap_keys),
        best_aam_heavy_pattern_matches_any_slap_modulo_score_preserving_symmetry=matching_best is not None,
        shared_best_example=matching_best,comparison_seconds=time.perf_counter()-start,
        scope='Saved terminal representatives only, not exhaustive AAM family extraction. '
              'Case 123 uses the separate cap-200 follow-up. Symmetry comparison projects mappings to heavy atoms, '
              'but endpoint graphs and event scores include explicit H. No stereo/energy/reference correctness claim.')
    args.run.mkdir(parents=True,exist_ok=True)
    save(args.run/f'{args.index}.json',result)


def report(args):
    rows=[json.loads((args.run/f'{i}.json').read_text()) for i in range(140)]
    comparison=Counter('aam_fewer' if r['aam_min_saved_events']<r['slap_min_representative_events'] else
        'slap_fewer' if r['aam_min_saved_events']>r['slap_min_representative_events'] else 'equal' for r in rows)
    summary=dict(cases=len(rows),event_comparison=dict(comparison),
        best_aam_heavy_pattern_matches_any_slap=sum(r['best_aam_heavy_pattern_matches_any_slap_modulo_score_preserving_symmetry'] for r in rows),
        exact_slap_heavy_map_seen_anywhere_in_aam=sum(r['any_slap_heavy_map_exactly_seen_in_aam_witnesses'] for r in rows),
        slap_candidates_with_noninvariant_h_scores=sum(not s['all_h_label_permutations_score_invariant'] for r in rows for s in r['slap']),
        analysis_seconds_sum=sum(r['comparison_seconds'] for r in rows),accuracy=None)
    save(args.run/'summary.json',summary);print(json.dumps(summary,indent=2))


def submit(args):
    args.run.mkdir(parents=True,exist_ok=True)
    shutil.copy2(__file__,args.run/'compare_driver.py')
    indices=[i for i in range(140) if not (args.run/f'{i}.json').exists()]
    command=['env','OPENBLAS_NUM_THREADS=1','OMP_NUM_THREADS=1',
        f'PYTHONPATH={Path(__file__).resolve().parent}',
        'timeout','--kill-after=5s','300',sys.executable,str(args.run/'compare_driver.py'),
        'compare','--source',str(args.source),'--cap200',str(args.cap200),'--run',str(args.run),'--index']
    options=['sbatch','--parsable','--partition=cpunodes','--exclude=bosque8',
        '--cpus-per-task=1','--mem=8G','--time=00:10:00',
        '--array='+','.join(map(str,indices))+'%16','--job-name=elem_compare',
        f'--output={args.run}/slurm_%A_%a.out','--wrap',
        shlex.join(command)+' "$SLURM_ARRAY_TASK_ID"']
    job=subprocess.check_output(options,text=True).strip()
    save(args.run/'submission.json',dict(job=job,command=options,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    print(job)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('compare','report','submit'))
    parser.add_argument('--source',type=Path,default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_feasibility_20260908'))
    parser.add_argument('--cap200',type=Path,default=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_case123_cap200_20260908'))
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--index',type=int)
    args=parser.parse_args();globals()[args.command](args)
