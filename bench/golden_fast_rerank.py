"""Cheap, exploratory reranking of cached Golden class representatives."""
import argparse
from collections import Counter
import hashlib
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from golden_policy_campaign import load_case,save,guarded
from rxn_core.artifacts import read_aam_checkpoint

# Common average bond energies from CDK 2.8 BondEnergies documentation, kJ/mol.
# These are generic bond-type costs, not molecule-specific BDEs or barriers.
ENERGY_SOURCE='https://cdk.github.io/cdk/2.8/docs/api/org/openscience/cdk/smsd/tools/BondEnergies.html'
ENERGIES={
    'H-H':(432,), 'B-H':(389,), 'C-H':(411,), 'H-Si':(318,), 'H-N':(386,),
    'H-P':(322,), 'H-O':(459,), 'H-S':(363,), 'F-H':(565,), 'Cl-H':(428,),
    'Br-H':(362,), 'H-I':(295,), 'B-B':(293,), 'B-O':(536,), 'B-F':(613,),
    'B-Cl':(456,), 'B-Br':(377,), 'C-C':(346,602,835), 'C-Si':(318,),
    'C-N':(305,615,887), 'C-P':(264,), 'C-O':(358,799,1072), 'B-C':(356,),
    'C-S':(272,573), 'C-F':(485,), 'C-Cl':(327,), 'Br-C':(285,), 'C-I':(213,),
    'Si-Si':(222,), 'N-Si':(355,), 'O-Si':(452,), 'S-Si':(293,),
    'F-Si':(565,), 'Cl-Si':(381,), 'Br-Si':(310,), 'I-Si':(234,),
    'N-N':(167,418,942), 'N-O':(201,607), 'F-N':(283,), 'Cl-N':(313,),
    'P-P':(201,), 'O-P':(335,544), 'P-S':(None,335), 'F-P':(490,),
    'Cl-P':(326,), 'Br-P':(264,), 'I-P':(184,), 'O-O':(142,494),
    'F-O':(190,), 'O-S':(None,522), 'S-S':(226,425), 'F-S':(284,),
    'Cl-S':(255,), 'F-F':(155,), 'Cl-Cl':(240,), 'Br-Br':(190,), 'I-I':(148,),
}


def bond_energy(elements,order):
    if order<=.2:return 0.
    pair='-'.join(sorted(elements));levels=ENERGIES.get(pair,())
    lo=int(np.floor(order));hi=int(np.ceil(order))
    if lo<1 or hi>len(levels) or levels[lo-1] is None or levels[hi-1] is None:
        raise KeyError((pair,float(order)))
    # Explicit aromatic interpolation convention; not an aromatic BDE claim.
    return levels[lo-1] if lo==hi else levels[lo-1]+(order-lo)*(levels[hi-1]-levels[lo-1])


def edit_features(problem,mapping):
    """Same event definition as rank_key; O(bonds) local union connectivity."""
    mapping=dict(mapping);inverse={p:r for r,p in mapping.items()}
    r,p=problem.reactant,problem.product;nr=problem.source_atom_count
    target_node=lambda a:inverse.get(a,nr+a)
    union=[];centers=set();edits=[];events=Counter(broken=0,formed=0,changed=0)
    for a,b in zip(*np.where(np.triu(r.wbo,1)>.2)):
        a,b=int(a),int(b);union.append((a,b));wr=float(r.wbo[a,b])
        if a not in mapping or b not in mapping:
            if a in mapping or b in mapping:
                events['broken']+=1;centers.update((a,b));edits.append(((r.elements[a],r.elements[b]),wr,0.))
        else:
            wp=float(p.wbo[mapping[a],mapping[b]])
            kind='broken' if wp<=.2 else 'changed' if abs(wr-wp)>.5 else None
            if kind:
                events[kind]+=1;centers.update((a,b));edits.append(((r.elements[a],r.elements[b]),wr,wp))
    for a,b in zip(*np.where(np.triu(p.wbo,1)>.2)):
        a,b=int(a),int(b);x,y=target_node(a),target_node(b);union.append((x,y))
        if a not in inverse or b not in inverse or r.wbo[inverse[a],inverse[b]]<=.2:
            events['formed']+=1;centers.update((x,y));edits.append(((p.elements[a],p.elements[b]),0.,float(p.wbo[a,b])))
    parents={a:a for a in centers}
    def root(a):
        while parents[a]!=a:parents[a]=parents[parents[a]];a=parents[a]
        return a
    for a,b in union:
        if a in centers and b in centers:parents[root(a)]=root(b)
    weighted=0.;missing=set()
    for elements,before,after in edits:
        try:weighted+=abs(bond_energy(elements,before)-bond_energy(elements,after))
        except KeyError as exc:missing.add(str(exc.args[0]))
    return dict(events=sum(events.values()),center_atoms=len(centers),
        center_components=len({root(a) for a in centers}),
        weighted_edits=None if missing else weighted/346.,missing_energies=sorted(missing))


VARIANTS=('baseline','compact_tie','compact_soft_0.5','compact_soft_1','compact_soft_2',
          'energy_tie','energy_primary','combined_tie','combined_soft')


def candidate_key(row,variant):
    f=row['features'];e=f['events'];c=max(0,f['center_components']-1);n=f['center_atoms'];w=f['weighted_edits']
    coverage=tuple(row['score'][:2]);tie=row['rank']
    if variant=='baseline':return coverage+(e,tie)
    if variant=='compact_tie':return coverage+(e,c,n,tie)
    if variant.startswith('compact_soft_'):return coverage+(e+float(variant.rsplit('_',1)[1])*c,n,e,tie)
    if variant=='energy_tie':return coverage+(e,w,tie)
    if variant=='energy_primary':return coverage+(w,e,tie)
    if variant=='combined_tie':return coverage+(e,c,w,n,tie)
    if variant=='combined_soft':return coverage+(e+c+.25*w,n,e,tie)
    raise ValueError(variant)


def init(args):
    args.run.mkdir(parents=True,exist_ok=False)
    for folder in ('src','bench'):
        shutil.copytree(folder,args.run/'engine'/folder,ignore=shutil.ignore_patterns('__pycache__'))
    manifest=json.loads((args.source/'manifest.json').read_text())
    save(args.run/'manifest.json',dict(source=str(args.source.resolve()),classes=str(args.classes.resolve()),
        indices=[r['index'] for r in manifest['records']],workers=args.workers,
        energy_source=ENERGY_SOURCE,variants=VARIANTS,
        convention='Unsigned bond-energy edit cost; aromatic order interpolated. No unsupported-case fallback. No core search or family enumeration.',
        source_hashes={str(p.relative_to(args.run/'engine')):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (args.run/'engine').rglob('*.py')}))


def extract(args):
    manifest=json.loads((args.run/'manifest.json').read_text());started=time.perf_counter()
    for index in manifest['indices'][args.slot::manifest['workers']]:
        out=args.run/str(index);out.mkdir(exist_ok=True)
        directory,plan=load_case(Path(manifest['source']),index)
        archive=directory/'cuts/aam.pkl.gz'
        if not archive.exists():save(out/'features.json',dict(index=index,status='incomplete_search'));continue
        start=time.perf_counter();aam=read_aam_checkpoint(archive);loaded=time.perf_counter()
        classes=json.loads((Path(manifest['classes'])/str(index)/'classes.json').read_text())
        for row in classes:
            mapping=plan.to_input_mapping(aam.graph.states[row['terminal']].mapping)
            row['features']=edit_features(plan.input_problem,mapping)
            assert row['features']['events']==row['score'][2]
        eligible=[r for r in classes if r['score'][:2]==classes[0]['score'][:2]] if classes else []
        energy_supported=all(r['features']['weighted_edits'] is not None for r in eligible)
        save(out/'features.json',dict(index=index,status='complete',classes=classes,
            energy_supported=energy_supported,archive_seconds=loaded-start,
            feature_seconds=time.perf_counter()-loaded,total_seconds=time.perf_counter()-start))
    print(json.dumps(dict(slot=args.slot,seconds=time.perf_counter()-started)))


def worker(args):
    out=args.run/f'worker_{args.slot}';out.mkdir(exist_ok=True)
    code=guarded([sys.executable,__file__,'extract','--run',str(args.run),'--slot',str(args.slot)],out/'log.txt',300)
    save(out/'status.json',dict(exit_code=code,job=os.environ.get('SLURM_JOB_ID'),finished=time.time()))


def report(args):
    manifest=json.loads((args.run/'manifest.json').read_text());rows=[];details=[]
    for index in manifest['indices']:
        path=args.run/str(index)/'features.json'
        if not path.exists():rows.append(dict(index=index,status='pending'));continue
        case=json.loads(path.read_text())
        if case['status']!='complete':rows.append(case);continue
        classes=case.pop('classes');eligible=[r for r in classes if r['score'][:2]==classes[0]['score'][:2]] if classes else []
        outcomes={}
        for variant in VARIANTS:
            uses_energy='energy' in variant or 'combined' in variant
            if uses_energy and not case['energy_supported']:
                outcomes[variant]=dict(status='unsupported_energy');continue
            ordered=sorted(eligible,key=lambda row:candidate_key(row,variant))
            if not ordered:outcomes[variant]=dict(status='empty',correct=False);continue
            best=ordered[0]
            outcomes[variant]=dict(status='scored',correct=best['reference_equivalent'],
                selected_original_rank=best['rank'],terminal=best['terminal'],features=best['features'],
                top5=any(r['reference_equivalent'] for r in ordered[:5]))
        case['outcomes']=outcomes;rows.append(case)
        details.append(dict(index=index,classes=classes))
    complete=[r for r in rows if r['status']=='complete'];total=len(rows);summary={}
    for variant in VARIANTS:
        scored=[r for r in complete if r['outcomes'][variant]['status']=='scored']
        wins=sum(r['outcomes'][variant]['correct'] for r in scored)
        baseline=sum(r['outcomes']['baseline']['correct'] for r in scored)
        summary[variant]=dict(scored=len(scored),correct=wins,all_records_percent=100*wins/total,
            scored_percent=100*wins/max(1,len(scored)),baseline_same_cases=baseline,
            gained=[r['index'] for r in scored if r['outcomes'][variant]['correct'] and not r['outcomes']['baseline']['correct']],
            lost=[r['index'] for r in scored if not r['outcomes'][variant]['correct'] and r['outcomes']['baseline']['correct']])
    result=dict(records=total,statuses=dict(Counter(r['status'] for r in rows)),variants=summary,
        archive_seconds=sum(r['archive_seconds'] for r in complete),feature_seconds=sum(r['feature_seconds'] for r in complete),
        candidates=sum(len(r['classes']) for r in details),
        note='Exploratory fixed-formula ablation on Golden; not held-out validation. Same cached class representative per method; no witness optimization within classes. Coverage priority unchanged. Unknown energies exclude the whole case from energy variants; compare on the identical supported subset.')
    args.output.mkdir(parents=True,exist_ok=True);save(args.output/'summary.json',result);save(args.output/'cases.json',rows)
    (args.output/'features.json.gz').write_bytes(gzip.compress(json.dumps(details).encode(),mtime=0))
    save(args.output/'manifest.json',manifest)
    if args.jobs:(args.output/'slurm_accounting.psv').write_text(subprocess.check_output(['sacct','-j',args.jobs,'-P','--format=JobID,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,Start,End'],text=True))
    print(json.dumps({k:({v:{a:b for a,b in d.items() if a not in ('gained','lost')} for v,d in value.items()} if k=='variants' else value) for k,value in result.items()},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['init','extract','worker','report'])
    p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path);p.add_argument('--classes',type=Path)
    p.add_argument('--slot',type=int);p.add_argument('--workers',type=int,default=32);p.add_argument('--output',type=Path);p.add_argument('--jobs')
    a=p.parse_args();dict(init=init,extract=extract,worker=worker,report=report)[a.mode](a)
