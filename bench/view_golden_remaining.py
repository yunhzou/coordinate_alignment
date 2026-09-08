"""Offline simultaneous reference/detected views from saved AAM archives only."""
import argparse
import colorsys
from collections import Counter
import heapq
import html
import json
from pathlib import Path
import time
import zipfile

from rdkit import Chem
from rdkit.Geometry import Point3D
from golden_policy_campaign import load_case,save
from view_golden_mapping import molecules,svg,mapping_fragments,PALETTE
from rxn_core.artifacts import read_aam_checkpoint,read_raw_cut,raw_cut_paths

INDICES=[7,19,590,603,780,833,845,865,867,871,986,1228,1285,1358,1377,1380,1475,1553,1568,1574,1793]


def packed(mol,hydrogens):
    """Translate disconnected 2D components into a wide, readable layout."""
    mol=Chem.Mol(mol) if hydrogens else Chem.RemoveHs(mol)
    conf=mol.GetConformer();parts=Chem.GetMolFrags(mol);boxes=[]
    for atoms in parts:
        xs=[conf.GetAtomPosition(i).x for i in atoms];ys=[conf.GetAtomPosition(i).y for i in atoms]
        boxes.append((min(xs),min(ys),max(max(xs)-min(xs),1.5),max(max(ys)-min(ys),1.5)))
    def shape(columns):
        rows=[boxes[i:i+columns] for i in range(0,len(boxes),columns)]
        width=max(sum(b[2] for b in row)+2.5*(len(row)-1) for row in rows)
        height=sum(max(b[3] for b in row) for row in rows)+2.5*(len(rows)-1)
        return max(width/1050,height/300)
    columns=min(range(1,len(parts)+1),key=shape);y=0
    for start in range(0,len(parts),columns):
        height=max(b[3] for b in boxes[start:start+columns]);x=0
        for atoms,box in zip(parts[start:start+columns],boxes[start:start+columns]):
            for atom in atoms:
                old=conf.GetAtomPosition(atom);conf.SetAtomPosition(atom,Point3D(old.x-box[0]+x,old.y-box[1]+y,0))
            x+=box[2]+2.5
        y+=height+2.5
    return mol


def draw_record(record,reactant,product,problem,owner):
    mapping=record['mapping']
    rcolors={r:f['color'] for f in record['regions'] for r in f['source']}
    pcolors={p:rcolors[r] for r,p in mapping.items()}
    bonds={tuple(sorted((r,s))) for r in mapping for s in mapping if r<s and problem.reactant.wbo[r,s]>.2 and problem.product.wbo[mapping[r],mapping[s]]>.2}
    pbonds={tuple(sorted((mapping[r],mapping[s]))) for r,s in bonds}
    record['drawings']={str(int(h)):dict(R=svg(packed(reactant,h),rcolors,True,1050,280,mapping=mapping,owner=owner,bonds=bonds,reference=record['context'] is None),
        P=svg(packed(product,h),pcolors,True,1050,300,mapping=mapping,owner=owner,side='P',bonds=pbonds,reference=record['context'] is None)) for h in (False,True)}


def ranker(problem,heavy_only=False):
    endpoints=(problem.reactant,problem.product)
    edges=[[(a,b,float(e.wbo[a,b])) for a in range(e.atom_count) for b in range(a+1,e.atom_count)
        if e.wbo[a,b]>.2 and (not heavy_only or e.elements[a]!='H' and e.elements[b]!='H')]
        for e in endpoints]
    def score(mapping):
        mapping=dict(mapping)
        if heavy_only:mapping={r:p for r,p in mapping.items() if endpoints[0].elements[r]!='H' and endpoints[1].elements[p]!='H'}
        inverse={p:r for r,p in mapping.items()};broken=formed=changed=0
        for a,b,w in edges[0]:
            if a not in mapping or b not in mapping:broken+=int(a in mapping or b in mapping)
            elif endpoints[1].wbo[mapping[a],mapping[b]]<=.2:broken+=1
            elif abs(w-endpoints[1].wbo[mapping[a],mapping[b]])>.5:changed+=1
        for a,b,w in edges[1]:
            if a not in inverse or b not in inverse or endpoints[0].wbo[inverse.get(a,0),inverse.get(b,0)]<=.2:formed+=1
        heavy=sum(endpoints[1].elements[p]!='H' for p in mapping.values())
        return (-heavy,-len(mapping),broken+formed+changed,tuple(sorted(mapping.items()))),dict(broken=broken,formed=formed,order_changed=changed)
    return score


def graph_for(directory):
    archive=directory/'cuts/aam.pkl.gz'
    if archive.exists():return read_aam_checkpoint(archive).graph,archive,'complete saved sweep'
    # Explicit incomplete-run display policy, never called a complete ranking.
    for archive in raw_cut_paths(directory/'cuts'):
        graph=read_raw_cut(archive)
        if graph.terminals:return graph,archive,'partial run: first completed cut with terminal states only'
    raise ValueError(f'No saved terminal states to visualize: {directory}')


def build(args):
    start=time.perf_counter();index=INDICES[args.slot]
    directory,plan=load_case(args.source/'random',index);problem=plan.input_problem
    reference=json.loads((directory/'reference.json').read_text())
    audit=next(json.loads(line) for line in args.audit.read_text().splitlines() if json.loads(line)['index']==index)
    reactant,product=molecules(audit['mapped_reaction'],problem)
    score=ranker(problem);heavy_score=ranker(problem,True);candidates=[];archives=[]
    for policy in ('random','distance'):
        directory,other=load_case(args.source/policy,index)
        graph,archive,scope=graph_for(directory)
        groups={}
        for terminal in graph.terminals:
            mapping=other.to_input_mapping(graph.states[terminal].mapping)
            key,events=score(mapping)
            signature=tuple((r,p) for r,p in sorted(mapping.items()) if problem.reactant.elements[r]!='H' and problem.product.elements[p]!='H')
            if signature not in groups or key<groups[signature][0]:groups[signature]=(key,terminal,mapping,events)
        for key,terminal,mapping,events in heapq.nsmallest(3,groups.values(),key=lambda x:x[0]):
            path=next(graph.paths(terminal));final_mapping=path.mapping
            steps=[]
            for eid in path.transitions:
                edge=graph.transitions[eid]
                if edge.match is None:continue
                group=edge.match['fragment']
                if other.reversed:group=[final_mapping[r] for r in group]
                steps.append(dict(fragment=list(group)))
            candidates.append(dict(label=policy,mapping=mapping,steps=steps,context={'cuts':path.context.cuts},
                terminal=terminal,archive=str(archive),scope=scope,key=key,events=events))
        archives.append(dict(policy=policy,path=str(archive),scope=scope,terminals=len(graph.terminals),
            distinct_raw_heavy=len(groups),capped=graph.capped,
            evaluation=json.loads((directory/'evaluation.json').read_text()).get('reference_recovery','unknown')))
        del graph,groups
    selected=[];seen=set()
    for record in sorted(candidates,key=lambda x:x['key']):
        signature=tuple((r,p) for r,p in sorted(record['mapping'].items()) if problem.reactant.elements[r]!='H' and problem.product.elements[p]!='H')
        if signature in seen:continue
        seen.add(signature);selected.append(record)
        if len(selected)==3:break
    records=[dict(label='Ground truth',mapping=dict(reference['mapping']),steps=[],context=None)]+selected
    colors={};components=Chem.GetMolFrags(reactant);owner={a:i for i,atoms in enumerate(components) for a in atoms}
    for record in records:
        record.pop('key',None)
        record['regions']=mapping_fragments(record,problem)
        for region in record['regions']:
            heavy=tuple(r for r in region['source'] if problem.reactant.elements[r]!='H')
            signature=('heavy',heavy) if heavy else ('H',tuple(region['source']))
            if signature not in colors:
                n=len(colors);rgb=colorsys.hsv_to_rgb((n*.61803398875)%1,.65,.8)
                colors[signature]=PALETTE[n] if n<len(PALETTE) else '#'+''.join(f'{round(x*255):02x}' for x in rgb)
            region['color']=colors[signature]
        mapping=record['mapping'];record['heavy_events']=heavy_score(mapping)[1]
        record['mapped_heavy']=sum(problem.product.elements[p]!='H' for p in mapping.values())
        record['mapped_total']=len(mapping)
        draw_record(record,reactant,product,problem,owner)
    payload=dict(index=index,records=records,archives=archives,
        p_elements=problem.product.elements,r_elements=problem.reactant.elements,
        target_heavy=sum(e!='H' for e in problem.product.elements),target_total=problem.target_atom_count,
        reaction=audit['mapped_reaction'],reference=dict(reference['mapping']),generation_seconds=time.perf_counter()-start)
    args.output.mkdir(parents=True,exist_ok=True);save(args.output/f'{index}.json',payload)
    print(json.dumps(dict(index=index,seconds=payload['generation_seconds'],records=len(records))),flush=True)


def redraw(args):
    for index in INDICES:
        payload=json.loads((args.output/f'{index}.json').read_text())
        _,plan=load_case(args.source/'random',index);problem=plan.input_problem
        reactant,product=molecules(payload['reaction'],problem)
        owner={a:i for i,atoms in enumerate(Chem.GetMolFrags(reactant)) for a in atoms}
        for record in payload['records']:
            record['mapping']={int(r):p for r,p in record['mapping'].items()}
            draw_record(record,reactant,product,problem,owner)
        save(args.output/f'{index}.json',payload)


def combine(args):
    cases=[json.loads((args.output/f'{i}.json').read_text()) for i in INDICES]
    template=Path(__file__).with_name('golden_remaining_template.html').read_text()
    data=json.dumps(cases,separators=(',',':')).replace('</',r'<\/')
    target=args.html_directory or args.output;target.mkdir(parents=True,exist_ok=True)
    (target/'viewer.html').write_text(template.replace('__DATA__',data))
    with zipfile.ZipFile(target/'viewer.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        archive.write(target/'viewer.html','viewer.html')
    save(target/'summary.json',[dict(index=c['index'],archives=c['archives'],generation_seconds=c['generation_seconds'],
        representatives=len(c['records'])-1) for c in cases])
    save(target/'mappings.json',[{**c,'records':[{k:v for k,v in record.items() if k!='drawings'} for record in c['records']]} for c in cases])
    print(str((target/'viewer.html').resolve()))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['build','combine','redraw'])
    p.add_argument('--source',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--slot',type=int)
    p.add_argument('--html-directory',type=Path)
    p.add_argument('--audit',type=Path,default=Path('data/aam_benchmarks/golden_original_20260906/audit.jsonl'))
    args=p.parse_args();globals()[args.mode](args)
