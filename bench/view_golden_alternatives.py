"""Browse every chemical class of saved representative mappings, without search."""
import argparse
from collections import Counter
from dataclasses import asdict
import html
import gzip
import json
from pathlib import Path
import time

import pynauty
from rdkit import Chem

from golden_evaluation import colored_graph,project,rank_key
from golden_policy_campaign import load_case,save
from view_golden_mapping import molecules,render_viewer
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.search_graph import SearchPath
from collect_golden_patterns import collect,first_paths,equivalence_from_molecules


def build(source,audit,index,output,*,direction=None,pattern_seconds=0,path_seconds=1,pattern_run=None):
    started=time.perf_counter()
    if direction is None:
        directory,plan=load_case(source,index)
        reference_path=directory/'reference.json'
    else:
        from golden_publication import plans
        pair,_=plans(source,index);plan=pair[direction]
        directory=source/'directions'/str(index)/direction
        reference_path=source/'inputs'/str(index)/'reference.json'
    load_started=time.perf_counter()
    aam=read_aam_checkpoint(directory/'cuts/aam.pkl.gz')
    archive_loading_seconds=time.perf_counter()-load_started
    first=first_paths(aam.graph)
    problem=plan.input_problem
    reference=json.loads(reference_path.read_text())
    expected=dict(reference['mapping']);features=reference['features']
    ref_certificate=pynauty.certificate(colored_graph(features,project(expected,features)))
    reactant,product=molecules(audit[index]['mapped_reaction'],problem)
    components=Chem.GetMolFrags(reactant)
    owner={atom:i for i,atoms in enumerate(components) for atom in atoms}
    groups={};heavy_relations=set()
    for terminal in aam.graph.terminals:
        mapping=plan.to_input_mapping(aam.graph.states[terminal].mapping)
        heavy=project(mapping,features);heavy_relations.add(tuple(sorted(heavy.items())))
        certificate=pynauty.certificate(colored_graph(features,heavy))
        key,events=rank_key(mapping,problem)
        item=groups.setdefault(certificate,dict(count=0,heavy=set(),all_atoms=set(),events=set(),best=None))
        item['count']+=1;item['heavy'].add(tuple(sorted(heavy.items())))
        item['all_atoms'].add(tuple(sorted(mapping.items())))
        item['events'].add(tuple(events.values()))
        if item['best'] is None or key<item['best'][0]:item['best']=(key,terminal,mapping,events)
    ordered=sorted(groups.items(),key=lambda pair:pair[1]['best'][0])
    def record(label,mapping,path=None,**extra):
        steps=[]
        if path is not None:
            for edge_id in path.transitions:
                edge=aam.graph.transitions[edge_id]
                if edge.match is None:continue
                fragment=edge.match['fragment']
                if plan.reversed:
                    physical=plan.to_search_mapping(mapping)
                    fragment=[physical[r] for r in fragment]
                steps.append(dict(transition=edge_id,seed=edge.seed,fragment=list(fragment)))
        context=None if path is None else dict(asdict(path.context),seed_side='p' if plan.reversed else 'r')
        return dict(label=label,mapping=mapping,context=context,steps=steps,
            heavy_donors=dict(Counter(owner[r] for r in mapping if problem.reactant.elements[r]!='H')),**extra)
    records=[record('Reference',expected)];classes=[]
    for rank,(certificate,item) in enumerate(ordered,1):
        key,terminal,mapping,events=item['best']
        path=first(terminal)
        correct=certificate==ref_certificate
        row=dict(rank=rank,terminal=terminal,reference_equivalent=correct,
            saved_terminals=item['count'],distinct_raw_heavy=len(item['heavy']),
            distinct_explicit_witnesses=len(item['all_atoms']),
            event_variants=[list(x) for x in sorted(item['events'])],
            heavy_mapped=-key[0],total_mapped=-key[1],**events)
        label=f'#{rank}'+(' — REFERENCE-EQUIVALENT' if correct else '')
        label+=f' | break {events["broken"]}, form {events["formed"]}, order {events["bond_order_changed"]}'
        label+=f' | {item["count"]} saved terminals'
        records.append(record(label,mapping,path,terminal=terminal,class_rank=rank))
        classes.append(row)
    out=output/f'case{index}';out.mkdir(parents=True,exist_ok=True)
    extraction=None
    if pattern_seconds>0 or pattern_run is not None:
        if pattern_run is not None:
            if pattern_seconds>0:raise ValueError('Choose saved collection or live extraction, not both')
            folder=pattern_run/'results'/str(index)/plan.direction
            summary=json.loads((folder/'summary.json').read_text())
            if summary['stage']!='finished':raise ValueError('Saved collection worker has not finished')
            from golden_publication import sha256
            if summary['archive_sha256']!=sha256(directory/'cuts/aam.pkl.gz'):
                raise ValueError('Saved patterns and viewer archive do not match')
            with gzip.open(folder/'patterns.jsonl.gz','rt') as stream:patterns=list(map(json.loads,stream))
            for p in patterns:
                p['origins']=[dict(terminal=p['terminal'],transitions=p['transitions'],actions=p['actions'],
                    mapping=p['mapping'],orientation=plan.direction)]
            equivalence=equivalence_from_molecules(plan,(reactant,product))
            extraction=dict(complete=summary['complete'],visited_paths=summary['visited_families'],
                            checkpoint=str(folder/'patterns.jsonl.gz'))
        else:
            state,patterns,equivalence=collect(directory/'cuts/aam.pkl.gz',aam,plan,(reactant,product),
                [item['best'][1] for _,item in ordered],out/'pattern_checkpoint.json',pattern_seconds,path_seconds)
            extraction=dict(complete=state['complete'],visited_paths=len(state['paths']),
                            checkpoint=str(out/'pattern_checkpoint.json'))
        # Preserve every old displayed witness. New patterns are deduplicated
        # by exact explicit-atom identity, never by score or heavy projection.
        existing={equivalence.key(plan.to_search_mapping(r['mapping'])):r for r in records[1:]}
        for p in patterns:
            mapping=plan.to_input_mapping(p['mapping'])
            if p['key'] in existing:
                existing[p['key']]['origins']=p['origins']
                existing[p['key']]['pattern_key']=p['key']
                continue
            key,events=rank_key(mapping,problem)
            correct=pynauty.certificate(colored_graph(features,project(mapping,features)))==ref_certificate
            path=SearchPath(aam.graph,p['terminal'],tuple(p['transitions']))
            classes.append(dict(rank=0,terminal=p['terminal'],reference_equivalent=correct,
                saved_terminals=len({o['terminal'] for o in p['origins']}),distinct_raw_heavy=1,
                distinct_explicit_witnesses=1,heavy_mapped=-key[0],total_mapped=-key[1],
                pattern_key=p['key'],origin='compressed_family',**events))
            records.append(record('Extracted pattern',mapping,path,terminal=p['terminal'],
                pattern_key=p['key'],origins=p['origins']))
            existing[p['key']]=records[-1]
        paired=sorted(zip(classes,records[1:]),key=lambda pair:rank_key(pair[1]['mapping'],problem)[0])
        classes=[];records=records[:1]
        for rank,(row,rec) in enumerate(paired,1):
            row['rank']=rank;rec['class_rank']=rank
            rec['label']=f'#{rank}'+(' — REFERENCE-EQUIVALENT' if row['reference_equivalent'] else '')
            rec['label']+=f' | break {row["broken"]}, form {row["formed"]}, order {row["bond_order_changed"]}'
            rec['label']+=' | '+row.get('origin','saved representative')
            classes.append(row);records.append(rec)
    evaluation=json.loads((directory/'evaluation.json').read_text())
    payload=dict(index=index,archive=str(directory/'cuts/aam.pkl.gz'),records=records,
        components=[list(c) for c in components],evaluation=evaluation,
        cap_stops=sum(s.reason=='capped' for s in aam.graph.stops))
    if extraction is not None:
        payload['pattern_extraction']=extraction
    render_started=time.perf_counter()
    page=render_viewer(payload,reactant,product,problem,aam.config)
    rendering_seconds=time.perf_counter()-render_started
    # The existing renderer gives linked R/P drawings, explicit H and actual
    # fragment colors. Add a class selector without changing old viewers.
    buttons=''.join(f'<tr class="class-row" data-choice="{c["rank"]}"><td>#{c["rank"]}</td>'
        f'<td>{"YES" if c["reference_equivalent"] else "—"}</td>'
        f'<td>{c["heavy_mapped"]}/{len(features[1]["heavy"])}</td>'
        f'<td>{c["broken"]}</td><td>{c["formed"]}</td><td>{c["bond_order_changed"]}</td>'
        f'<td>{c["saved_terminals"]}</td><td>{c["distinct_raw_heavy"]}</td>'
        f'<td>{c["distinct_explicit_witnesses"]}</td></tr>' for c in classes)
    max_events=max(sum(c[k] for k in ('broken','formed','bond_order_changed')) for c in classes)
    plot=[]
    for c in classes:
        x=65+690*(c['rank']-1)/max(1,len(classes)-1)
        events=sum(c[k] for k in ('broken','formed','bond_order_changed'))
        y=210-165*events/max(1,max_events)
        plot.append(f'<g class="class-point" data-choice="{c["rank"]}" tabindex="0" role="button">'
            f'<title>#{c["rank"]}: {events} bond events; {c["saved_terminals"]} terminals</title>'
            f'<circle cx="{x}" cy="{y}" r="7" fill="{"#009e73" if c["reference_equivalent"] else "#0072b2"}"/>'
            f'<text x="{x}" y="{y-12}" text-anchor="middle" font-size="10">{c["rank"]}</text></g>')
    overview=f'''<div class="overview"><h2>Alternative atom assignments</h2>
<p>{len(aam.graph.terminals):,} saved terminal states → {len(heavy_relations)} distinct raw heavy mappings → <b>{len(classes)} chemical-equivalence classes</b>.</p>
<p>All classes are shown, ordered by the existing reference-blind ranker. Each displays its best-ranked saved witness and one actual fragment path. This does not enumerate alternatives inside compressed families. Green means this displayed representative is reference-equivalent; blue does not prove its whole family lacks the reference.</p>
<p>Identical atom assignments may recur across seeds, cuts and histories. Raw-heavy variants within one class are chemically symmetry-equivalent. Additional explicit-atom variants change H assignments. Bond-event scores include H; graph-only ranking does not establish chemical plausibility.</p>
<svg viewBox="0 0 820 260" style="max-width:1000px;background:white"><path d="M55 25V220H775" fill="none" stroke="black"/>
<text x="10" y="18">Bond events (0–{max_events})</text><text x="360" y="248">Class rank — click a point to inspect</text>{''.join(plot)}</svg>
<details><summary>All class counts and scores — click a row</summary><table><tr><th>Rank</th><th>Reference</th><th>Heavy P</th><th>Break</th><th>Form</th><th>Order</th><th>Terminals</th><th>Raw heavy</th><th>Explicit witnesses</th></tr>{buttons}</table></details></div>'''
    page=page.replace('<div class="controls">',overview+'<div class="controls">',1)
    if extraction is not None:
        status='COMPLETE within the saved archive' if extraction['complete'] else 'INCOMPLETE — additional patterns remain unexamined'
        page=page.replace('This does not enumerate alternatives inside compressed families.',
            'Additional strict, explicit-atom patterns extracted from compressed families are also shown; equal scores are not merged.')
        page=page.replace('<h2>Alternative atom assignments</h2>',
            f'<h2>Alternative atom assignments</h2><p style="color:#b45309;font-weight:bold">Pattern extraction: {status}. '
            f'{extraction["visited_paths"]} paths/families visited. Original compressed archive retained; extraction can resume. '
            'Ranks are provisional until extraction completes. Reference labels are checked only after extraction.</p>')
        page=page.replace(f'<b>{len(classes)} chemical-equivalence classes</b>',
            f'<b>{len(classes)} displayed candidates (preserved representatives + extracted patterns)</b>')
    page+='''<style>.class-row,.class-point{cursor:pointer}.class-row:hover{background:#daeafa}.class-point:hover circle{stroke:#e69f00;stroke-width:4}.overview{border:1px solid #ccd5dd;padding:18px}.class-active{background:#d1f0df}</style>
<script>document.querySelectorAll('[data-choice]').forEach(el=>{function choose(){mode.selectedIndex=Number(el.dataset.choice);render();document.querySelectorAll('.class-active').forEach(x=>x.classList.remove('class-active'));el.classList.add('class-active');document.querySelector('.controls').scrollIntoView({behavior:'smooth',block:'start'});}el.onclick=choose;el.onkeydown=e=>{if(e.key==='Enter')choose();};});</script>'''
    write_started=time.perf_counter()
    save(out/'mapping.json',payload);save(out/'classes.json',classes)
    (out/'viewer.html').write_text(page)
    artifact_writing_seconds=time.perf_counter()-write_started
    summary=dict(index=index,terminals=len(aam.graph.terminals),raw_heavy=len(heavy_relations),classes=len(classes),
        representative_reference_ranks=[c['rank'] for c in classes if c['reference_equivalent']],
        generation_seconds=time.perf_counter()-started,html_bytes=len(page.encode()),
        archive_loading_seconds=archive_loading_seconds,rendering_seconds=rendering_seconds,
        artifact_writing_seconds=artifact_writing_seconds)
    if extraction is not None:summary['pattern_extraction']=payload['pattern_extraction']
    save(out/'summary.json',summary)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--audit',type=Path,required=True)
    parser.add_argument('--indices',type=int,nargs='+',required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--direction',choices=('R_to_P','P_to_R'),help='Use publication archive layout')
    parser.add_argument('--pattern-seconds',type=float,default=0,help='Resumable extraction budget per case; no new AAM search')
    parser.add_argument('--path-seconds',type=float,default=1,help='Soft query budget per path per pass; unresolved paths remain in checkpoint')
    parser.add_argument('--pattern-run',type=Path,help='Render an existing benchmark collection without solving again')
    args=parser.parse_args()
    audit={r['index']:r for r in map(json.loads,args.audit.read_text().splitlines())}
    summaries=[build(args.source,audit,index,args.output,direction=args.direction,
        pattern_seconds=args.pattern_seconds,path_seconds=args.path_seconds,pattern_run=args.pattern_run) for index in args.indices]
    save(args.output/'summary.json',summaries)
    links=''.join(f'<li><a href="case{s["index"]}/viewer.html">Case {s["index"]}: {s["classes"]} classes</a> — reference representative rank {html.escape(str(s["representative_reference_ranks"]))}</li>' for s in summaries)
    (args.output/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>AAM alternatives</title><h1>Inspect saved AAM alternatives</h1><p>No new searches or chirality filtering. Each viewer is standalone and works offline.</p><ul>'+links+'</ul>')
    print(json.dumps(summaries,indent=2))
