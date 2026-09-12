"""Independent witness audit and report for the fixed-budget 140-case experiment."""
import sys,json,time,hashlib,collections,shutil,statistics
from pathlib import Path
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');OUT=ROOT/'outputs/competition140_budget128';CODE=ROOT/'work/aam-event-improved'
sys.path[:0]=[str(CODE/'src'),str(CODE/'bench')]
from rxn_core.artifacts import read_aam_checkpoint
from rxn_core.event_patterns import SignedEventIndex
from rxn_core.family_scoring import validate_representative
from rxn_core.family_query import query_path
from rxn_core.frag import build_graph
from rxn_core.matcher import _nauty_orbits
from rxn_core.fragment import match_fragment,FragmentMatchConfig,FragmentMatchContext
from rxn_core import viewers
from fragment_competition_optimized import save
read=lambda p:json.loads(Path(p).read_text());hashes={}
def load(p):
 p=Path(p);hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return read(p)
base=load(ROOT/'outputs/holdout140_decoded_results/results.json')['cases'];execution=load(OUT/'execution.json');dec_execution=load(OUT/'decode_execution.json');manifest=load(OUT/'manifest.json')
assert all(x['status']=='passed' for e in (execution,dec_execution) for x in e['jobs'])
assert len(execution['jobs'])==len(dec_execution['jobs'])==140
stats={k:dict(baseline_classes=0,union_classes=0,total_classes=0,baseline_complete_cases=0,union_complete_cases=0,missing=[]) for k in ('native_slap','slap_sweep')}
rows=[];documents=[];proofs=[];counts=collections.Counter();dec_counts=collections.Counter();start=time.perf_counter();cpu=time.process_time()
for br in base:
 c=br['case'];folder=OUT/f'case{c}';b=load(br['source']);d=load(folder/'result.json');f=load(folder/'full_decode.json')
 assert f['complete_saved_windows'];assert d['operation_budget']==128 and d['config']['iso_tolerance']==1 and d['config']['branch_limit']==2000
 assert d['archive_sha256']==manifest['inputs'][str(c)]['sha256']
 counts.update(d['counts']);dec_counts.update(f['counts']);old=set(b['patterns']);union=old|d['patterns'].keys()|f['patterns'].keys();window=br['max_events']
 new_window=set(f['patterns'])-old;comparisons={};refs=set()
 for name,st in stats.items():
  wanted=set(b['slap'][name]['minimum_ids']);refs|=wanted
  st['baseline_classes']+=len(wanted&old);st['union_classes']+=len(wanted&union);st['total_classes']+=len(wanted);st['baseline_complete_cases']+=wanted<=old;st['union_complete_cases']+=wanted<=union
  missing=sorted(wanted-union)
  if missing:st['missing'].append(dict(case=c,ids=missing))
  comparisons[name]=dict(baseline=len(wanted&old),union=len(wanted&union),total=len(wanted),new=sorted((wanted-old)&union),missing=missing)
 a=read_aam_checkpoint(Path(d['archive']));idx=SignedEventIndex(a.problem);archives={};records=[]
 def get_archive(path):
  if path not in archives:archives[path]=read_aam_checkpoint(Path(path));hashes[path]=hashlib.sha256(Path(path).read_bytes()).hexdigest()
  return archives[path]
 def add(p,label,provenance):
  checked=idx.describe(p['mapping']);assert checked['id']==p['id'] and checked['total']==p['total']
  m=p['mapping'];events=[]
  for sign,pairs in checked['events'].items():
   for r,t in pairs:
    before=float(idx.r[r,t]);after=float(idx.p[m[r],m[t]])
    kind=('weakened' if after>.2 else 'broken') if sign=='broken' else ('strengthened' if before>.2 else 'formed')
    events.append(dict(kind=kind,r=[r,t],p=[m[r],m[t]],wbo=[before,after]))
  records.append(dict(name=label,mapping=m,events=events,counts=dict(collections.Counter(e['kind'] for e in events),total=len(events)),pattern=p['id'],provenance=provenance))
 # Audit all new classes in the fully decoded window, including non-minimum ones.
 for key in sorted(new_window,key=lambda k:(f['patterns'][k]['total'],k)):
  p=f['patterns'][key];aa=get_archive(p['archive']);path=next(x for x in aa.graph.paths(p['terminal']) if list(x.transitions)==p['transitions']);validate_representative(path,aa.problem)
  mapping=dict(enumerate(p['mapping']));assert all(mapping[r]==v for r,v in aa.config.anchors)
  status,certificate=query_path(path,aa.problem,mapping,source_atoms=tuple(range(idx.n)),timeout_ms=5000)
  assert status=='recovered',(c,key,status,certificate)
  assert dict(certificate['mapping'])==mapping
  proof=dict(case=c,id=key,total=p['total'],archive=p['archive'],terminal=p['terminal'],exact_membership=True,anchors_preserved=True,checked_fragment_edges=certificate['checked_fragment_edges'])
  # Replay the B collision for every newly recovered reference minimum.
  if key in refs:
   rep=d['patterns'].get(key)
   if rep is not None:
    prov=rep['provenance'];rgraph,tgraph=[build_graph(ep.elements,ep.wbo,bond_cut=.2) for ep in (aa.problem.reactant,aa.problem.product)]
    sub=rgraph.subgraph(set(prov['B_atoms'])|{prov['contact']}).copy()
    proposal=match_fragment(sub,tgraph,seed=prov['seed'],context=FragmentMatchContext(target_orbits=_nauty_orbits(tgraph,wbo_tol=1)),config=FragmentMatchConfig(iso_tolerance=1,branch_limit=2000))
    assert not proposal.capped and any(dict(x)==dict(prov['B_placement']) for x in proposal.matches)
    witness=rep['mapping'];parent=prov['parent']
    assert idx.describe(witness)['id']==key
    witness_archive=get_archive(prov['archive'])
    assert dict(witness_archive.graph.states[prov['terminal']].mapping)==dict(enumerate(witness))
    for witness_path in witness_archive.graph.paths(prov['terminal']):validate_representative(witness_path,witness_archive.problem)
    if isinstance(parent,int):parentmap=dict(a.graph.states[parent].mapping)
    else:
     ap,tid=parent.rsplit(':',1);parentmap=dict(get_archive(ap).graph.states[int(tid)].mapping)
    assert all(parentmap[r]==v and witness[r]==v for r,v in prov['retained_A_pairs'])
    assert all(witness[r]==v for r,v in prov['B_placement'])
    proof.update(collision_replayed=True,unaffected_A_preserved=True,first_repair=int(Path(prov['archive']).name.split('takeover')[1].split('.')[0]),depth=prov['depth'])
  proofs.append(proof)
  label='Competition · new reference-minimum class' if key in refs else 'Competition · additional event class'
  if key not in d['patterns']:label+=' · decoded from symmetry'
  add(p,label,dict(archive=p['archive'],terminal=p['terminal']))
 for key in sorted(old,key=lambda k:(b['patterns'][k]['total'],k)):
  if key in refs or b['patterns'][key]['total']==br['best_found']:add(b['patterns'][key],'Baseline AAM',dict(source=br['source']))
 for name in ('slap_sweep','native_slap'):
  ref=b['slap'][name]
  for key in ref['minimum_ids']:
   label=('SLAP sweep' if name=='slap_sweep' else 'Native SLAP')+(' · still missing' if key not in union else ' · covered')
   add(ref['patterns'][key],label,dict(source=br['source'],method=name))
 raw=load(ROOT/f'work/full140_inputs/{c}/input.json')
 note=f'Uniform 140-case experiment: saved forward AAM at iso_tol=1, one seed order, cap 2,000; add up to 128 B-priority repairs and two collision levels. Baseline retained by union. All accepted repair families decoded through {window} events (same window as baseline), raw 0.5 ordinary / 0.3 metal changes, all H. Rejected repair graphs and omitted priority branches are excluded; this is not exhaustive search. New classes use no SLAP guidance during search.'
 doc=viewers.comparison_document(dict(index=c,name=b.get('name',str(c)),endpoints=[raw['reactant'],raw['product']],records=records,note=note))
 if doc['mechanisms']:doc['default_mech_id']=doc['mechanisms'][0]['id']
 documents.append(doc)
 best=min([br['best_found']]+[p['total'] for p in f['patterns'].values()])
 rows.append(dict(case=c,name=b.get('name'),comparisons=comparisons,new_window_classes=sorted(new_window),baseline_best=br['best_found'],union_best=best,window=window,complete_saved_windows=True,counts=d['counts'],pending_states=d['pending'],search_wall_seconds=d['search_wall_seconds'],search_cpu_seconds=d['search_cpu_seconds'],reported_total_cpu_seconds=d['total_cpu_seconds'],decode_cpu_seconds=f['cpu_seconds'],decode_wall_seconds=f['wall_seconds']))
script_names=('fragment_competition_optimized.py','run_competition140.py','decode_competition140.py','report_competition140.py','check_competition_optimization.py','test_competition_optimized.py','trace_invalid_competition.py','probe_invalid_repair.py')
(OUT/'scripts').mkdir(exist_ok=True)
for name in script_names:
 target=OUT/'scripts'/name;shutil.copy2(ROOT/'work'/name,target);hashes[str(target)]=hashlib.sha256(target.read_bytes()).hexdigest()
assert hashes[str(OUT/'scripts/fragment_competition_optimized.py')]==manifest['sha256']
parity=load(ROOT/'outputs/competition_optimization/pilot320/parity.json');budgets=[]
for budget in (64,128,320):
 folder=ROOT/f'outputs/competition_optimization/pilot{budget}';runs=[load(folder/f'case{c}/result.json') for c in (6,11,64,101)]
 budgets.append(dict(budget=budget,new_sweep_minimum_classes=sum(len(r['newly_recovered']) for r in runs),reported_cpu_seconds=sum(r['cpu_seconds'] for r in runs),wall_seconds=load(folder/'execution.json')['wall_seconds']))
report=dict(comparisons=stats,cases=rows,counts=dict(counts),decode_counts=dict(dec_counts),execution=execution,decode_execution=dec_execution,search_reported_cpu_seconds=sum(r['reported_total_cpu_seconds'] for r in rows),decode_reported_cpu_seconds=sum(r['decode_cpu_seconds'] for r in rows),new_window_class_count=sum(len(r['new_window_classes']) for r in rows),new_window_cases=[r['case'] for r in rows if r['new_window_classes']],best_count_improved_cases=[r['case'] for r in rows if r['union_best']<r['baseline_best']],proofs=proofs,proof_audit_wall_seconds=time.perf_counter()-start,proof_audit_cpu_seconds=time.process_time()-cpu,pilot_parity=parity,pilot_budgets=budgets,source_sha256=hashes,scope='Additional competition over saved baseline archives, not a fresh end-to-end rerun. Full decoding means all accepted saved paths within the existing per-case windows. Fixed budget chosen on four known failures, so this dataset is developmental, not independent validation.')
save(OUT/'summary.json',report);save(OUT/'viewer_data.json',documents)
page=viewers.collection_html(documents,title='B-priority competition · all 140 cases · 128 repairs')
page=page.replace('let showAtomIndices = false;','let showAtomIndices = true;').replace(json.dumps('id="showAtomIndices">')[1:-1],json.dumps('id="showAtomIndices" checked>')[1:-1]);assert viewers.stylesheet() in page
(OUT/'viewer.html').write_text(page)
lines=['# B-priority competition: optimized 140-case experiment','',
'All 140 cases completed with four workers and five-minute per-case watchdogs. The added competition pass improves saved SLAP-sweep minimum-class coverage from **163/169 to 167/169**, and cases with all such classes from **136/140 to 139/140**. Only case 64 still has missing reference alternatives. This measures coverage of SLAP output, not annotated mapping accuracy.','',
'| Reference | Baseline classes | Baseline + competition | Cases with every reference minimum class |','|---|---:|---:|---:|']
for name,s in stats.items():lines.append(f"| {name} | {s['baseline_classes']}/{s['total_classes']} | {s['union_classes']}/{s['total_classes']} | {s['baseline_complete_cases']} → {s['union_complete_cases']} / 140 |")
lines+=['','## Improvements and controlled pilot','',
'The implementation batches terminal scoring, avoids scoring duplicate mapping vectors, caches dependent-component release calculations, uses a deque, and replaces repeated full progress serialization with compact throttled checkpoints. At 320 repairs, the four-case pilot preserves exactly the old parents, operation/rejection counts, output classes, representative mappings and remaining queue length. Recorded CPU fell from 47.25 to 15.36 seconds (3.08×); batch wall time fell from 33.91 to 9.88 seconds. These are selected-case timings, not a measured whole-dataset speedup.','',
'| Fixed repair budget | Recovered missing sweep classes (four-case pilot) | Recorded CPU seconds | Batch wall seconds |','|---:|---:|---:|---:|']
for r in budgets:lines.append(f"| {r['budget']} | {r['new_sweep_minimum_classes']}/6 | {r['reported_cpu_seconds']:.2f} | {r['wall_seconds']:.2f} |")
lines+=['','The 128-repair budget was fixed before the full 140-case run. Its choice used the four known gap cases, so the larger run is development-set validation, not an independent test set. No SLAP label or witness is read during proposal search.','',
'## Larger-run cost','',
f"Competition: **{execution['wall_seconds']:.2f} s wall**, {report['search_reported_cpu_seconds']:.2f} recorded CPU seconds including initial archive loading; peak combined sampled worker RSS **{execution['peak_total_mib']/1024:.2f} GiB**.",
f"Full accepted-family decoding: **{dec_execution['wall_seconds']:.2f} s wall**, {report['decode_reported_cpu_seconds']:.2f} recorded CPU seconds; peak combined sampled worker RSS {dec_execution['peak_total_mib']/1024:.2f} GiB.",
f"Sequential pass wall total: **{execution['wall_seconds']+dec_execution['wall_seconds']:.2f} s**. Recorded worker CPU total: **{(report['search_reported_cpu_seconds']+report['decode_reported_cpu_seconds'])/60:.2f} minutes**. All processes passed; no watchdog or memory-limit stops. Longest competition case was case 35 at 30.79 s including startup, closely followed by case 52.",
'',
'These are additional costs after the existing one-seed, cap-2,000 baseline search. The older successful baseline searches recorded approximately 2.75 CPU minutes, so this extra competition search (4.21 CPU minutes) is substantial relative to matching alone, despite modest wall time. The existing baseline decoder was not rerun. CPU counters exclude interpreter imports before measurement, final output serialization after measurement, parent orchestration, and killed-process costs (none here); development and this separate proof/report audit are additional. RSS is sampled at 0.5 s intervals and may miss short peaks.','',
'## Coverage and complete decoding','',
'| Case | Sweep baseline | Sweep with competition | Native-SLAP baseline | Native-SLAP with competition |','|---:|---:|---:|---:|---:|']
for r in rows:
 if any(v['new'] for v in r['comparisons'].values()):
  s=r['comparisons']['slap_sweep'];n=r['comparisons']['native_slap'];lines.append(f"| {r['case']} | {s['baseline']}/{s['total']} | {s['union']}/{s['total']} | {n['baseline']}/{n['total']} | {n['union']}/{n['total']} |")
lines+=['',
f"Every one of **{dec_counts['paths']:,} accepted full paths** was decoded completely in its baseline event window; {dec_counts['solver_queries']} solver queries sufficed, with most paths handled by exact invariance or reachable-image lower-bound certificates. There was no class-count cap. All 140 windows completed. No permutation-group or full-bijection enumeration was introduced.",
'',
f"The pass adds **{report['new_window_class_count']} event classes across {len(report['new_window_cases'])} cases** within those windows. Six classes in case 129 require decoding the compressed families and are absent from the repair representatives. These extra classes do not lower any case's previously found minimum event count. All {len(proofs)} new window classes passed independent exact family-membership queries, anchor checks, and raw event-class recomputation; each newly recovered reference-minimum representative also passed a replay of its local B collision and preservation of unaffected A assignments.",
'',
'Case 64 still misses two SLAP-sweep minimum classes at four events and five native-SLAP classes at seven events. Those are different reference minima, so their missing counts should not be combined. Full saved-family decoding does not prove that the pending priority search or omitted parent partitions cannot recover them.','',
'## Search limits and kernel finding','',
f"The run made {counts['growth_calls']:,} native proposal calls and {counts['completion_calls']:,} repair attempts, with cap 2,000 per native call, up to eight baseline parent partitions, and two collision levels. {sum(bool(r['pending_states']) for r in rows)} cases retained pending queued states; 13 native repair calls hit their branch caps. Baseline classes are retained by construction, so zero lost baseline classes is a property of union, not a test of replacing the original search.",
'',
f"Independent validation rejected **{counts['rejected_repair_graphs']} whole repair graphs** (out of {counts['completion_calls']:,} attempts), containing {counts['invalid_paths']} invalid paths across {sum(bool(r['counts'].get('rejected_repair_graphs')) for r in rows)} cases. They were neither expanded nor credited; otherwise-valid paths in the same rejected graphs were also excluded. Thus accepted-witness claims are sound, but this is still an experimental integration, not a production-ready complete search.",
'',
'A traced case-64 rejection from the 320-repair pilot reproduces in both Python growth and the native engine. Before extending atom 17, source atom 16 retains target pool {12,14} while anchored source atom 1 maps to target 3. The current representative uses 16→12. Extending 17 admits 16→14, making previously preserved source bond (1,16), WBO 1.10987949312885, map to target bond (3,14), WBO 0. The exact final family is UNSAT. This localizes the problem to retained support in a compressed pool during growth, rather than merely export or the later event thresholds. The safe repair needs to retain the earlier bond-support constraint when refining the pool; changing the event threshold alone would not repair it. No broad kernel change was made on the strength of this single trace.',
'',
'## Artifacts','',
'- `summary.json`: all per-case coverage, timing, proof checks, configuration, and source hashes.',
'- `viewer.html`: all 140 cases using the existing white R/P viewer style. New competition classes appear first, with baseline and both SLAP references available.',
'- `case*/full_decode.json`: complete accepted-path certificates in the same per-case windows as the baseline.',
'- `case*/takeover*.pkl.gz`: accepted compressed repair archives; `rejected*.json` records excluded graphs.',
'- `scripts/`: exact experiment, decoder, audit, and diagnostic sources. The production repository was not modified by this experiment.',
'- `../competition_optimization/invalid_growth_trace.json`: diagnostic Python growth trace; the companion probe script checks the corresponding family with the exact solver.','']
(OUT/'README.md').write_text('\n'.join(lines))
print(json.dumps({k:report[k] for k in ('comparisons','new_window_class_count','new_window_cases','best_count_improved_cases','proof_audit_wall_seconds','proof_audit_cpu_seconds')},indent=2));print('verified',len(proofs),'new window classes; viewer cases',len(documents))
