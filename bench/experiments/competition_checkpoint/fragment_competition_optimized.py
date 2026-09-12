"""Reference-blind, bounded local B-priority takeover experiment.

Start from selected saved terminal witnesses, retain their compressed partition,
let one existing island grow with priority, release conflicting source/target
assignments, and regrow holes with all unaffected assignments anchored. Native
fragment placements are used directly; no atom-permutation expansion. Only the
new placement representatives are screened in this pilot; no negative family
coverage claims are made. Original A-priority results remain in the union.
"""
import sys,os,json,time,hashlib,collections
from pathlib import Path
from dataclasses import asdict
from concurrent.futures import ThreadPoolExecutor
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');CODE=ROOT/'work/aam-event-improved';OUT=ROOT/'outputs/fragment_competition_experiment'
# Fixed validated collision policy. Budgets are independent of reference results.
ORIGINAL_CONTACTS=UNMAPPED_COMPLETION=RELEASE_BOUNDARIES=False
VERIFIED_PILOT=DEPENDENT_REGROWTH=RELOCATE_B=LOCAL_COLLISIONS=True
OUT=Path(os.environ.get('COMPETITION_OUTPUT',str(ROOT/'outputs/competition_optimization/pilot320')))
BUDGET=int(os.environ.get('COMPETITION_BUDGET','320'))
WORKERS=int(os.environ.get('COMPETITION_WORKERS','2'))
CASES=tuple(map(int,os.environ.get('COMPETITION_CASES','64,11,101,6').split(',')))
sys.path[:0]=[str(CODE/'src'),str(CODE/'bench')]
def save(p,d):
 p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');tmp.replace(p)
def takeover_plan(old,owners,owner,new,*,allow_b_relocation=False):
 """Exact pairs that survive B's priority; reject malformed proposals."""
 old=dict(old);new=dict(new);owners=dict(owners)
 if len(set(new.values()))!=len(new):raise ValueError('noninjective takeover')
 base={r:p for r,p in old.items() if owners[r]==owner}
 if not allow_b_relocation and not all(new.get(r,p)==p for r,p in base.items()):raise ValueError('B prefix was changed')
 bmap={} if allow_b_relocation else dict(base);bmap.update(new);used=set(bmap.values())
 retained={r:p for r,p in old.items() if r not in bmap and p not in used}
 anchors=dict(retained);anchors.update(bmap)
 holes=sorted(set(old)-anchors.keys())
 eaten=sorted(r for r in bmap if owners[r]!=owner)
 displaced=sorted(r for r in old if r not in bmap and old[r] in used)
 changed=sorted(r for r in bmap if bmap[r]!=old[r])
 return dict(anchors=anchors,retained=retained,bmap=bmap,holes=holes,eaten=eaten,displaced=displaced,changed=changed,
             touched_owners=sorted({owners[r] for r in eaten+displaced}))
def child(case):
 import numpy as np
 from rxn_core import AAMProblem,MolecularEndpoint,AAMSearchConfig
 from rxn_core.artifacts import read_aam_checkpoint,write_aam_checkpoint
 from rxn_core.domain import AAMResult,AAMSearchMetrics
 from rxn_core.frag import build_graph
 from rxn_core.fragment import match_fragment,FragmentMatchConfig,FragmentMatchContext
 from rxn_core.matcher import _nauty_orbits
 from rxn_core.native_search import find_islands_native
 from rxn_core.search_symmetry import finalize_graph_symmetry
 folder=OUT/f'case{case}';folder.mkdir(parents=True,exist_ok=True)
 archive=ROOT/f'outputs/holdout140_cap2000_seed1/case{case}/cuts/aam.pkl.gz'
 if case==25:archive=ROOT/'outputs/holdout140_case25_retry/case25/cuts/aam.pkl.gz'
 total_start=time.perf_counter();total_cpu=time.process_time()
 a=read_aam_checkpoint(archive);problem=a.problem;cfg=a.config
 from rxn_core.event_patterns import SignedEventIndex
 idx=SignedEventIndex(problem);start=time.perf_counter();cpu=time.process_time();deadline=start+270
 # Search-side selection uses only our terminal scores/partitions. Saved SLAP
 # labels/results are not read until the search loop has ended.
 terminals=[];batch=[]
 def flush_parents():
  if not batch:return
  scores=idx.counts([v for t,v,part in batch]).sum(axis=1)
  terminals.extend((int(score),t,part) for (t,v,part),score in zip(batch,scores));batch.clear()
 for t in a.graph.terminals:
  st=a.graph.states[t]
  if len(st.mapping)!=idx.n:continue
  groups=collections.defaultdict(list)
  for atom,g in st.islands:groups[g].append(atom)
  if len(groups)<2:continue
  partition=tuple(sorted(tuple(v) for v in groups.values()))
  m=dict(st.mapping);batch.append((t,[m[r] for r in range(idx.n)],partition))
  if len(batch)>=64:flush_parents()
 flush_parents()
 terminals.sort();selected=[];partitions=set()
 for score,t,part in terminals:
  if part in partitions:continue
  partitions.add(part);selected.append(t)
  if len(selected)>=(8 if LOCAL_COLLISIONS else 24):break
 target=build_graph(problem.product.elements,problem.product.wbo,bond_cut=cfg.graph_floor)
 po=_nauty_orbits(target,wbo_tol=cfg.iso_tolerance)
 sources={};offers_cache={};completion_cache={};queue=collections.deque();seen_states=set();patterns={};offers=[];counts=collections.Counter();scored_vectors=set();release_cache={};checkpoint_at=start
 for t in selected:
  state=a.graph.states[t];context=a.graph.contexts[state.context]
  queue.append(dict(mapping=dict(state.mapping),owners=dict(state.islands),cuts=tuple(context.cuts),order=tuple(context.seed_order),parent=t,depth=0))
 def prepared(cuts):
  if cuts not in sources:
   g=build_graph(problem.reactant.elements,problem.reactant.wbo,bond_cut=cfg.graph_floor);g.remove_edges_from(cuts);sources[cuts]=g
  return sources[cuts]
 def score_store(mapping,proof):
  vector=tuple(mapping[r] for r in range(idx.n))
  if vector in scored_vectors:counts['duplicate_scoring_avoided']+=1;return
  scored_vectors.add(vector);p=idx.describe(vector)
  if p['id'] not in patterns:patterns[p['id']]=dict(p,provenance=proof)
 def checkpoint(force=False):
  nonlocal checkpoint_at
  now=time.perf_counter()
  if not force and now-checkpoint_at<2:return
  checkpoint_at=now
  save(folder/'progress.json',dict(case=case,counts=dict(counts),selected_parents=selected,pattern_count=len(patterns),pending=len(queue),wall_seconds=now-start))
 # Depth is limited to two successive local priority changes. A fixed global
 # operation budget bounds alternative generation independently of reference hits.
 while queue and counts['growth_calls']<BUDGET and counts['completion_calls']<BUDGET and time.perf_counter()<deadline:
  state=queue.popleft();old=state['mapping'];owners=state['owners'];cuts=state['cuts'];order=state['order'];source=prepared(cuts);offer_source=prepared(()) if ORIGINAL_CONTACTS or LOCAL_COLLISIONS else source
  statekey=(tuple(sorted(old.items())),tuple(sorted(owners.items())),cuts)
  if statekey in seen_states:continue
  seen_states.add(statekey);groups=collections.defaultdict(list)
  for r,g in owners.items():groups[g].append(r)
  # Visit smaller established fragments first: this is the challenger B.
  tasks=[]
  for owner,atoms in sorted(groups.items(),key=lambda item:(len(item[1]),min(item[1]))):
   if len(atoms)==idx.n:continue
   boundary=[r for r in atoms if any(owners[s]!=owner for s in offer_source.neighbors(r))]
   if not boundary:continue
   # One deterministic contact atom per B; no reference/event selection.
   seed=min(boundary,key=lambda r:order.index(r) if r in order else r)
   if LOCAL_COLLISIONS:
    for contact in sorted({n for r in boundary for n in offer_source.neighbors(r) if owners[n]!=owner}):
     contact_seed=min(r for r in boundary if offer_source.has_edge(r,contact))
     tasks.append((owner,atoms,contact_seed,contact))
   else:tasks.append((owner,atoms,seed,None))
  for owner,atoms,seed,contact in tasks:
   if counts['growth_calls']>=BUDGET or counts['completion_calls']>=BUDGET or time.perf_counter()>=deadline:break
   base={r:old[r] for r in atoms};key=(cuts,tuple(sorted(base.items())),seed,contact)
   if key in offers_cache:result=offers_cache[key];counts['growth_cache_hits']+=1
   else:
    counts['growth_calls']+=1
    proposal_source=offer_source.subgraph(set(atoms)|{contact}).copy() if LOCAL_COLLISIONS else offer_source
    result=match_fragment(proposal_source,target,seed=seed,context=FragmentMatchContext({} if RELOCATE_B else base,{} if RELOCATE_B else {r:0 for r in atoms},(),target_orbits=po),config=FragmentMatchConfig(graph_floor=cfg.graph_floor,iso_tolerance=cfg.iso_tolerance,branch_limit=2000,allow_mapped_seed=not RELOCATE_B))
    offers_cache[key]=result
   if result.capped:counts['growth_capped']+=1;continue
   for placement in result.matches:
    plan=takeover_plan(old,owners,owner,dict(placement),allow_b_relocation=RELOCATE_B);counts['offers']+=1
    if RELOCATE_B and contact not in placement:counts['no_contact_reached']+=1;continue
    if not plan['touched_owners'] or not plan['changed']:counts['no_changed_takeover']+=1;continue
    if counts['completion_calls']>=BUDGET or time.perf_counter()>=deadline:break
    if DEPENDENT_REGROWTH:
     import networkx as nx
     affected=set(plan['changed'])|set(plan['displaced']);release=set()
     for side,(graph0,centers,inverse) in enumerate(((prepared(()),affected,{r:r for r in old}),(target,{old[r] for r in affected}|{plan['bmap'][r] for r in affected if r in plan['bmap']},{p:r for r,p in old.items()}))):
      component_key=(side,frozenset(centers))
      if component_key not in release_cache:
       released_nodes=set()
       for component in nx.connected_components(nx.subgraph_view(graph0,filter_node=lambda n:n not in centers)):
        if len(component)<=8 and any(v in centers for n in component for v in graph0.neighbors(n)):released_nodes.update(component)
       release_cache[component_key]=released_nodes
      else:counts['dependent_cache_hits']+=1
      release.update(inverse[n] for n in release_cache[component_key] if n in inverse)
     release-=plan['bmap'].keys()
     plan['released_dependents']=sorted(release&plan['retained'].keys())
     for r in release:plan['anchors'].pop(r,None);plan['retained'].pop(r,None)
     plan['holes']=sorted(set(old)-plan['anchors'].keys())
    anchors=plan['anchors'];completion_key=(cuts,tuple(sorted(anchors.items())))
    if completion_key in completion_cache:counts['completion_cache_hits']+=1;continue
    counts['completion_calls']+=1
    released_edges=tuple(sorted((min(r,t),max(r,t)) for r,t in source.edges if r in anchors and t in anchors and (problem.product.wbo[anchors[r],anchors[t]]<cfg.graph_floor or abs(problem.reactant.wbo[r,t]-problem.product.wbo[anchors[r],anchors[t]])>cfg.iso_tolerance+1e-9))) if RELEASE_BOUNDARIES else ()
    effective_cuts=tuple(sorted(set(cuts)|set(released_edges)))
    completion_source=prepared(effective_cuts)
    # New graph starts from revised constraints. Keeping the ancestor A path
    # would incorrectly reimpose the very assignments/bonds B has displaced.
    completion_order=tuple(r for r in order if r not in anchors) if UNMAPPED_COMPLETION else order
    graph=find_islands_native(completion_source,target,completion_order,graph_floor=cfg.graph_floor,iso_tol=cfg.iso_tolerance,max_branches=2000,p_orbits=po,cuts=effective_cuts,anchor_map=anchors)
    graph,_=finalize_graph_symmetry(graph,target,iso_tolerance=cfg.iso_tolerance)
    if RELEASE_BOUNDARIES or VERIFIED_PILOT:
     from rxn_core.family_scoring import validate_representative
     invalid=[]
     for path in graph.paths():
      if len(path.mapping)!=idx.n:continue
      try:validate_representative(path,problem);counts['validated_paths']+=1
      except AssertionError:
       if not VERIFIED_PILOT:raise
       invalid.append(path.terminal)
     if invalid:
      counts['rejected_repair_graphs']+=1;counts['invalid_paths']+=len(invalid)
      save(folder/f'rejected{counts["completion_calls"]}.json',dict(terminal_ids=invalid,anchors=sorted(anchors.items()),cuts=effective_cuts,reason='preserved_fragment_constraints_fail_independent_validation',scope='Rejected experimental repair, never counted or expanded'))
      completion_cache[completion_key]=True
      continue
    completion_cache[completion_key]=True
    slot=counts['completion_calls'];out_archive=folder/f'takeover{slot}.pkl.gz'
    repaircfg=AAMSearchConfig(**dict(asdict(cfg),anchors=tuple(sorted(anchors.items()))))
    out_aam=AAMResult(problem,repaircfg,graph,AAMSearchMetrics.from_record({},0));write_aam_checkpoint(out_aam,out_archive)
    proof=dict(parent=state['parent'],depth=state['depth']+1,seed=seed,contact=contact,B_atoms=atoms,B_prefix=sorted(base.items()),B_placement=sorted(dict(placement).items()),cuts=cuts,released_preservation_edges=released_edges,effective_cuts=effective_cuts,
               eaten=plan['eaten'],displaced=plan['displaced'],holes=plan['holes'],released_dependents=plan.get('released_dependents',[]),changed=plan['changed'],retained_A_pairs=sorted(plan['retained'].items()),archive=str(out_archive),capped=graph.capped)
    counts['completion_capped']+=int(graph.capped);full=0
    for t in graph.terminals:
     m=dict(graph.states[t].mapping)
     if len(m)!=idx.n:continue
     assert all(m[r]==p for r,p in anchors.items());assert len(set(m.values()))==idx.n
     assert all(problem.reactant.elements[r]==problem.product.elements[p] for r,p in m.items())
     full+=1;counts['full_witnesses']+=1;score_store(m,dict(proof,terminal=t))
     if state['depth']<1 and len(queue)<512:
      # B is the new dominant island. Preserve A's residual partition; each
      # regrown hole starts as a distinct island for the next encounter.
      newowners={r:owners[r] for r in plan['retained']};fresh=max(owners.values())+1
      newowners.update({r:fresh for r in plan['bmap']})
      for j,r in enumerate(plan['holes'],fresh+1):newowners[r]=j
      queue.append(dict(mapping=m,owners=newowners,cuts=effective_cuts,order=order,parent=str(out_archive)+f':{t}',depth=state['depth']+1))
    offers.append(dict(proof,full_witnesses=full));checkpoint()
 checkpoint(force=True)
 search_wall=time.perf_counter()-start;search_cpu=time.process_time()-cpu
 # Evaluation begins only now.
 br=next(r for r in json.loads((ROOT/'outputs/holdout140_decoded_results/results.json').read_text())['cases'] if r['case']==case)
 baseline=json.loads(Path(br['source']).read_text());wanted=set(baseline['slap']['slap_sweep']['minimum_ids']);oldids=set(baseline['patterns']);union=oldids|patterns.keys()
 result=dict(case=case,original_contacts=ORIGINAL_CONTACTS,local_collisions=LOCAL_COLLISIONS,relocating_challenger=RELOCATE_B,dependent_regrowth=DEPENDENT_REGROWTH,released_boundaries=RELEASE_BOUNDARIES,unmapped_completion=UNMAPPED_COMPLETION,verified_pilot=VERIFIED_PILOT,operation_budget=BUDGET,config=asdict(cfg),selected_parents=selected,parent_selection=f'up to {8 if LOCAL_COLLISIONS else 24} lowest-event distinct source partitions; no SLAP input',counts=dict(counts),patterns=patterns,offers=offers,baseline_recovered=sorted(wanted&oldids),recovered=sorted(wanted&union),newly_recovered=sorted((wanted-oldids)&patterns.keys()),missing=sorted(wanted-union),total=len(wanted),new_event_classes=sorted(patterns.keys()-oldids),pending=len(queue),wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,total_wall_seconds=time.perf_counter()-total_start,total_cpu_seconds=time.process_time()-total_cpu,search_wall_seconds=search_wall,search_cpu_seconds=search_cpu,archive=str(archive),archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),scope='Observed terminal witnesses from bounded local takeover; baseline preserved by union. No negative completeness claim for takeover families or omitted parents/contacts. Search reads no SLAP labels.')
 save(folder/'result.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ('patterns','offers','new_event_classes','selected_parents')}),flush=True)
def main():
 from run_event_campaign import MemoryGuard,bounded_process
 guard=MemoryGuard(3072,8192,6144);env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1')
 def run(case):
  folder=OUT/f'case{case}';folder.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
  if (folder/'result.json').exists() and (folder/'execution.json').exists():
   previous=json.loads((folder/'execution.json').read_text())
   if previous['status']=='passed':return previous
  with (folder/'run.log').open('w') as f:status,peak=bounded_process([sys.executable,str(Path(__file__).resolve()),str(case)],env,f,300,guard)
  d=dict(case=case,status=status,peak_mib=peak,wall_seconds=time.perf_counter()-start);save(folder/'execution.json',d);print(json.dumps(d),flush=True);return d
 start=time.perf_counter()
 with ThreadPoolExecutor(max_workers=WORKERS) as pool:rows=list(pool.map(run,CASES))
 save(OUT/'execution.json',dict(jobs=rows,wall_seconds=time.perf_counter()-start,peak_total_mib=guard.peak_total_kib/1024,workers=WORKERS,watchdog=300,operation_budget=BUDGET,cases=CASES))
if __name__=='__main__':
 if len(sys.argv)>1:child(int(sys.argv[1]))
 else:main()
