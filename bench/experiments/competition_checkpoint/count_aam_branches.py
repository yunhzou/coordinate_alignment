"""Count exactly the literal relation keys used by AAMResult.branches.

Cache fragment keys and intern them to count without retaining every SearchPath.
No symmetry or atom-mapping expansion.
"""
import sys,os,json,time,collections,gc
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');RUN=ROOT/'outputs/competition140_budget128';OUT=RUN/'branch_counts';CODE=ROOT/'work/aam-event-improved'
sys.path[:0]=[str(CODE/'src'),str(CODE/'bench')]
from fragment_competition_optimized import save

def child(case):
 from rxn_core.artifacts import read_aam_checkpoint
 from rxn_core.search_graph import frozen_value
 start=time.perf_counter();interner={};base=set();repairs=set();counts=collections.Counter()
 def collect(file,keys,label):
  a=read_aam_checkpoint(file);graph=a.graph;edgekeys={};local=set()
  for path in graph.paths():
   counts[label+'_paths']+=1
   fragments=[]
   for eid in path.transitions:
    edge=graph.transitions[eid]
    if edge.match is None:continue
    if eid not in edgekeys:
     symmetry={k:v for k,v in edge.match['symmetry'].items() if k not in {'multiplicity','automorph_group_source'}}
     f=(tuple(edge.match['fragment']),edge.preserved_bonds,frozen_value(symmetry),frozen_value(edge.match['deferred_edges']))
     edgekeys[eid]=interner.setdefault(f,len(interner))
    fragments.append(edgekeys[eid])
   key=(tuple(sorted(graph.states[path.terminal].mapping)),tuple(fragments));keys.add(key);local.add(key)
  counts[label+'_sum_local_branches']+=len(local)
  # Check against the public API on small archives, where no large memory cost arises.
  if len(graph.terminals)<=10 and label=='baseline':assert len(local)==len(a.branches)
  return len(graph.terminals)
 input_path=ROOT/f'outputs/holdout140_cap2000_seed1/case{case}/cuts/aam.pkl.gz' if case!=25 else ROOT/'outputs/holdout140_case25_retry/case25/cuts/aam.pkl.gz'
 terminals=collect(input_path,base,'baseline');gc.collect()
 for file in (RUN/f'case{case}').glob('takeover*.pkl.gz'):collect(file,repairs,'repair')
 d=dict(case=case,baseline_branches=len(base),repair_branches=len(repairs),combined_branches=len(base|repairs),repair_new_branches=len(repairs-base),baseline_terminals=terminals,counts=dict(counts),wall_seconds=time.perf_counter()-start)
 save(OUT/f'case{case}.json',d)

def main():
 from run_event_campaign import MemoryGuard,bounded_process
 OUT.mkdir(exist_ok=True);guard=MemoryGuard(3072,8192,6144);start=time.perf_counter();env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 def run(c):
  with (OUT/f'case{c}.log').open('w') as log:status,peak=bounded_process([sys.executable,str(Path(__file__).resolve()),str(c)],env,log,300,guard)
  if status!='passed':return dict(case=c,status=status,peak_mib=peak)
  return dict(json.loads((OUT/f'case{c}.json').read_text()),status=status,peak_mib=peak)
 with ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(run,range(140)))
 counts={k:sum(r.get(k,0) for r in rows) for k in ('baseline_branches','repair_branches','combined_branches','repair_new_branches','baseline_terminals')}
 result=dict(cases=rows,**counts,complete=all(r['status']=='passed' for r in rows),wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024,definition='Exact AAMResult.branches relation-key equality, evaluated within each reaction across baseline and all accepted repairs. No decoded mappings or event classes used.')
 save(OUT/'summary.json',result);print(json.dumps({k:v for k,v in result.items() if k!='cases'}),flush=True)
if __name__=='__main__':
 if len(sys.argv)>1:child(int(sys.argv[1]))
 else:main()
