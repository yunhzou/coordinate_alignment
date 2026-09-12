"""Count saved representatives; do not enumerate their symmetry families."""
import sys,os,json,time,collections
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');RUN=ROOT/'outputs/competition140_budget128';OUT=RUN/'path_counts';CODE=ROOT/'work/aam-event-improved'
sys.path[:0]=[str(CODE/'src'),str(CODE/'bench')]
from fragment_competition_optimized import save

def child(case):
 from rxn_core.artifacts import read_aam_checkpoint
 from rxn_core.event_patterns import SignedEventIndex
 start=time.perf_counter();row=next(r for r in json.loads((RUN/'summary.json').read_text())['cases'] if r['case']==case);minimum=row['union_best']
 hist=collections.Counter();unique=collections.Counter();seen=set();idx=None;batch=[];archive_count=0
 def flush():
  if not batch:return
  scores=idx.counts(batch).sum(axis=1)
  for vector,score in zip(batch,scores):
   delta=int(score)-minimum;hist[delta]+=1
   if vector not in seen:seen.add(vector);unique[delta]+=1
  batch.clear()
 for file in (RUN/f'case{case}').glob('takeover*.pkl.gz'):
  a=read_aam_checkpoint(file);archive_count+=1
  if idx is None:idx=SignedEventIndex(a.problem)
  for path in a.graph.paths():
   m=path.mapping
   if len(m)!=idx.n:continue
   batch.append(tuple(m[r] for r in range(idx.n)))
   if len(batch)>=64:flush()
 flush();expected=json.loads((RUN/f'case{case}/full_decode.json').read_text())['counts'].get('paths',0)
 assert sum(hist.values())==expected,(case,sum(hist.values()),expected)
 d=dict(case=case,minimum=minimum,archives=archive_count,paths=sum(hist.values()),distinct_mappings=len(seen),by_excess=dict(sorted(hist.items())),distinct_mappings_by_excess=dict(sorted(unique.items())),wall_seconds=time.perf_counter()-start)
 save(OUT/f'case{case}.json',d)

def main():
 from run_event_campaign import MemoryGuard,bounded_process
 OUT.mkdir(exist_ok=True);guard=MemoryGuard(2048,6144,6144);start=time.perf_counter();env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
 def run(c):
  with (OUT/f'case{c}.log').open('w') as log:status,peak=bounded_process([sys.executable,str(Path(__file__).resolve()),str(c)],env,log,300,guard)
  assert status=='passed',(c,status)
  return json.loads((OUT/f'case{c}.json').read_text())
 with ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(run,range(140)))
 hist=collections.Counter();unique=collections.Counter()
 for r in rows:
  hist.update({int(k):v for k,v in r['by_excess'].items()});unique.update({int(k):v for k,v in r['distinct_mappings_by_excess'].items()})
 result=dict(cases=rows,paths=sum(hist.values()),distinct_mappings=sum(unique.values()),by_excess=dict(sorted(hist.items())),distinct_mappings_by_excess=dict(sorted(unique.items())),window_paths=sum(hist[k] for k in range(7)),window_distinct_mappings=sum(unique[k] for k in range(7)),wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024,scope='Scores of saved path representatives only, relative to each case\'s best known AAM event count. Unique mappings deduplicated within each case across all repairs. Equal representatives need not have identical compressed families. This is not an exhaustive family count in the expanded +6 window.')
 save(OUT/'summary.json',result);print(json.dumps({k:v for k,v in result.items() if k!='cases'}),flush=True)
if __name__=='__main__':
 if len(sys.argv)>1:child(int(sys.argv[1]))
 else:main()
