"""Decode every accepted repair family in the same event window as the baseline."""
import sys,os,json,time,collections
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');CODE=ROOT/'work/aam-event-improved';OUT=ROOT/'outputs/competition140_budget128'
sys.path[:0]=[str(CODE/'src'),str(CODE/'bench')]
from fragment_competition_optimized import save

def child(case):
 from rxn_core.artifacts import read_aam_checkpoint
 from rxn_core.event_patterns import SignedEventIndex,extract_path_events
 start=time.perf_counter();cpu=time.process_time();folder=OUT/f'case{case}'
 baseline=next(r for r in json.loads((ROOT/'outputs/holdout140_decoded_results/results.json').read_text())['cases'] if r['case']==case)
 window=int(baseline['max_events']);patterns={};rows=[];counts=collections.Counter();unfinished=[];idx=None
 archives=sorted(folder.glob('takeover*.pkl.gz'),key=lambda p:int(p.name.split('takeover')[1].split('.')[0]))
 last_saved=start
 for archive in archives:
  if time.perf_counter()-start>270:break
  a=read_aam_checkpoint(archive)
  if idx is None:idx=SignedEventIndex(a.problem)
  count=0;complete=True
  for path in a.graph.paths():
   if len(path.mapping)!=idx.n:continue
   if time.perf_counter()-start>270:complete=False;break
   def record(p):
    patterns.setdefault(p['id'],dict(p,archive=str(archive),terminal=path.terminal,transitions=list(path.transitions)))
   r=extract_path_events(path,a.problem,idx,max_events=window,seconds=min(5,max(.01,270-(time.perf_counter()-start))),max_patterns=None,on_pattern=record)
   counts['paths']+=1;counts['solver_queries']+=r['solver_queries'];counts[r['reason']]+=1;count+=1;complete &= r['complete']
   if not r['complete']:unfinished.append(dict(archive=str(archive),terminal=path.terminal,transitions=list(path.transitions),reason=r['reason']))
   if time.perf_counter()-last_saved>2:
    save(folder/'decode_progress.json',dict(case=case,max_events=window,counts=dict(counts),patterns=patterns,completed_archives=len(rows),unfinished=unfinished));last_saved=time.perf_counter()
  rows.append(dict(archive=str(archive),paths=count,complete=complete))
 complete=len(rows)==len(archives) and all(r['complete'] for r in rows)
 d=dict(case=case,max_events=window,complete_saved_windows=complete,archives=rows,expected_archives=len(archives),patterns=patterns,counts=dict(counts),unfinished=unfinished,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,scope='All full paths from accepted repair archives, within the baseline event window. Does not cover rejected repairs, capped/omitted branches, or pending priority states.')
 save(folder/'full_decode.json',d);print(json.dumps({k:v for k,v in d.items() if k not in ('patterns','archives','unfinished')}),flush=True)

def main():
 from run_event_campaign import MemoryGuard,bounded_process
 guard=MemoryGuard(3072,8192,6144);env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1');start=time.perf_counter()
 def run(c):
  folder=OUT/f'case{c}'
  if (folder/'full_decode.json').exists() and (folder/'decode_execution.json').exists():return json.loads((folder/'decode_execution.json').read_text())
  t=time.perf_counter()
  with (folder/'decode.log').open('w') as log:status,peak=bounded_process([sys.executable,str(Path(__file__).resolve()),str(c)],env,log,300,guard)
  row=dict(case=c,status=status,peak_mib=peak,wall_seconds=time.perf_counter()-t);save(folder/'decode_execution.json',row);print(json.dumps(row),flush=True);return row
 with ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(run,range(140)))
 save(OUT/'decode_execution.json',dict(jobs=rows,wall_seconds=time.perf_counter()-start,peak_total_mib=guard.peak_total_kib/1024,workers=4,watchdog=300))
if __name__=='__main__':
 if len(sys.argv)>1:child(int(sys.argv[1]))
 else:main()
