import sys,json,gzip,time,os,hashlib
from pathlib import Path
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');CODE=ROOT/'work/aam-event-improved';OUT=ROOT/'outputs/final_fragment_dedup/lossless_pilot/case129'
sys.path[:0]=[str(CODE/'src'),str(CODE/'bench')]
def save(path, value):
 path.parent.mkdir(parents=True, exist_ok=True)
 temporary = path.with_suffix(path.suffix + '.tmp')
 temporary.write_text(json.dumps(value, indent=2) + '\n')
 temporary.replace(path)

def child():
 from rxn_core.artifacts import read_aam_checkpoint
 from rxn_core.final_branches import FinalFamily
 from rxn_core.event_patterns import SignedEventIndex,extract_path_events
 start=time.perf_counter();cpu=time.process_time();r=json.loads((OUT/'result.json').read_text())
 with gzip.open(OUT/'final_branches.json.gz','rt') as f:d=json.load(f)
 archive=ROOT/'outputs/holdout140_cap2000_seed1/case129/cuts/aam.pkl.gz';problem=read_aam_checkpoint(archive).problem;idx=SignedEventIndex(problem);results=[]
 for missing in r['unfinished']:
  fid=missing['family'];family=FinalFamily.from_record(d['families'][fid])
  q=extract_path_events(family.as_path(problem),problem,idx,max_events=r['window'],seconds=20,max_patterns=None)
  for p in q['patterns']:r['patterns'].setdefault(p['id'],dict(p,family=fid))
  results.append(dict(family=fid,**{k:v for k,v in q.items() if k!='patterns'}))
 complete=all(x['complete'] for x in results)
 report=dict(complete=complete,results=results,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu)
 save(OUT/'resume.json',report)
 if complete:
  r['initial_unfinished']=r['unfinished'];r['unfinished']=[];r['complete_saved_window']=True;r['resume']=report
  save(OUT/'result.json',r)
 print(json.dumps(report),flush=True)
if __name__=='__main__':
 if '--child' in sys.argv:child()
 else:
  from run_event_campaign import MemoryGuard,bounded_process
  with (OUT/'resume.log').open('w') as log:status,peak=bounded_process([sys.executable,str(Path(__file__).resolve()),'--child'],dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),log,120,MemoryGuard(2048,4096,6144))
  print(status,peak)
