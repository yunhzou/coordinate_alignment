"""Bounded saved-result final-fragment reconstruction and event validation."""
import sys,os,json,gzip,time,collections,hashlib,gc
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');CODE=ROOT/'work/aam-event-improved';BASE=ROOT/'outputs/competition140_budget128';OUT=Path(os.environ.get('FINAL_DEDUP_OUTPUT',str(ROOT/'outputs/final_fragment_dedup/pilot')))
sys.path[:0]=[str(CODE/'src'),str(CODE/'bench')]
def save(path, value):
 path.parent.mkdir(parents=True, exist_ok=True)
 temporary = path.with_suffix(path.suffix + '.tmp')
 temporary.write_text(json.dumps(value, indent=2) + '\n')
 temporary.replace(path)
CASES=tuple(map(int,os.environ.get('FINAL_DEDUP_CASES','6,11,59,64,101,129').split(',')))
MODE=os.environ.get('FINAL_DEDUP_MODE','intrinsic')
DECODE=os.environ.get('FINAL_DEDUP_DECODE','1')=='1'
WORKERS=int(os.environ.get('FINAL_DEDUP_WORKERS','2'))

def child(case):
 from rxn_core.artifacts import read_aam_checkpoint
 from rxn_core.final_branches import FinalBranchCatalogue
 from rebuild_intrinsic import rebuild_intrinsic_catalogue
 from rxn_core.event_patterns import SignedEventIndex,extract_path_events
 from rxn_core.family_scoring import validate_representative
 start=time.perf_counter();cpu=time.process_time();folder=OUT/f'case{case}';folder.mkdir(parents=True,exist_ok=True)
 row=next(r for r in json.loads((BASE/'summary.json').read_text())['cases'] if r['case']==case)
 previous=next(r for r in json.loads((BASE/'branch_counts/summary.json').read_text())['cases'] if r['case']==case)
 inp=ROOT/f'outputs/holdout140_cap2000_seed1/case{case}/cuts/aam.pkl.gz' if case!=25 else ROOT/'outputs/holdout140_case25_retry/case25/cuts/aam.pkl.gz'
 a=read_aam_checkpoint(inp);problem=a.problem;catalogue=FinalBranchCatalogue(problem);catalogue.add_graph(a.graph,str(inp));baseline_branches=len(catalogue.branches);del a;gc.collect()
 archives=[inp]
 for archive in (BASE/f'case{case}').glob('takeover*.pkl.gz'):
  if time.perf_counter()-start>260:raise TimeoutError('catalogue budget')
  a=read_aam_checkpoint(archive);catalogue.add_graph(a.graph,str(archive));archives.append(archive);del a
 assert catalogue.path_count==previous['counts']['baseline_paths']+previous['counts']['repair_paths']
 saved_stats=dict(input_paths=catalogue.path_count,old_literal_branches=previous['combined_branches'],baseline_pair_combinations=baseline_branches,pair_combinations=len(catalogue.branches),flat_saved_families=len(catalogue.families))
 t=time.perf_counter()
 if MODE=='lossless':
  rebuilt=catalogue;rebuilt.rebuild_all_symmetries()
 else:rebuilt=rebuild_intrinsic_catalogue(catalogue)
 rebuild_seconds=time.perf_counter()-t
 metrics=dict(case=case,mode=MODE,**saved_stats,final_branches=rebuilt.branch_count,independent_pair_branches=len(rebuilt.branches),coupled_branches=len(rebuilt.coupled_branches),rebuilt_families=len(rebuilt.families),intrinsic_fragment_groups=len(rebuilt.symmetries),nontrivial_intrinsic_groups=sum(bool(g['generators']) for g in rebuilt.symmetries.values()),reconstruction=rebuilt.reconstruction,rebuild_seconds=rebuild_seconds,predecode_wall_seconds=time.perf_counter()-start)
 save(folder/'catalogue_summary.json',metrics)
 with gzip.open(folder/'final_branches.json.gz','wt') as f:json.dump(rebuilt.to_record(),f,separators=(',',':'))
 # Validate that every new orbit starts from a valid saved representative.
 validation_started=time.perf_counter()
 for family in rebuilt.families:family.validate_representative(problem)
 metrics['all_rebuilt_representatives_valid']=True;metrics['representative_validation_seconds']=time.perf_counter()-validation_started
 patterns={};reasons=collections.Counter();unfinished=[];done=0;last=time.perf_counter();window=row['window'];idx=SignedEventIndex(problem)
 if DECODE:
  for fid,family in enumerate(rebuilt.families):
   if time.perf_counter()-start>270:break
   path=family.as_path(problem)
   def record(p):patterns.setdefault(p['id'],dict(p,family=fid))
   d=extract_path_events(path,problem,idx,max_events=window,seconds=min(3,max(.01,270-(time.perf_counter()-start))),max_patterns=None,on_pattern=record)
   done+=1;reasons[d['reason']]+=1;reasons['solver_queries']+=d['solver_queries']
   if not d['complete']:unfinished.append(dict(family=fid,reason=d['reason']))
   if time.perf_counter()-last>2:
    save(folder/'decode_progress.json',dict(**metrics,decoded_families=done,patterns=patterns,unfinished=unfinished,reasons=dict(reasons)));last=time.perf_counter()
 original=next(r for r in json.loads((ROOT/'outputs/holdout140_decoded_results/results.json').read_text())['cases'] if r['case']==case)
 base=json.loads(Path(original['source']).read_text());old=set(base['patterns'])|set(json.loads((BASE/f'case{case}/full_decode.json').read_text())['patterns'])
 old={k for k in old if (base['patterns'].get(k) or json.loads((BASE/f'case{case}/full_decode.json').read_text())['patterns'][k])['total']<=window}
 complete=DECODE and done==len(rebuilt.families) and not unfinished
 metrics.update(decode_requested=DECODE,complete_saved_window=complete,decoded_families=done,decode_reasons=dict(reasons),unfinished=unfinished,window=window,patterns=patterns,previous_window_classes=len(old),retained_previous_classes=len(old&patterns.keys()),unrecovered_previous_classes=sorted(old-patterns.keys()) if DECODE else None,additional_window_classes=sorted(patterns.keys()-old),wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,source_sha256=hashlib.sha256((CODE/'src/rxn_core/final_branches.py').read_bytes()).hexdigest())
 save(folder/'result.json',metrics);print(json.dumps({k:v for k,v in metrics.items() if k not in ('patterns','unrecovered_previous_classes','additional_window_classes','unfinished')}),flush=True)

def main():
 from run_event_campaign import MemoryGuard,bounded_process
 OUT.mkdir(parents=True,exist_ok=True);guard=MemoryGuard(3072,8192,6144);start=time.perf_counter();env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
 def run(c):
  folder=OUT/f'case{c}';folder.mkdir(exist_ok=True)
  with (folder/'run.log').open('w') as log:status,peak=bounded_process([sys.executable,str(Path(__file__).resolve()),str(c)],env,log,300,guard)
  row=dict(case=c,status=status,peak_mib=peak);save(folder/'execution.json',row);print(json.dumps(row),flush=True);return row
 with ThreadPoolExecutor(max_workers=WORKERS) as pool:rows=list(pool.map(run,CASES))
 save(OUT/'execution.json',dict(jobs=rows,wall_seconds=time.perf_counter()-start,peak_mib=guard.peak_total_kib/1024,workers=WORKERS,watchdog=300))
if __name__=='__main__':
 if len(sys.argv)>1:child(int(sys.argv[1]))
 else:main()
