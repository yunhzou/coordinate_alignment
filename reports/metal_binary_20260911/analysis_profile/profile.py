from pathlib import Path
import collections, functools, json, sys, time
import numpy as np
run=Path('/project/yunhengzou/coordinate_alignment/aam_benchmarks/metal_binary_20260911')
sys.path.insert(0,str(run))
import metal_binary_events as events
from rxn_core.artifacts import read_aam_checkpoint
index=int(sys.argv[1]);variant='original';folder=run/f'runs/holdout/{variant}/results/holdout/{index}/R_to_P/original'
raw=json.loads((run/f'raw/holdout/{index}/input.json').read_text());expected=json.loads((folder/'raw_evaluation.json').read_text())
start=time.process_time();aam=read_aam_checkpoint(folder/'cuts/aam.pkl.gz');loading=time.process_time()-start
profile=collections.defaultdict(lambda:dict(cpu=0.,calls=0))
def timed(name,fn):
 @functools.wraps(fn)
 def wrapper(*a,**k):
  t=time.process_time()
  try:return fn(*a,**k)
  finally:profile[name]['cpu']+=time.process_time()-t;profile[name]['calls']+=1
 return wrapper
for name in ('scalar_events','colored_graph'):
 setattr(events,name,timed(name,getattr(events,name)))
events.pynauty.certificate=timed('nauty_certificate',events.pynauty.certificate)
events.DeltaPatterns.counts=timed('vector_event_counts',events.DeltaPatterns.counts)
events.DeltaPatterns.describe=timed('describe_inclusive',events.DeltaPatterns.describe)
features=functools.cached_property(timed('score_response_features',events.DeltaPatterns.features.func));features.__set_name__(events.DeltaPatterns,'features');events.DeltaPatterns.features=features
start=time.process_time();canonical=events.DeltaPatterns(raw)
t=time.process_time();vectors=sorted({tuple(dict(aam.graph.states[t].mapping)[a] for a in range(canonical.n)) for t in aam.graph.terminals if len(aam.graph.states[t].mapping)==canonical.n});extraction=time.process_time()-t
best=None;patterns={};histogram=collections.Counter();target=None
if index in (11,64,101):target=canonical.describe(json.loads((run/'historical_targets.json').read_text())[str(index)]['pattern']['mapping'])
for offset in range(0,len(vectors),128):
 batch=vectors[offset:offset+128];scores=canonical.counts(batch).sum(axis=1)
 for vector,score in zip(batch,scores):
  score=int(score);histogram[score]+=1
  if best is None or score<best:best=score;patterns={}
  if score==best or target is not None and score==target['total']:
   pattern=canonical.describe(vector)
   if score==best:patterns.setdefault(pattern['id'],pattern)
assert {str(k):v for k,v in histogram.items()}==expected['event_histogram']
assert set(patterns)==set(expected['minimum_patterns'])
result=dict(index=index,variant=variant,atoms=canonical.n,representatives=len(vectors),loading_cpu=loading,analysis_cpu=time.process_time()-start,extraction_cpu=extraction,profile=dict(profile),frozen_analysis_cpu=expected['analysis_cpu'],verified_same_histogram_and_minimum_patterns=True)
output=run/f'analysis_profile/{index}.json';output.parent.mkdir(exist_ok=True);output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
