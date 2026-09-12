"""Upgrade and validate flat exports without materializing an entire JSON archive."""
import sys,json,gzip,time,os
from pathlib import Path
ROOT=Path('/Users/yunhengz/Documents/Codex/2026-09-11/aa');CODE=ROOT/'work/aam-event-improved';RUN=ROOT/'outputs/final_fragment_dedup/lossless140'
sys.path[:0]=[str(CODE/'src'),str(CODE/'bench')]

class Reader:
 def __init__(self,stream): self.stream=stream;self.buffer='';self.decoder=json.JSONDecoder()
 def whitespace(self):
  self.buffer=self.buffer.lstrip()
  while not self.buffer:
   self.buffer=self.stream.read(65536).lstrip()
   if not self.buffer:raise EOFError('Unexpected end of JSON')
 def peek(self):self.whitespace();return self.buffer[0]
 def take(self,char):
  assert self.peek()==char,(char,self.buffer[:30]);self.buffer=self.buffer[1:]
 def value(self):
  self.whitespace()
  while True:
   try: value,end=self.decoder.raw_decode(self.buffer)
   except json.JSONDecodeError:
    chunk=self.stream.read(65536)
    if not chunk:raise
    self.buffer+=chunk;continue
   # A numeric token may end at the buffer boundary.
   if end==len(self.buffer):
    chunk=self.stream.read(65536)
    if chunk:self.buffer+=chunk;continue
   self.buffer=self.buffer[end:];return value

def child():
 from rxn_core import AAMProblem,MolecularEndpoint
 from rxn_core.final_branches import FinalBranchCatalogue,FinalFamily
 start=time.perf_counter();rows=[]
 for c in range(140):
  raw=json.loads((ROOT/f'work/full140_inputs/{c}/input.json').read_text());problem=AAMProblem(*(MolecularEndpoint(**raw[s]) for s in ('reactant','product')))
  fingerprint=FinalBranchCatalogue(problem).problem_sha256
  archive=RUN/f'case{c}/final_branches.json.gz';temporary=archive.with_suffix('.gz.tmp')
  summary=json.loads(archive.with_name('result.json').read_text());counts={};paths=0;input_paths=None
  with gzip.open(archive,'rt') as source,gzip.open(temporary,'wt') as target:
   reader=Reader(source);reader.take('{');target.write('{"problem_sha256":'+json.dumps(fingerprint))
   while reader.peek()!='}':
    key=reader.value();reader.take(':')
    if key=='problem_sha256':assert reader.value()==fingerprint
    else:
     target.write(','+json.dumps(key)+':')
     if reader.peek()=='[':
      reader.take('[');target.write('[');count=0
      while reader.peek()!=']':
       value=reader.value()
       if key=='families':
        family=FinalFamily.from_record(value);family.validate_representative(problem)
        assert json.dumps(family.to_record(),sort_keys=True)==json.dumps(value,sort_keys=True)
        paths+=len(family.provenance)
       elif key=='branches':assert all(0<=i<summary['flat_saved_families'] for i in value['families'])
       if count:target.write(',')
       json.dump(value,target,separators=(',',':'));count+=1
       if reader.peek()==',':reader.take(',')
       else:break
      reader.take(']');target.write(']');counts[key]=count
     else:
      value=reader.value();json.dump(value,target,separators=(',',':'))
      if key=='input_paths':input_paths=value
      if key=='schema':assert value==FinalBranchCatalogue.schema
    if reader.peek()==',':reader.take(',')
    else:break
   reader.take('}');target.write('}')
  assert counts['branches']==summary['final_branches'] and counts['families']==summary['flat_saved_families']
  assert not counts.get('coupled_branches') and paths==input_paths
  temporary.replace(archive)
  rows.append(dict(case=c,branches=counts['branches'],families=counts['families'],all_provenance_retained=True,family_roundtrip_validated=True,problem_sha256=fingerprint))
  print(c,flush=True)
 (RUN/'export_validation.json').write_text(json.dumps(dict(cases=rows,wall_seconds=time.perf_counter()-start),indent=2)+'\n')
 print('Validated all 140 standalone final catalogues')
if __name__=='__main__':
 if '--child' in sys.argv:child()
 else:
  from run_event_campaign import MemoryGuard,bounded_process
  with (RUN/'export_validation.log').open('w') as log:status,peak=bounded_process([sys.executable,str(Path(__file__).resolve()),'--child'],dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'),log,300,MemoryGuard(3072,6144,6144))
  (RUN/'export_validation_execution.json').write_text(json.dumps(dict(status=status,peak_mib=peak)))
  print(status,peak)
