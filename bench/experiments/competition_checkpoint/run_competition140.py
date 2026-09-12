import os,json,hashlib
from pathlib import Path
import fragment_competition_optimized as m
m.BUDGET=128;m.WORKERS=4;m.CASES=tuple(range(140));m.OUT=m.ROOT/'outputs/competition140_budget128'
m.OUT.mkdir(parents=True,exist_ok=True)
m.os.environ.update(COMPETITION_BUDGET='128',COMPETITION_WORKERS='4',COMPETITION_CASES=','.join(map(str,m.CASES)),COMPETITION_OUTPUT=str(m.OUT))
inputs={}
for c in m.CASES:
 p=m.ROOT/f'outputs/holdout140_cap2000_seed1/case{c}/cuts/aam.pkl.gz' if c!=25 else m.ROOT/'outputs/holdout140_case25_retry/case25/cuts/aam.pkl.gz'
 assert p.is_file(),p
 inputs[str(c)]=dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
script=Path(m.__file__)
m.save(m.OUT/'manifest.json',dict(budget=128,workers=4,hard_case_seconds=300,soft_search_seconds=270,per_worker_mib=3072,total_workers_mib=8192,system_reserve_mib=6144,source=str(script),sha256=hashlib.sha256(script.read_bytes()).hexdigest(),inputs=inputs,selection='Budget selected on the four known gap cases; subsequent evaluation covers all 140 of the same dataset. This is not independent held-out validation.'))
m.main()
