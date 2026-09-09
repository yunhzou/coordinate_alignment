"""Native SLAP alternatives with symmetry controls and blind order/direction diversity."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random
import time

from golden_competitors import save, signatures


def prepare(args):
    from rdkit import Chem
    args.run.mkdir(exist_ok=False, parents=True)
    source = Path('data/aam_benchmarks/golden_original_20260906')
    rows = [json.loads(s) for s in (source/'audit.jsonl').read_text().splitlines()]
    inputs, plans = [], []
    for row in rows:
        sides = row['input_reaction'].split('>>')
        rng = random.Random(args.seed + row['index'])
        for order in range(args.orders):
            if order == 0:
                shuffled = sides
            else:
                shuffled = []
                for side in sides:
                    mol = Chem.MolFromSmiles(side)
                    indices = list(range(mol.GetNumAtoms()))
                    rng.shuffle(indices)
                    shuffled.append(Chem.MolToSmiles(Chem.RenumberAtoms(mol,indices), canonical=False))
            for reverse in (False,True):
                serial = len(inputs)
                reaction = '>>'.join(shuffled[::-1] if reverse else shuffled)
                inputs.append(dict(index=serial,input_reaction=reaction))
                plans.append(dict(index=serial,case=row['index'],order=order,reverse=reverse))
    (args.run/'inputs.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in inputs))
    save(args.run/'plans.json',plans)
    save(args.run/'manifest.json',dict(dataset=str(source.resolve()),orders=args.orders,
        seed=args.seed,cases=len(rows),attempts_per_mode=len(inputs),
        modes=[f"slap_{'all_' if args.symmetry=='all' else ''}{m}" for m in ('binary','weighted')],
        break_sym=args.symmetry,add_Hs=True,
        watchdog_seconds=300,output_limit=None,
        input_sha256=hashlib.sha256((args.run/'inputs.jsonl').read_bytes()).hexdigest(),
        interpretation='Native minimum-cost alternatives only, no core changes or reference-guided constraints. '
                       'Input-order and direction diversity is an external search-budget experiment, not native top-K.'))


def evaluate(args):
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.warning')
    manifest=json.loads((args.run/'manifest.json').read_text())
    rows={r['index']:r for r in map(json.loads,(Path(manifest['dataset'])/'audit.jsonl').read_text().splitlines())}
    plans=json.loads((args.run/'plans.json').read_text())
    outputs=[]
    start=time.perf_counter()
    destination=args.run/f'evaluation_{args.method}_{args.shard}.json'
    fingerprint=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    previous=json.loads(destination.read_text()) if destination.exists() else {}
    cached={c['case']:c for c in previous.get('cases',[])
            if previous.get('evaluator_sha256')==fingerprint
            and all(a['status']!='missing' for a in c['attempts'])}
    for case in range(args.shard,manifest['cases'],args.shards):
        if case in cached:
            outputs.append(cached[case]);continue
        expected_endpoints,expected=signatures(rows[case]['mapped_reaction'])
        seen=set();attempts=[]
        for plan in plans[case*2*manifest['orders']:(case+1)*2*manifest['orders']]:
            path=args.run/args.method/f"{plan['index']}.json"
            if not path.exists():
                attempts.append(dict(**plan,status='missing'));continue
            raw=json.loads(path.read_text())
            hits=[];invalid=[]
            for i,candidate in enumerate(raw.get('candidates',[])):
                reaction=candidate['mapped_rxn']
                if plan['reverse']:
                    left,agents,right=reaction.split('>')
                    reaction='>'.join((right,agents,left))
                try:
                    endpoints,actual=signatures(reaction)
                    if endpoints!=expected_endpoints:
                        raise ValueError('Endpoint chemistry changed')
                    seen.add(actual)
                    if actual==expected:hits.append(i)
                except Exception as error:
                    invalid.append(dict(candidate=i,error=str(error)))
            attempts.append(dict(**plan,status=raw['status'],hits=hits,
                candidates=len(raw.get('candidates',[])),invalid=invalid,
                request_seconds=raw.get('request_seconds'),
                mapping_seconds=raw.get('mapping_seconds'),mapping_cpu_seconds=raw.get('mapping_cpu_seconds')))
        outputs.append(dict(case=case,attempts=attempts,unique_heavy_mappings=len(seen),
                            recovered=any(a.get('hits') for a in attempts)))
    save(destination,dict(cases=outputs,evaluator_sha256=fingerprint,
        cached_cases=len(cached),evaluation_seconds=time.perf_counter()-start))


def summarize(args):
    manifest=json.loads((args.run/'manifest.json').read_text());summary={}
    for method in manifest['modes']:
        cases=[c for f in sorted(args.run.glob(f'evaluation_{method}_*.json'))
               for c in json.loads(f.read_text())['cases']]
        assert len({c['case'] for c in cases})==manifest['cases']
        attempts=[a for c in cases for a in c['attempts']]
        summary[method]=dict(cases=len(cases),statuses=dict(Counter(a['status'] for a in attempts)),
            returned_candidates=sum(a.get('candidates',0) for a in attempts),
            unique_heavy_mappings=sum(c['unique_heavy_mappings'] for c in cases),
            recovered=sum(c['recovered'] for c in cases),
            recovery_by_orders={str(n):sum(any(a.get('hits') and a['order']<n for a in c['attempts'])
                for c in cases) for n in sorted({min(x,manifest['orders']) for x in (1,2,5,manifest['orders'])})},
            cpu_hours=sum(a.get('mapping_cpu_seconds') or 0 for a in attempts)/3600,
            all_request_seconds_sum=sum(a.get('request_seconds') or 0 for a in attempts),
            timing_note='Mapping CPU excludes killed workers without a final CPU reading; watchdog costs remain in request times and Slurm accounting.',
            mapping_seconds_sum=sum(a.get('mapping_seconds') or 0 for a in attempts),
            unresolved_cases=[c['case'] for c in cases if not c['recovered']])
    summary['union_recovered']=manifest['cases']-len(set.intersection(
        *(set(summary[m]['unresolved_cases']) for m in manifest['modes'])))
    save(args.run/'summary.json',summary)
    print(json.dumps(summary,indent=2))


def retry(args):
    from golden_competitors import worker
    plan=json.loads((args.run/args.retry_plan).read_text())[args.shard]
    args.method=plan['method'];args.indices=str(plan['index'])
    args.shards=1;args.shard=0
    worker(args)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=('prepare','evaluate','summarize','retry'))
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--orders',type=int,default=10)
    p.add_argument('--seed',type=int,default=20260908)
    p.add_argument('--method',choices=('slap_all_binary','slap_all_weighted','slap_binary','slap_weighted'))
    p.add_argument('--symmetry',choices=('heavy','all'),default='heavy')
    p.add_argument('--shard',type=int,default=0)
    p.add_argument('--shards',type=int,default=1)
    p.add_argument('--timeout',type=int,default=300)
    p.add_argument('--retry-plan',default='retry_plan.json')
    a=p.parse_args();{'prepare':prepare,'evaluate':evaluate,'summarize':summarize,'retry':retry}[a.command](a)


if __name__=='__main__':main()
