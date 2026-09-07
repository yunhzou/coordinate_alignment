"""Summarize diagnosis, preserving the discovered domain-verification limitation."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from golden_policy_campaign import save


def main(args):
    manifest=json.loads((args.run/'manifest.json').read_text());rows=[]
    for job in manifest['jobs']:
        folder=args.run/f"{job['index']}_{job['variant']}"
        def read(name):
            path=folder/name
            return json.loads(path.read_text()) if path.exists() else {}
        diagnosis=read('diagnosis.json');partial=read('partial_witness.json');search=read('search_status.json')
        positive=diagnosis.get('reference_recovery')=='recovered' or partial.get('reference_recovery')=='recovered'
        rows.append(dict(index=job['index'],variant=job['variant'],verified_positive=positive,
            scope='partial cut' if partial.get('reference_recovery')=='recovered' else 'full checkpoint',
            generator_only_outcome=diagnosis.get('reference_recovery','interrupted'),
            source_eligible=diagnosis.get('source_eligible'),
            search_exit=search.get('exit_code'),capped=diagnosis.get('capped'),
            search_seconds=read('search.json').get('seconds'),
            warning='Negative generator-only verification is NOT proof of absence from assignment domains',
            directory=str(folder)))
    singletons=[json.loads(p.read_text()) for p in args.run.glob('singletons/*/result.json')]
    summary=dict(original_remaining_records=len({j['index'] for j in manifest['jobs']}),
        variants=len(rows),opposite_recovered=[r['index'] for r in rows if r['variant']=='opposite' and r['verified_positive']],
        cap2000_recovered=[r['index'] for r in rows if r['variant']=='cap2000' and r['verified_positive']],
        cap2000_watchdogs=[r['index'] for r in rows if r['variant']=='cap2000' and r['search_exit']=='timeout'],
        saved_domain_recovered=[r['index'] for r in singletons if r['reference_recovery']=='recovered'],
        singleton_checks=len(singletons),
        finding='Evaluator consumes target_generators but omits symmetry_domains; singleton domains can span multiple conditioned automorphism orbits.',
        warning='Do not merge these ablations into the fixed benchmark. Positive witnesses are valid; generator-only negatives require domain-aware reverification.',
        source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in map(Path,[
            'bench/golden_remaining.py','bench/golden_singleton_probe.py','bench/golden_evaluation.py'])})
    args.output.mkdir(parents=True,exist_ok=True)
    save(args.output/'summary.json',summary);save(args.output/'manifest.json',manifest)
    save(args.output/'singletons.json',singletons)
    with (args.output/'cases.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);main(p.parse_args())
