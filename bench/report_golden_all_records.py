"""Combine certified fixed-search results under exact relation scoring."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import subprocess
from golden_policy_campaign import save


def main(args):
    rows=[];witnesses=[]
    manifest=json.loads((args.source/'manifest.json').read_text())
    for record in manifest['records']:
        index=record['index'];original=args.source/str(index)/'evaluation.json'
        base=json.loads(original.read_text());path=original
        for run in (args.corrected,args.one_sided):
            candidate=run/str(index)/'evaluation.json'
            if candidate.exists():path=candidate
        result=json.loads(path.read_text())
        if not record['reference_annotation_complete'] and path==original:
            raise ValueError(f'One-sided case {index} has not been strictly rescored')
        top1=result.get('top1_correct',base.get('top1_correct'))
        rows.append(dict(index=index,one_sided=not record['reference_annotation_complete'],
            recovery=result['reference_recovery'],top1=top1,evaluation=str(path),
            seconds=result.get('total_seconds',result.get('evaluation_seconds'))))
        if result['reference_recovery']=='recovered' and path!=original:
            witnesses.append(dict(index=index,**result))
    def metrics(items):
        return dict(total=len(items),outcomes=dict(Counter(r['recovery'] for r in items)),
            recovered=sum(r['recovery']=='recovered' for r in items),
            coverage_percent=100*sum(r['recovery']=='recovered' for r in items)/len(items),
            top1_correct=sum(r['top1'] is True for r in items),
            top1_unknown=sum(r['top1'] is None for r in items),
            top1_percent=100*sum(r['top1'] is True for r in items)/len(items))
    opposite=set(json.loads(args.diagnostics.read_text())['opposite_recovered'])
    combined={r['index'] for r in rows if r['recovery']=='recovered'}|opposite
    summary=dict(all_records=metrics(rows),one_sided=metrics([r for r in rows if r['one_sided']]),
        fully_paired_product=metrics([r for r in rows if not r['one_sided']]),
        combined_search_recovered=len(combined),combined_search_percent=100*len(combined)/len(rows),
        not_recovered=[r['index'] for r in rows if r['recovery']=='not_recovered'],
        unknown=[r['index'] for r in rows if r['recovery']=='unknown'],
        semantics='Exact heavy-atom relation and unmatched status, modulo endpoint chemical symmetry. Explicit H in AAM; no H-identity accuracy claim.',
        caveat='Combined-search coverage is a union with prior opposite-direction diagnostics, not fixed-policy accuracy. This is our evaluator, not a reproduction of the original CGR scoring implementation.')
    args.output.mkdir(parents=True,exist_ok=True)
    save(args.output/'summary.json',summary);save(args.output/'witnesses.json',witnesses)
    save(args.output/'manifest.json',dict(source=str(args.source),corrected=str(args.corrected),
        one_sided=json.loads((args.one_sided/'manifest.json').read_text())))
    with (args.output/'cases.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    if args.jobs:
        accounting=subprocess.check_output(['sacct','-j',args.jobs,'-P',
            '--format=JobID,State,ElapsedRaw,AllocCPUS,TotalCPU,MaxRSS,Start,End'],text=True)
        (args.output/'slurm_accounting.psv').write_text(accounting)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('source','corrected','one-sided','diagnostics','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--jobs')
    main(parser.parse_args())
