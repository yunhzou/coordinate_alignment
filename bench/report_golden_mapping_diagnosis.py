"""Export saved mapping diagnostics without remapping or repeating verification."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import re


def read(path):return json.loads(path.read_text())


def csv_write(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def export(run,output):
    output.mkdir(parents=True,exist_ok=True)
    diagnoses=[]
    for p in sorted(run.glob('*/diagnosis.json')):
        d=read(p);audit=read(p.parent/'reference_audit.json')
        # The first raw-order audit retained atom-map labels during RDKit
        # stereochemistry assignment. Use the independently normalized audit.
        d['original_reference_verified']=audit['original_order_reference_matches_saved']
        diagnoses.append(d)
    assert len(diagnoses)==137
    probes=[]
    for folder in sorted(list(run.glob('probe*'))+list(run.glob('empty*'))):
        if not folder.is_dir() or not (folder/'evaluation.json').exists():continue
        e=read(folder/'evaluation.json');s=read(folder/'search.json')
        design=read(folder/'diagnostic_design.json') if (folder/'diagnostic_design.json').exists() else {}
        probes.append(dict(experiment=folder.name,reference_recovery=e['reference_recovery'],
            terminals=e['candidate_terminals'],capped=e['capped'],seconds=s['seconds'],
            reference_directed=bool(design.get('reference_directed')),archive=str(folder/'cuts/aam.pkl.gz')))
    verification=[]
    for folder in sorted(run.glob('unknown_*')):
        if not folder.is_dir():continue
        index=folder.name.split('_')[1];states={}
        row=dict(index=int(index),original='unknown')
        for label,prefix in [('filtered','unknown'),('finite','finite'),('long_query','longquery'),
                             ('orbit_query','orbitquery')]:
            path=run/f'{prefix}_{index}'/'verification.json'
            if path.exists():
                data=read(path);states[label]=data['reference_recovery']
                row[label]=data['reference_recovery'];row[label+'_seconds']=data['evaluation_seconds']
                row[label+'_evidence']=str(path)
        assert not ('recovered' in states.values() and 'not_recovered' in states.values())
        row['resolved_status']=('recovered' if 'recovered' in states.values() else
            'not_recovered' if 'not_recovered' in states.values() else 'unknown')
        verification.append(row)
    summary=dict(run=str(run),audited_nonempty_misses=len(diagnoses),
        uncapped=sum(not x['cap_stages'] for x in diagnoses),
        capped=sum(bool(x['cap_stages']) for x in diagnoses),
        incorrect_reference_conversions=sum(not x['original_reference_verified'] for x in diagnoses),
        cases_with_non_endpoint_symmetry_actions=sum(x['non_endpoint_automorphisms']>0 for x in diagnoses),
        cases_with_unfinalized_terminal_edges=sum(x['unfinalized_terminal_edges']>0 for x in diagnoses),
        no_reference_source_selection=sum(x['reference_source_orbit_selection_matches']==0 for x in diagnoses),
        no_reference_orbit_pairing=sum(x['reference_atom_orbit_pairing_matches']==0 for x in diagnoses),
        original_unknowns=len(verification),verification_resolution=dict(Counter(x['resolved_status'] for x in verification)),
        unresolved_verification_indices=[x['index'] for x in verification if x['resolved_status']=='unknown'],
        probes=probes)
    sampling=set()
    for probe in probes:
        match=re.fullmatch(r'(?:probe|empty)(\d+)_seeds(\d+)(?:_cap100|_tol1)?',probe['experiment'])
        if match and int(match[2])>3 and probe['reference_recovery']=='recovered' and not probe['reference_directed']:
            sampling.add(int(match[1]))
    summary['confirmed_blind_higher_seed_recoveries']=sorted(sampling)
    resolved={x['index']:x['resolved_status'] for x in verification}
    source_report=Path(read(run/'manifest.json')['source_report'])
    with source_report.open() as f:original=list(csv.DictReader(f))
    fixed=[]
    for row in original:
        if row['reference_complete']!='True':continue
        index=int(row['index']);previous=row['reference_recovery']
        fixed.append(dict(index=index,original_reference_recovery=previous,
            verified_reference_recovery=resolved.get(index,previous),
            original_top1_correct=row['top1_correct'],archive=row['archive']))
    summary['unchanged_search_complete_reference']=dict(Counter(x['verified_reference_recovery'] for x in fixed))
    csv_write(output/'unchanged_search_reference_status.csv',fixed)
    csv_write(output/'misses137.csv',[{k:v for k,v in x.items() if not isinstance(v,(dict,list))} for x in diagnoses])
    csv_write(output/'probes.csv',probes);csv_write(output/'verification23.csv',verification)
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='probes'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();export(args.run,args.output)
