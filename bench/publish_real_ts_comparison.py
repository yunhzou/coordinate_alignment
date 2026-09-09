"""Package saved real-TS comparisons and an offline preferred-style viewer."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

from compare_real_ts_mappings import METHODS, save
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events


def scalar_events(problem, mapping):
    """Independent all-pairs check of a complete physical mapping."""
    counts = dict(broken=0, formed=0, order_changed=0)
    for a in range(problem.atom_count):
        for b in range(a + 1, problem.atom_count):
            left = float(problem.reactant.wbo[a, b])
            right = float(problem.product.wbo[mapping[a], mapping[b]])
            if left > .2 and right <= .2:
                counts['broken'] += 1
            elif left <= .2 and right > .2:
                counts['formed'] += 1
            elif left > .2 and right > .2 and abs(left - right) > .5:
                counts['order_changed'] += 1
    return dict(**counts, total=sum(counts.values()))


def main(args):
    args.output.mkdir(parents=True,exist_ok=False)
    analysis=json.loads((args.run/'analysis.json').read_text())
    totals={m:dict(search_cpu=0.,symmetry_cpu=0.,compute_cpu=0.,all_H_scoring_cpu=0.,
                   capped_cuts=0,heavy_witness_classes=0) for m in METHODS}
    table=[];checked=0;score_rows={}
    for row in analysis['cases']:
        raw=json.loads((args.run/f"inputs/{row['index']}.json").read_text())
        problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
        entry=dict(index=row['index'],name=row['name'],direction=row['direction'],methods={})
        for method in METHODS:
            value=row['methods'][method]
            scored=json.loads((args.run/f"results/{row['index']}/{method}/full_atom_scores.json").read_text())
            s=scored['summary'];score_rows[(row['index'],method)]=s
            assert s['minimum_proven_over_saved_families'] and not s['inconsistent_native_units']
            for u in scored['families']:
                if u.get('upper') is None:continue
                mapping=dict(u['mapping'])
                assert len(mapping)==len(set(mapping.values()))==problem.atom_count
                assert all(problem.reactant.elements[a]==problem.product.elements[b] for a,b in mapping.items())
                assert scalar_events(problem,mapping)==bond_events(problem,mapping)==u['events']
                checked+=1
            entry['methods'][method]=dict(full_H_events=s['upper'],lower=s['lower'],
                full_H_minimum_proven=s['minimum_proven_over_saved_families'],
                heavy_witness_classes=value['heavy_classes'],capped_cuts=value['capped_units'],
                cut_or_order_count=value['total_units'],compute_cpu=value['compute_cpu'],
                all_H_scoring_cpu=s['postprocess_cpu'],best=s['best'])
            t=totals[method];t['search_cpu']+=value['phases']['search']['cpu']
            t['symmetry_cpu']+=value['phases'].get('symmetry',{}).get('cpu',0.)
            t['compute_cpu']+=value['compute_cpu'];t['all_H_scoring_cpu']+=s['postprocess_cpu']
            t['capped_cuts']+=value['capped_units'];t['heavy_witness_classes']+=value['heavy_classes']
        table.append(entry)
    save(args.output/'summary.json',dict(rows=table,totals=totals,rescored_full_atom_witnesses=checked,
        reference_scope='Historical heavy-core assignments, not curated chemical truth.',
        diversity_scope='AAM heavy witness classes are lower bounds on retained compressed alternatives.'))
    # Two methods toggle over the same original endpoint geometry, not side-by-side
    # independently reoriented illustrations. Everything is generated from saved maps.
    displays=[]
    labels=dict(aam_weighted='AAM continuous',aam_binary='AAM binary',slap_binary='SLAP binary',
                slap_weighted='SLAP raw WBO',slap_float_cost='SLAP float-cost diagnostic')
    def append(index,left,left_map,right,right_map,note):
        raw=json.loads((args.run/f'inputs/{index}.json').read_text())
        endpoints=[raw[k] for k in ('reactant','product')]
        problem=AAMProblem(*(MolecularEndpoint(**d) for d in endpoints)); records=[]
        for name,mapping,provenance in ((left,left_map['mapping'],left_map),(right,right_map['mapping'],right_map)):
            mapping=dict(mapping) if isinstance(mapping[0],list) else dict(enumerate(mapping))
            counts=bond_events(problem,mapping);events=[]
            for a in range(problem.atom_count):
                for b in range(a+1,problem.atom_count):
                    w=float(problem.reactant.wbo[a,b]);v=float(problem.product.wbo[mapping[a],mapping[b]])
                    kind=('broken' if w>.2 and v<=.2 else 'formed' if w<=.2 and v>.2 else
                          'order_changed' if w>.2 and v>.2 and abs(w-v)>.5 else None)
                    if kind:events.append(dict(kind=kind,r=[a,b],p=[mapping[a],mapping[b]],wbo=[w,v]))
            assert len(events)==counts['total']
            records.append(dict(name=name,mapping=sorted(mapping.items()),counts=counts,events=events,provenance=provenance))
        displays.append(dict(index=index,name=f"{raw['name']} · {raw['direction']} · {left} / {right}",
                             endpoints=endpoints,records=records,source=str(args.run),note=note))
    extra=analysis['cases'][4]['methods']['aam_weighted']['historical_core']['2']['example']
    append(4,'AAM extra oxygen-role assignment',extra,'SLAP binary best',score_rows[(4,'slap_binary')]['best'],
           'Historical heavy core #2 appears only on the AAM side for R→TS; higher event cost, NOT ground truth.')
    for row in table:
        for method in METHODS[1:]:
            append(row['index'],labels['aam_weighted'],score_rows[(row['index'],'aam_weighted')]['best'],
                   labels[method],score_rows[(row['index'],method)]['best'],
                   'Best all-H event score within saved families. Core assignments are not independently validated.')
    template=Path(__file__).with_name('elementary_comparison_viewer.html').read_text()
    replacements={
        '<b>R</b><span>Original reactants</span>':'<b>Source</b><span>R or P endpoint, as named above</span>',
        '<b>P</b><span>Product · mapped source identities</span>':'<b>TS</b><span>Cached reference TS geometry · mapped source identities</span>',
        'Our AAM target':'Choice A target','SLAP target':'Choice B target',
        'AAM → p${mapping(0)[r]}, SLAP → p${mapping(1)[r]}':'${current.records[0].name} → p${mapping(0)[r]}, ${current.records[1].name} → p${mapping(1)[r]}',
        '`Our AAM · ${current.records[0].counts.total} events`':'`${current.records[0].name} · ${current.records[0].counts.total} events`',
        '`SLAP · ${current.records[1].counts.total} events`':'`${current.records[1].name} · ${current.records[1].counts.total} events`',
        'trace R → P':'trace source → TS',
        'R atoms → P atoms':'Source atoms → TS atoms',
        'WBO R → P':'WBO source → TS',
        'Both choices are actual saved mappings, not ground truth.':'Both choices are saved full-atom mappings, not ground truth. Here r denotes the source endpoint atom index and p denotes the TS atom index.',
    }
    for before,after in replacements.items():
        assert before in template;template=template.replace(before,after)
    library=(Path(__file__).resolve().parents[1]/'src/rxn_core/static/3Dmol-min.js').read_text()
    html=template.replace('__LIBRARY__',library).replace('__DATA__',json.dumps(displays))
    (args.output/'viewer.html').write_text(html)
    save(args.output/'viewer_data.json',displays)
    for name in ('analysis.json','membership.json','core_equivalence.json','manifest.json','references.json',
                 'submission.json','relocation.json','scoring_submission.json'):
        shutil.copy2(args.run/name,args.output/name)
    jobs=[json.loads((args.run/name).read_text())['job'] for name in ('submission.json','relocation.json','scoring_submission.json')]
    accounting=subprocess.check_output(['sacct','-j',','.join(jobs),'-P',
        '--format=JobID,State,ExitCode,ElapsedRaw,TotalCPU,AllocCPUS,MaxRSS,NodeList'],text=True)
    (args.output/'slurm_accounting.psv').write_text(accounting)
    shutil.copy2(args.slap_license,args.output/'SLAP_LICENSE')
    with tarfile.open(args.output/'saved_outputs.tar.gz','x:gz') as archive:
        for folder in ('inputs','results','status'):
            archive.add(args.run/folder,arcname=folder)
        archive.add(args.run/'slap_core.py',arcname='slap_core.py')
        archive.add(args.run/'engine/bench/compare_real_ts_mappings.py',arcname='search_driver.py')
        archive.add(args.run/'score_real_ts_families.py',arcname='score_driver.py')
    print(json.dumps(dict(rescored=checked,viewer=str((args.output/'viewer.html').resolve()),totals=totals),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--slap-license',type=Path,required=True)
    main(p.parse_args())
