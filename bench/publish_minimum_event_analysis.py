"""Build an audited minimum-event catalogue, figures, and an offline viewer."""

from rxn_core.viewers import style_document
import argparse
from collections import Counter
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
import tarfile

import numpy as np
from holdout_minimum_events import AAM,EventPatterns,read,save,sha
from rxn_core import AAMProblem
from rxn_core.domain import MolecularEndpoint
from rxn_core.family_scoring import bond_events

METHODS=('aam','native_slap','slap_sweep')


def comparison(left,right,membership,complete,same_score):
    left,right=set(left),set(right)
    if not same_score:
        return dict(aam=len(left),slap=len(right),shared=[],aam_only=[],
                    slap_only_proven=[],unresolved=[],identical=None,
                    observed_sets_equal=None,certified_difference=False,
                    comparable=False,reason='Different benchmark minimum scores')
    missing=right-left
    excluded={k for k in missing if membership[k]['aam'].startswith('excluded')}
    unresolved=missing-excluded
    return dict(aam=len(left),slap=len(right),shared=sorted(left & right),
        aam_only=sorted(left-right),slap_only_proven=sorted(excluded),unresolved=sorted(unresolved),
        identical=left==right and complete,observed_sets_equal=left==right,
        certified_difference=bool(left-right or excluded),comparable=True)


def collect(run):
    rows,hashes=[],{}
    checked=0
    for index in range(140):
        snapshot=read(run/f'snapshots/{index}.json')
        family_path=run/f'families_enriched/{index}.json'
        if not family_path.exists():family_path=run/f'families/{index}.json'
        family=read(family_path)
        enumeration_path=run/f'enumeration_extended/{index}.json'
        if not enumeration_path.exists():enumeration_path=run/f'enumeration/{index}.json'
        native=read(enumeration_path)
        aam=read(run/f'aam_enumeration/{index}.json')
        assert 'finished' in read(run/f'status/aam_enumeration_{index}.json')
        assert all(v['complete'] for v in native['methods'].values()),index
        raw=read(AAM/f'inputs/{index}/input.json');canonical=EventPatterns(raw)
        problem=AAMProblem(*(MolecularEndpoint(**raw[k]) for k in ('reactant','product')))
        sets={'aam':dict(aam['patterns']),**{m:dict(native['methods'][m]['patterns']) for m in METHODS[1:]}}
        for key,value in snapshot['methods']['aam']['patterns'].items():sets['aam'].setdefault(key,value)
        for key,value in family['methods']['aam'].items():
            if value['status']=='represented':sets['aam'].setdefault(key,value['witness'])
        minima={m:snapshot['methods'][m]['minimum'] for m in METHODS}
        patterns={}
        for method,values in sets.items():
            for key,value in values.items():
                actual=canonical.describe(value['mapping'])
                assert actual['id']==key and actual['total']==minima[method]
                scalar=bond_events(problem,dict(enumerate(actual['mapping'])))
                c=actual['counts']
                assert scalar==dict(broken=c['broken'],formed=c['formed'],
                    order_changed=c['strengthened']+c['weakened'],total=actual['total'])
                checked+=1
                patterns.setdefault(key,actual)
        membership={}
        for key,pattern in patterns.items():
            state={}
            for method in METHODS:
                if key in sets[method]:state[method]='represented'
                elif pattern['total']!=minima[method]:state[method]='different_score'
                elif method!='aam' or aam['complete']:state[method]='excluded_by_complete_enumeration'
                elif family['methods']['aam'].get(key,{}).get('status')=='excluded_from_saved_families':
                    state[method]='excluded_by_family_queries'
                else:state[method]='unresolved'
            membership[key]=state
        coords=np.asarray(raw['reactant']['coordinates'],dtype=float)
        centered=coords-coords.mean(axis=0)
        _,_,vt=np.linalg.svd(centered,full_matrices=False)
        xy=centered@vt[:2].T
        for column in range(2):
            if xy[np.argmax(np.abs(xy[:,column])),column]<0:xy[:,column]*=-1
        w=np.asarray(raw['reactant']['wbo'])
        bonds=[(int(a),int(b),float(w[a,b])) for a,b in zip(*np.where(np.triu(w,1)>.2))]
        rows.append(dict(index=index,name=raw['name'],minima=minima,patterns=patterns,membership=membership,
            pattern_ids={m:sorted(v) for m,v in sets.items()},
            catalogue_complete=dict(aam=aam['complete'],native_slap=True,slap_sweep=True),
            comparisons={m:comparison(sets['aam'],sets[m],membership,aam['complete'],
                                      minima['aam']==minima[m]) for m in METHODS[1:]},
            elements=raw['reactant']['elements'],coordinates_2d=xy.tolist(),bonds=bonds,
            aam_representative_patterns=len(snapshot['methods']['aam']['patterns']),
            aam_full_terminals_scanned=snapshot['methods']['aam']['full_terminals_scanned'],
            aam_enumeration_metrics=aam['metrics'],family_query_metrics=family['query_metrics']))
        for p in [run/f'snapshots/{index}.json',family_path,enumeration_path,run/f'aam_enumeration/{index}.json']:
            hashes[str(p.relative_to(run))]=sha(p)
    summaries={}
    for method in METHODS[1:]:
        ties=[r for r in rows if r['minima']['aam']==r['minima'][method]]
        comp=[r['comparisons'][method] for r in ties]
        summaries[method]=dict(tied_cases=len(ties),
            verified_aam_patterns=sum(c['aam'] for c in comp),slap_patterns=sum(c['slap'] for c in comp),
            shared_patterns=sum(len(c['shared']) for c in comp),
            aam_only_patterns=sum(len(c['aam_only']) for c in comp),
            slap_only_proven_patterns=sum(len(c['slap_only_proven']) for c in comp),
            unresolved_memberships=sum(len(c['unresolved']) for c in comp),
            complete_equal_catalogue_cases=sum(c['identical'] for c in comp),
            equal_observed_catalogue_cases=sum(c['observed_sets_equal'] for c in comp),
            certified_different_cases=sum(c['certified_difference'] for c in comp),
            aam_multiple_cases=sum(c['aam']>1 for c in comp),slap_multiple_cases=sum(c['slap']>1 for c in comp),
            both_multiple_cases=sum(c['aam']>1 and c['slap']>1 for c in comp),
            difference_cases=[r['index'] for r in ties if r['comparisons'][method]['certified_difference']],
            incomplete_aam_cases=[r['index'] for r in ties if not r['catalogue_complete']['aam']])
    summary=dict(cases=140,comparisons=summaries,
        total_verified_patterns={m:sum(len(r['pattern_ids'][m]) for r in rows) for m in METHODS},
        complete_catalogue_cases={m:sum(r['catalogue_complete'][m] for r in rows) for m in METHODS},
        full_aam_terminals_scanned=sum(r['aam_full_terminals_scanned'] for r in rows),
        independently_rescored_witnesses=checked,
        scope='Patterns at each method\'s recorded benchmark minimum, with explicit H. '
              'AAM catalogues are lower bounds where full family enumeration is incomplete. '
              'Individual negative membership claims require complete family queries or complete enumeration. '
              'Distinct bond-edit patterns are not experimentally validated mechanisms.')
    return rows,summary,hashes


def publish(args):
    rows,summary,hashes=collect(args.run)
    args.destination.mkdir(parents=True,exist_ok=False)
    for name,value in [('summary.json',summary),('per_case.json',rows),('analysis_hashes.json',hashes)]:
        save(args.destination/name,value)
    with (args.destination/'per_case.csv').open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n')
        writer.writerow(['index','name','aam_score','native_slap_score','slap_sweep_score','aam_patterns',
            'native_slap_patterns','slap_sweep_patterns','aam_enumeration_complete','sweep_score_tied',
            'shared_with_sweep','aam_only_vs_sweep','sweep_only_proven','unresolved_sweep_membership'])
        for row in rows:
            c=row['comparisons']['slap_sweep']
            writer.writerow([row['index'],row['name'],*(row['minima'][m] for m in METHODS),
                *(len(row['pattern_ids'][m]) for m in METHODS),row['catalogue_complete']['aam'],c['comparable'],
                len(c['shared']),len(c['aam_only']),len(c['slap_only_proven']),len(c['unresolved'])])
    for folder in ('snapshots','families','families_enriched','enumeration','enumeration_extended','aam_enumeration'):
        with gzip.open(args.destination/f'{folder}.json.gz','wt') as stream:
            json.dump({p.stem:read(p) for p in sorted((args.run/folder).glob('*.json'))},stream)
    with gzip.open(args.destination/'worker_status.json.gz','wt') as stream:
        json.dump({p.stem:read(p) for p in sorted((args.run/'status').glob('*.json'))},stream)
    with tarfile.open(args.destination/'analysis_sources.tar.gz','w:gz') as archive:
        for path in sorted(args.run.glob('*.py')):
            archive.add(path,arcname=f'frozen/{path.name}')
        for path in sorted(Path(__file__).parent.glob('*minimum*')):
            if path.suffix in ('.py','.html'):
                archive.add(path,arcname=f'published/{path.name}')
        for path in sorted((args.run/'validation').glob('*.json')):
            archive.add(path,arcname=f'validation/{path.name}')
    for name in ['manifest.json','snapshot_submission.json','snapshot_relocation.json','family_submission.json',
                 'enumeration_submission.json','enumeration_relocation.json','aam_enumeration_submission.json',
                 'aam_enumeration_relocation.json','analysis_provenance.json','slurm_accounting.tsv']:
        shutil.copy2(args.run/name,args.destination/name)
    shutil.copy2(args.run/'validation/finite_checks.json',args.destination/'finite_checks.json')
    report(rows,summary,args.run,args.destination)
    viewer(rows,summary,args.destination)
    figure(rows,summary,args.destination)
    print(json.dumps(summary,indent=2),flush=True)


def report(rows,summary,run,destination):
    c=summary['comparisons']['slap_sweep'];n=summary['comparisons']['native_slap']
    interesting=[11,25,64,101,129,132]
    table='\n'.join(f"| {i} | {rows[i]['name']} | {rows[i]['minima']['aam']} | {'≥' if not rows[i]['catalogue_complete']['aam'] else ''}{len(rows[i]['pattern_ids']['aam'])} | {len(rows[i]['pattern_ids']['slap_sweep'])} | {len(rows[i]['comparisons']['slap_sweep']['shared'])} | {len(rows[i]['comparisons']['slap_sweep']['slap_only_proven'])} |"
        for i in interesting)
    (destination/'README.md').write_text(f"""# Minimum-event alternatives on the 140 XYZ/WBO cases

[Open the offline bond-edit viewer](viewer.html), or download it and open locally.
`comparison.pdf`, `comparison.svg` and `comparison.png` show the aggregate result;
`per_case.csv` and `per_case.json` retain every case and its membership status.

## AAM versus SLAP with the single-edge sweep

At the **{c['tied_cases']} tied benchmark scores**, the verified AAM catalogue contains
**{c['verified_aam_patterns']}** distinct bond-edit patterns; the complete saved
SLAP-sweep minimum-score label-family catalogue contains **{c['slap_patterns']}**.
There are **{c['shared_patterns']} verified shared patterns**, **{c['aam_only_patterns']}
AAM-only patterns**, **{c['slap_only_proven_patterns']} patterns proved present only
in SLAP's saved families**, and **{c['unresolved_memberships']} unresolved cross-method memberships**.
Counts for an incompletely enumerated AAM catalogue are lower bounds.

Both methods have more than one verified minimum-score pattern in
**{c['both_multiple_cases']}** tied cases. Their observed pattern sets agree in
{c['equal_observed_catalogue_cases']} cases; complete enumeration certifies equal
catalogues in {c['complete_equal_catalogue_cases']} cases. A difference is certified
in {c['certified_different_cases']} tied cases: {c['difference_cases']}.

| Case | Name | Common minimum score | Verified AAM patterns | SLAP-sweep patterns | Shared | SLAP-only proved |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
{table}

The original one-representative comparison hid alternatives in both methods.
In case 129, for example, the terminal/optimized representatives gave 12 AAM
patterns and seven SLAP-sweep patterns. Family enumeration and membership
queries establish at least 16 shared patterns in the saved outputs. A smaller representative
list does not establish that the algorithm missed those alternatives.

Against **native SLAP without the sweep**, {n['tied_cases']} scores tie. The
catalogues contain {n['shared_patterns']} verified shared patterns,
{n['aam_only_patterns']} AAM-only patterns, {n['slap_only_proven_patterns']}
SLAP-only patterns and {n['unresolved_memberships']} unresolved memberships on
those tied cases. The native and swept protocols remain separate.

## What constitutes an alternative

An alternative is the full set of **broken, formed, strengthened and weakened
bonds**, pulled back to original reactant atom indices. Order changes use the
same 0.5 WBO threshold and bond presence uses the same 0.2 floor as the benchmark.
Explicit H is included. Multiple atom maps with the same edits count once.
All edits are jointly canonicalized on the full reactant graph, with exact signed
event-response colors. Symmetry-equivalent atom relabelings are merged; different
affected bonds or event types remain distinct. Geometry and stereochemistry do
not distinguish patterns in this analysis. These are candidate bond-edit
mechanisms, not experimentally established reaction pathways.

## Completeness and checks

The analysis scanned **{summary['full_aam_terminals_scanned']:,} full AAM terminal
witnesses** from the frozen cap-1000, one-seed, tolerance-1.0 archives. It then
queried the recorded correlated path families and enumerated event patterns at
the **previously recorded benchmark minimum score**. It does not replace that
score with an independently optimized global minimum over all implicit mappings.

SLAP enumeration is complete for all 140 saved label-family catalogues. It blocks
whole event sets and removes only safe source/target twin-atom permutations,
preserving every event class. AAM enumeration is complete for
{summary['complete_catalogue_cases']['aam']}/140 cases. Incomplete tied cases are
{c['incomplete_aam_cases']}; their AAM counts carry a lower-bound marker in the
viewer. An incomplete enumeration is never used to assert absence. The separate
cross-family membership checks can still prove exclusion of an individual pattern
by scanning every applicable saved family.

Every published mapping was independently rescored with the scalar bond-event
routine ({summary['independently_rescored_witnesses']:,} method/witness checks).
Direct enumeration of 13,872 complete assignments in ten small real SLAP families
exactly matches the symbolic event-pattern sets. `finite_checks.json` records
these tests and event-type/symmetry checks. All analysis is offline; none of the
mapping searches or reported mapping timings was rerun.

Full analysis and frozen sources: `{run}`.
Compressed snapshots, membership witnesses, enumeration records, input/source
hashes, and job submissions accompany this report. The fixed manuscript baseline
is separate from this cap-1000 follow-up.

## Reproduction

The array submission records retain the exact environment, commands and budgets.
`analysis_sources.tar.gz` includes the frozen worker versions, published analysis
scripts and targeted-reproduction checks. Mapping inputs and large AAM archives
are identified by the preceding benchmark reports and their hash manifests.

The targeted case-25 SLAP retry and case-129 enrichment can be reproduced with
`bench/complete_minimum_event_checks.py retry_slap --index 25` and
`bench/complete_minimum_event_checks.py enrich_aam --index 129`, supplying
`--run` and a fresh `--destination`. Both wrappers were checked against the saved
results. The default retry uses a 30-second solver-query limit.

To rebuild this report, run `bench/publish_minimum_event_analysis.py --run`
with the source directory above and a fresh `--destination`, using the frozen
source paths in the submission records and the manuscript plotting dependencies.
`slurm_accounting.tsv` records analysis workers, including cancelled allocations
that never started. These postprocessing costs do not alter the mapping-speed
comparison in the preceding benchmark report.
""")


def figure(rows,summary,destination):
    os.environ.setdefault('MPLCONFIGDIR',str(destination/'.mpl-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none',
        'axes.spines.top':False,'axes.spines.right':False})
    tied=[r for r in rows if r['minima']['aam']==r['minima']['slap_sweep']]
    counts=Counter((len(r['pattern_ids']['slap_sweep']),len(r['pattern_ids']['aam'])) for r in tied)
    fig,axes=plt.subplots(1,2,figsize=(9.5,4),layout='constrained')
    limit=max(max(p) for p in counts)+1
    axes[0].plot([0,limit],[0,limit],ls='--',color='#ABB6BA',lw=1)
    for (x,y),n in counts.items():
        axes[0].scatter(x,y,s=35+12*np.sqrt(n),color='#007F68' if y==x else '#CC6B26',alpha=.8)
        if n>1:axes[0].annotate(str(n),(x,y),xytext=(6,5),textcoords='offset points',fontsize=8)
    axes[0].set(xlabel='SLAP sweep: minimum-event patterns',ylabel='AAM: verified minimum-event patterns',
        xlim=(0,limit),ylim=(0,limit),title='a  Alternative counts on 138 score ties')
    c=summary['comparisons']['slap_sweep']
    labels=['Shared','AAM only','SLAP sweep only','Unresolved membership']
    values=[c['shared_patterns'],c['aam_only_patterns'],c['slap_only_proven_patterns'],c['unresolved_memberships']]
    axes[1].barh(range(4),values,color=['#007F68','#3269A8','#CC6B26','#ABB6BA'])
    for y,value in enumerate(values):axes[1].text(value+2,y,str(value),va='center')
    axes[1].set(yticks=range(4),yticklabels=labels,xlim=(0,max(values)*1.13),
        xlabel='Distinct case / bond-edit patterns',title='b  Membership in saved output families')
    axes[1].invert_yaxis()
    for ax in axes:ax.grid(alpha=.12);ax.set_axisbelow(True)
    fig.get_layout_engine().set(rect=(0,.17,1,.83))
    fig.text(.02,.075,'Full-H bond edits; joint symmetry normalization | AAM: one seed, cap 1000 | Native XYZ SLAP + all single-edge cuts',fontsize=8)
    fig.text(.02,.025,f"AAM enumeration complete in {summary['complete_catalogue_cases']['aam']}/140 cases; remaining AAM counts are lower bounds. See per-case records.",fontsize=8)
    for ext in ('png','pdf','svg'):fig.savefig(destination/f'comparison.{ext}',dpi=200,facecolor='white')
    svg=destination/'comparison.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig);shutil.rmtree(destination/'.mpl-cache',ignore_errors=True)


def viewer(rows,summary,destination):
    payload=json.dumps(dict(cases=rows,summary=summary),separators=(',',':')).replace('</','<\\/')
    template=Path(__file__).with_name('minimum_events_viewer.html').read_text()
    template=style_document(template, layout='minimum')
    (destination/'viewer.html').write_text(template.replace('__PAYLOAD__',payload))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--destination',type=Path,required=True)
    publish(parser.parse_args())
