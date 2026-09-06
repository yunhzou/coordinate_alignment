"""Prepare the historical eight targets against one bank and submit budgeted cases."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import subprocess


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--cpu-budget',type=int,default=2048)
    args=parser.parse_args()
    old=Path('data/retro_runs/eight_target_merged_fast_delivery_with_inventory')
    catalog=Path('data/mcule/merged_fast_delivery_with_inventory.csv.gz').resolve()
    with gzip.open(catalog,'rt') as stream:
        rows={r['Bank ID']:r for r in csv.DictReader(stream)}
    targets=[line.split('\t') for line in (old/'targets.tsv').read_text().splitlines()]
    budget=args.cpu_budget//len(targets)
    submissions=[]
    for case_id,smiles in targets:
        directory=args.root/case_id
        directory.mkdir(parents=True,exist_ok=False)
        report=json.loads((old/case_id/'results.json').read_text())
        reference=dict(reactants=[dict(id=i,smiles=rows[i]['SMILES']) for i in report['expected_ids']])
        run=directory/'run'
        if case_id.startswith('t05_'):
            reference=json.loads(Path('docs/example_runs/t05_ground_truth.json').read_text())
            run=Path('/project/yunhengzou/coordinate_alignment/retro_runs/beta_full_t05_20260906')
            reference['scan_summary']=str(run/'assemblies_viewer.json')
        manifest=dict(reference,id=case_id,target_smiles=smiles,run=str(run),catalog=str(catalog))
        path=directory/'case.json'
        path.write_text(json.dumps(manifest,indent=2)+'\n')
        job=subprocess.check_output(['sbatch','--parsable',
            '--exclude=bosque38,bosque39,bosque48,bosque56,bosque61,bosque75',
            f'--output={directory}/coordinator.out',f'--error={directory}/coordinator.err',
            'hpc/beta_case.sbatch',str(path),str(budget)],text=True).strip().split(';')[0]
        submissions.append(dict(case=case_id,job=job,cpu_budget=budget,manifest=str(path)))
        (args.root/'submissions.json').write_text(json.dumps(submissions,indent=2)+'\n')
        print(json.dumps(submissions[-1]),flush=True)


if __name__=='__main__':
    main()
