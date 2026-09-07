"""Freeze a completion engine and launch bounded, memory-budgeted node pools."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--version',required=True)
    parser.add_argument('--nodes',nargs='+',required=True)
    parser.add_argument('--workers',type=int,required=True)
    parser.add_argument('--cut-workers',type=int,default=1)
    parser.add_argument('--cpus',type=int,required=True)
    parser.add_argument('--memory',required=True)
    parser.add_argument('--dependency')
    args=parser.parse_args();repo=Path(__file__).resolve().parents[1]
    engine=args.run/('engine_'+args.version)
    shutil.copytree(repo/'src',engine/'src',ignore=shutil.ignore_patterns('__pycache__'))
    (engine/'bench').mkdir()
    for name in ('complete_golden_campaign.py','golden_campaign.py','golden_evaluation.py','report_golden_campaign.py'):
        shutil.copy2(repo/'bench'/name,engine/'bench'/name)
    manifest=dict(revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        files={str(p.relative_to(engine)):hashlib.sha256(p.read_bytes()).hexdigest()
               for p in engine.rglob('*') if p.is_file()},
        nodes=args.nodes,workers_per_node=args.workers,cut_workers=args.cut_workers,
        cpus_per_node=args.cpus,memory_per_node=args.memory,jobs=[])
    for shard,node in enumerate(args.nodes):
        command=['env',f'PYTHONPATH={engine}/src','RXN_CORE_NATIVE=1','OPENBLAS_NUM_THREADS=1',
            'OMP_NUM_THREADS=1','MKL_NUM_THREADS=1',sys.executable,str(engine/'bench/complete_golden_campaign.py'),
            'pool','--run',str(args.run),'--shard',str(shard),'--shards',str(len(args.nodes)),
            '--workers',str(args.workers),'--cut-workers',str(args.cut_workers)]
        job=subprocess.check_output(['sbatch','--parsable','--partition=cpunodes,cpunodes_nia',
            f'--nodelist={node}',f'--cpus-per-task={args.cpus}',f'--mem={args.memory}','--time=00:30:00',
            f'--job-name=golden_finish_{args.version}',f'--chdir={repo}',
            f'--output={args.run}/pool_{args.version}_{shard}.out',f'--error={args.run}/pool_{args.version}_{shard}.err',
            *([f'--dependency={args.dependency}'] if args.dependency else []),
            '--wrap',shlex.join(command)],text=True).strip().split(';')[0]
        manifest['jobs'].append(job);print(job,flush=True)
        (engine/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':main()
