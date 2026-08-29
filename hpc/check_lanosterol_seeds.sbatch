#!/bin/bash
#SBATCH --job-name=lan_seed_check
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --array=0-2
#SBATCH --output=data/retro_runs/lanosterol_seed_check/%A_%a.out
#SBATCH --error=data/retro_runs/lanosterol_seed_check/%A_%a.err

set -euo pipefail
repo_dir="/h/399/yunhengzou/coordinate_alignment"
seeds=(6 12 24)
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src:$repo_dir/tools"
"$repo_dir/.venv/bin/python" tools/check_explicit_h_pair.py \
  --source-smiles 'CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C=C(\C)CC/C=C(\C)CC[C@@H]1OC1(C)C' \
  --target-smiles 'C[C@H](CCC=C(C)C)[C@H]1CC[C@@]2([C@@]1(CCC3=C2CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C' \
  --n-seeds "${seeds[$SLURM_ARRAY_TASK_ID]}" \
  --branch-limit 100
