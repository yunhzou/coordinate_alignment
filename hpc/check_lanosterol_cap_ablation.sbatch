#!/bin/bash
#SBATCH --job-name=lan_cap_ablate
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=24G
#SBATCH --time=00:08:00
#SBATCH --array=0-3
#SBATCH --output=data/retro_runs/lanosterol_cap_ablation/%A_%a.out
#SBATCH --error=data/retro_runs/lanosterol_cap_ablation/%A_%a.err

set -euo pipefail
repo_dir="/h/399/yunhengzou/coordinate_alignment"
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src:$repo_dir/tools"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
caps=(50 100 200 500)
cap="${caps[$SLURM_ARRAY_TASK_ID]}"
run_id="${SLURM_ARRAY_JOB_ID}_cap${cap}"
"$repo_dir/.venv/bin/python" tools/check_explicit_h_cut_sweep.py \
  --source-smiles 'CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C=C(\C)CC/C=C(\C)CC[C@@H]1OC1(C)C' \
  --target-smiles 'C[C@H](CCC=C(C)C)[C@H]1CC[C@@]2([C@@]1(CCC3=C2CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C' \
  --only-cut 13-14 \
  --only-cut 14-15 \
  --only-cut 16-17 \
  --only-cut 17-19 \
  --n-seeds 10 \
  --branch-limit "$cap" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --intermediate-dir "data/retro_runs/lanosterol_cap_ablation/intermediate_${run_id}" \
  --best-result "data/retro_runs/lanosterol_cap_ablation/best_${run_id}.json" \
  --all-results "data/retro_runs/lanosterol_cap_ablation/aam_pool_${run_id}.json.gz"
