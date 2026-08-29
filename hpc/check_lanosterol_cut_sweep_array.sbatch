#!/bin/bash
#SBATCH --job-name=lan_cut_array
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G
#SBATCH --time=00:12:00
#SBATCH --array=0-9
#SBATCH --output=data/retro_runs/lanosterol_cut_sweep_array/%A_%a.out
#SBATCH --error=data/retro_runs/lanosterol_cut_sweep_array/%A_%a.err

set -euo pipefail
repo_dir="/h/399/yunhengzou/coordinate_alignment"
cd "$repo_dir"
mkdir -p data/retro_runs/lanosterol_cut_sweep_array
export PYTHONPATH="$repo_dir/src:$repo_dir/tools"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run_id="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
"$repo_dir/.venv/bin/python" tools/check_explicit_h_cut_sweep.py \
  --source-smiles 'CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C=C(\C)CC/C=C(\C)CC[C@@H]1OC1(C)C' \
  --target-smiles 'C[C@H](CCC=C(C)C)[C@H]1CC[C@@]2([C@@]1(CCC3=C2CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C' \
  --cut-shard-index "$SLURM_ARRAY_TASK_ID" \
  --cut-shard-count 10 \
  --n-seeds 10 \
  --branch-limit 100 \
  --workers "$SLURM_CPUS_PER_TASK" \
  --intermediate-dir "data/retro_runs/lanosterol_cut_sweep_array/intermediate_${run_id}" \
  --best-result "data/retro_runs/lanosterol_cut_sweep_array/best_${run_id}.json" \
  --all-results "data/retro_runs/lanosterol_cut_sweep_array/aam_pool_${run_id}.json.gz"
