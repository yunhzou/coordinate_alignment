#!/bin/bash
#SBATCH --job-name=lan_win_cuts
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --output=data/retro_runs/lanosterol_cut_sweep/%j_winner.out
#SBATCH --error=data/retro_runs/lanosterol_cut_sweep/%j_winner.err

set -euo pipefail
repo_dir="/h/399/yunhengzou/coordinate_alignment"
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src:$repo_dir/tools"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
seed_count="${AAM_SEEDS:-10}"
run_label="${AAM_LABEL:-winning_cut}"
"$repo_dir/.venv/bin/python" tools/check_explicit_h_cut_sweep.py \
  --source-smiles 'CC(C)=CCC/C(C)=C/CC/C(C)=C/CC/C=C(\C)CC/C=C(\C)CC[C@@H]1OC1(C)C' \
  --target-smiles 'C[C@H](CCC=C(C)C)[C@H]1CC[C@@]2([C@@]1(CCC3=C2CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C' \
  --only-cut 13-14 \
  --n-seeds "$seed_count" \
  --branch-limit 100 \
  --workers "$SLURM_CPUS_PER_TASK" \
  --intermediate-dir "data/retro_runs/lanosterol_cut_sweep/intermediate_${SLURM_JOB_ID}" \
  --best-result "data/retro_runs/lanosterol_cut_sweep/${run_label}_best_mapping.json" \
  --all-results "data/retro_runs/lanosterol_cut_sweep/${run_label}_aam_pool.json.gz"
