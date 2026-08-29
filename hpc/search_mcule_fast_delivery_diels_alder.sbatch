#!/bin/bash
#SBATCH --job-name=mcule_diels_alder
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --output=data/retro_runs/mcule_fast_class2_diels_alder/logs/%A_%a.out
#SBATCH --error=data/retro_runs/mcule_fast_class2_diels_alder/logs/%A_%a.err

set -euo pipefail
repo_dir="/h/399/yunhengzou/coordinate_alignment"
run_dir="$repo_dir/data/retro_runs/mcule_fast_class2_diels_alder"
mkdir -p "$run_dir/logs" "$run_dir/parts"
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"$repo_dir/.venv/bin/python" tools/search_mcule_retro.py \
  --target-smiles 'CC12C=CC(C)(O1)C1C(=O)OC(=O)C12' \
  --catalog data/mcule/mcule_fast_delivery_bb_cls_2_260820_1H2y9f.csv.gz \
  --catalog-format csv \
  --output "$run_dir/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --shard-index "$SLURM_ARRAY_TASK_ID" \
  --shard-count "$SLURM_ARRAY_TASK_COUNT" \
  --batch-size 16 \
  --minimum-fragment-size 5 \
  --iso-tolerance 1.0 \
  --branch-limit 100 \
  --candidate-limit 100 \
  --maximum-boundary-bonds 1 \
  --maximum-leftover-fragments 1
