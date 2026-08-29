#!/bin/bash
#SBATCH --job-name=paa_full_cover
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-9
#SBATCH --output=data/retro_runs/phenylacetic_acid/full_fragment_cover/logs/%A_%a.out
#SBATCH --error=data/retro_runs/phenylacetic_acid/full_fragment_cover/logs/%A_%a.err

set -euo pipefail

repo_dir="/h/399/yunhengzou/coordinate_alignment"
run_dir="$repo_dir/data/retro_runs/phenylacetic_acid/full_fragment_cover"
mkdir -p "$run_dir/logs" "$run_dir/parts"
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"$repo_dir/.venv/bin/python" tools/search_mcule_retro.py \
  --target-smiles 'O=C(O)Cc1ccccc1' \
  --catalog data/mcule/mcule_fast_delivery_bb_cls_2_260820_augmented.csv.gz \
  --catalog-format csv \
  --output "$run_dir/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --shard-index "$SLURM_ARRAY_TASK_ID" \
  --shard-count "$SLURM_ARRAY_TASK_COUNT" \
  --batch-size 16 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100
