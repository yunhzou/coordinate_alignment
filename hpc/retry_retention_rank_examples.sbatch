#!/bin/bash
#SBATCH --job-name=retro_retention
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-31
#SBATCH --output=data/retro_runs/retention_rank_retry/logs/%A_%a.out
#SBATCH --error=data/retro_runs/retention_rank_retry/logs/%A_%a.err

set -euo pipefail
repo_dir="/h/399/yunhengzou/coordinate_alignment"
root_dir="$repo_dir/data/retro_runs/retention_rank_retry"
shards_per_target=8
target_index=$((SLURM_ARRAY_TASK_ID / shards_per_target))
shard_index=$((SLURM_ARRAY_TASK_ID % shards_per_target))

targets=(
  'Clc1ccccc1-c1ccccc1'
  'CC12C=CC(C)(O1)C1C(=O)OC(=O)C12'
  'CCOC(=O)c1ccc(-c2ccc(C(F)(F)F)cc2)cc1'
  'N1(C(N2C=CC=N2)N3N=CC=C3)N=CC=C1'
)

run_dir="$root_dir/target_$((target_index + 1))"
mkdir -p "$run_dir/parts"
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"$repo_dir/.venv/bin/python" tools/search_mcule_retro.py \
  --target-smiles "${targets[$target_index]}" \
  --catalog data/mcule/mcule_fast_delivery_bb_cls_2_260820_deduplicated.csv.gz \
  --catalog-format csv \
  --output "$run_dir/parts/part_${shard_index}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --shard-index "$shard_index" \
  --shard-count "$shards_per_target" \
  --batch-size 8 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100
