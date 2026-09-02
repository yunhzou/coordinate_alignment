#!/bin/bash
#SBATCH --job-name=retro_fast_inventory
#SBATCH --partition=cpunodes_nia,cpunodes
#SBATCH --array=0-223
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=data/retro_runs/eight_target_merged_fast_delivery_with_inventory/logs/%A_%a.out
#SBATCH --error=data/retro_runs/eight_target_merged_fast_delivery_with_inventory/logs/%A_%a.err

set -euo pipefail

repo_root=/h/399/yunhengzou/coordinate_alignment
run_root="$repo_root/data/retro_runs/eight_target_merged_fast_delivery_with_inventory"
shard_count=28
target_index=$((SLURM_ARRAY_TASK_ID / shard_count))
shard_index=$((SLURM_ARRAY_TASK_ID % shard_count))
target_line=$((target_index + 1))
target_id=$(sed -n "${target_line}p" "$run_root/targets.tsv" | cut -f1)
target_smiles=$(sed -n "${target_line}p" "$run_root/targets.tsv" | cut -f2-)
target_root="$run_root/$target_id"

mkdir -p "$target_root/parts"
cd "$repo_root"
export PYTHONPATH="$repo_root/src"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

srun "$repo_root/.venv/bin/python" tools/search_mcule_retro.py \
  --target-smiles "$target_smiles" \
  --catalog data/mcule/merged_fast_delivery_with_inventory.csv.gz \
  --catalog-format csv \
  --id-column 'Bank ID' \
  --output "$target_root/parts/part_${shard_index}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --scheduling adaptive \
  --shard-index "$shard_index" \
  --shard-count "$shard_count" \
  --batch-size 8 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100 \
  --seed-mode orbit_representatives \
  --save-all-results
