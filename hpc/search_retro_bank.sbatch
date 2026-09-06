#!/bin/bash
#SBATCH --job-name=t05_bank_tol1
#SBATCH --partition=cpunodes_nia,cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=80G
#SBATCH --time=00:10:00

set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
run_root=$1
target_smiles=$2
shard_count=$3
export PYTHONPATH=src RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun timeout --kill-after=10s 580s .venv/bin/python tools/search_mcule_retro.py \
  --target-smiles "$target_smiles" \
  --catalog data/mcule/merged_fast_delivery_with_inventory.csv.gz \
  --catalog-format csv --id-column 'Bank ID' \
  --output "$run_root/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --checkpoint-directory "$run_root/checkpoints" \
  --workers "$SLURM_CPUS_PER_TASK" --scheduling adaptive \
  --shard-index "$SLURM_ARRAY_TASK_ID" --shard-count "$shard_count" --batch-size 1 \
  --minimum-fragment-size 1 --branch-limit 100 --seed-mode all --save-all-results
