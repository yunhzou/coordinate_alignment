#!/bin/bash
#SBATCH --job-name=retro_regression
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --time=00:10:00
#SBATCH --array=0-27
#SBATCH --output=data/retro_runs/principled_bian_full_bank_20260905/logs/%A_%a.out
#SBATCH --error=data/retro_runs/principled_bian_full_bank_20260905/logs/%A_%a.err

set -euo pipefail
repo_root=/h/399/yunhengzou/coordinate_alignment
run_root=$1
target_smiles='CC(C)C1=CC=CC(C(C)C)=C1/N=C2/C(C3=C4C(C=CC=C42)=CC=C3)=N/C5=C(C(C)C)C=CC=C5C(C)C'
cd "$repo_root"
export PYTHONPATH="$repo_root/src"
export RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Same target, bank, sharding, workers and seeds as the 6.7655 s baseline.
# No candidate cap: measure the current workflow, not the old truncation.
# Ten minutes is an operational watchdog; timed-out scans are incomplete.
srun "$repo_root/.venv/bin/python" tools/search_mcule_retro.py \
  --target-smiles "$target_smiles" \
  --catalog data/inventory/processed/inventory_structure_bank.csv.gz \
  --catalog-format csv --id-column 'Inventory ID' \
  --output "$run_root/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" --scheduling adaptive \
  --shard-index "$SLURM_ARRAY_TASK_ID" --shard-count 28 --batch-size 1 \
  --minimum-fragment-size 1 --iso-tolerance 0.5 --branch-limit 100 \
  --seed-mode orbit_representatives --save-all-results
