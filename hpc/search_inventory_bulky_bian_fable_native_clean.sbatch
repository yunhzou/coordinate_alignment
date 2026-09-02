#!/bin/bash
#SBATCH --job-name=bian_fable_clean
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --array=0-27
#SBATCH --output=data/retro_runs/bulky_bian_fable_symmetry_orbit/logs/%A_%a.out
#SBATCH --error=data/retro_runs/bulky_bian_fable_symmetry_orbit/logs/%A_%a.err

set -euo pipefail

REPO_ROOT=/h/399/yunhengzou/coordinate_alignment
PYTHON_BIN="$REPO_ROOT/.venv/bin/python"
RUN_ROOT="$REPO_ROOT/data/retro_runs/bulky_bian_fable_symmetry_orbit"
TARGET_SMILES="CC(C)C1=CC=CC(C(C)C)=C1/N=C2/C(C3=C4C(C=CC=C42)=CC=C3)=N/C5=C(C(C)C)C=CC=C5C(C)C"
SHARD_COUNT=28

cd "$REPO_ROOT"

PYTHONPATH=src RXN_CORE_NATIVE=1 srun "$PYTHON_BIN" tools/search_mcule_retro.py \
  --target-smiles "$TARGET_SMILES" \
  --catalog data/inventory/processed/inventory_structure_bank.csv.gz \
  --catalog-format csv \
  --id-column "Inventory ID" \
  --output "$RUN_ROOT/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --scheduling precursor \
  --shard-index "$SLURM_ARRAY_TASK_ID" \
  --shard-count "$SHARD_COUNT" \
  --batch-size 1 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100 \
  --seed-mode all \
  --save-all-results
