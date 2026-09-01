#!/bin/bash
#SBATCH --job-name=vancomycin_inventory
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=160G
#SBATCH --time=04:00:00
#SBATCH --array=0-3
#SBATCH --output=data/retro_runs/vancomycin_inventory/logs/%A_%a.out
#SBATCH --error=data/retro_runs/vancomycin_inventory/logs/%A_%a.err

set -euo pipefail

REPO_ROOT=/h/399/yunhengzou/coordinate_alignment
PYTHON_BIN=/h/399/yunhengzou/micromamba/envs/el-agente/bin/python
RUN_ROOT="$REPO_ROOT/data/retro_runs/vancomycin_inventory"
TARGET_SMILES=$(cut -f1 "$REPO_ROOT/docs/example_runs/vancomycin_pubchem_14969.smi")

mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/parts"
cd "$REPO_ROOT"

PYTHONPATH=src srun "$PYTHON_BIN" tools/search_mcule_retro.py \
  --target-smiles "$TARGET_SMILES" \
  --catalog data/inventory/processed/inventory_structure_bank.csv.gz \
  --catalog-format csv \
  --id-column "Inventory ID" \
  --output "$RUN_ROOT/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --workers 48 \
  --shard-index "$SLURM_ARRAY_TASK_ID" \
  --shard-count "$SLURM_ARRAY_TASK_COUNT" \
  --batch-size 4 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100 \
  --save-all-results
