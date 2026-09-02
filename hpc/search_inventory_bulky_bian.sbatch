#!/bin/bash
#SBATCH --job-name=bian_inventory
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --array=0-27
#SBATCH --output=data/retro_runs/bulky_bian_inventory_rough/logs/%A_%a.out
#SBATCH --error=data/retro_runs/bulky_bian_inventory_rough/logs/%A_%a.err

set -euo pipefail

REPO_ROOT=/h/399/yunhengzou/coordinate_alignment
PYTHON_BIN=/h/399/yunhengzou/micromamba/envs/el-agente/bin/python
RUN_ROOT="$REPO_ROOT/data/retro_runs/bulky_bian_inventory_rough"
TARGET_SMILES=$(cut -d ' ' -f1 "$REPO_ROOT/docs/example_runs/bulky_bian_target.smi")

mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/parts"
cd "$REPO_ROOT"

PYTHONPATH=src srun "$PYTHON_BIN" tools/search_mcule_retro.py \
  --target-smiles "$TARGET_SMILES" \
  --catalog data/inventory/processed/inventory_structure_bank.csv.gz \
  --catalog-format csv \
  --id-column "Inventory ID" \
  --output "$RUN_ROOT/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --shard-index "$SLURM_ARRAY_TASK_ID" \
  --shard-count "$SLURM_ARRAY_TASK_COUNT" \
  --batch-size 1 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100 \
  --seed-mode fragment_cover \
  --rough-retention-threshold 0.5 \
  --save-all-results
