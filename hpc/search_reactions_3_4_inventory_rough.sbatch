#!/bin/bash
#SBATCH --job-name=retro_r3_r4_rough
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --array=0-55
#SBATCH --output=data/retro_runs/reactions_3_4_inventory_rough/logs/%A_%a.out
#SBATCH --error=data/retro_runs/reactions_3_4_inventory_rough/logs/%A_%a.err

set -euo pipefail

REPO_ROOT=/h/399/yunhengzou/coordinate_alignment
PYTHON_BIN=/h/399/yunhengzou/micromamba/envs/el-agente/bin/python
RUN_ROOT="$REPO_ROOT/data/retro_runs/reactions_3_4_inventory_rough"
SHARDS_PER_TARGET=28
TARGET_INDEX=$((SLURM_ARRAY_TASK_ID / SHARDS_PER_TARGET))
SHARD_INDEX=$((SLURM_ARRAY_TASK_ID % SHARDS_PER_TARGET))

TARGETS=(
  'CC(C=C(C)C=C1C)=C1/N=C(N(C2CCCCC2)C3CCCCC3)/NC4=C(C)C=C(C)C=C4C'
  'CP(C=P(NC1=C(C)C=CC=C1C)(C)C)(C)=NC2=C(C)C=CC=C2C'
)

TARGET_NAME=$((TARGET_INDEX == 0 ? 3 : 4))
TARGET_ROOT="$RUN_ROOT/reaction_${TARGET_NAME}"
mkdir -p "$TARGET_ROOT/parts"
cd "$REPO_ROOT"

PYTHONPATH=src srun "$PYTHON_BIN" tools/search_mcule_retro.py \
  --target-smiles "${TARGETS[$TARGET_INDEX]}" \
  --catalog data/inventory/processed/inventory_structure_bank.csv.gz \
  --catalog-format csv \
  --id-column "Inventory ID" \
  --output "$TARGET_ROOT/parts/part_${SHARD_INDEX}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --shard-index "$SHARD_INDEX" \
  --shard-count "$SHARDS_PER_TARGET" \
  --batch-size 1 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100 \
  --seed-mode fragment_cover \
  --rough-retention-threshold 0.5 \
  --save-all-results
