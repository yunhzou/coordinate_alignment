#!/bin/bash
#SBATCH --job-name=retro_norv_gap
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=120G
#SBATCH --time=00:10:00
#SBATCH --array=0-7
#SBATCH --output=data/retro_runs/%x_%A_%a.out
#SBATCH --error=data/retro_runs/%x_%A_%a.err

set -euo pipefail

REPO_ROOT=/h/399/yunhengzou/coordinate_alignment
PYTHON_BIN=/h/399/yunhengzou/micromamba/envs/el-agente/bin/python
TARGET_SMILES=$(cut -f1 "$REPO_ROOT/docs/example_runs/n_demethylvancomycin_mcule_3152543364.smi")
RUN_ROOT="$REPO_ROOT/data/retro_runs/$RUN_NAME"

mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/parts"
cd "$REPO_ROOT"

PYTHONPATH=src srun "$PYTHON_BIN" tools/search_mcule_retro.py \
  --target-smiles "$TARGET_SMILES" \
  --catalog "$CATALOG" \
  --catalog-format csv \
  --id-column "Mcule ID" \
  --output "$RUN_ROOT/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --workers 48 \
  --shard-index "$SLURM_ARRAY_TASK_ID" \
  --shard-count "$SLURM_ARRAY_TASK_COUNT" \
  --batch-size 1 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100 \
  --seed-limit 10 \
  --target-region-report "$REGION_REPORT" \
  --target-region-field best_partial_uncovered_target_atoms \
  --save-all-results
