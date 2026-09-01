#!/bin/bash
#SBATCH --job-name=vanco_methyl_donors
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=data/retro_runs/vancomycin_methyl_donors/search.out
#SBATCH --error=data/retro_runs/vancomycin_methyl_donors/search.err

set -euo pipefail

REPO_ROOT=/h/399/yunhengzou/coordinate_alignment
PYTHON_BIN=/h/399/yunhengzou/micromamba/envs/el-agente/bin/python
RUN_ROOT="$REPO_ROOT/data/retro_runs/vancomycin_methyl_donors"
TARGET_SMILES=$(cut -f1 "$REPO_ROOT/docs/example_runs/vancomycin_pubchem_14969.smi")

mkdir -p "$RUN_ROOT"
cp "$REPO_ROOT/docs/example_runs/vancomycin_methyl_donors.csv" \
  "$RUN_ROOT/donors.csv"
gzip -f "$RUN_ROOT/donors.csv"
cd "$REPO_ROOT"

PYTHONPATH=src srun "$PYTHON_BIN" tools/search_mcule_retro.py \
  --target-smiles "$TARGET_SMILES" \
  --catalog "$RUN_ROOT/donors.csv.gz" \
  --catalog-format csv \
  --id-column "Mcule ID" \
  --output "$RUN_ROOT/part_0.jsonl.gz" \
  --workers 3 \
  --batch-size 1 \
  --minimum-fragment-size 1 \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100 \
  --seed-limit 10 \
  --save-all-results
