#!/bin/bash
#SBATCH --job-name=paa_augmented
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=0-9
#SBATCH --output=data/retro_runs/phenylacetic_acid/augmented_full/logs/%A_%a.out
#SBATCH --error=data/retro_runs/phenylacetic_acid/augmented_full/logs/%A_%a.err

set -euo pipefail

repo_dir="/h/399/yunhengzou/coordinate_alignment"
run_dir="$repo_dir/data/retro_runs/phenylacetic_acid/augmented_full"
target_smiles='O=C(O)Cc1ccccc1'
coverage_fraction='0.5'
mkdir -p "$run_dir/logs" "$run_dir/parts"
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

minimum_fragment_size=$("$repo_dir/.venv/bin/python" -c \
  "import math; from rdkit import Chem; print(math.ceil(float('$coverage_fraction') * Chem.AddHs(Chem.MolFromSmiles('$target_smiles')).GetNumAtoms()))")

"$repo_dir/.venv/bin/python" tools/search_mcule_retro.py \
  --target-smiles "$target_smiles" \
  --catalog data/mcule/mcule_fast_delivery_bb_cls_2_260820_augmented.csv.gz \
  --catalog-format csv \
  --output "$run_dir/parts/part_${SLURM_ARRAY_TASK_ID}.jsonl.gz" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --shard-index "$SLURM_ARRAY_TASK_ID" \
  --shard-count "$SLURM_ARRAY_TASK_COUNT" \
  --batch-size 16 \
  --minimum-fragment-size "$minimum_fragment_size" \
  --minimum-target-coverage-fraction "$coverage_fraction" \
  --iso-tolerance 0.5 \
  --branch-limit 100 \
  --candidate-limit 100
