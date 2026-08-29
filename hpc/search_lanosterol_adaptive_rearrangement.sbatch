#!/bin/bash
#SBATCH --job-name=lano_adapt_aam
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=data/retro_runs/lanosterol_current_workflow/rearrangement_tier/logs/%j.out
#SBATCH --error=data/retro_runs/lanosterol_current_workflow/rearrangement_tier/logs/%j.err

set -euo pipefail
repo_dir="/h/399/yunhengzou/coordinate_alignment"
run_dir="$repo_dir/data/retro_runs/lanosterol_current_workflow/rearrangement_tier"
mkdir -p "$run_dir/logs" "$run_dir/parts"
cd "$repo_dir"
export PYTHONPATH="$repo_dir/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

"$repo_dir/.venv/bin/python" tools/search_mcule_rearrangement.py \
  --target-smiles 'C[C@H](CCC=C(C)C)[C@H]1CC[C@@]2([C@@]1(CCC3=C2CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C' \
  --catalog data/retro_runs/lanosterol_current_workflow/exact_composition_catalog.csv.gz \
  --id-column 'compound mcule ID' \
  --output "$run_dir/parts/part_0.jsonl.gz" \
  --workers 30 \
  --batch-size 1 \
  --minimum-excess-atoms 0 \
  --maximum-excess-atoms 0 \
  --maximum-events 30 \
  --branch-limit 100 \
  --no-chirality
