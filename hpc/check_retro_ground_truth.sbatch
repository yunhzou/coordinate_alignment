#!/bin/bash
#SBATCH --job-name=t05_correctness
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --time=00:10:00

set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
export PYTHONPATH=src RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source_args=()
if [[ "$1" == detect ]]; then
  source_args=(--source-index "$SLURM_ARRAY_TASK_ID")
fi
srun timeout --kill-after=10s 580s .venv/bin/python bench/check_retro_ground_truth.py \
  "$1" --case "$2" --output-dir "$3" --workers "$SLURM_CPUS_PER_TASK" "${source_args[@]}" "${@:4}"
