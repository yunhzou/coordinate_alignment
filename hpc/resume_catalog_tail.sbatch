#!/bin/bash
#SBATCH --job-name=retro_tail_complete
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --time=01:00:00
#SBATCH --array=0-17

set -euo pipefail
repo_root=/h/399/yunhengzou/coordinate_alignment
cd "$repo_root"
export PYTHONPATH="$repo_root/src"
export RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun "$repo_root/.venv/bin/python" bench/resume_catalog_source.py \
  --prior-run "$1" --output-dir "$2" \
  --index "$SLURM_ARRAY_TASK_ID" --workers "$SLURM_CPUS_PER_TASK"
