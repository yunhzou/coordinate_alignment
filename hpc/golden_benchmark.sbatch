#!/bin/bash
#SBATCH --job-name=golden_aam
#SBATCH --partition=cpunodes_nia,cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
export PYTHONPATH="$1/engine/src" RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python "$1/engine/bench/golden_campaign.py" worker --run "$1" --slot "$((SLURM_ARRAY_TASK_ID + $2))"
