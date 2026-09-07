#!/bin/bash
#SBATCH --job-name=golden_policy
#SBATCH --partition=cpunodes_nia,cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:10:00
set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
export PYTHONPATH="$1/engine/src:$1/engine/bench" RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python "$1/engine/bench/golden_policy_campaign.py" worker --run "$1" --slot "$((SLURM_ARRAY_TASK_ID + $2))"
