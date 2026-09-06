#!/bin/bash
#SBATCH --job-name=beta_index
#SBATCH --partition=cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:10:00
set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
export PYTHONPATH=src RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun timeout --kill-after=10s 580s .venv/bin/python bench/index_beta_query.py \
    --query "$1" --shard "$SLURM_ARRAY_TASK_ID"
