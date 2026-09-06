#!/bin/bash
#SBATCH --job-name=beta_connected
#SBATCH --partition=cpunodes_nia,cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=80G
#SBATCH --time=00:10:00
set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
export PYTHONPATH=src RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun timeout --kill-after=10s 580s .venv/bin/python bench/scan_beta_catalog.py \
    --query "$1" --shard "$SLURM_ARRAY_TASK_ID" --workers "$SLURM_CPUS_PER_TASK"
