#!/bin/bash
#SBATCH --job-name=aam_memory_diagnosis
#SBATCH --partition=cpunodes_nia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:10:00

set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
export RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
srun timeout --kill-after=10s 580s .venv/bin/python bench/aam_memory_boundary.py "$@"
