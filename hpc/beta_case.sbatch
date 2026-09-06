#!/bin/bash
#SBATCH --job-name=beta_case
#SBATCH --partition=cpunodes_nia,cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=08:00:00
set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
export PYTHONPATH=src RXN_CORE_NATIVE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python bench/run_beta_case.py --case "$1" --cpu-budget "$2"
