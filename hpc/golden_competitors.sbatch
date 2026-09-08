#!/bin/bash
#SBATCH --job-name=gold_competitor
#SBATCH --partition=cpunodes_nia,cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=02:00:00
set -euo pipefail
cd /h/399/yunhengzou/coordinate_alignment
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=""
exec "$1/bin/python" "${5:-$2/golden_competitors.py}" worker --run "$2" --method "$3" --shard "${SLURM_ARRAY_TASK_ID}" --shards "$4"
