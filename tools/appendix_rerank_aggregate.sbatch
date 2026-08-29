#!/usr/bin/env bash
#SBATCH --job-name=bgcp_appendix_aggregate
#SBATCH --partition=cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/h/399/yunhengzou/appendix_final/bgcp_rerank_latest/logs/%x-%j.out
#SBATCH --error=/h/399/yunhengzou/appendix_final/bgcp_rerank_latest/logs/%x-%j.err

set -euo pipefail

cd /h/399/yunhengzou/coordinate_alignment

export PYTHONPATH=/h/399/yunhengzou/coordinate_alignment/src
export RXN_CORE_PROJECT=/h/399/yunhengzou/coordinate_alignment
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

/h/399/yunhengzou/micromamba/envs/el-agente/bin/python \
  tools/appendix_rerank.py aggregate \
  --run-root /h/399/yunhengzou/appendix_final/bgcp_rerank_latest
