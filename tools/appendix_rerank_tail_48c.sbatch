#!/usr/bin/env bash
#SBATCH --job-name=bgcp_appendix_tail48
#SBATCH --partition=cpunodes,cpunodes_7d
#SBATCH --array=66-154%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=170G
#SBATCH --time=3-00:00:00
#SBATCH --output=/h/399/yunhengzou/appendix_final/bgcp_rerank_latest/logs/%x-%A_%a.out
#SBATCH --error=/h/399/yunhengzou/appendix_final/bgcp_rerank_latest/logs/%x-%A_%a.err

set -euo pipefail

cd /h/399/yunhengzou/coordinate_alignment
module load xtb/6.7.1-dev

export PYTHONPATH=/h/399/yunhengzou/coordinate_alignment/src
export RXN_CORE_PROJECT=/h/399/yunhengzou/coordinate_alignment
export TSDISCO_SLURM_CPUS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

/h/399/yunhengzou/micromamba/envs/el-agente/bin/python \
  tools/appendix_rerank.py run-step \
  --run-root /h/399/yunhengzou/appendix_final/bgcp_rerank_latest \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --xtb-mode auto \
  --xtb-omp-threads auto \
  --xtb-max-threads 8 \
  --xtb-workers auto \
  --fragment-separation 5.0 \
  --save-alignment-files
