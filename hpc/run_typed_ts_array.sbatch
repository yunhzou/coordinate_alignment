#!/usr/bin/env bash
#SBATCH --partition=cpunodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=360G
#SBATCH --time=24:00:00

set -euo pipefail

: "${RXN_CORE_ROOT:?set RXN_CORE_ROOT to an immutable source snapshot}"
: "${RXN_CORE_PYTHON:?set RXN_CORE_PYTHON to the project interpreter}"
: "${TS_MANIFEST:?set TS_MANIFEST to one tier manifest}"
: "${TS_OUTPUT_ROOT:?set TS_OUTPUT_ROOT}"
: "${TS_WORK_ROOT:?set TS_WORK_ROOT}"
: "${TS_BENCHMARK_ROOT:?set TS_BENCHMARK_ROOT}"

export PYTHONPATH="$RXN_CORE_ROOT/src"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONUNBUFFERED=1

"$RXN_CORE_PYTHON" \
  "$RXN_CORE_ROOT/tools/run_typed_ts_case.py" \
  --manifest "$TS_MANIFEST" \
  --work-root "$TS_WORK_ROOT" \
  --benchmark-root "$TS_BENCHMARK_ROOT" \
  --output-root "$TS_OUTPUT_ROOT" \
  --workers "$SLURM_CPUS_PER_TASK" \
  --post-workers "$SLURM_CPUS_PER_TASK"
