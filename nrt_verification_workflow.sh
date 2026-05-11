#!/bin/bash
# nrt_verification_workflow.sh
#
# Run the full BGCP verification pipeline (R↔P + R↔GT + 20× R↔IG cut_sweeps
# per step, score, rank) for all 155 step directories on the NRT cluster's
# cpu_short partition. Strategy:
#
#   - Slurm job array: ONE array task per step.
#   - Each task allocates N_CPUS cores and runs the inner-parallel pipeline:
#       python build_bgcp_views_v2.py --steps <step> --inner-workers $N_CPUS
#     so cut_sweep inside that step uses all N_CPUS cores.
#   - Slurm scheduler runs as many tasks concurrently as the partition allows;
#     the rest queue.
#
# CPU budget: with N_CPUS=5 per task and ~155 concurrent slots target, total
# core demand peaks at 5 * concurrency. Stay under the user's 880-core
# allotment by capping array concurrency at 880 / N_CPUS ≈ 176 (with N_CPUS=5).
# In practice on cpu_short most steps finish fast, so wall time is dominated
# by the few large-molecule steps (e.g. pr14 ~15 min with 14 inner workers).
#
# Usage:
#   bash nrt_verification_workflow.sh                      # run all steps
#   bash nrt_verification_workflow.sh --dry-run            # print the sbatch cmd, don't submit
#   N_CPUS=14 bash nrt_verification_workflow.sh            # 14 cores per step
#   STEPS_FILE=my_steps.txt bash nrt_verification_workflow.sh
#
# Output:
#   $LOG_DIR/step_<N>.{out,err}   per-task stdout/stderr
#   $PROJECT/out/bgcp_views/<step>/view.html
#   $PROJECT/out/bgcp_alignment_eval_v2.json   per-step record (one per task)
#
# Note: each task writes its own slim eval record. After all tasks complete,
# you may want to concatenate them into a single JSON; see the postprocess
# stub at the bottom of this script.

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJECT=${PROJECT:-/lustre/fsw/portfolios/nvr/users/yunhengz/Projects/coordinate_alignment}
WORK_DIR=$PROJECT/appendix_perparation/xtb_frequency_calculations
OUT_DIR=$PROJECT/out
LOG_DIR=$OUT_DIR/slurm_logs
PARTITION=${PARTITION:-cpu_short}
ACCOUNT=${ACCOUNT:-nvr_lpr_agentic}  # Slurm account (NRT requires --account)
N_CPUS=${N_CPUS:-5}                  # cores per step (inner-workers)
TIME_LIMIT=${TIME_LIMIT:-03:30:00}   # under cpu_short 4:00:00 cap
ARRAY_CONCURRENCY=${ARRAY_CONCURRENCY:-176}   # max concurrent tasks (876 / 5)
JOB_NAME=${JOB_NAME:-bgcp_align}

# Optional file listing steps to run, one per line.
# Defaults to "all step dirs in $WORK_DIR".
STEPS_FILE=${STEPS_FILE:-$OUT_DIR/_nrt_all_steps.txt}

DRY_RUN=0
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    *) echo "unknown arg: $arg" >&2; exit 2 ;;
  esac
done

# ---------------------------------------------------------------------------
# Pre-flight
# ---------------------------------------------------------------------------
mkdir -p "$OUT_DIR" "$LOG_DIR" "$OUT_DIR/bgcp_views"

if [ ! -d "$WORK_DIR" ]; then
  echo "ERROR: $WORK_DIR not found. Did you scp the xtb cache?" >&2
  exit 1
fi

# Build / refresh the step list if STEPS_FILE wasn't provided
if [ ! -f "$STEPS_FILE" ] || [ "$STEPS_FILE" = "$OUT_DIR/_nrt_all_steps.txt" ]; then
  find "$WORK_DIR" -mindepth 1 -maxdepth 1 -type d \
    | xargs -n1 basename | sort > "$STEPS_FILE"
fi
N_STEPS=$(wc -l < "$STEPS_FILE")
echo "Submitting $N_STEPS steps (file: $STEPS_FILE)"
echo "  Partition:    $PARTITION"
echo "  Account:      $ACCOUNT"
echo "  CPUs/step:    $N_CPUS"
echo "  Concurrency:  $ARRAY_CONCURRENCY"
echo "  Time limit:   $TIME_LIMIT"
echo "  Log dir:      $LOG_DIR"

# ---------------------------------------------------------------------------
# Submit Slurm job array. Each task reads its step name from $STEPS_FILE
# using SLURM_ARRAY_TASK_ID as the 1-based line number.
# ---------------------------------------------------------------------------
SBATCH_WRAP=$(cat <<'WRAP'
set -e
STEPS_FILE="__STEPS_FILE__"
PROJECT="__PROJECT__"
N_CPUS="__N_CPUS__"

STEP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$STEPS_FILE")
if [ -z "$STEP" ]; then
  echo "ERROR: empty step at task ${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

echo "[task ${SLURM_ARRAY_TASK_ID}] step=$STEP  cpus=$N_CPUS  host=$(hostname)  start=$(date -Is)"
cd "$PROJECT"

python -u build_bgcp_views_v2.py \
  --steps "$STEP" \
  --inner-workers "$N_CPUS"

echo "[task ${SLURM_ARRAY_TASK_ID}] step=$STEP  end=$(date -Is)"
WRAP
)
SBATCH_WRAP=${SBATCH_WRAP//__STEPS_FILE__/$STEPS_FILE}
SBATCH_WRAP=${SBATCH_WRAP//__PROJECT__/$PROJECT}
SBATCH_WRAP=${SBATCH_WRAP//__N_CPUS__/$N_CPUS}

CMD=(sbatch
  --partition="$PARTITION"
  --account="$ACCOUNT"
  --array="1-${N_STEPS}%${ARRAY_CONCURRENCY}"
  --cpus-per-task="$N_CPUS"
  --time="$TIME_LIMIT"
  --job-name="$JOB_NAME"
  --output="$LOG_DIR/step_%a.out"
  --error="$LOG_DIR/step_%a.err"
  --wrap "$SBATCH_WRAP"
)

if [ "$DRY_RUN" = "1" ]; then
  echo
  echo "DRY RUN — would submit:"
  printf '%q ' "${CMD[@]}"
  echo
  exit 0
fi

echo
echo "Submitting array job..."
"${CMD[@]}"
echo
echo "Submitted. Monitor with:"
echo "  squeue -u $USER -n $JOB_NAME"
echo "  tail -F $LOG_DIR/step_<N>.out"
echo
echo "When all tasks complete, concatenate per-step JSONs into one eval file:"
cat <<'POST'
  python - <<'PY'
  import json, glob, pathlib
  records = []
  for p in sorted(glob.glob("out/bgcp_views/*/*_eval_v2_slim.json")):
      records.append(json.loads(pathlib.Path(p).read_text()))
  pathlib.Path("out/bgcp_alignment_eval_v2.json").write_text(json.dumps(records))
  print(f"merged {len(records)} step records")
  PY
POST
