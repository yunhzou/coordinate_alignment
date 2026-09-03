#!/bin/bash
set -euo pipefail

target_id=$1
shift
run_root=data/retro_runs/eight_target_merged_fast_delivery_with_inventory
target_smiles=$(awk -F '\t' -v id="$target_id" \
  '$1 == id {print substr($0, index($0, $2))}' "$run_root/targets.tsv")
target_root="$run_root/$target_id"
ground_truth_status=$(awk -F '\t' -v id="$target_id" \
  '$1 == id {print $2}' "$run_root/ground_truth.tsv")
ground_truth_note=$(awk -F '\t' -v id="$target_id" \
  '$1 == id {print $3}' "$run_root/ground_truth.tsv")

merge_start=$(date +%s.%N)
.venv/bin/python tools/merge_retro_catalog.py \
  --parts "$target_root/parts" \
  --target-smiles "$target_smiles" \
  --output "$target_root/results.json" \
  --assembly-limit 20 \
  --maximum-precursors 6 \
  --index-workers "${SLURM_CPUS_PER_TASK:-1}" \
  "$@" \
  > "$target_root/merge.log" 2>&1
merge_end=$(date +%s.%N)

.venv/bin/python tools/build_retro_db_viewer.py \
  --results "$target_root/results.json" \
  --output "$target_root/viewer.html" \
  --top 20 \
  --title "$target_id" \
  --ground-truth-status "$ground_truth_status" \
  --ground-truth-note "$ground_truth_note" \
  > "$target_root/viewer.log" 2>&1
viewer_end=$(date +%s.%N)

awk -v merge_start="$merge_start" -v merge_end="$merge_end" \
    -v viewer_end="$viewer_end" \
    'BEGIN {printf "merge_seconds\t%.6f\nviewer_seconds\t%.6f\n", merge_end-merge_start, viewer_end-merge_end}' \
    > "$target_root/postprocess_time.tsv"
