# Example 5: current-workflow correctness check — FAILED

This is a three-entry, bank-verified known-ingredient test, **not a blind
full-bank scan**. No mappings, target regions, or reaction-specific matching
rules were supplied. Core AAM, detection, and assembly code were not changed.

The exact bank rows and raw SMILES were verified against
`data/mcule/merged_fast_delivery_with_inventory.csv.gz`.

| Known raw ingredient | Bank ID | Detection seconds | Saved candidates | Cap hit |
| --- | --- | ---: | ---: | --- |
| Acetylacetone | INVENTORY-000283 | 0.958 | 402 | Yes |
| 2,6-Diisopropylaniline | INVENTORY-000436 | 0.352 | 11 | No |
| 2-Amino-6-cyanobenzothiazole | MCULE-7889053722 | 1.372 | 1,475 | Yes |

Settings: explicit H, tolerance 0.5, branch cap 100, all source seeds,
no sweep, no candidate cap. Three simultaneous Slurm tasks with 48 CPUs each
(432170_0–2). Detection plus checkpoint/encoding took 0.54–2.46 seconds inside
the scanners; Slurm task durations including startup were 6–9 seconds.

Production assembly (432173) returned zero complete covers. The index took
33.73 seconds; the complete merge/viewer Slurm task took 68 seconds. The
reported 0.005-second ranked-stream traversal excludes index persistence and
decision-graph preparation, so it is not the total post-processing time.

## Correctness finding

The union of **all indexed correlated occupations**, including alternative
target placements, covers only 54/58 atoms. P0, P15, P16, and P17 (all carbon)
have no support. Thus no number of repeated copies can complete this target
from these saved detections. This is not a ranking exclusion or bank absence.

For visual inspection, a separate diagnostic maximizes coverage using exactly
one copy of each known ingredient. It selects complete saved occupation
relations, not atomwise mixtures; tied coverage prefers less overlap. Its best
coverage is 53/58, missing P0, P1, P15, P16, and P17. This diagnostic is **not**
inserted into the recommendation stream or labeled as a successful assembly.

Acetylacetone's 15 initial searches did not hit the cap; their recorded carbon
placements are outside the expected central backbone. Its cap hit is later in
augmentation. Therefore the initial placement gap cannot be attributed simply
to an initial branch-cap cutoff. Growth commits extensions shared by surviving
candidates; the detector then anchors those returned fragments for augmentation.
This check does not establish whether that behavior is a historical regression.

The earlier passing `progressive_fragment_matching` t05 unit test rematches
against residual target graphs. It is a different workflow and does not prove
this independent-per-R production scan recovers the ground truth.

## Saved artifacts

Cluster directory:
`/h/399/yunhengzou/coordinate_alignment/data/retro_runs/t05_current_correctness_20260905/`

- `viewer.html`: three R panels and one combined P; five missing atoms labeled.
- `results.json`: unchanged failed recommendation result plus explicitly
  separate `diagnostic_assembly` with original detection references.
- `correctness_audit.json`: indexed support and original per-seed placements.
- `checkpoints/`: complete typed AAM/detection results for all three sources.
- `parts/`: full v7 archives; `results.occupations.json`: saved occupation index.
- `provenance.json`: original merged-bank rows.

Reproduction: `bench/check_retro_ground_truth.py` with
`docs/example_runs/t05_ground_truth.json`; Slurm wrapper:
`hpc/check_retro_ground_truth.sbatch`. Diagnostic/viewer rebuilding reuses the
saved index and never reruns AAM.
