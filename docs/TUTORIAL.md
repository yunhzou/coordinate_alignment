# rxn_core Tutorial

This tutorial shows the three resumable stages and the composed full pipeline.
The same command is installed as both `rxn-core` and `rxn-core-pipeline`; new
workflows should use `rxn-core`.

## 0. Setup

Install in a conda environment:

```bash
conda create -n rxn-core python=3.12 -y
conda activate rxn-core
conda install -c conda-forge xtb -y
python -m pip install -e .
```

Use an existing xtb cache when possible:

```bash
export BGCP_WORK=/path/to/xtb_frequency_calculations
```

Expected cached step layout:

```text
<BGCP_WORK>/<step>/
  R/                    reactant-complex xyz + wbo
  P/                    product-complex xyz + wbo
  sp_iter<N>/           optional IG xyz + wbo
  hess_iter<N>/         optional IG g98.out
  sp_groundtruth/       optional GT xyz + wbo
  hess_groundtruth/     optional GT g98.out
```

Use `--xtb-mode cache-only` to guarantee the run does not launch xtb. Use
`--xtb-mode auto` to fill missing `wbo` or `g98.out` files from available XYZ
inputs.

## 1. Stage 1: R-P Alignment And Mechanisms

Stage 1 discovers mechanism-dependent R-P alignments.

```bash
rxn-core --stage rp --steps <step> \
  --workers 8 \
  --parallel-mode inner \
  --inner-workers 8 \
  --xtb-mode cache-only
```

Main output:

```text
out/bgcp_stages/<step>/rp_stage.json
out/bgcp_views/<step>/view.html
out/bgcp_views/<step>/_eval_slim.json
```

`rp_stage.json` contains:

```text
mechanisms[]
  id
  mapping_RP
  broken_bonds_R
  formed_bonds_R
  formed_bonds_P
  core_atoms
  product_xyz_in_R
```

For NEB/path setup, also write clean aligned coordinate files:

```bash
rxn-core --stage rp --steps <step> \
  --save-alignment-files \
  --xtb-mode cache-only
```

This adds:

```text
out/bgcp_alignments/<step>/
  manifest.json
  R.xyz
  P_original.xyz
  mechanisms/
    mechanism_001/
      R.xyz
      P_aligned.xyz
      neb_endpoints.xyz
      mapping_R_to_P.csv
      mechanism.json
```

`P_aligned.xyz` is product geometry reindexed into the R atom order for that
mechanism. No Kabsch/spatial fitting is applied.

## 2. Stage 2: Verify TS / IG / GT Against Mechanisms

Stage 2 resumes from `rp_stage.json`, maps each GT/IG/TS target through each
selected mechanism, and scores the selected imaginary mode.

Verify all mechanisms:

```bash
rxn-core --stage ts --steps <step> \
  --include-gt \
  --workers 8 \
  --parallel-mode inner \
  --inner-workers 8 \
  --xtb-mode cache-only
```

Verify only one mechanism:

```bash
rxn-core --stage ts --steps <step> \
  --mechanism 2 \
  --include-gt \
  --xtb-mode cache-only
```

Main output:

```text
out/bgcp_stages/<step>/ts_stage.json
out/bgcp_views/<step>/view.html
out/bgcp_views/<step>/_eval_slim.json
```

`ts_stage.json` contains per-mechanism scores:

```text
mechanisms[]
  gt
  igs[]
    S
    beta
    rho
    kappa
    freq
    n_imag
    core_map
    core_sources
    xyz
    picked_disp
    xyz_in_R
    picked_disp_R
```

To export the selected best-S TS materializations:

```bash
rxn-core --stage ts --steps <step> \
  --mechanism 2 \
  --include-gt \
  --save-alignment-files \
  --xtb-mode cache-only
```

This adds:

```text
out/bgcp_alignments/<step>/ts_alignments/
  manifest.json
  mechanisms/
    mechanism_002/
      gt_GT/
        TS_native.xyz
        TS_core_aligned_R_frame.xyz
        picked_mode_R_frame.xyz
        score.json
      ig_iter1/
        TS_native.xyz
        TS_core_aligned_R_frame.xyz
        picked_mode_R_frame.xyz
        score.json
```

Stage 2 scores core mappings, not full spectator bijections.
`TS_core_aligned_R_frame.xyz` replaces mapped core atoms with selected target
coordinates and leaves spectators at the R endpoint. `TS_native.xyz` is the
original target atom order.

## 3. Stage 3: Regenerate Views

Stage 3 is presentation-only. It reloads existing stage artifacts and rewrites
the HTML/eval output.

```bash
rxn-core --stage view --steps <step> --include-gt
```

Use this after editing viewer code or after moving stage artifacts.

Output:

```text
out/bgcp_views/<step>/view.html
out/bgcp_views/<step>/_eval_slim.json
```

The HTML has a `Download` button that packages R/P, optional GT, IGs,
per-mechanism score metadata, and `viewer_data.json`.

## 4. Full Pipeline

The full pipeline composes Stage 1, Stage 2, and Stage 3 in one run:

```bash
rxn-core --stage full --steps <step> \
  --include-gt \
  --workers 8 \
  --parallel-mode inner \
  --inner-workers 8 \
  --xtb-mode cache-only
```

Full run with all coordinate exports:

```bash
rxn-core --stage full --steps <step> \
  --include-gt \
  --save-alignment-files \
  --workers 8 \
  --parallel-mode inner \
  --inner-workers 8 \
  --xtb-mode cache-only
```

This writes:

```text
out/bgcp_stages/<step>/rp_stage.json
out/bgcp_stages/<step>/ts_stage.json
out/bgcp_views/<step>/view.html
out/bgcp_views/<step>/_eval_slim.json
out/bgcp_alignments/<step>/...
```

## 5. Python API

Run the same stages in Python:

```python
from rxn_core.pipeline import (
    load_step_inputs,
    load_ts_targets,
    run_rp_stage,
    run_ts_stage,
    write_view_stage,
    write_rp_alignment_files,
    write_ts_alignment_files,
)

step = "<step>"
inputs = load_step_inputs(step)

rp = run_rp_stage(inputs, inner_workers=8)
write_rp_alignment_files(inputs, rp)

targets = load_ts_targets(inputs, include_gt=True)
ts = run_ts_stage(inputs, rp, targets, mechanism_ids=[2], inner_workers=8)
write_ts_alignment_files(inputs, ts)

view = write_view_stage(inputs, rp, ts, include_gt=True)
print(view["view_html"])
```

For in-memory molecules instead of a BGCP cache:

```python
from rxn_core.pipeline import (
    step_inputs_from_arrays,
    ts_target_from_arrays,
    run_rp_stage,
    run_ts_stage,
)

inputs = step_inputs_from_arrays(
    "my_step",
    elR, xyzR, wboR,
    elP, xyzP, wboP,
)
rp = run_rp_stage(inputs)

target = ts_target_from_arrays(
    "ig", "guess1",
    elT, xyzT, wboT,
    freqs, modes,
)
ts = run_ts_stage(inputs, rp, [target])
```

## 6. Useful Runtime Knobs

Common controls:

```text
--workers                 total CPU budget in auto mode
--parallel-mode           auto | outer | inner
--inner-workers           per-step inner workers for R-P sweep and TS matching
--include-gt              score sp_groundtruth/hess_groundtruth
--mechanism ID            restrict Stage 2 to selected mechanism IDs
--save-alignment-files    write clean coordinate exports
--xtb-mode cache-only     never run xtb
--xtb-mode auto           fill missing cache files with xtb
--charge                  charge used only for xtb cache fill
--multiplicity            spin multiplicity used only for xtb cache fill
```

Important hypothesis knobs:

```text
--iso-tol
--dwbo-threshold
--symmetry-wbo-tol
--w-rxn
--w-core
--imag-pen
```

The defaults are documented in `README.md`.
