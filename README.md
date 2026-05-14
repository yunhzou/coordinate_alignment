# rxn_core

Symmetry-aware WBO graph alignment and BGCP transition-state ranking.

The package discovers R-P mechanisms with sweep cut, aligns IG transition
states plus optional GT through mechanism-local core atoms, scores normal
modes, and writes a self-contained HTML view per step.

## Install

In a conda environment:

```bash
conda create -n rxn-core python=3.12 -y
conda activate rxn-core
conda install -c conda-forge xtb -y
python -m pip install -e .
```

`xtb` is only needed when `BGCP_XTB_MODE=auto` fills missing `wbo` or
`g98.out` cache files. Fully cached runs can use `--xtb-mode cache-only`.

Verify the installed import and command line:

```bash
python -c "import rxn_core; print(rxn_core.__file__)"
rxn-core --help
rxn-core-pipeline --help
```

## Usage

For a stage-by-stage walkthrough, open
[`docs/TUTORIAL.ipynb`](docs/TUTORIAL.ipynb).

```bash
# Direct XYZ mode: no benchmark step schema required
rxn-core --stage rp --name my_reaction \
  --reactant-xyz R.xyz --product-xyz P.xyz \
  --workdir work/my_reaction \
  --charge 0 --multiplicity 1 --xtb-mode auto

# Direct Stage 2 target verification after direct Stage 1 has written rp_stage.json
rxn-core --stage ts --name my_reaction \
  --reactant-xyz R.xyz --product-xyz P.xyz \
  --target-xyz guess_1.xyz --target-label guess_1 --target-kind ig \
  --workdir work/my_reaction --save-alignment-files

# Use the default cache root: data/xtb_frequency_calculations
rxn-core-pipeline --steps pr7.V.dodh_ts910 --workers 8

# Use an external cache root
BGCP_WORK=/path/to/xtb_frequency_calculations \
  rxn-core-pipeline --steps pr7.V.dodh_ts910 --workers 64

# Forbid xtb execution and require all caches to already exist
rxn-core-pipeline --steps pr7.V.dodh_ts910 --xtb-mode cache-only

# Include ground-truth TS scoring when GT cache directories are available
rxn-core-pipeline --steps pr7.V.dodh_ts910 --include-gt

# When auto-filling missing xtb caches for an open-shell/charged system
rxn-core-pipeline --steps pr15.example --charge 0 --multiplicity 4
```

`rxn-core` and `rxn-core-pipeline` are the same command.  The shorter
`rxn-core` name is intended for new workflows.

### Real Example

This direct-XYZ example uses `pr1.tempo_ts3` from the appendix benchmark. That
step has one complete reactant-complex XYZ and one complete product-complex XYZ,
so it can be passed directly to Stage 1.

```python
from pathlib import Path
import rxn_core.pipeline as rxnp

benchmark_root = Path("~/Downloads/appendix_final 2/benchmark").expanduser()
step_dir = benchmark_root / "pr1.tempo_ts3"

name = "pr1.tempo_ts3"
reactant_xyz = step_dir / "reactants/reactant_01_reactant_01_5.xyz"
product_xyz = step_dir / "products/product_01_product_01_6.xyz"
workdir = f"work/{name}"

target_specs = [
    {
        "kind": "ig",
        "label": "iter1",
        "xyz": step_dir / "initial_guess/pr1.tempo_ts3_benchmark_plain_iter1_87ea3b8f.xyz",
    },
    {
        "kind": "gt",
        "label": "GT",
        "xyz": step_dir / "groundtruth/ts_groundtruth_01_reference_ts_01_TS3.xyz",
    },
]

result = rxnp.process_xyz_stage(
    name,
    reactant_xyz,
    product_xyz,
    workdir=workdir,
    stage="full",
    target_specs=target_specs,
    charge=0,
    multiplicity=1,
    xtb_mode="auto",
    inner_workers=8,
    save_alignment_files=True,
)
```

For raw benchmark steps with multiple separate reactant/product fragments, first
assemble each endpoint into one reaction-complex XYZ before calling Stage 1.

The tutorial notebook has been executed on this example. Its stored artifacts
are under `docs/example_runs/pr1.tempo_ts3/`, including the generated
`view.html`, Stage 1/2 JSON, xtb caches, and aligned coordinate exports.
The same example is also prepared as a cache-only step folder at
`docs/example_runs/pr1.tempo_ts3/prepared_steps/pr1.tempo_ts3/` with `R/`,
`P/`, `sp_iter1/`, `hess_iter1/`, `sp_groundtruth/`, and
`hess_groundtruth/`.

### Staged Workflows

The pipeline is split into three resumable stages.

```bash
# Stage 1: R-P alignment and mechanism discovery only
rxn-core --stage rp --steps pr7.V.dodh_ts910 --workers 8

# Stage 1 plus clean mechanism-specific aligned R/P files for NEB/path setup
rxn-core --stage rp --steps pr7.V.dodh_ts910 --save-alignment-files

# Stage 2 plus selected best-S TS core-aligned files
rxn-core --stage ts --steps pr7.V.dodh_ts910 \
  --mechanism 2 --include-gt --save-alignment-files

# Stage 2: verify GT/IG/TS targets from the saved Stage 1 mechanisms
rxn-core --stage ts --steps pr7.V.dodh_ts910 --mechanism 2 --include-gt

# Stage 3: regenerate only the HTML/eval view from saved stage artifacts
rxn-core --stage view --steps pr7.V.dodh_ts910 --include-gt

# Full pipeline: compose rp + ts + view in one run
rxn-core --stage full --steps pr7.V.dodh_ts910 --workers 8
```

The same pieces are importable:

```python
import rxn_core.pipeline as rxnp

# Preferred Stage 1 API: arbitrary endpoint XYZ files plus charge/multiplicity.
inputs = rxnp.alignment_inputs_from_xyz(
    "R.xyz", "P.xyz",
    workdir="work/my_reaction",
    name="my_reaction",
    charge=0,
    multiplicity=1,
    xtb_mode="auto",
)
rp = rxnp.run_rp_stage(inputs, inner_workers=8)
rxnp.write_rp_alignment_files(inputs, rp)

# Preferred Stage 2 API: arbitrary TS/IG/GT XYZ files plus the same molecular
# charge/multiplicity.  The explicit cache directories are optional; they are
# shown here to make the data ownership clear.
targets = [
    rxnp.ts_target_from_xyz(
        "ig", "guess_1", "guess_1.xyz",
        sp_workdir="work/my_reaction/guess_1_sp",
        hess_workdir="work/my_reaction/guess_1_hess",
        charge=0,
        multiplicity=1,
        xtb_mode="auto",
    ),
]
ts = rxnp.run_ts_stage(inputs, rp, targets, mechanism_ids=[2], inner_workers=8)
rxnp.write_ts_alignment_files(inputs, ts)
view = rxnp.write_view_stage(inputs, rp, ts, include_gt=True)
```

`load_step_inputs(...)` and `load_ts_targets(...)` are benchmark adapters over
the same file-based API. For callers that already have WBO matrices and
normal modes in memory, `step_inputs_from_arrays(...)`,
`ts_target_from_arrays(...)`, and `discover_mechanisms_from_arrays(...)`
avoid filesystem cache loading entirely. The resulting mechanisms include
`mapping_RP` and `product_xyz_in_R`, so each mechanism has its own aligned
product coordinate frame.

## Inputs

The principal API input is:

| stage | required molecule data | cache-fill inputs |
|---|---|---|
| Stage 1 R-P mechanism discovery | R endpoint XYZ, P endpoint XYZ, charge, multiplicity | per-endpoint cache directory containing or receiving `wbo` |
| Stage 2 TS/IG/GT verification | TS/IG/GT XYZ, charge, multiplicity, Stage 1 mechanisms | one single-point cache directory containing or receiving `wbo`; one Hessian cache directory containing or receiving `g98.out` |
| Stage 3 view/export | Stage 1 result, optional Stage 2 result, original loaded molecules | no new chemistry calculation |

The benchmark step schema is a convenience wrapper. By default it looks in:

```text
data/xtb_frequency_calculations/<step>/
  R/                    reactant xyz + wbo
  P/                    product xyz + wbo
  sp_iter<N>/           IG xyz + wbo
  hess_iter<N>/         IG g98.out
  sp_groundtruth/       optional GT TS xyz + wbo
  hess_groundtruth/     optional GT g98.out
```

Use `BGCP_WORK=/path/to/xtb_frequency_calculations` when the cache lives
outside the repo.

Required inputs for IG ranking are:

| path | required | contents | used for |
|---|---:|---|---|
| `R/` | yes | one reactant-complex XYZ and `wbo` | R-P mechanism discovery and R-frame indexing |
| `P/` | yes | one product-complex XYZ and `wbo` | R-P mechanism discovery |
| `sp_iter<N>/` | yes for IG ranking | one IG TS XYZ and `wbo` | TS/IG core matching |
| `hess_iter<N>/` | yes for IG ranking | `g98.out` plus optional/copyable XYZ | imaginary-mode scoring |
| `sp_groundtruth/` | optional | GT TS XYZ and `wbo` | GT score/view when `--include-gt` is set |
| `hess_groundtruth/` | optional | GT `g98.out` plus optional/copyable XYZ | GT mode scoring when `--include-gt` is set |

Each endpoint directory is one complete molecule/complex graph. For
multi-reactant or multi-product cases, the fragments must already be present
in the single XYZ under `R/`, `P/`, `sp_groundtruth/`, or `sp_iter<N>/`; the
pipeline does not merge separate per-fragment xtb outputs. The atom order only
has to be internally consistent within each XYZ/WBO pair. R, P, GT, and IG
files do not need to share the same atom order because alignment computes the
cross-index mapping, but they must contain the same element composition.

Without `sp_iter<N>/` and `hess_iter<N>/`, the pipeline can still discover R-P
mechanisms and optionally score GT, but it is not validating or ranking initial
guesses.

GT is optional and disabled by default. Pass `--include-gt` or set
`BGCP_INCLUDE_GT=1` to load and score `sp_groundtruth/` plus
`hess_groundtruth/`.

By default `rxn-core-pipeline` runs in `BGCP_XTB_MODE=auto`: if an expected
`wbo` or `g98.out` file is missing and an XYZ is available, it checks for
`xtb` on `PATH` and fills the missing cache with `xtb --sp` or `xtb --hess`.
Use `BGCP_XTB_MODE=cache-only` or `--xtb-mode cache-only` to fail fast instead.
Each xtb subprocess uses `OMP_NUM_THREADS` capped by `BGCP_XTB_MAX_THREADS=8`
by default. The cache-fill input uses molecular `charge` and spin
`multiplicity`; the xtb adapter converts multiplicity to
`--uhf=multiplicity-1` internally.

## Outputs

```text
out/bgcp_views/<step>/view.html
out/bgcp_views/<step>/_eval_slim.json
out/bgcp_alignment_eval.json
out/bgcp_stages/<step>/rp_stage.json
out/bgcp_stages/<step>/ts_stage.json
out/bgcp_alignments/<step>/manifest.json
out/bgcp_alignments/<step>/mechanisms/mechanism_<id>/R.xyz
out/bgcp_alignments/<step>/mechanisms/mechanism_<id>/P_aligned.xyz
out/bgcp_alignments/<step>/mechanisms/mechanism_<id>/neb_endpoints.xyz
out/bgcp_alignments/<step>/ts_alignments/manifest.json
out/bgcp_alignments/<step>/ts_alignments/mechanisms/mechanism_<id>/<target>/TS_native.xyz
out/bgcp_alignments/<step>/ts_alignments/mechanisms/mechanism_<id>/<target>/TS_core_aligned_R_frame.xyz
out/bgcp_alignments/<step>/ts_alignments/mechanisms/mechanism_<id>/<target>/picked_mode_R_frame.xyz
```

Each `view.html` has a step-level `Download` button. The downloaded archive
contains `R.xyz`, `P.xyz`, unique `IG/<label>.xyz` files, optional
`GT/GT.xyz`, a root `mechanism.json` manifest with per-mechanism IG/GT scores
and score decomposition, one `mechanisms/mechanism_<id>.json` file per
mechanism, and `viewer_data.json` with the full data used by the HTML view.

When `--save-alignment-files` is set, Stage 1 also writes a clean aligned
coordinate package under `out/bgcp_alignments/<step>/`. Each mechanism
directory is self-contained and includes `R.xyz`, `P_aligned.xyz`, a
two-frame `neb_endpoints.xyz`, `mapping_R_to_P.csv`, and `mechanism.json`.
`P_aligned.xyz` is product geometry reindexed into the R atom order; no
Kabsch or spatial fitting is applied.

When Stage 2 is run with `--save-alignment-files`, the same output root gets
`ts_alignments/`. Each scored GT/IG/TS target has native target coordinates,
the selected best-S core-aligned R-frame materialization, the picked mode in
R-frame extended XYZ, and `score.json`. The Stage 2 aligned file is core-only:
mapped core atoms use the target coordinates, while spectator atoms remain at
the reactant endpoint because spectator bijections are intentionally not
enumerated.

Generated caches, views, and paper artifacts are intentionally not part of the
main repository.

## Pipeline

The public full pipeline is:

```text
process_step(...)
  run_rp_stage(...)
  run_ts_stage(...)
  write_view_stage(...)
```

1. `run_rp_stage(...)` / `cut_sweep(...)` runs R-P mechanism discovery:
   no-cut plus one-edge R cuts above `BGCP_CUT_FLOOR`.
2. Mechanisms are deduped by symmetry-canonical broken/formed bond changes.
3. The reactive core is the atoms touching any broken or formed bond.
4. `run_ts_stage(...)` / `ts_core_pool(...)` runs endpoint-to-TS matching from
   both R and P. P-side
   core mappings are pulled back through the R-P mechanism witness, then merged
   with R-side mappings in R-core indexing.
5. Each GT/IG candidate mapping is scored on the selected imaginary mode:

```text
S = beta * (1 + W_RXN * rho) * (1 + W_CORE * kappa) / n_imag^IMAG_PEN
```

where `beta` is bond-axis overlap, `rho` is reaction-coordinate overlap, and
`kappa` is the core-mode fraction. `W_RXN`, `W_CORE`, and `IMAG_PEN` are
exposed as hypothesis knobs with defaults preserved.

## Repository Layout

```text
src/rxn_core/
  alignment/              molecule-level WBO alignment and sweep cut
  growth/                 connected-fragment growth and trace events
  matcher/                symmetry-compressed candidate matching
  chemistry_computations/ xyz, xtb cache, and frame helpers
  frag.py                 WBO graph and bond-change classification
  modes.py                Hessian parsing and mode-score features
  pipeline.py             BGCP full-view pipeline and CLI entry point

tests/                    focused unit tests
docs/TUTORIAL.ipynb       stage-by-stage usage tutorial notebook
docs/ARCHITECTURE.md      module boundary notes
ALGORITHM.md              algorithm details
```

## Runtime Inputs

These controls select files, outputs, cache-fill behavior, and parallelism.

| CLI | environment | default | meaning |
|---|---|---:|---|
| none | `RXN_CORE_PROJECT` | package root | base path used to resolve default data/output paths |
| none | `BGCP_WORK` | `data/xtb_frequency_calculations` | xtb cache root containing step directories |
| none | `BGCP_OUT_ROOT` | `out/bgcp_views` | per-step HTML/eval output root |
| `--stage-root` | `BGCP_STAGE_ROOT` | `out/bgcp_stages` | resumable `rp_stage.json` and `ts_stage.json` artifact root |
| `--alignment-out-root` | `BGCP_ALIGNMENT_OUT_ROOT` | `out/bgcp_alignments` | output root for optional clean Stage 1 aligned-coordinate exports |
| none | `BGCP_EVAL_JSON` | `out/bgcp_alignment_eval.json` | merged JSON summary output |
| `--stage` | `BGCP_STAGE` | `full` | run `rp`, `ts`, `view`, or composed `full` stage |
| `--mechanism` | none | all | restrict Stage 2 verification to one mechanism id; repeat for multiple ids |
| `--save-alignment-files` | `BGCP_SAVE_ALIGNMENT_FILES` | `0` | write mechanism-specific aligned R/P files during Stage 1 or full runs |
| `--steps` | none | all steps | explicit cached step names to process |
| `--limit` | none | none | process first N step directories after sorting |
| `--include-gt` | `BGCP_INCLUDE_GT` | `0` | load and score optional GT cache directories |
| `--xtb-mode` | `BGCP_XTB_MODE` | `auto` | `auto` fills missing xtb caches; `cache-only` never runs xtb |
| `--xtb-omp-threads` | `BGCP_XTB_OMP_THREADS` | `auto` | requested OMP threads for each xtb molecule |
| `--xtb-max-threads` | `BGCP_XTB_MAX_THREADS` | `8` | hard cap on OMP threads for each xtb molecule |
| `--charge` | `BGCP_CHARGE` | `0` | molecular charge used only when auto-filling missing xtb caches |
| `--multiplicity` | `BGCP_MULTIPLICITY` | `1` | spin multiplicity for auto xtb cache-fill; converted to `--uhf=multiplicity-1` |
| `--workers` | none | `os.cpu_count()-1` | total CPU budget in auto mode, or outer workers in outer mode |
| `--parallel-mode` | `BGCP_PARALLEL_MODE` | `auto` | `auto`, `outer`, or `inner` scheduling |
| `--inner-workers` | none | `0` | explicit inner workers per step; `0` lets the mode choose |
| `--auto-inner-workers` | `BGCP_AUTO_INNER_WORKERS` | `8` | target inner workers per concurrent step in auto mode |
| none | `BGCP_CUTSWEEP_CHUNKSIZE` | `1` | multiprocessing chunk size for R-P cut-sweep work units |
| none | `BGCP_TIMING` | `0` | print per-target timing diagnostics |

## Hypothesis Defaults

These parameters define the algorithmic hypotheses used for mechanism
discovery and TS/IG scoring. Defaults are conservative values used by the
current benchmark workflow; expose them when testing sensitivity.

| CLI | environment/API name | default | rationale |
|---|---|---:|---|
| none | `graph_floor` | `0.2` | Active WBO graph edge floor used by fragment growth and TS core-edge preservation; it is fixed in the BGCP pipeline and available in lower-level APIs. |
| none | `BGCP_CUT_FLOOR` | `0.2` | R-P mechanism discovery sweeps cuts over every R edge at or above this WBO; `0.2` includes weak but chemically relevant WBO graph edges while excluding near-zero pairs. |
| `--iso-tol` | `BGCP_ISO_TOL` / `iso_tol` | `1.0` | Active R-side growth edges must have a P-side active edge with `abs(WBO_R-WBO_P) <= iso_tol`; the loose default tolerates endpoint/TS distortion while bond-change ranking decides the mechanism. |
| `--dwbo-threshold` | `BGCP_DWBO_THRESHOLD` / `dwbo_threshold` | `0.5` | Broken/forming events require `abs(delta WBO) >= 0.5`; smaller WBO differences are treated as spectator variation. |
| `--symmetry-wbo-tol` | `BGCP_SYMMETRY_WBO_TOL` / `symmetry_wbo_tol` | `0.2` | Nauty orbit detection buckets WBO values within this tolerance; this collapses xtb/noise-level symmetry without changing exact active-edge validity checks. |
| none | `BGCP_VIEW_MAX_BRANCHES` / `max_branches` | `100` | Per-cut branch cap for R-P sweep; cuts that exceed it are discarded as pathological branch multipliers. |
| none | `BGCP_TS_CORE_EDGE_FLOOR` | `0.2` | Minimum target WBO for preserving a core edge during TS/IG core matching; mirrors the active graph floor. |
| none | `BGCP_TS_CORE_MAX_CANDIDATES` | `20000` | Cap on mechanism-local TS/IG core mappings to prevent runaway core enumeration. |
| none | `BGCP_SYMMETRY_REPAIR` | `1` | Enables local reshuffling inside product symmetry orbits after R-P matching to remove witness-choice artifacts. |
| none | `BGCP_SYMMETRY_REPAIR_MIN_CHANGES` | `5` | Only run symmetry repair when the initial witness has at least this many changed bonds; avoids unnecessary local search on already-clean mappings. |
| none | `BGCP_SYMMETRY_REPAIR_MAX_EVALS` | `20000` | Evaluation cap for the local symmetry-repair search. |
| `--w-rxn` | `BGCP_W_RXN` | `1.0` | Weight on reaction-coordinate overlap `rho` in the final TS/IG score. |
| `--w-core` | `BGCP_W_CORE` | `0.2` | Weight on core-mode fraction `kappa`; lower than `W_RXN` so localized core motion helps without dominating bond-axis overlap. |
| `--imag-pen` | `BGCP_IMAG_PEN` | `0.3` | Penalty exponent for multiple imaginary modes; soft penalty because IG Hessians may not be optimized TSs. |
| none | `N_SEEDS_PER_RUN` | `3` | Fixed seed-order count per cut-sweep run; cut diversity plus three seeds has been enough for the benchmark while keeping runtime bounded. |

CLI options:

```text
--workers              total CPU budget in auto mode
--parallel-mode        auto | outer | inner
--inner-workers        explicit per-step inner worker count
--auto-inner-workers   target inner workers per concurrent step in auto mode
--iso-tol              WBO tolerance for active graph matching
--dwbo-threshold       WBO delta for broken/formed bond classification
--symmetry-wbo-tol     WBO tolerance for symmetry-orbit bucketing
--w-rxn                reaction-coordinate score weight
--w-core               core-mode score weight
--imag-pen             imaginary-mode count penalty exponent
--xtb-mode             auto | cache-only
--xtb-omp-threads      requested OMP_NUM_THREADS per xtb molecule
--xtb-max-threads      hard cap on OMP_NUM_THREADS per xtb molecule
--charge               molecular charge for auto xtb cache-fill
--multiplicity         spin multiplicity for auto xtb cache-fill
--include-gt           load and score optional GT cache directories
--steps                explicit cached step names
--limit                first N cached steps
```

In `auto` or `inner` mode, the per-step inner worker pool is used for both
expensive phases: R-P cut sweep and independent TS/IG endpoint core matching.
If missing caches trigger xtb, each individual xtb subprocess gets
`OMP_NUM_THREADS=min(requested, BGCP_XTB_MAX_THREADS)`, with the default cap at
8 per molecule.

## Public API

Use top-level imports for stable pieces:

```python
from rxn_core import align_from_arrays, cut_sweep, ts_core_pool
from rxn_core import build_graph, classify_bonds
from rxn_core import parse_g98_modes, bond_overlap_per_mode
```

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for module boundaries and
[ALGORITHM.md](ALGORITHM.md) for the matching algorithm.
