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
SMILES/CXSMILES formal-WBO mode uses RDKit; install it with
`python -m pip install -e '.[smiles]'` or from conda-forge in the same
environment. The local webapp can persist sessions in MongoDB/GridFS; install
that optional dependency with `python -m pip install -e '.[sessions]'`.
The optional internal-coordinate path/viewer tooling uses ASE's IDPP
implementation; install it with `python -m pip install -e '.[path]'`.

Verify the installed import and command line:

```bash
python -c "import rxn_core; print(rxn_core.__file__)"
rxn-core --help
rxn-core-pipeline --help
rxn-core-webapp --help
```

## Usage

For a stage-by-stage walkthrough, open
[`docs/TUTORIAL.ipynb`](docs/TUTORIAL.ipynb).
For candidate electronic descriptors to retain for mechanism analysis, see
[`docs/REACTIVITY_DESCRIPTORS.md`](docs/REACTIVITY_DESCRIPTORS.md).

```bash
# Direct XYZ mode: no benchmark step schema required
rxn-core --stage rp --name my_reaction \
  --reactant-xyz R.xyz --product-xyz P.xyz \
  --workdir work/my_reaction \
  --charge 0 --multiplicity 1 --xtb-mode auto

# Direct R-P AAM with hard atom anchors. Here R atom 13 must map to P atom 9.
rxn-core --stage rp --name my_reaction_anchor_13_9 \
  --reactant-xyz R.xyz --product-xyz P.xyz \
  --workdir work/my_reaction_anchor_13_9 \
  --anchor 13:9

# Direct SMILES/CXSMILES mode: formal bond orders, no xtb.
rxn-core --stage rp --name acid_deprotonation \
  --reactant-smiles '[O:1][H:2]' \
  --product-smiles '[O-:1].[H+:2]'

# Local interactive AAM/subgraph workbench
rxn-core-webapp --port 8765 --open

# Optional local MongoDB for saved webapp sessions
docker compose -f docker-compose.mongo.yml up -d
rxn-core-webapp --port 8765 \
  --mongo-uri mongodb://localhost:27018 \
  --mongo-db rxn_core --open

# Post-Stage-1 collective validation after Stage 1 has written rp_stage.json
rxn-core --stage post-rp --name my_reaction \
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

### SMILES / CXSMILES Formal-WBO Mode

For graph-only examples, `rxn-core` can build the R/P endpoint graphs directly
from SMILES or CXSMILES. This mode does not run xtb. It parses the written
molecular graph with RDKit, keeps explicitly written hydrogens, and sets the
WBO matrix to formal bond orders: single `1.0`, double `2.0`, triple `3.0`,
aromatic `1.5`. By default RDKit atom hydrogen counts are materialized as
separate H nodes in the AAM graph, including hydrogens implied by bracket atoms
such as `[CH]` and `[CH2]`. Pass `--smiles-preserve-explicit-only` only when
you need the older graph containing just atoms explicitly present in the parsed
SMILES.

Generated coordinates are planar RDKit depictions for the viewer only.
CXSMILES atom-map labels are written as source metadata, but they are not hard
AAM constraints unless you also pass `--anchor R:P`.

Programmatic use:

```python
from rxn_core import smiles_to_formal_wbo, smiles_inputs_from_strings
from rxn_core import run_rp_stage

endpoint = smiles_to_formal_wbo("[CH3:1][O:2][H:3]")
inputs = smiles_inputs_from_strings(
    "[O:1][H:2]",
    "[O-:1].[H+:2]",
    name="acid_deprotonation",
)
rp = run_rp_stage(inputs)
```

### Subgraph Matching And Anchors

The matcher can also be used as a standalone weighted-subgraph search. This is
the same symmetry-compressed fragment-growth engine used by AAM, but exposed
with replaceable node compatibility rules. The default node rule is same
element; the edge rule still uses WBO active-edge matching with `--iso-tol`.

Weighted graph JSON inputs have this shape:

```json
{
  "nodes": [
    {"element": "C", "features": {"outer_shell": 4}},
    {"element": "O", "features": {"outer_shell": 6}}
  ],
  "weights": [[0.0, 1.1], [1.1, 0.0]],
  "weight_name": "wbo",
  "coords": [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]
}
```

Run standalone subgraph matching from the CLI:

```bash
rxn-core --subgraph-query-json query.json \
  --subgraph-target-json target.json \
  --subgraph-node-policy outer_shell \
  --subgraph-anchor 0:12 \
  --subgraph-output subgraph_matches.json
```

`--subgraph-node-policy` may be repeated to build a multi-field node key.
Fields can live either as top-level node attributes or inside `features`.
`--subgraph-anchor q:t` hard-locks query node `q` to target node `t`; repeated
anchors must be one-to-one. Anchored atoms are preloaded as locked single-atom
islands, but they can still seed growth so their local environment is checked.

The same anchor syntax is available for direct R-P AAM:

```bash
rxn-core --stage rp --name anchored_rp \
  --reactant-xyz R.xyz --product-xyz P.xyz \
  --anchor 13:9
```

For direct R-P mode, `--anchor R:P` is enforced through every cut-sweep island
search. If an anchor is incompatible with the available chemistry graph, the
anchored run returns no matching mechanism instead of silently choosing a
different atom. `--subgraph-anchor` is accepted as an alias in direct R-P mode
for anchor maps exported by the picker.

For an interactive local picker, open:

```bash
open tools/anchor_picker.html
```

The picker renders R and P with the same 3Dmol stick/sphere style as the main
viewer. Click one R atom and one P atom to add an anchor pair, then copy the
JSON `anchor_map` or CLI flags.

For a browser UI that runs the actual R-P matcher and subgraph matcher, launch:

```bash
rxn-core-webapp --port 8765 --open
```

The workbench has two tabs. The R-P tab accepts SMILES/CXSMILES or XYZ
endpoints, previews atoms for anchor picking, runs R-P AAM with the selected
anchors, and embeds the same generated `view.html` used by the pipeline. XYZ
mode uses xTB for endpoint WBO generation. The subgraph tab accepts drawn
element/bond-order graphs, SMILES, XYZ, or WeightedGraph JSON for query and
target matching.

Saved webapp sessions use MongoDB document storage for the UI state and GridFS
for generated artifacts such as `view.html`, source files, and `rp_stage.json`.
If MongoDB is unavailable, the workbench still runs but the Save/Load controls
report the session-store error. A local Mongo instance can be started with:

```bash
docker compose -f docker-compose.mongo.yml up -d
```

The compose file binds the container to host port `27018` by default to avoid
conflicting with an existing local Mongo on `27017`. Override with
`RXN_CORE_MONGO_PORT=27017 docker compose -f docker-compose.mongo.yml up -d`
if you want the standard host port.

Programmatic use:

```python
from rxn_core import WeightedGraph, match_weighted_subgraph

query = WeightedGraph(
    nodes=[
        {"element": "C", "features": {"outer_shell": 4}},
        {"element": "O", "features": {"outer_shell": 6}},
    ],
    weights=[[0.0, 1.1], [1.1, 0.0]],
)
target = WeightedGraph(
    nodes=[
        {"element": "C", "features": {"outer_shell": 4}},
        {"element": "O", "features": {"outer_shell": 6}},
        {"element": "O", "features": {"outer_shell": 6}},
    ],
    weights=[
        [0.0, 1.1, 0.0],
        [1.1, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ],
)
matches = match_weighted_subgraph(
    query,
    target,
    node_policy="outer_shell",
    anchor_map={0: 0},
)
```

### Real Example

This direct-XYZ example uses `pr1.tempo_ts3` from the stored docs example. The
repo includes the source XYZ files and cache files under
`docs/example_runs/pr1.tempo_ts3/`, so it is self-contained.

```python
from pathlib import Path
import rxn_core.pipeline as rxnp

name = "pr1.tempo_ts3"
example_root = Path("docs/example_runs") / name
prepared_step = example_root / "prepared_steps" / name

reactant_xyz = prepared_step / "R/reactant_01_reactant_01_5.xyz"
product_xyz = prepared_step / "P/product_01_product_01_6.xyz"
workdir = example_root / "work"

target_specs = [
    {
        "kind": "ig",
        "label": "iter1",
        "xyz": prepared_step / "sp_iter1/pr1.tempo_ts3_benchmark_plain_iter1_87ea3b8f.xyz",
    },
    {
        "kind": "gt",
        "label": "GT",
        "xyz": prepared_step / "sp_groundtruth/ts_groundtruth_01_reference_ts_01_TS3.xyz",
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
    xtb_mode="cache-only",
    inner_workers=8,
    save_alignment_files=True,
)
```

If Stage 1 has already run and IGs were generated from its mechanisms, resume
from the saved `rp_stage.json` instead of rerunning R-P discovery:

```python
post_rp = rxnp.process_xyz_stage(
    name,
    reactant_xyz,
    product_xyz,
    workdir=workdir,
    stage="post-rp",
    target_specs=target_specs,
    charge=0,
    multiplicity=1,
    xtb_mode="cache-only",
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

# Stage 1 plus clean mechanism-specific aligned R/P files for path setup
rxn-core --stage rp --steps pr7.V.dodh_ts910 --save-alignment-files

# Post-Stage-1 collective validation plus selected best-S TS core-aligned files
rxn-core --stage post-rp --steps pr7.V.dodh_ts910 \
  --mechanism 2 --include-gt --save-alignment-files

# Post-Stage-1: verify GT/IG/TS targets from saved Stage 1 mechanisms
rxn-core --stage post-rp --steps pr7.V.dodh_ts910 --mechanism 2 --include-gt

# Stage 3: regenerate only the HTML/eval view from saved stage artifacts
rxn-core --stage view --steps pr7.V.dodh_ts910 --include-gt

# Full pipeline: compose rp + ts + view in one run
rxn-core --stage full --steps pr7.V.dodh_ts910 --workers 8

# Full validation/view/export, but reuse existing rp_stage.json
rxn-core --stage full --resume-rp --steps pr7.V.dodh_ts910 --workers 8
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
by default. Target cache filling uses a separate `BGCP_XTB_WORKERS` pool, so
multi-threaded xtb jobs do not share the alignment worker count directly. The
cache-fill input uses molecular `charge` and spin
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
out/bgcp_alignments/<step>/mechanisms/mechanism_<id>/path_endpoints.xyz
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
two-frame `path_endpoints.xyz`, `mapping_R_to_P.csv`, and `mechanism.json`.
`P_aligned.xyz` is product geometry reindexed into the R atom order; no
spatial fitting is applied.

### Native index-chirality-preserving AAM

`--index-chirality preserve` makes the final Stage 1 atom assignment preserve
the R-index orientation of every defined tetrahedral frame that the recorded
AAM symmetry choices can change. It is disabled by default; `off` is a true
no-op.

The selected mapping is still a full, element-preserving R-to-P bijection from
the assignment family already encoded by AAM. The method can materialize an
implicit odd/even choice from a recorded closed symmetry block, but it never
discovers a new automorphism, invents a new atom choice, or changes the
mechanism to repair geometry.

This runs natively while Stage 1 finalizes each mechanism, before
`rp_stage.json` is written. `_SymCand` growth, cut sweeping, graph matching,
and mechanism deduplication keep their existing logic; the new selector
consumes their serialized witnesses, alternates, and symmetry blocks.

#### Run one XYZ pair

Use real 3D endpoint geometries. This command computes missing WBO caches with
xtb, runs Stage 1 with the native constraint, and writes both the certificate
and R/P files:

```bash
BGCP_OUT_ROOT=out/native_chirality/views \
rxn-core --stage rp --name my_case \
  --reactant-xyz R.xyz --product-xyz P.xyz \
  --workdir work/my_case \
  --charge 0 --multiplicity 1 --xtb-mode auto \
  --index-chirality preserve --workers 8 \
  --stage-root out/native_chirality/stages \
  --alignment-out-root out/native_chirality/alignments \
  --save-alignment-files
```

The important outputs are:

```text
out/native_chirality/stages/my_case/rp_stage.json
out/native_chirality/views/my_case/view.html
out/native_chirality/alignments/my_case/manifest.json
out/native_chirality/alignments/my_case/mechanisms/mechanism_<NNN>/R.xyz
out/native_chirality/alignments/my_case/mechanisms/mechanism_<NNN>/P_aligned.xyz
out/native_chirality/alignments/my_case/mechanisms/mechanism_<NNN>/path_endpoints.xyz
out/native_chirality/alignments/my_case/mechanisms/mechanism_<NNN>/mapping_R_to_P.csv
```

`P_aligned.xyz` is P reindexed into the selected R atom order. Stage 1 does
not translate, rotate, Kabsch-fit, reflect, or otherwise modify either
endpoint geometry.

For a fully cached run, place the WBO files at:

```text
work/my_case/endpoints/R/wbo
work/my_case/endpoints/P/wbo
```

Their atom order must match the corresponding input XYZ files. Then use the
same command with `--xtb-mode cache-only`. R and P must have the same atom
count and elemental composition.

Standard `--steps` mode uses a different historical layout:

```text
<BGCP_WORK>/<case>/R/<endpoint>.xyz
<BGCP_WORK>/<case>/R/wbo
<BGCP_WORK>/<case>/P/<endpoint>.xyz
<BGCP_WORK>/<case>/P/wbo
```

Run multiple cases from that layout with:

```bash
BGCP_WORK=/path/to/work \
BGCP_OUT_ROOT=out/native_chirality/views \
rxn-core --stage rp --steps TS_01 TS_04 \
  --index-chirality preserve \
  --workers 40 --parallel-mode auto \
  --stage-root out/native_chirality/stages \
  --alignment-out-root out/native_chirality/alignments \
  --save-alignment-files --xtb-mode cache-only
```

`--workers` is the total CPU budget. For one direct XYZ case it is used for
the inner cut-sweep work; for multiple standard cases, `--parallel-mode auto`
divides it between case-level and inner work.

Do not use `--stage view`, `--stage post-rp`, or `--resume-rp` to convert an
old `rp_stage.json` made with `index_chirality=off`: those modes reuse the
stored mapping. Run Stage 1 again with `preserve`. Planar coordinates generated
from SMILES are useful for graph examples but generally cannot certify 3D
tetrahedral endpoint orientation.

#### Intended algorithm in plain language

1. Find each R atom with exactly four active graph neighbors whose four mapped
   neighbors remain attached to the mapped P center.
2. Record the orientation sign of those four neighbor labels in R and P.
3. Read possible assignments only from the existing AAM object: the source
   mapping, complete branch witnesses, correlated nested alternates, and
   eligible closed symmetry blocks.
4. Freeze anchors, `exact_fixed` atoms, and the exact broken/formed-bond
   signature. Solve the remaining closed-block odd/even choices together as
   XOR equations over `GF(2)`.
5. Keep only authorized mappings with zero mismatches on every defined,
   switchable frame. Select the one with the fewest changes from the source
   witness, then use canonical mapping order as the deterministic tie-break.

This is a bounded parity finalization on top of AAM metadata. It is not a
greedy swap repair and does not enumerate all permutations of a symmetry
group. The base AAM search can still be the expensive part, but this
finalization has no factorial candidate expansion.

#### Frame definition

For the four R-index-ordered neighbors `n0 < n1 < n2 < n3`, the frame sign is
the affine orientation

```text
det(X[n1] - X[n0], X[n2] - X[n0], X[n3] - X[n0])
```

of the four neighbor points. The center atom is deliberately not the
determinant origin. An odd relabeling of the four neighbors therefore reverses
the sign, while an even relabeling preserves it.

This is R-index label parity, not a CIP `R`/`S` assignment. Three-coordinate
centers and centers with more than four R graph neighbors do not define this
four-label constraint.

The determinant is evaluated in long-double arithmetic. A frame is undefined
only when it is exactly degenerate or its sign cannot be distinguished from
floating-point roundoff using a deterministic, scale-aware error bound derived
from machine epsilon. A planar four-neighbor set is one example. Undefined
frames are neutral diagnostics.

There is no volume tolerance and no chirality tuning hyperparameter. The
existing active-graph threshold still determines graph adjacency; the
chirality method introduces no geometric cutoff.

#### Assignment authority

Complete branch witnesses and nested alternates remain atomic. An alternate
is patched onto its complete owner as one correlated mapping and is never
split into independent per-atom image choices.

One implicit parity route may use only untagged, closed, equal-size, pairwise
disjoint `_SymBlock`s that coexist in one nested fragment symmetry state. Its
base is that fragment's complete owner witness or one atomic alternate from
the same fragment. Blocks from different historical fragments are never
composed. The base must occupy each block's exact P pool; anchors and
`exact_fixed` R atoms remain frozen.

A block contributes one parity bit to a frame only when the whole block lies
inside that frame's four-neighbor shell. A block containing the center or only
partially intersecting the shell is not treated as an independent parity
choice. Bit `0` keeps the base representative; bit `1` applies a deterministic
canonical odd transposition. Gaussian elimination solves all frame equations
together over `GF(2)`, with unconstrained free bits fixed to zero.

Every explicit or parity-materialized candidate must retain the source
mapping's exact R-index broken/formed event signature. Whether a frame touches
a reaction-event atom is recorded for diagnosis but does not disable the
constraint.

A defined frame is hard only when it contains an atom that the recorded AAM
choices can switch. A frame with no switchable atom cannot be repaired by
assignment selection, so a mismatch there is reported under
`immutable_frames`; it is not hidden and does not authorize a new mapping.

#### Verify the Stage 1 certificate

The authoritative selected mapping is
`mechanisms[*].mapping_RP` in `rp_stage.json`. Its audit is under
`mechanisms[*].branch_symmetry.index_chirality`. For example:

```python
import json
from pathlib import Path

from rxn_core.alignment.index_chirality import mapping_sha256

stage = json.loads(Path(
    "out/native_chirality/stages/my_case/rp_stage.json"
).read_text())
for mechanism in stage["mechanisms"]:
    branch = mechanism["branch_symmetry"]
    audit = branch["index_chirality"]
    assert branch["active_assignment_family"] == (
        "index_chirality.native_symcand/v3"
    )
    assert audit["schema_version"] == "rxn_core.index_chirality/v3"
    assert audit["policy"] == "preserve"
    assert audit["status"] == "applied"
    assert audit["selected_index_chirality_violation_count"] == 0
    assert audit["selected_mapping_sha256"] == mapping_sha256(
        mechanism["mapping_RP"]
    )
    print(
        mechanism["id"],
        audit["source_index_chirality_violation_count"],
        audit["selected_index_chirality_violation_count"],
        audit["immutable_source_mismatch_count"],
        audit["undefined_frame_count"],
    )
```

A successful v3 selector application has:

```text
active_assignment_family = index_chirality.native_symcand/v3
schema_version = rxn_core.index_chirality/v3
policy = preserve
status = applied
selected_index_chirality_violation_count = 0
```

That is the hard guarantee for every defined frame that AAM can switch. If
your acceptance rule is that every defined frame must match, also require
`immutable_source_mismatch_count = 0`. If it additionally requires every
eligible four-neighbor frame to have a defined sign, require
`undefined_frame_count = 0` as well.

`mapping_changes` records each R index changed from the initial AAM witness.
`frames` contains the hard switchable-frame certificate. `candidate_search`
records the atomic seeds, fragment-local parity variables, `GF(2)` equations,
solutions, and evaluated candidates. `invariants` records the authority and
mechanism-preservation checks.

Always read `immutable_source_mismatch_count` separately: a nonzero value means AAM
encoded no legal choice capable of changing those reported frames. Likewise,
`undefined_frame_count` records degenerate or numerically indeterminate frames
whose sign was not constrained. Therefore a selected violation count of zero
certifies every defined AAM-switchable frame, not arbitrary stereochemistry
outside the encoded assignment family.

If no authorized mapping satisfies all hard frames, core finalization raises
`IndexChiralityConflict`. The direct CLI prints `<name>: ERROR: ...`, records
the error in its evaluation JSON, and does not write a new successful Stage 1
artifact. Check the printed/evaluation result rather than relying only on the
shell exit status.

#### Reproducible cached-case runner

The repository helper used for selected benchmark reruns expects:

```text
<work-root>/<case>/endpoints/R/reactant_combined.xyz
<work-root>/<case>/endpoints/R/wbo
<work-root>/<case>/endpoints/P/product_combined.xyz
<work-root>/<case>/endpoints/P/wbo
```

If a preferred combined XYZ name is absent, its endpoint directory must
contain exactly one XYZ file. Run one case with:

```bash
python tools/run_native_index_chirality_case.py \
  --case TS_01 \
  --work-root /path/to/work \
  --output-root /path/to/native_v3 \
  --workers 8 \
  --policy preserve
```

It atomically writes:

```text
/path/to/native_v3/runs/TS_01/rp_stage.json
/path/to/native_v3/runs/TS_01/run_summary.json
/path/to/native_v3/runs/TS_01/alignment/
/path/to/native_v3/runs/TS_01/views/TS_01/view.html
```

The helper propagates a conflict without publishing a partial case and refuses
to overwrite an existing run directory.

### Optional endpoint fitting, interpolation audit, and viewer

These helpers live under `tools/`, outside the installable `rxn_core` package.
They are packaging and path-diagnostic tools, not part of AAM and not an NEB
optimizer. For a native v3 result, the data flow is:

```text
native Stage 1 mechanism.mapping_RP
  -> freeze that exact mapping
  -> reindex P into R row order
  -> fit on fixed mapped atoms and apply that one proper transform to all P
  -> generate a signed IDPP internal-coordinate band
  -> write R.xyz, P_final.xyz, the band, and an offline viewer
```

The proper global fit has `det(rotation) = +1` and changes coordinates only.
It cannot change the atom mapping, reflect P, or rotate fragments
independently. For native v3 records, the downstream candidate family is a
singleton and the recorded affine frames are reused exactly. Legacy non-native
records retain a separate compatibility path that can read correlated choices
from the selected AAM witness.

The fit requires at least three non-collinear fixed mapped atoms. It prefers
non-core heavy fixed atoms, then any non-core fixed atoms, then all fixed
atoms; hydrogen anchors receive weight `0.1`. The resulting one proper
rotation and translation are applied to every P atom.

To reproduce the complete multi-case package, first run each case with the
cached-case helper above. Then build an indexed source archive. The viewer
template must be the historical all-case coordinate-alignment archive, not a
per-case view produced by the helper. Its root must contain `index.json` and
the original `viewer.html` with embedded `const DATA`, the requested cases,
and the standard remote 3Dmol script tag. The final packager replaces that tag
with the local JavaScript supplied below:

```bash
python tools/build_native_chirality_source_archive.py \
  --run-root /path/to/native_v3 \
  --viewer-template-root /path/to/original_viewer_archive \
  --output-root /path/to/native_v3/source_archive \
  --case TS_01 \
  --case TS_04
```

Freeze each native mapping and apply the proper global fit:

```bash
python tools/run_neb_orientation_batch.py \
  --source-root /path/to/native_v3/source_archive \
  --wbo-work-root /path/to/work \
  --output-root /path/to/native_v3/orientation
```

Each intermediate mechanism directory contains `R.xyz`,
`P_neb_ordered.xyz`, `neb_endpoints.xyz`, `mapping_R_to_P.csv`, and
`neb_orientation.json`. The `neb` names are historical; this command does not
run an NEB calculation.

Finally build the portable package:

```bash
python tools/build_neb_orientation_deliverable.py \
  --source-root /path/to/native_v3/source_archive \
  --orientation-root /path/to/native_v3/orientation \
  --output-root /path/to/native_v3/deliverable \
  --3dmol-js /path/to/3Dmol-min.js \
  --3dmol-license /path/to/3Dmol-LICENSE \
  --archive /path/to/native_v3/deliverable.tar.gz \
  --interpolation-images 21
```

Use fresh destinations for every command in this packaging chain. The source
archive builder, orientation batch, deliverable builder, and archive writer
refuse to overwrite existing outputs.

The package has one self-contained `viewer.html` for all cases and, for every
case/mechanism:

```text
cases/<case>/mechanisms/mechanism_<NNN>/R.xyz
cases/<case>/mechanisms/mechanism_<NNN>/P_final.xyz
cases/<case>/mechanisms/mechanism_<NNN>/metadata.json
cases/<case>/mechanisms/mechanism_<NNN>/internal_coordinate_interpolation.xyz
cases/<case>/mechanisms/mechanism_<NNN>/interpolation_report.json
```

`R.xyz` and `P_final.xyz` have identical element rows in R-index order.
`P_final.xyz` uses the exact native `mapping_RP` plus one proper global rigid
fit.

The animation is not a Cartesian endpoint interpolation. Its internal
coordinates are all mapped atom-pair distances, with targets interpolated
between R and `P_final` and reconstructed as one IDPP band. Pair distances
alone cannot distinguish mirror images, so each native hard tetrahedron is
placed deterministically on its recorded R-sign branch by choosing the local
mirror embedding that best retains the IDPP distance targets. This changes
neither the mapping nor the six pair distances inside that tetrahedron.

The viewer plays the exact generated images embedded in the HTML; it does not
recompute `R + t(P_final - R)` in JavaScript. The report applies the exact
signed-volume cubic test to every adjacent image segment and analytically
checks every atom pair's closest approach on that piecewise path. The separate
close-approach thresholds remain a geometry-review heuristic, not part of the
native mapping constraint.

This is still an initializer/animation, not an optimized NEB trajectory or a
minimum-energy path. A failed path audit can coexist with a valid endpoint
mapping. The diagnostic package is still written so the failing images and
report can be inspected; a failure means the path should not be used as an NEB
initializer until it is improved. It must not widen the AAM mapping family.

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
4. `run_ts_stage(...)` runs no-cut endpoint-to-TS realignment from both R and
   P using the same symmetry-aware fragment-growth matcher as R-P alignment.
   P-side core mappings are pulled back through the R-P mechanism witness,
   then merged with R-side mappings in R-core indexing. Compressed branch
   degeneracy is expanded only for mechanism core atoms before scoring.
5. Each GT/IG candidate mapping is scored on the selected imaginary mode:

```text
S = beta * wbo_progress^WBO_PROGRESS_POWER
```

`beta` is the normal-mode overlap with the broken/formed bond-axis vector,
where each event is weighted by its detected R-P `abs(delta WBO)`.
`wbo_progress` is the event-weighted TS WBO progress in the same direction:
forming events require `WBO_TS > WBO_R`, while broken events require
`WBO_TS < WBO_R`. Multiple imaginary modes are not penalized because IG
Hessians are not optimized transition states.

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
tools/anchor_picker.html  local R/P anchor-picking viewer
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
| `--anchor` | none | none | hard R:P anchor for direct R-P AAM; can be repeated |
| `--reactant-smiles` / `--product-smiles` | none | none | direct R-P formal-bond-order input from SMILES/CXSMILES; skips xtb |
| `--smiles-expand-hydrogens` | none | on | materialize SMILES atom hydrogen counts as explicit H atoms before formal-WBO graph construction |
| `--smiles-preserve-explicit-only` | none | off | keep only atoms explicitly present in the parsed SMILES graph; atom hydrogen counts remain implicit |
| `--subgraph-query-json` | none | none | standalone weighted-subgraph query graph JSON |
| `--subgraph-target-json` | none | none | standalone weighted-subgraph target graph JSON |
| `--subgraph-node-policy` | none | same element | node attribute/feature field used for standalone subgraph compatibility; repeat for multi-field keys |
| `--subgraph-anchor` | none | none | hard query:target anchor for standalone subgraph matching; also accepted as a direct R-P anchor alias |
| `--subgraph-output` | none | `BGCP_EVAL_JSON` | output JSON path for standalone subgraph matching |
| `--save-alignment-files` | `BGCP_SAVE_ALIGNMENT_FILES` | `0` | write mechanism-specific aligned R/P files during Stage 1 or full runs |
| `--steps` | none | all steps | explicit cached step names to process |
| `--limit` | none | none | process first N step directories after sorting |
| `--include-gt` | `BGCP_INCLUDE_GT` | `0` | load and score optional GT cache directories |
| `--xtb-mode` | `BGCP_XTB_MODE` | `auto` | `auto` fills missing xtb caches; `cache-only` never runs xtb |
| `--xtb-omp-threads` | `BGCP_XTB_OMP_THREADS` | `auto` | requested OMP threads for each xtb molecule |
| `--xtb-max-threads` | `BGCP_XTB_MAX_THREADS` | `8` | hard cap on OMP threads for each xtb molecule |
| `--xtb-workers` | `BGCP_XTB_WORKERS` | `auto` | concurrent xtb target-cache jobs; `auto` uses available CPUs divided by xtb threads, capped by inner workers |
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
| `--index-chirality` | `BGCP_INDEX_CHIRALITY` / `index_chirality` | `off` | `preserve` enforces exact affine tetrahedral index parity over atomic AAM witnesses/alternates and coherent fragment-local closed-block parity equations; `off` leaves assignment selection unchanged. |
| none | `BGCP_VIEW_MAX_BRANCHES` / `max_branches` | `100` | Per-cut branch cap for R-P sweep; cuts that exceed it are discarded as pathological branch multipliers. |
| none | `BGCP_TS_ALIGN_GRAPH_FLOOR` / `graph_floor` | `0.2` | Active WBO graph edge floor for no-cut R/P-to-TS fragment growth. |
| none | `BGCP_TS_ALIGN_MAX_CORE_MAPS` / `max_core_maps` | `20000` | Cap on mechanism-local TS/IG core maps after expanding compressed core degeneracy. |
| none | `BGCP_PREFER_ENDPOINT_CONSENSUS` / `prefer_endpoint_consensus` | `1` | Prefer the highest-S exact core map recovered from both R-side and P-side endpoint matching; if no consensus map exists, fall back to the highest-S endpoint-union map. |
| `--event-weight-power` | `BGCP_EVENT_WEIGHT_POWER` | `1.0` | Exponent on each detected event's R-P `abs(delta WBO)` when building the weighted bond-motion vector for `beta`. |
| `--wbo-progress-power` | `BGCP_WBO_PROGRESS_POWER` | `1.0` | Exponent on the TS WBO progress factor in the final TS/IG score. |
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
--index-chirality      off | preserve
--event-weight-power   exponent on R-P event delta-WBO weights
--wbo-progress-power   exponent on TS WBO progress factor
--xtb-mode             auto | cache-only
--xtb-omp-threads      requested OMP_NUM_THREADS per xtb molecule
--xtb-max-threads      hard cap on OMP_NUM_THREADS per xtb molecule
--charge               molecular charge for auto xtb cache-fill
--multiplicity         spin multiplicity for auto xtb cache-fill
--include-gt           load and score optional GT cache directories
--anchor               hard R:P direct AAM anchor; repeatable
--reactant-smiles      direct R endpoint SMILES/CXSMILES formal-WBO input
--product-smiles       direct P endpoint SMILES/CXSMILES formal-WBO input
--smiles-expand-hydrogens
                       materialize SMILES atom hydrogen counts as explicit H atoms
                       (default)
--smiles-preserve-explicit-only
                       keep only atoms explicitly present in the parsed SMILES
                       graph
--subgraph-query-json  standalone weighted-subgraph query graph JSON
--subgraph-target-json standalone weighted-subgraph target graph JSON
--subgraph-node-policy node field for standalone subgraph compatibility
--subgraph-anchor      hard query:target subgraph anchor; direct AAM alias
--subgraph-output      standalone subgraph output JSON
--steps                explicit cached step names
--limit                first N cached steps
```

In `auto` or `inner` mode, the per-step inner worker pool is used for both
expensive phases: R-P cut sweep and independent TS/IG endpoint core matching.
If missing caches trigger xtb, each individual xtb subprocess gets
`OMP_NUM_THREADS=min(requested, BGCP_XTB_MAX_THREADS)`, with the default cap at
8 per molecule. Target cache filling uses `BGCP_XTB_WORKERS` instead of the
alignment worker pool directly, which lets runs use fewer concurrent xtb
processes with more OMP threads per process.

## Public API

Use top-level imports for stable pieces:

```python
from rxn_core import align_from_arrays, cut_sweep, run_cut_sweep_chunk
from rxn_core import build_graph, classify_bonds
from rxn_core import WeightedGraph, WeightedNode, match_weighted_subgraph
from rxn_core import parse_g98_modes, bond_overlap_per_mode
```

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for module boundaries and
[ALGORITHM.md](ALGORITHM.md) for the matching algorithm.
