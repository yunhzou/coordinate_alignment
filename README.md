# rxn_core

Symmetry-aware WBO atom mapping, analytical R/P alignment, and
mechanism-local transition-state analysis.

## Design

The computational API consists of immutable typed stages:

```text
search_aam
    -> compile_mapping_families
        -> select_rp_mappings
            -> analyze_transition_state
```

AAM is the authoritative source of mapping information. Its result retains
mechanism classes, every unique completed branch, fragment hierarchy,
assignment domains, exact fixed roles, target automorphism generators, cuts,
provenance, and search metrics. R/P and TS processing consume those objects;
they do not reconstruct an alternative AAM model from serialized records.

See [ALGORITHM.md](ALGORITHM.md) for the complete data model and algorithms.

## Install

```bash
python -m pip install -e .
```

Install xTB separately only when endpoint WBO matrices must be computed. The
typed core API accepts already materialized coordinates and WBO matrices and
does not invoke xTB.

## Python API

```python
from rxn_core import (
    AAMProblem,
    AAMSearchConfig,
    MolecularEndpoint,
    TransitionStateTarget,
    VibrationalModes,
    analyze_transition_state,
    compile_mapping_families,
    search_aam,
    select_rp_mappings,
)

reactant = MolecularEndpoint(elements_R, xyz_R, wbo_R, label="R")
product = MolecularEndpoint(elements_P, xyz_P, wbo_P, label="P")
problem = AAMProblem(reactant, product, name="reaction")

aam = search_aam(problem, AAMSearchConfig(), workers=8)
families = compile_mapping_families(
    aam, workers=8, minimum_events_only=True)
rp = select_rp_mappings(families)

target = TransitionStateTarget(
    MolecularEndpoint(elements_TS, xyz_TS, wbo_TS, label="TS"),
    VibrationalModes(frequencies, normal_modes),
)
ts = analyze_transition_state(rp, target)
```

For R/P only, `align_reaction(problem, workers=8)` is the convenience
composition of the first three stages. The individual stages remain available
when callers need to inspect, cache, audit, or transform AAM information.

Serialization, CLI workflows, and self-contained HTML views are typed
artifact adapters, not computational data models.

## CLI

NPZ endpoint files contain `elements`, `coordinates`, and `wbo` arrays:

```bash
rxn-core \
  --stage rp \
  --reactant-npz R.npz \
  --product-npz P.npz \
  --workers 48 \
  --post-workers 8 \
  --output alignment
```

Existing xTB cache directories containing one XYZ and a `wbo` file can be
used directly:

```bash
rxn-core --reactant-cache cache/R --product-cache cache/P \
  --workers 48 --output alignment
```

The output contains `rp.json`, one `R.xyz` and `P_aligned.xyz` pair per
mechanism, a reusable `reaction.json`, and a self-contained `view.html`.
The TS verifier/scorer can then be entered independently without R/P search:

```bash
rxn-core --stage ts \
  --reactant-npz R.npz --product-npz P.npz \
  --reaction-json alignment/reaction.json \
  --target-npz guess_1.npz --target-npz guess_2.npz \
  --output ts_scores
```

Use `--stage full` to compose both stages in one process. Target NPZ files
contain `elements`, `coordinates`, `wbo`, `frequencies`, and `modes`.

The executed [tutorial notebook](docs/TUTORIAL.ipynb) demonstrates the staged
Python API, CPU controls, AAM inspection, R/P alignment, and TS scoring.

## Result hierarchy

```text
AAMResult
`- AAMMechanism[]
   `- AAMBranch[]
      |- AtomBijection representative
      |- AAMHierarchy / FragmentMatch[]
      |- exact symmetry domains and generators
      `- provenance and branch counts

AnalyticalAAMResult
`- maximal exact mapping families per mechanism

RPResult
`- chirality-valid, minimum fixed-mapping-RMSD mapping per mechanism

TSResult
`- R->TS and P->TS CoreAAMResult plus exact scored core tuples
```

## Tests

```bash
.venv/bin/pytest -q
```

The suite includes a non-empty TS integration case that performs endpoint
AAM, analytical-family compilation, R/P selection, two partial core searches,
endpoint-consensus merging, and imaginary-mode scoring.
