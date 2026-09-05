# Native AAM growth engine

`rxn_core._engine` is an optional compiled extension (pybind11, C++17) that
executes fragment growth and its inner symmetry matcher. Python retains the
AAM public objects, search orchestration, and result hierarchy.

| Native operation | Python boundary | Current use |
|---|---|---|
| `_engine.grow_island` | `growth.island.grow_island` through `growth.native` | Active for supported calls: default element policy, exact orbit map, no trace events, and compatible graph data |
| `_native.paired_mapping_invariant` | `fragment_matching.detection._paired_mapping_invariant` | Active when both orbit maps provide structural zero buckets |
| `_engine.freeze_analytical` | Potential replacement for `alignment.sweep._freeze_analytical` | Exported by the extension; not called by the current Python pipeline |
| `_engine.repair_group` | Potential replacement for the per-group search in `alignment.branch.symmetry_repair_mapping` | Exported; not called by the current Python pipeline |
| `_engine.AutGraph` | Potential replacement for standalone `pynauty.autgrp` calls | Exported; current `matcher.canonical` still calls pynauty |

The C++ growth implementation ports the growth loop, candidate extension,
compressed candidate state, and candidate deduplication together. Their Python
implementations remain in `growth/` and `matcher/` as the reference engine.
Nauty is also used internally by the C++ growth engine.

## Python abstraction and native boundary

The public call remains `search_aam(AAMProblem, AAMSearchConfig) -> AAMResult`.
The following is structural pseudocode; scheduling and classification remain
Python operations:

```text
search_aam(problem, config)                         # Python
    schedule cut/seed searches                    # alignment.sweep / branch
    for each fragment-growth request:
        grow_island(graph_R, graph_P, seed, ...)   # Python entry point
            prepare/cache native graph views     # growth.native
            raw = _engine.grow_island(...)        # C++ growth + matcher
            restore original atom indices
            return list[_IsoResult]              # Python result objects
    merge branches and classify mechanisms
    attach exact fragment groups
    return AAMResult                              # Python domain object
```

`_IsoResult` contains the mapping, retained fragment, deferred edges, and
compressed symmetry state. Native candidate state lives in C++ while the
kernel runs; the Python bridge reconstructs these result fields when it
returns. It also converts a native branch-cap report into
`IslandBranchLimitExceeded`.

The completed Python output remains:

```text
AAMResult
    mechanisms: tuple[AAMMechanism]
        branches: tuple[AAMBranch]
            representative: AtomBijection
            hierarchy: AAMHierarchy
            mapping_family, target_group, cuts, path_provenance, ...
    metrics: AAMSearchMetrics
```

`AtomBijection` holds a representative assignment; the hierarchy and group
fields carry the family information. Using C++ does not remove those fields.
The public objects are defined in `src/rxn_core/domain.py` and
`src/rxn_core/alignment/post_aam.py`.

Fragment detection is a separate consumer of the same `grow_island` entry
point. Its augmentation and detection orchestration remain Python, as do
retrosynthesis assembly, ranking, and visualization. The separate `_native`
extension accelerates the mapping invariant used during detection.

## Build

```bash
.venv/bin/python native/build_engine.py
```

This compiles `native/src/*.cpp` and the vendored nauty sources in
`native/nauty/` and copies `_engine.cpython-*.so` into `src/rxn_core/`.
Requirements: a C/C++ compiler, `pybind11` and `setuptools` in the
environment. The `_engine` binary is not committed; when unavailable, growth
uses the Python implementation. The separate `_native` extension is built by
`setup.py` during package installation and is imported by fragment detection.

## Switches

- The native growth engine is used automatically when it is built.
- `RXN_CORE_NATIVE=0` selects Python fragment growth. It does not disable the
  separate `_native` mapping-invariant extension or pynauty.

## Catalog validation

The recorded explicit-hydrogen BIAN inventory benchmark searched all 1,919 precursors
with 28 shards and 48 CPUs per shard. The integrated native engine completed
in 73.17 seconds versus 78.32 seconds for the previous exact implementation.
All 1,919 serialized chemistry records were identical. These are historical
measurements, not timings for subsequent workflow changes.

## Verification

- `tests/test_native_engine.py`: every `grow_island` call of a full search on
  the tetraphenylmethane benchmark through both engines, plus equality of the
  typed results.
- `tests/test_native.py`: mapping-invariant kernel checks.

## Third-party code

`native/nauty/` contains the nauty 2.8.8 sources (Brendan McKay and Adolfo
Piperno, Apache License 2.0; see `native/nauty/LICENSE-2.0.txt` and
`native/nauty/COPYRIGHT`) as shipped inside the pynauty source distribution,
so the compiled engine uses the same nauty version as the pynauty package the
Python engine calls.
