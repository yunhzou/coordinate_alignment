# Native AAM growth engine

`rxn_core._engine` is an optional compiled extension (pybind11, C++17) that
replaces the Python growth engine (`grow_island`, `_extend_sym_cands`,
`_dedup_sym_cands`, the `_SymCand` state) and a few exact kernels used by
the sweep:

| Symbol | Replaces | Dispatch |
|---|---|---|
| `grow_island` | `rxn_core.growth.island.grow_island` | `rxn_core.growth.native` (default element policy, exact orbit map, no trace events) |
| `freeze_analytical` | `rxn_core.alignment.sweep._freeze_analytical` | bound at import |
| `repair_group` | per-group search in `symmetry_repair_mapping` | `rxn_core.alignment.branch` |
| `AutGraph` | `pynauty.autgrp` in `atom_generators` | `rxn_core.matcher.canonical` |

Every kernel is output-identical to the Python code it replaces; the Python
code stays in the tree as the reference and as the fallback.

## Build

```bash
.venv/bin/python native/build_engine.py
```

This compiles `native/src/*.cpp` and the vendored nauty sources in
`native/nauty/` and copies `_engine.cpython-*.so` into `src/rxn_core/`.
Requirements: a C/C++ compiler, `pybind11` and `setuptools` in the
environment.  The extension is not committed; without it everything runs in
Python.

## Switches

- `RXN_CORE_NATIVE=1` enables the native growth engine; Python remains the
  default when the variable is unset.
- `RXN_CORE_VERIFY_REPAIR=1` runs the repair kernel and the Python search side
  by side and asserts they agree.
- `RXN_CORE_VERIFY_ROLES=1` asserts the incremental role cache of the Python
  engine against the from-scratch definition.

## Verification

- `tests/test_native_engine.py`: every `grow_island` call of a full search on
  the tetraphenylmethane benchmark through both engines, plus equality of the
  typed results.
- `tests/test_freeze_native.py`, `tests/test_repair_native.py`,
  `tests/test_autgrp_native.py`: differential tests of the kernels.
- `bench/record_grow_calls.py <case> <out.pkl>` then
  `bench/compare_grow_calls.py <out.pkl>`: record every growth call of a
  Python run and replay it through the engine.
- `bench/replay_harness.py record/compare`: end-to-end pool identity.

## Third-party code

`native/nauty/` contains the nauty 2.8.8 sources (Brendan McKay and Adolfo
Piperno, Apache License 2.0; see `native/nauty/LICENSE-2.0.txt` and
`native/nauty/COPYRIGHT`) as shipped inside the pynauty source distribution,
so the compiled engine uses the same nauty version as the pynauty package the
Python engine calls.
