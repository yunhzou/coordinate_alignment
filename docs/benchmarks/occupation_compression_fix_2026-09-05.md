# Exact competitor compression: implemented and verified

Both occupation-expansion timeouts now finish with their original settings.
The change is in the retro occupation consumer, not native AAM island growth,
seed selection, branch admission, or chemical matching rules.

## Implementation

`occupation_orbit` accepts an optional observable-atom set. Augmented detection
supplies the real target atoms; other callers retain full-output semantics.

1. Join atom indices connected by any recorded permutation generator.
2. Keep every atom in a component intersecting the observable target distinct.
3. Give the remaining, permanently unobservable atoms one label in occupation
   keys. Preserve multiplicities, fragment classes, attachments, and bonds.
4. Walk the recorded stages using these keys, retaining a complete injective
   witness and permutation action for each state.
5. Reuse state/generator pairs evaluated in previous stages.

The unobservable complement is invariant under every recorded generator.
Consequently key projection commutes with their actions: two equivalent keys
cannot later produce different observable keys under the same action. Closing
over all stages is conservative; it may keep more competitors distinct than
strict chronological reachability requires, but cannot lose a target choice.

Only equality of permanently unobservable image labels changes. The groups and
raw search paths are unchanged, and contain the competitor alternatives. No
independent component products, random witnesses, hydrogen removal, truncation,
new search caps, or fallback paths were introduced.

## Full affected-precursor reruns

Same BIAN target, explicit H, tolerance 0.5, branch cap 100, no sweep, no new
candidate/seed limit; one 48-CPU allocation per source. Job 432164, tasks 2 and 4.

| Source | Previously | Detection | Typed checkpoint | Complete saved result |
| --- | --- | ---: | ---: | ---: |
| INVENTORY-000538 | Exceeded ten-minute watchdog | 39.91 s | 14.06 s | 123.58 s |
| INVENTORY-000594 | Exceeded ten-minute watchdog | 58.58 s | 0.93 s | 65.40 s |

The final column is elapsed time inside the source worker, including checkpoint,
archive construction, encoding, and compression. Slurm allocation times were
2:11 and 1:12 including startup/termination. Serialization is now most of the
remaining time for 538; it is not being counted as AAM.

The saved results contain 17,551 and 3,388 candidates respectively. Both retain
`status=capped` and `complete=false`: completing the computation is not a claim
that a branch-limited search exhausts every matching possibility. Their gzip,
JSON, row IDs, and source SMILES were verified against the bank.

Combined bank progress is now **1,918/1,919** saved. The separate native-growth
timeout `INVENTORY-000620` was not rerun or changed by this fix. Existing bank
records were not regenerated.

## Correctness and measured expansion improvement

- All **1,472 captured calls** were replayed against the pre-fix native logic
  preserved in the diagnostic extension. Every observed fragment relation and
  **exact source-to-target atom-pair set** is identical.
- Explicit states: **814,748 → 16,537**, with competitor choices retained in
  the original compressed groups.
- Sum of per-call elapsed timings: **1,949.86 → 15.09 seconds**, about **129×**
  for this occupation stage. These are not CPU accounting times or full-bank
  wall times. The verification suite used 48 workers and took 66.91 s overall,
  including old-logic replay, comparisons, and saved output states.
- Four isolated real-call comparisons also verify that omitting the observable
  set preserves every original full state, witness, action, and output order.
- 350 randomized observed-quotient comparisons against an exhaustive Python
  oracle, including actions crossing the target/competitor boundary.
- A specific two-stage test ensures competitors that can later enter different
  target positions are **not** merged prematurely.
- The full 256-test suite passed, followed by four new explicit-H molecular
  integration tests: methane/ethane, CO2/acid, water/acid, and alcohol matching.
  Those compare full expansion with compression, including candidate identities,
  cap diagnostics, and complete unchanged search graphs.

## Reproduction and saved results

Rebuild the native bookkeeping extension after updating source:

```sh
.venv/bin/python setup.py build_ext --inplace
```

`bench/verify_observed_occupations.py` checks isolated calls against the pre-fix
`_occupation_probe` extension. `bench/verify_occupation_suite.py` validates an
entire capture directory in parallel. Both save full and compressed outputs.

Cluster-relative artifacts:

- `data/retro_runs/occupation_diagnosis_20260905/fixed_verification/`: isolated
  comparisons and complete output states.
- `data/retro_runs/occupation_diagnosis_20260905/fixed_suite/`: all 1,472 input
  references, per-case timings/comparisons, and saved full/compressed states.
- `data/retro_runs/full_bank_core_fixed_20260905/occupation_fixed/`: the two
  completed source archives, full typed detection checkpoints, logs, summaries,
  and validated record hashes.
- `data/retro_runs/full_bank_core_fixed_20260905/combined_progress.json`: combined
  bank progress; the original frozen resume list remains untouched.

Compact timings are in
[`occupation_compression_fix_2026-09-05.json`](occupation_compression_fix_2026-09-05.json).
