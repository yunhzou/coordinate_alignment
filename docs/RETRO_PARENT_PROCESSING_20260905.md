# Parent-side processing after the AAM memory repair

## Scope

This follows `fa198e2`/`b891676`. The core AAM/C++ search, matching rules,
explicit hydrogens, cap100, no-sweep policy, and occupations are unchanged.
The baseline is the immediately preceding exact workflow, **not** the old
approximate recommendation scanner.

## Measured and applied

Replayed the saved 26,487-occupation family at
`data/retro_runs/aam_memory_boundary_20260905/shared/2.post.pkl`, with one CPU
for each measurement. No AAM was rerun for these comparisons.

| Operation | Before | After |
| --- | ---: | ---: |
| Build archive | 14.161 s | 4.572 s |
| JSON encoding | 8.352 s | 5.712 s |
| Compression | 3.603 s | 1.543 s |
| Total above | 26.116 s | 11.827 s |
| Peak RSS through encoding/compression | 1,343.7 MB | 647.4 MB |
| Compressed bytes | 1,318,398 | 5,443,069 |

The faster gzip setting trades file size for throughput; it is lossless. All
135,756,715 decompressed bytes are identical, SHA-256:
`62450997bb25b08055e814c43daf8aef463b8454e43d6a56332fe4abbaea72e6`.
The measured stage timings use `bench/retro_archive_profile.py`. The catalog
writer additionally writes each record into a gzip stream rather than joining
an entire batch into one large string.

Changes:

- Archive construction shares encoded arrays/path records and directly reuses
  immutable state fields. It no longer recursively copies every state mapping
  with `dataclasses.asdict`. JSON schema v7 and its list-valued fields stay the
  same. Encoded arrays may share storage; treat archives as read-only snapshots.
  Mutating an encoded array cannot mutate its typed input tuple.
- Merge provenance accumulates in lists, then freezes once per retained
  candidate. It previously constructed/replaced candidates for every discovery
  and repeatedly concatenated growing provenance tuples. A saved-family replay
  reduced merge time from 0.838 s to 0.512 s, with all 1,177 retained candidates
  and all derivations equal, in the same order.
- The writer omits circular-reference checks because this schema represents
  graph links using IDs, not cyclic Python containers. It uses gzip level1 and
  preserves arbitrarily large integer symmetry counts.
- `_detect_one` and `_record_detection` separate search from persistence.
  `_search_one` retains its external behavior. `bench/retro_tail_stages.py` saves
  typed detection and archive checkpoints, records individual stage times, and
  emits internal stack traces each minute. It never silently substitutes a
  partial result or repeats a missing search during an archive replay.

Native RapidJSON was evaluated, but its measured JSON time (7.77 s) did not
justify a production dependency. No native JSON dependency was added.

## Full precursor evidence

`INVENTORY-000400`, same target and bank row as the prior OOM/timeouts:

- Job432070: detection 316.585 s; complete typed-result checkpoint 52.335 s.
  This search started before the merge-list change, so it does not measure
  that change's full-precursor benefit.
- Job432074: archive construction 119.371 s; archive checkpoint 39.865 s.
  This reused the typed checkpoint and performed no AAM.

The checkpoints are roughly 1 GB each, at
`data/retro_runs/aam_memory_boundary_20260905/staged_tail/2/detection.pkl`
and `archive.pkl`. These are trusted local Python intermediates, not inputs
to accept from unknown sources. Separate checkpoint/load times must not be
confused with algorithm time or hidden inside fresh-run speed claims.

All jobs have 10-minute allocations and a 580-second process watchdog. No
four-source rerun or full-bank speed claim follows from the single-family
measurements. The full precursor's stages are explicitly timed separately.

## Validation

`pytest tests -q`: **219 passed**. Regression coverage includes immutable-state
reuse, encoded-array sharing without input mutation, one freeze per candidate,
preserved derivation order/count, separate search/archive orchestration, and
lossless gzip output with shared arrays and a 101-digit symmetry count.

Evidence directories under `data/retro_runs/aam_memory_boundary_20260905/`:
`archive_before`, `archive_profile`, `archive_shared`, `archive_native`,
`archive_final`, `merge_profile`, and `staged_tail/2`.
