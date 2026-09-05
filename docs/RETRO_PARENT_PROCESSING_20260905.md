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
| Build archive | 14.161 s | 4.530 s |
| JSON encoding | 8.352 s | 1.248 s |
| Compression | 3.603 s | 1.441 s |
| Total above | 26.116 s | 7.219 s |
| Peak RSS through encoding/compression | 1,343.7 MB | 648.4 MB |
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
- The writer uses `msgspec`'s native JSON encoder and gzip level1. It preserves
  the existing JSON schema and integer symmetry counts, including the tested
  101-digit value. Unicode is emitted as UTF-8 rather than ASCII escapes; the
  decoded value is unchanged. No native codec fallback or witness sampling is
  involved. Install the catalog tools with `pip install -e '.[catalog]'`.
  Core AAM's Python objects and C++ engine have no dependency on this writer.
- `_detect_one` and `_record_detection` separate search from persistence.
  `_search_one` retains its external behavior. `bench/retro_tail_stages.py` saves
  typed detection and archive checkpoints, records individual stage times, and
  emits internal stack traces each minute. It never silently substitutes a
  partial result or repeats a missing search during an archive replay.

Native RapidJSON was evaluated, but its measured JSON time (7.77 s) did not
justify a production dependency. A subsequent `msgspec` test reduced this stage
to 1.248 s and was adopted. Both evaluated native encoders preserve the saved
family's bytes; production uses only `msgspec`. The standard-library-only
intermediate improvement measured 11.827 s for the three stages together.

## Full precursor evidence

`INVENTORY-000400`, same target and bank row as the prior OOM/timeouts:

- Job432070: detection 316.585 s; complete typed-result checkpoint 52.335 s.
  This search started before the merge-list change, so it does not measure
  that change's full-precursor benefit.
- Job432074: archive construction 119.371 s; archive checkpoint 39.865 s.
  This reused the typed checkpoint and performed no AAM.
- Job432075: standard-library JSON/gzip-level1 encoding 189.800 s, excluding
  its 38.638 s checkpoint load. Its full 40,062-candidate archive was saved.
- Job432079: native JSON/gzip-level1 encoding **69.304 s**, excluding its
  29.497 s checkpoint load. Both encoding jobs used one CPU on bosque67.
  Peak RSS including the loaded archive fell from 13,665.8 MB to 10,116.0 MB.
  The result is 130,923,317 compressed bytes (about 125 MiB), with all 40,062
  candidates and the original capped/incomplete flags retained.

The measured detection + archive + native encoding stage sum is **505.3 s**
(8.42 min), excluding checkpoint saves/loads and about 1.4 s final file write.
This is a sum of separately timed stages, not a newly timed end-to-end scan;
in particular the detection measurement predates the merge-list improvement.

The full native-written result is available at:
`/h/399/yunhengzou/coordinate_alignment/data/retro_runs/aam_memory_boundary_20260905/staged_native_tail/2/result.jsonl.gz`.

Streaming gzip integrity/checksum jobs432080 and432081 confirmed the entire
standard-library and native-written JSON streams are byte-for-byte identical:
SHA-256 `a40badac8912dfd0b16079158e3a7ea2c8a3cf6f3b232a7013fe645b4de86d49`.

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
`archive_final`, `archive_msgspec`, `merge_profile`, `staged_tail/2`, and
`staged_native_tail/2`.
