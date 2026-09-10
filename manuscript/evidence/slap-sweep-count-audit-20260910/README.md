# SLAP sweep counting audit

The standalone SLAP sweep count is **1,795/1,851 = 96.9746%**, or **97.0%** to
one decimal place. Fresh re-scoring of **572,753 saved candidate outputs** in
**352,210 committed call records** reproduced every saved case outcome and
direction/mode hit list. No discrepancy in the recovery numerator was found.
This audit used the frozen strict evaluator from the original campaign; it did
not introduce a second equivalence metric or rerun the mapper.

## Recovery accounting

| Output union | Recovered cases | Percentage |
| --- | ---: | ---: |
| Uncut, both directions and both modes | 1,661 | 89.74% |
| Binary mode, both directions, uncut plus cuts | 1,785 | 96.43% |
| Weighted mode, both directions, uncut plus cuts | 1,762 | 95.19% |
| Binary and weighted sweep union | 1,795 | 96.97% |
| Sweep plus the separate older expanded-SLAP experiment | 1,797 | 97.08% |

Binary and weighted recovery sets overlap on **1,752** cases. Binary alone adds
33 cases, weighted alone adds 10, so `1,785 + 1,762 - 1,752 = 1,795`.
The uncut set is contained in the union; the sweep adds 134 cases. In this
archive, the cut outputs alone also recover all 1,795 cases. The protocol still
includes the uncut call and must be described consistently.

The older expanded experiment adds **1071 and 1595**, which the standalone sweep
does not recover. The 1,797 result must not be substituted for the 1,795 ablation
result: it has an additional input-ordering/search budget.

Each case is counted once if a single saved candidate matches the reference
relation. Atom correspondences from different candidates are not pooled into a
synthetic success. The metric preserves the original endpoint chemistry and
compares the complete heavy-atom relation modulo endpoint chemical symmetry,
including unmatched heavy atoms. Explicit H participates in search, but this is
not full-H reference accuracy or top-1 prediction accuracy.

## Completed misses and unresolved cases

The **56 standalone nonrecoveries** divide into:

- **52 fully swept, valid cases without a reference witness** under the tested
  policy. They are completed misses for this budget and policy, not proof that
  every possible SLAP configuration must fail.
- **4 unresolved cases** with unfinished variants: **1503, 1665, 1723, 1776**.
  A timeout does not establish that further computation would recover the
  reference. These cases remain in the full denominator and are not credited
  as successes.

For the separate 1,797 union, the corresponding split is **50 completed misses
plus the same 4 unresolved cases**. Confusing these two unions would produce
an incorrect completion breakdown.

All identifiers below are zero-based. A cell gives committed calls/planned calls,
including the uncut call.

| Case | R-to-P binary | R-to-P weighted | P-to-R binary | P-to-R weighted |
| --- | ---: | ---: | ---: | ---: |
| 1503 | 151/151 | 151/151 | **1/42** | **1/42** |
| 1665 | **9/144** | **6/144** | **6/70** | **8/70** |
| 1723 | 145/145 | 145/145 | **0/59** | **0/59** |
| 1776 | 106/106 | 106/106 | **6/68** | 68/68 |

The bold cells are **nine** unfinished variants, with **661** calls not committed.
Their task records have exit code 124 and approximately 300 seconds elapsed.
The frozen batch driver implements this status by killing the process group
when its 300-second watchdog expires. The watchdog is for a complete
case/direction/mode variant across its cuts, including startup and persistence;
it is not a fresh 300 seconds for every cut.

The statement that only P-to-R needs attention is incomplete: **case 1665 also
timed out in both R-to-P modes**. Case 1723 produced no committed P-to-R uncut
call. Case 1776's P-to-R weighted variant completed normally.

Across the entire campaign there were **65** watchdog-limited variants affecting
**36** cases, not just four cases. Of those 36 cases, **32 already have a valid
reference witness from another completed output** and legitimately count as
recovered. The other four are the unresolved nonrecoveries above.

Of 355,104 planned calls, 352,210 were recorded, leaving 2,894 unfinished calls.
Two recorded calls raised upstream mapping exceptions, both on case 836;
another output recovers that case. All **572,753 returned candidates** passed
the frozen format and endpoint checks. There are 1,814 cases with all planned
calls complete and no mapping/validation error.

## Rerun handoff

[unresolved_variants.csv](unresolved_variants.csv) records the nine relevant task
slots, modes, completion counts, and the first uncommitted ordinal/cut. No
mapping search was launched by this audit. Completing only those nine variants
could raise the standalone count to at most **1,799/1,851 = 97.19%**; that is an
upper bound on the possible gain, not a prediction.

A targeted follow-up should use a new output directory and retain the original
300-second results. Keep the pinned SLAP source, input ordering, explicit-H
handling, modes, cut enumeration and evaluator unchanged, and state the revised
budget. Preserve every returned alternative and job status. Select the unfinished
variants in advance and run the prescribed cuts without stopping based on a
reference hit. Report the resulting union as an extended-budget follow-up.

Finishing these four cases resolves their recovery uncertainty but does not
complete the whole campaign: 56 other unfinished variants belong to already
recovered cases. A claim about total sweep completion or complete campaign CPU
would require accounting for those variants too. Recorded mapping CPU omits
interrupted in-flight work, even for a case recovered elsewhere.

Useful follow-up evidence consists of the exact command/configuration and source
revision, saved raw outputs, task exit/status records, per-case witness locations,
an evaluator hash, and CPU totals with their case set and timing scope. For the
140 coordinate cases, the documented absence of annotated reference mappings
still limits claims to feasibility, scores, and alternative coverage.

## Three-seed timing cross-check

The quoted AAM figures are correct **on 1,837 mutually completed Golden cases**:

| Timing set | AAM, three seeds | AAM, ten seeds | SLAP sweep |
| --- | ---: | ---: | ---: |
| AAM-only common set: 1,837 cases | 59,466.8 CPU-s | 176,323.7 CPU-s | — |
| Paper's shared comparison: 1,807 cases | 53,790.9 CPU-s | 158,774.6 CPU-s | 31,831.6 CPU-s |

Three-seed AAM's full 1,851-case completed Golden total is 70,570.7 CPU-seconds.
Therefore 59k is a paired-subset total, not the full three-seed workload cost.
The CPU instrumentation also differs between AAM and SLAP, as documented in
the manuscript; these values are not controlled wall-time speedups.

Both seed settings verify **1,836/1,851**, but the verified sets differ. Three
seeds newly verify 1786, 1799 and 1823; ten seeds uniquely verify 850, 1314 and
1568. Say "the same number of verified references," rather than "identical
coverage" or "three seeds reproduce every ten-seed result."

## Evidence and reproduction

- [summary.json](summary.json): re-scored recovery, completion and timing totals;
  original dataset/evaluator/audit-script hashes; recovered case indices.
- [cases.csv](cases.csv): all 1,851 case-level outcomes and completion checks.
- [case_checks.json.gz](case_checks.json.gz): journal hashes, exact per-variant
  coverage/status records and a checksum-verified native witness frame per
  recovered case. Frames were hash-checked; the audit does not claim to have
  decoded and independently re-exported every native LAP record.
- [related-count-checks.json](related-count-checks.json): separate expanded-SLAP
  union and AAM timing/case-set checks, with their source hashes.
- Audit implementation: `bench/audit_slap_sweep_counts.py` in the repository.

The recheck verifies contiguous unique cut ordinals and actual cut identities
against the saved input graphs, task/status counts, reference-file identity,
all candidate signatures, every saved hit list, and one native archive frame
hash for each of the 1,795 recovered cases. All assertions passed.

From the repository root, with the original full archive and Python dependencies
available, reproduce into a **new** output directory:

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../.venv/bin/python bench/audit_slap_sweep_counts.py \
  --run /project/yunhengzou/coordinate_alignment/aam_benchmarks/slap_edge_sweep_full_20260910 \
  --output manuscript/build/slap-count-audit-repeat --workers 8 \
  2> manuscript/build/slap-count-audit-repeat-stderr.log
```

This re-scores existing outputs using the frozen evaluator; it does not change
the frozen campaign or establish a different definition of chemical correctness.
The compact manuscript bundle contains the audit evidence, while reproduction
also needs the full repository and original archives.
