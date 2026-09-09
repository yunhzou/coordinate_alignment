# Learning from SLAP: proposal-directed AAM search

Exploratory experiment on the 140 native XYZ/WBO elementary steps. **No production
AAM code changed.** This is a faster experimental workflow, not a claim of equal
compressed-family coverage or chemical accuracy.

## Fresh workflow result

The selected workflow was rerun from scratch on all 140 cases, after the ablations:

- **140/140** reproduce the full workflow's best saved bond-event score.
- **102.70 CPU-seconds total**, averaging **0.734 CPU-seconds per reaction**.
- **1.128 seconds mean elapsed**, maximum **6.757 seconds**, one CPU per case.
- Includes proposal generation, within-pool H refinement, AAM matching, symmetry
  finalization, WBO scoring and internal graph/proposal artifact writes. Excludes
  initial module imports/input loading and final summary writing. Queue excluded.
- All 1,293 emitted concrete witnesses independently rescored; assignments checked.

The preceding full bidirectional benchmark used **9,562 CPU-seconds** at cap100
(another 34 seconds for its separately reported case123 cap200 follow-up).
The observed CPU ratio is approximately **93×**. Hardware, process layout and
timing scope differ: this is not a microbenchmark of the same core operation.
The pilot preserves observed best scores, **not proven equality of alternatives**.

The reference here is the previous best saved AAM score (including its explicitly
labeled cap200 follow-up for case123), not an exact oracle or chemical ground truth.
The fresh workflow itself uses cap100 throughout, including case123.

## What we learned from SLAP

SLAP's inexpensive sequential linear assignments provide a global proposal before
expensive fragment search. The useful transfer is **where to allocate search**,
not merely porting more functions to C++. Its unchanged core was invoked directly
on binary connectivity derived from the same WBO matrices (edge presence >0.2).
No SMILES conversion or distance-based XYZ bond perception was used.
[Pinned SLAP source](https://github.com/shin1koda/slap-mapper/blob/ea248fd9494f52f4865193e87a98cc92c62b5f9e/src/slapmapper/core.py).

SLAP's binary objective is only a proposal. Actual candidate scores always use
original WBOs, including H, with order-change threshold0.5. AAM matching tolerance
remains1.0. Native SLAP H assignment pools are refined with the existing bounded
symbolic scorer; no H bijection enumeration. All pool bounds resolved in the pilot.

The frozen full-run profile explains the opportunity: 6,737 CPU-seconds in
cut/seed search (~70%), 2,128 in symmetry finalization (~22%). Explicit merge
phases total117 seconds (~1%; worker merges are also included in the search
phase, so do not sum these as disjoint categories). Optimizing merge alone
cannot close this gap.

## Selected experimental algorithm

1. Produce SLAP mapping proposals and score their H-refined representatives.
2. Run one ordinary, uncut AAM seed in each direction. Retain its compressed graph.
3. Take the best SLAP proposal and each available directional AAM best witness.
4. Identify their proposed broken/order-changed edges on each source endpoint.
   Reverse-direction processing also exposes bonds formed in the original direction.
5. Deduplicate these edge choices and run one ordinary AAM seed per selected cut.
6. Keep the candidate union, original WBO scores and all compressed cut graphs.

No atom pair is hard-anchored from SLAP. Cut selection sees neither the previous
full-search answer nor a reference mapping. The fresh driver reads the full-run
score only after its search and timing are finished. Single-cut searches reuse
the existing deterministic cut-seed policy and cap100. No new matching rule,
per-molecule exception or production fallback path was added.

This does not discover every cut the full sweep explores. Wrong or incomplete
proposals can miss useful alternatives. This experiment cannot certify completeness.

## Ablations

Both directions in all AAM rows; unchanged compressed matching machinery.
These CPU numbers are **component sums**, excluding import/setup, persistence,
benchmark certificate evaluation and local cut-list construction. Use the fresh
102.70-second measurement above for the selected workflow's broader cost.

| Configuration | Matches full best score | Measured component CPU, all140 |
|---|---:|---:|
| Direct SLAP on WBO connectivity, H-refined | 137/140 | 7.44 s |
| Uncut AAM, one random seed | 133/140 | 11.61 s |
| Uncut AAM, ten random seeds | 138/140 | 71.92 s |
| SLAP-ordered uncut AAM, one seed | 134/140 | 18.41 s |
| SLAP + uncut random-one-seed candidate union | 139/140 | 19.06 s |
| Above + proposed changed-edge cuts | **140/140** | **50.98 s** |
| Above + all edges incident to proposed event atoms | 140/140 | 147.60 s |

Changed-edge refinement selected874 directional cuts; incident-edge refinement
selected3144. Both used one seed per cut. The component-only ratio (~188×)
is **not** the fresh-workflow ratio (~93×).

### Ideas that did not work as hoped

- Softly ordering seeds away from SLAP's predicted event atoms was inferior to
  retaining independent random-seed and SLAP proposals. Do not adopt this rule.
- Agreement is not an optimality certificate: in case64 both cheap methods give
  six events, while full AAM and targeted-cut refinement find five. Refining only
  when the cheap methods disagree would miss this case. A retrospective replay
  of that trigger selects eight cases, adds629 CPU-seconds, and still misses64.
- Increasing uncut seeds from one to ten costs substantially more without
  recovering case64. Directing cut exploration is the more effective lever here.

## Alternative mappings: unresolved trade-off

The uncut one-seed AAM witnesses include126/140 tracked baseline-best heavy
patterns and80 tracked extra patterns; ten seeds give138 and103 respectively.
These are **positive witness-level lower bounds**, not exhaustive membership
tests of the compressed output families. A missing representative is not proof
the family lacks it. The extra-pattern probe uses the cap100 archive (including
case123); it is not a complete catalog of minimum-score alternatives.

Targeted cuts retain additional graphs, but their full alternative-pattern
coverage has not yet been measured. Matching140 best scores cannot replace that
test. The workflow is exploratory, and was designed while examining this same
holdout: **Golden/reference recovery and independent cases must be tested with
the policy frozen before promoting it or making paper accuracy claims.**

## Reproducibility and execution

- `bench/slap_guided_pilot.py`: ablations, proposal-driven cuts, scoring/audits.
- `bench/slap_local_workflow.py`: fresh selected workflow; no baseline answer read
  before search completion.
- `tests/test_slap_guided_pilot.py`: seed ordering retains all atoms and does not
  mutate assignments. Four relevant unit tests passed; local-cut audit checked
  4,068 witnesses, and fresh-workflow audit checked1,293.
- SLAP core pin above; isolated SciPy1.18.1 installed under the pilot directory.
  The initial pilot lacked SciPy and failed before matching; queued work was
  cancelled, dependency isolated, smoke test passed, and missing cases resubmitted.
  Those failed attempts are not included in reported successful compute totals.
- Tasks use one CPU,8GB,300-second process watchdog and10-minute Slurm limit.
  Only unstarted CONFIGURING tasks were relocated; completed work was not rerun.
  Full compressed graphs and native LAP alternatives were saved.

Compact evidence and frozen drivers are in `evidence.tar.gz`. Full archives:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/slap_guided_pilot_20260909
/project/yunhengzou/coordinate_alignment/aam_benchmarks/slap_local_fresh_20260909
```
