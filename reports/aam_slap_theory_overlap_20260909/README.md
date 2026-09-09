# AAM versus SLAP: search and retained alternatives

**Tolerance follow-up:** the [full tolerance1 rerun](../elementary140_tol1_20260909/README.md)
recovers case135's 5-event pattern. Results below describe the original
tolerance0.5 archives, not a general inability of AAM to find that mapping.

[Offline 3D mapping comparison viewer](viewer.html): cases 135, 59 and 64.
Opens on case 135 (AAM 7 events; SLAP 5), using the compact dark AAM layout:
two large original-geometry R/P panels, linked rotation, AAM/SLAP selection,
explicit hydrogens, event overlays and collapsible per-bond WBO tables.
Source-identity colors are optional; element colors are the default.
These are saved concrete witnesses, not ground truth.

Our AAM is not simply a weaker SLAP. Their saved solution sets differ in both
directions. SLAP is a strong, inexpensive approximate mapper; our experiments
show additional low-edit alternatives and better recovery on some controlled
inputs, not universal superiority or completeness.

## Theory: what differs

| Aspect | Our current AAM | SLAP 1.0.0 |
|---|---|---|
| Search | Sequential conditional fragment matching; seeds and sweep cuts diversify paths | WL-like refinement and sequential linear assignment, with symmetry breaking |
| Restriction | Earlier fragments constrain later fragments; finite seeds/cuts and branch caps can miss mappings | Approximate refinement, assignment choices and pruning can miss globally best mappings |
| Results | Branch/path families with correlated symmetry actions | Multiple label-group results, plus internal assignment alternatives; not just one mapping |
| Retention | Alternative fragment paths, including different edit scores | Removes results above the minimum native cost reached by its search |
| Symmetry | Explicit correlated actions within stored families | Label groups, internal LAP alternatives, and WL-signature-based output deduplication |
| Guarantee | No guarantee of all chemically valid or minimum-edit mappings | No guarantee of all chemically valid or globally minimum-edit mappings |

SLAP describes itself as **approximate** graph matching. Solving each linear
assignment efficiently does not make the full graph-matching problem exact.
Its source perturbs assignments to obtain alternatives and stops on repeated
assignments or increased cost; this is not exhaustive permutation enumeration.
Its isomorphic-result filter uses a WL signature, which is not a general exact
graph-isomorphism certificate. This is a theoretical limitation, not evidence
that a particular holdout result was incorrectly deduplicated.
[SLAP source, pinned revision](https://github.com/shin1koda/slap-mapper/blob/ea248fd9494f52f4865193e87a98cc92c62b5f9e/src/slapmapper/core.py),
[authors' preprint abstract](https://www.cambridge.org/engage/chemrxiv/article-details/691d9f5ea10c9f5ca1852be4).

Compression alone is not a scientific advantage if the represented mappings are
the same. Any finite family can, in principle, be listed as individual mappings.
The useful questions are which distinct low-edit possibilities survive, whether
correlations remain valid, and their computational cost. Neither implementation
has been proved a superset of the other.

The native SLAP XYZ adapter uses distance-derived binary connectivity; our
holdout uses WBO graphs. SLAP's core also supports weighted graphs: the adapter
difference must not be portrayed as a fundamental binary-only limitation.

The later [continuous-weight audit](../ts_weight_information_20260909/README.md)
identifies integer assignment-cost storage in that pinned implementation. Weighted
discrete bond orders and direct fractional WBO input must be distinguished.

## Saved elementary holdout: 140 cases

No reference atom mappings exist here. These are structural comparisons, not
chemical accuracy. All explicit H atoms participate in feasibility and bond-event
scoring. Heavy mapping relations are compared modulo endpoint symmetries that
preserve the shared score (bond presence above 0.2; order-change threshold 0.5).
All saved SLAP heavy labels are singleton. We do not claim equality of the full
hydrogen assignment spaces; native `lap_sols` were not saved by the earlier adapter.

- 159 of 176 SLAP heavy mapping classes are represented in our saved families.
- 11 classes are excluded from those families, in cases **59, 64 and 135**.
- Six classes remain unresolved in cases **4, 6, 7, 27 and 31** under the bounded
  query budget. Unknown is not absence.
- All SLAP classes are represented for 132 of 140 cases.
- In **97 cases**, an AAM heavy mapping pattern absent from saved SLAP outputs
  ties our lowest observed bond-event count. This rises to 99 cases within +1,
  132 within +2, and 137 within +5 events.

Thus the difference is not merely duplicate atom labels or exclusively high-edit
extra outputs. It still does not establish that these alternatives are chemically
correct. Nor does absence from one saved run mean unreachable with another policy.

The [preceding score comparison](../elementary140_output_comparison_20260908/README.md)
found equal best observed event counts in 133 cases, fewer for ours in six, and
fewer for SLAP in one (135). SLAP hydrogen assignments were refined within saved
label pools before that comparison. Our minimum is over saved witnesses, not a
proof of the optimum over every compressed family.

### Overlap procedure and verification

Stream both directions' saved terminal witnesses, compare exact colored-graph
certificates, and query compressed families for still-missing SLAP classes.
An entire path can be excluded without expanding it when every allowed action
fixes the heavy relation or is an endpoint score-preserving automorphism: its
representative already determines its equivalence class. Other paths use the
existing symbolic family query. No production bijection enumeration is added.

The 140 one-CPU analysis tasks used up to 32 concurrent slots, with a 180-second
soft query budget and a 300-second process watchdog. Summed analysis wall time
was 1159.65 seconds, not mapper search time or a CPU benchmark. Five unstarted
Slurm tasks were relocated; completed tasks were not rerun. Case 123 uses the
previously saved cap-200 follow-up; other cases retain cap 100.

An independent audit checked all **298 positive/extra mapping witnesses**, full
element-preserving assignments, symmetry certificates and extra-witness scores.
Negative family exclusions depend on the recorded path-action analysis; the
witness audit is not a second independent proof of those exclusions.

## Identical-input exact oracle control

**Internal diagnostic only; excluded from paper evidence by the user's decision
on 2026-09-09.** The test and all saved results are retained for debugging. It is
not an independently established benchmark and must not be used to support the
paper's comparative performance claims.

100 reproducible connected six-vertex binary graph pairs (all vertices C,
maximum degree four, RNG 20260909) were supplied directly to both engines.
These are abstract graphs, **not chemical reactions**. Exact edge-edit minima
were obtained from all 720 permutations solely for this tiny testing oracle.
Neither search was given the oracle mapping.

| Search policy | Exact minimum recovered |
|---|---:|
| SLAP binary, one original-order run | 90/100 |
| SLAP binary, 10 input orderings, both directions (2,000 calls) | 92/100 |
| Unchanged AAM, seed 10, cap 100, tolerance 0.5, sweep cuts | 100/100 |

The expanded SLAP misses are 13, 38, 39, 43, 47, 62, 71 and 91.
Every saved mapping and event count, including expanded SLAP, was independently
rescored; the exact oracle was recomputed. This removes endpoint perception as
the explanation for these differences. It is not an equal-time comparison,
an asymptotic speed result, or evidence of universal AAM optimality.

## Conclusion for the paper

Defensible: the methods retain different alternatives; ours demonstrably retains
additional low-edit heavy mapping patterns, while SLAP also finds patterns absent
from our saved search. SLAP's efficiency is a real strength. Chemical superiority
on this holdout requires expert validation, and matched-budget comparisons remain
necessary for a joint accuracy/speed claim.

Not defensible: SLAP returns only one bijection; compressed storage alone is
more expressive; our search covers everything; or these reference-free tests
prove chemical accuracy.

## Reproduction and artifacts

Scripts: `bench/elementary_family_overlap.py`, `bench/aam_small_graph_oracle.py`,
`bench/audit_aam_slap_overlap.py`. Only benchmark scripts were added; no core AAM
or SLAP algorithm was changed. Compact evidence is in `saved_evidence.tar.gz`.

Full cluster paths:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary_family_overlap_20260909
/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_slap_small_graph_oracle_20260909
/project/yunhengzou/coordinate_alignment/aam_benchmarks/elementary140_output_comparison_20260908
```

Full AAM checkpoints remain in these runs/their recorded source directories.
Compact evidence includes per-case overlap results, oracle inputs and mappings,
AAM summaries, submitted driver snapshots and Slurm submission records.
