# Output information and continuous weights for TS structures

## Conclusions

- Our output provides reusable matched-fragment structure, conditional correlated
  symmetry actions and search provenance. These are computational information,
  not extra experimental evidence about a reaction mechanism.
- Continuous weights preserve real TS information that binary connectivity loses.
  Our separate TS stage uses it with vibrational modes. Better TS mapping accuracy
  has **not** been demonstrated by this audit.
- The synthetic exact oracle remains saved but is excluded from paper evidence,
  as requested on 2026-09-09. No new oracle or full benchmark was run here.

## Output: core AAM versus downstream interpretation

| Information | Our saved interface | Use |
|---|---|---|
| Matched pieces | `FragmentTransition.match`, `preserved_bonds`, `FragmentMatch.r_atoms` and representative pairs | Query connected pieces; reuse them in partial matching and retro assembly |
| Conditional alternatives | Fragment domains, fixed atoms, multiplicity, finalized `target_generators` | Apply correlated motions, rather than independently shuffling atoms |
| Dependencies | `AAMSearchGraph` states/transitions; `SearchPath.realize` | Replay choices while transporting subsequent decisions into the correct frame |
| Search scope | Seed orders, cuts, tolerance, limits, terminals and `SearchStop` | Explain caps and preserve reproducibility; not a completeness guarantee |
| Physical input | `AAMResult.problem` retains coordinates, full continuous WBO matrices and explicit H | Re-score stored possibilities without reconstructing chemistry from a displayed witness |
| Optional interpretation | Mechanism grouping, RP geometry/chirality selection, TS analysis | Apply additional evidence separately from the core fragment search |

The DAG is a **computational search history**, not a chemical trajectory. Fragment
boundaries depend on search choices and tolerance; they are not proven intermediates
or uniquely conserved chemical units. Finalized groups describe recorded families,
not completeness across every seed/cut. Marginal domains alone do not authorize
independent permutations.

SLAP's native results contain refined graph pairs, costs and internal assignment
alternatives (`lap_sols`, group-pair/fingerprint information). It is not merely one
mapping or independent atom choices. Its schema does not expose our conditional
fragment DAG/group-action interface; it retains only its lowest reached native-cost
results. This does not prove our mapping set contains theirs.
[Pinned source](https://github.com/shin1koda/slap-mapper/blob/ea248fd9494f52f4865193e87a98cc92c62b5f9e/src/slapmapper/core.py).

The old native XYZ wrapper did not serialize every internal assignment object.
The comparison above uses the full source-level result, not that reduced archive.
Code inspected: `src/rxn_core/domain.py`, `search_graph.py`,
`alignment/post_aam.py`, `fragment.py`, and `ts.py`.

## Actual cached TS: vanadium case 114

These values come from the first archived mapping interpretation of
`pr7.V.dodh_ts910`, checked against cached R and reference-TS WBO arrays. Its saved
mapping supplies the R/P/TS atom correspondence. A reference TS geometry does not
make that algorithm-generated correspondence independently curated ground truth.

| Mapped bond | R WBO | TS WBO | P WBO | Binary R / TS / P at floor 0.2 |
|---|---:|---:|---:|---|
| Breaking V–O, R atoms 8–15 | 1.216286 | 0.599188 | 0 | 1 / 1 / 0 |
| Breaking O–H, R atoms 12–22 | 0.909423 | 0.397667 | 0 | 1 / 1 / 0 |
| Forming V–O, R atoms 8–12 | 0 | 0.548400 | 1.146039 | 0 / 1 / 1 |
| Forming O–H, R atoms 15–22 | 0 | 0.406296 | 0.920104 | 0 / 1 / 1 |

Binary TS connectivity calls all four present. Continuous values distinguish
partly weakened/developed contacts from endpoint bond strengths. That is an
information advantage over binarization, not a measured accuracy advantage.

Our growth checks floating-point edge differences within `iso_tolerance`, but
contacts below the connectivity floor are not active edges. At tolerance 1.0,
many positive weights intentionally remain compatible. Symmetry also uses
tolerance-based weight classes. The search therefore does not exploit every
fractional distinction or automatically prefer the closest WBO. Tightening
tolerance may exclude legitimate endpoint-to-TS matches.

The separate `analyze_transition_state` stage already uses per-event WBO progress
between endpoints, weighted event bond directions, and alignment with imaginary
vibrational modes. It returns core assignments, endpoint provenance, mode/frequency,
overlap, progress and per-event terms. That is a structural score, **not** an energy
barrier, probability, IRC coordinate or proof of endpoint connectivity. Another
mapper's assignments could also feed such a scorer; this is not exclusive to AAM.

## Fractional-cost diagnostic

The pinned weighted core accepts floating weights but builds an integer assignment
cost matrix, truncating neighborhood costs from direct WBO input. Weighted SMILES
supplies scaled discrete orders; native XYZ supplies binary adjacency. This is an
implementation limitation for continuous input, not inherent to linear assignment.
[Cost builder](https://github.com/shin1koda/slap-mapper/blob/ea248fd9494f52f4865193e87a98cc92c62b5f9e/src/slapmapper/core.py#L595-L625).

Unmodified initial-cost-builder calls on three actual cached R/reference-TS pairs:

| Case | Atoms | Initial matrix cells | Positive costs cast to zero |
|---|---:|---:|---:|
| `pr1.tempo_ts1` | 26 | 238 | 168 |
| `pr16.carbocation_ts5` | 18 | 162 | 122 |
| `pr7.V.dodh_ts910` | 28 | 228 | 98 |

For example, the first case's iodine neighborhood difference 0.852671 becomes
zero. All 628 cells lose fractional parts; 388 positive costs become zero. These
include tiny numerical differences, not necessarily chemically important ones.
Later refinement retains other information; these counts are not failure rates.

The corrected numerical audit used 0.051 CPU seconds after imports; it is not a
mapper timing benchmark. `audit_v2/` saves inputs, complete initial matrices,
atom orders, hashes and archived event terms. The first attempt met `gt=null` in a
historical mechanism. Unscored mechanism IDs are now explicitly recorded, not
assigned invented events; see `initial_attempt.txt`.

## What remains to demonstrate

Freeze real R/TS inputs and budgets; compare binary versus continuous input within
our engine, then the same WBO graphs with a separately labeled, validated
continuous-cost competitor. Measure expert-checked core mapping recovery and
distinct plausible alternatives versus CPU time. Include mode/IRC evidence where
available. No current evidence isolates the accuracy benefit of weights alone.

The 140-case mapper benchmark evaluated **R-to-P**, despite TS names in case IDs.
It is not an R-to-TS/P-to-TS accuracy benchmark. Neither mapper core was modified
here. **75 existing tests passed in 1.76 seconds**, covering weighted subgraph
matching, search graphs, hierarchy transport and the typed AAM/TS API.

## Reproduce the numerical audit only

```bash
env PYTHONPATH=/project/yunhengzou/coordinate_alignment/aam_benchmarks/slap_guided_pilot_20260909/dependencies:src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --kill-after=5s 60 .venv/bin/python bench/audit_ts_weight_information.py --source /h/399/yunhengzou/appendix_final/bgcp_rerank_latest --slap-core /project/yunhengzou/coordinate_alignment/aam_benchmarks/competitor_env_20260908/lib/python3.10/site-packages/slapmapper/core.py --output NEW_DIRECTORY --cases pr1.tempo_ts1 pr16.carbocation_ts5 pr7.V.dodh_ts910
```

Full report path:
`/h/399/yunhengzou/coordinate_alignment/reports/ts_weight_information_20260909/README.md`
