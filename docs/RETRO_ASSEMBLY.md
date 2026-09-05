# Fragment detection and geometric assembly

Given a target and a precursor bank, recommend sets of building blocks whose
matched fragments jointly support every target atom, including explicit H.
This establishes a geometric foundation, **not chemical feasibility, a reaction
mechanism, unique atom ownership, or a multistep synthesis route**.

```text
R bank + P target
    -> shared fragment matcher
    -> anchored AAM against P target + residual competitors
    -> saved fragment relations and correlated occupations
    -> exact coverage decision graph
    -> certified set ranking and construction-pattern grouping
    -> assembled viewer
```

## Detection owns matching

`fragment_matching.detection` calls `match_fragment` for its initial searches.
For each discovered family, `augmentation` copies the residual source graph
beside P, cuts the initial boundary/deferred bonds, and calls the existing
`find_islands` scheduler with the initial assignment anchored. That scheduler
uses the shared, native-enabled fragment matcher. There is no sweep here.

Singletons and hydrogen atoms go through the same AAM continuation as other
residuals. There are no hand-built singleton witnesses, retained-neighbour
attachment tests, independent component-witness products, or maximum-ownership
filters. A competitor is an actual graph component with a recorded assignment.

Each detection retains its search graphs and derivations. Projection carries:

- Source and target fragment atom sets and preserved source bonds.
- The actual mapping witness and the action relating it to the saved path.
- Proven equivalence classes of source fragment units.
- Leftover/competitor assignments and AAM cap diagnostics.

Conditional generator groups act in their recorded order; different branches
are alternatives, not a union of their generators. Projection distinguishes
joint fragment partitions and carried bonds, not just their union of atoms.
Equivalent source units are unordered within their proven class. Ten identical
singleton H units therefore do not expand into 10! atom assignments.

The occupation walk transports compact integer images. Hierarchies and exact
generators are transported only for surviving final occupations.

## Assembly owns combinations, not matching

`retrosynthesis.assembly.AssemblyProblem` indexes complete correlated occupation
relations. Identical occupation relations form a substitution pool; distinct
fragment partitions remain separate slots even if they cover identical atoms.
Repeated use of the same precursor at different occupations is allowed.

`CoverageDecisionGraph` shares equal `(next slot, covered atoms)` suffix
problems. Include/exclude edges represent the complete slot-set search. It
allows overlap and redundant distinct occupations. There is no assumed two-
or three-reactant limit, beam, marginal-coverage cutoff, per-pool shortlist, or
special treatment of a known answer. It is not an atom-bijection search.

Complete covers satisfy only the target-union condition. Overlap stays as a
support relation; it is never assigned to whichever reactant appeared first.
The viewer colors shared support grey. A target connection unsupported by a
matched fragment is a construction connection, **not a proven bond-formation
event**. Source cuts likewise are not a balanced reaction edit count.

## Ranking and exact stopping

Ranking compares full supplier sets in this order:

1. Number of distinct precursor structures, not total copies.
2. Symmetry-adjusted explicit-atom retention.
3. Direct explicit-atom retention, followed by deterministic supplier IDs.

Direct retention counts the target union once and counts every input copy in
the denominator. For symmetry adjustment, each input's atom cost is divided
by the maximum number of **disjoint whole retained-fragment copies** supported
by source automorphisms. This packing is solved exactly; independent atom-orbit
capacities are not treated as realizable copies. The adjusted set score divides
target atom count by the sum of these effective input costs.

Best-first traversal uses optimistic suffix input-cost bounds and the number
of distinct structures already chosen. It cannot overestimate a completion's
rank. A popped complete assembly is therefore the next ranked result over the
saved occupation index. Unit tests compare this order against exhaustive
enumeration. Exact colored incidence-graph certificates group construction
patterns up to target symmetry; fragment and source-copy roles remain distinct.

Display limits request the best patterns and supplier sets. They do not prune
the index or impose search budgets. The merge stops when that requested prefix
is certified, or exhausts the graph. `--exhaustive` writes every assembly.
Exact cover/packing remain combinatorial; no constant runtime is promised.

## Persistence and completeness

Detection schema is now `rxn_core.fragment_detection/v4`. Older implicit-
augmentation records are not silently promoted to the new semantics: rerun
detection into a new directory. The precursor inventory itself is unchanged.
The bank scanner saves all search records by default, including capped and
no-match searches. Completeness also checks shard-wide diagnostics, so a dropped
or filtered record cannot silently turn an incomplete scan into a complete one.

The merge writes, in order:

- `*.occupations.json`: the lossless index and source shard locations.
- `*.decisions.json`: the coverage decision graph.
- `*.assemblies.jsonl`: every evaluated complete assembly in rank order.
- `*.checkpoint.json`: phase/progress, including unfinished status.
- The report and standalone assembled HTML.

Ground truth is checked against this same ranked stream. There is no separate
progressive rematching or injected answer. Absence from a certified prefix does
not prove non-assemblability; absence of a saved detection does not prove bank
absence.

`recommendations_certified` and `assembly_complete` are distinct: the former can
be true without exhausting all worse assemblies. `detection_complete` separately
reports incomplete AAM evidence. Certification is conditional on that evidence,
not a claim of globally optimal chemistry or exhaustive unconstrained AAM.

AAM retains its default branch cap of 100 with explicit cap records. Detection
has no default candidate cap. Explicit caller-supplied matching constraints or
budgets remain possible in the reusable matching API; the recommendation CLI
adds none. The C++ growth/search code is unchanged. The only core representation
adjustment extends generator frames when relabeling augmentation atom indices.

Reproducible saved smoke experiment:

```bash
.venv/bin/python bench/principled_retro.py \
  --output-dir data/retro_runs/principled_retro_smoke
```

It uses a two-entry **test bank**, not a blind inventory scan. Re-running that
command replays its saved detections without repeating AAM.
