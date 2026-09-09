# Concrete AAM / SLAP comparison on real TS inputs

**Result: no best-event-score advantage for AAM in these six tests. There are
additional retained assignments, and continuous weights sometimes greatly
reduce capped searches, but also lose one historical alternative.** This is an
exploratory comparison on three real cached structures, not curated TS accuracy.

Open [the offline viewer](viewer.html). It uses actual saved full-atom mappings,
original XYZ coordinates and WBO matrices. Switch the two method buttons over
the same source and TS geometries; enable atom labels or click an event/atom.
Identity colors represent source atom identities, not matched fragments.
The first selection exposes the concrete oxygen-role contrast described below.

Full viewer path:

```text
/h/399/yunhengzou/coordinate_alignment/reports/real_ts_comparison_20260909/viewer.html
```

## 1. Best bond-event counts, including every explicit hydrogen

Every method is scored on the **same original continuous WBO matrices**: an edge
exists above 0.2; a retained edge changes order when its WBO difference exceeds
0.5. Broken, formed and order-changed edges each count once. Scores concern
endpoint-to-TS structural differences, not activation energies or complete
reaction mechanisms.

| Real input | Atoms | AAM continuous | AAM binary | SLAP binary | SLAP raw WBO | SLAP float-cost diagnostic |
|---|---:|---:|---:|---:|---:|---:|
| `pr1.tempo_ts1`, R → TS | 26 | 1 | 1 | 1 | 1 | 7 |
| `pr1.tempo_ts1`, P → TS | 26 | 1 | 1 | 1 | 1 | 11 |
| `pr16.carbocation_ts5`, R → TS | 18 | 1 | 1 | 1 | 5 | 11 |
| `pr16.carbocation_ts5`, P → TS | 18 | 4 | 4 | 4 | 8 | 13 |
| `pr7.V.dodh_ts910`, R → TS | 28 | 4 | 4 | 4 | 4 | 5 |
| `pr7.V.dodh_ts910`, P → TS | 28 | 3 | 3 | 3 | 3 | 3 |

The minimum is certified over each method's **saved family union** in all 30
comparisons. This is not a global optimality certificate over unseen search
branches. Hydrogen assignments are optimized inside the retained correlated
families; the table does not penalize an arbitrary displayed H witness.

SLAP binary and AAM tie in all six rows. Supplying raw WBOs to this pinned SLAP
implementation does not improve them. Its LAP matrix uses integer costs; the
last column replaces only that matrix with floating costs. All other upstream
search, tie and perturbation behavior remains unchanged. This diagnostic is
**not a validated float-native SLAP redesign**, and its worse results do not
establish an inherent limitation of LAP-based mapping.

## 2. A concrete additional assignment, not just a richer container

For vanadium **R → TS**, source oxygen atoms 12 and 15 have these alternatives
in the saved outputs (indices are zero-based and unchanged between methods):

| Oxygen roles | AAM continuous | SLAP binary |
|---|---|---|
| O12 → TS17; O15 → TS9 | Retained | Retained |
| O12 → TS9; O15 → TS17 | Retained | Not retained in these ten runs |

The second row corresponds to historical heavy core #2, including C11 → TS16.
Its absence from SLAP is checked modulo common score-preserving endpoint
symmetries, not by mistaking a symmetric relabeling for a new assignment.
The first viewer selection shows an actual AAM witness with this oxygen-role
swap (5 full-H events), versus a SLAP binary minimum witness (4 events).

Important qualifications:

- These historical core assignments came from previous algorithm outputs,
  **not independently curated chemical ground truth**. Additional output is not
  automatically a correct additional mechanism.
- For **P → TS**, SLAP retains both historical heavy-core alternatives. Thus
  this R → TS difference does not prove that only AAM's complete two-endpoint
  workflow can discover the second mechanism.
- Literal full witnesses may differ in peripheral atoms/H; the displayed extra
  AAM witness is illustrative, not claimed optimal subject to that core.

Across the six pairs we observed **59** continuous-AAM versus **18** binary-SLAP
heavy mapping classes, modulo the same score-preserving endpoint symmetries.
AAM's number is a positive witness lower bound, not enumeration of its complete
compressed space. SLAP heavy labels were singleton and every returned native
candidate was retained. Of the AAM witness classes, 42 were absent from the
saved SLAP binary outputs; 17 were shared. One SLAP class absent from AAM's
terminal witnesses has not been excluded from AAM's full compressed families.

Most observed extra AAM alternatives have higher event cost. Among witnesses
obtained by the all-H family minimization, we found **no additional minimum-score
AAM heavy class absent from SLAP binary**. This is not exhaustive exclusion of
all tied mappings inside the families.

## 3. Do continuous edges actually help?

The controlled comparison is **the same AAM with continuous versus binary
inputs**, using the same seed policy, cap and source cut list.

| Measurement | AAM continuous | AAM binary |
|---|---:|---:|
| Vanadium P → TS, capped cut contexts | **0 / 33** | **33 / 33** |
| Vanadium P → TS, best full-H events | 3 | 3 |
| Vanadium P → TS, historical heavy cores retained | 2 / 2 | 2 / 2 |
| Vanadium R → TS, capped cut contexts | 0 / 33 | 4 / 33 |
| Vanadium R → TS, observed heavy classes | 7 | 27 |
| Vanadium R → TS, mapping CPU seconds | 0.735 | 1.358 |
| Carbocation R → TS, historical heavy core #2 | **Absent** | **Present** |

This demonstrates useful branch discrimination on the vanadium inputs, **not a
universal accuracy improvement**. Continuous carbocation R → TS had 5 capped
cut contexts, versus zero for binary; this experiment does not isolate whether
its missing alternative is caused by the cap, growth ordering or weight
constraints. It is absent from the saved correlated families, not merely from
their sampled representatives.

Historical-core query totals: continuous AAM **7/8**, binary AAM **8/8**, binary
SLAP **5/8**, raw-WBO SLAP **3/8**, float-cost diagnostic **3/8**. These are eight
queries across four endpoint/TS pairs from two reactions, **not eight independent
reactions or an accuracy percentage**. The TEMPO case has no historical core
reference. See [the symmetry-normalized queries](core_equivalence.json).

## 4. Actual CPU cost, with scoring separate

Totals over all six input pairs; each mapping task used one CPU.

| Configuration | Mapping compute, CPU s | Subsequent all-H family scoring, CPU s |
|---|---:|---:|
| AAM continuous | 4.064 | 35.974 |
| AAM binary | 4.824 | 46.863 |
| SLAP binary | 0.609 | 2.285 |
| SLAP raw WBO | 1.384 | 4.632 |
| SLAP float-cost diagnostic | 0.422 | 2.072 |

Here binary SLAP is about **6.7× cheaper in mapping compute** than continuous
AAM; including this family-minimization evaluation, about **13.8× cheaper**.
Different search policies retain different amounts of output. These are pilot
costs, not an equal-output-work benchmark or a general speed ratio.

Mapping compute includes preparation, search, symmetry finalization and witness
collection; import/input load and persistence are excluded. Scoring CPU includes
loading saved family files, encoding, solver checks and witness realization;
it excludes writing the final score JSON. Ground-reference queries and viewer
generation are separate evaluation work, not included in these mapping/scoring
columns. Per-phase CPU/wall records and Slurm accounting are saved.

Search tasks completed in 2–6 s including process startup/I/O. Mapping array
455891 ran on bosque7; two tasks that never started remained CONFIGURING, were
canceled after approximately 184 s and relocated to bosque11 as array 455928.
All 30 actual mapping tasks succeeded without mapper timeouts. The scheduling
wait is excluded from compute totals. Scoring array 455951 completed all 30
tasks on bosque11. This is not a strict same-node paired timing experiment.

## Protocol, saved evidence and reproduction

- Inputs: original cached R, P and reference TS XYZ/WBO files. No SMILES
  conversion, new electronic-structure computation or synthetic oracle.
- Six independent searches are R → TS and P → TS. They are **not** the prior
  R → P benchmark or reciprocal TS → endpoint search.
- AAM: production native fragment growth, 10 seed orders per cut, root seed 42,
  cap 100, tolerance 1.0, uncut plus every source-edge cut, explicit H.
- SLAP: pinned upstream core; ten input orders including the original, seeded
  by 20260909 plus pair index; heavy-atom symmetry breaking. Binary input is
  constructed before graph labeling, avoiding WBO leakage into binary labels.
  Raw-WBO and float-cost runs use the same raw weighted adjacency.
- Search never reads the historical mapping references. They are separate
  evaluator-only files. No core production AAM changes were made for this pilot.
- All-H scoring uses AAM's correlated path constraints. SLAP scoring keeps
  final label domains **and all returned LAP fingerprint-count correlations**;
  it does not treat H atoms as independently interchangeable. All 30 saved
  unions reached closed minimum bounds; no inconsistent native units.
- Hard mapping watchdog 300 s, Slurm allocation 10 min. Scoring has a 180 s
  watchdog, 120 s soft worker budget and 0.5 s soft per-family solve budget.
  Arrays used at most 30 concurrent single-CPU tasks.

[Summary and full minimum witnesses](summary.json), [class/cap/phase analysis](analysis.json),
[manifest](manifest.json), [original core references](references.json),
[Slurm accounting](slurm_accounting.psv), [saved native outputs](saved_outputs.tar.gz).
The archive contains inputs, every full native AAM cut graph and SLAP output
(including LAP internals), full-H scores, logs, pinned SLAP core and frozen
search/scoring drivers. The accompanying [SLAP license](SLAP_LICENSE) applies
to that third-party source. Git stores the current analysis/publishing scripts.

Full run, including frozen AAM source and native binary:

```text
/project/yunhengzou/coordinate_alignment/aam_benchmarks/real_ts_compare_20260909_EwGzgT/run
```

Search code pin: `848e44bf889c250d179a95b654387b7e4be16370`; driver, native and
input hashes are in the manifest/inputs. SLAP core SHA256:
`62964b12866e691698d594bfdc65f84ae714fb049df4267a5f5964b5c7e5fa3c`.
See [the pinned upstream core](https://github.com/shin1koda/slap-mapper/blob/ea248fd9494f52f4865193e87a98cc92c62b5f9e/src/slapmapper/core.py).

[Validation](validation.json): 64 related tests passed in 2.44 s, including new
checks that SLAP fingerprint correlations and their returned alternatives survive
the scoring adapter. An independent scalar all-pairs calculation verified 5,581
full-H witnesses and all 50 viewer records. Both inline scripts passed JavaScript
syntax checks. A live-browser visual inspection was not performed.

Use `bench/compare_real_ts_mappings.py --help` for preparation/submission and
saved-result query commands. `bench/score_real_ts_families.py` works only on saved
outputs; `bench/publish_real_ts_comparison.py` builds the viewer without mapping
reruns. Always choose a new directory for a new mapping experiment.

**Conclusion:** the demonstrated additional information is specific alternative
atom-role assignments and search-completeness/cap information. Continuous WBO
changes branch discrimination substantially in some TS inputs. Neither result
yet establishes better chemical TS accuracy; curated validation is still needed.
