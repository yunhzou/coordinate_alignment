---
name: pathway-discovery
description: Product-blind reaction pathway discovery from XYZ intermediates using xTB-derived descriptors, WBO changes, geometry, and AAM identity tracking. Use when proposing, ranking, explaining, or validating possible bond formation, bond cleavage, coupled bond-change events, native atom-index viewers, or reaction-path progress without hard-coded empirical mechanism templates.
---

# Pathway Discovery

## Principle

Derive proposals from measurable descriptors of the current structure. Do not
rank a proposed next step using known child structures, product labels, or
mechanism templates. Use downstream structures only for verification after the
proposal is made.

Keep charge, multiplicity, method, and atom-index frame explicit. For all xTB
descriptor jobs in one comparison, use the same charge/multiplicity convention.

## Data Sources

For each intermediate, compute or load:

- Native XYZ atom order and Cartesian coordinates.
- xTB WBO matrix.
- xTB partial charges.
- xTB condensed Fukui indices: `f_plus`, `f_minus`, and optionally `f_zero`.
- Optional molecule-level descriptors such as HOMO-LUMO gap or global softness.
- Optional AAM maps for verification or path identity tracking.

These are xTB descriptor measurements, not NBO measurements. Do not call them
NBO unless a separate NBO-capable backend is actually run and reported.

## Atom And Pair Metrics

Build explicit tables before producing any ranked call.

| scope | metric | definition |
|---|---|---|
| atom | native index | 1-based atom index in the current XYZ |
| atom | element | chemical element label |
| atom | `q` | xTB partial charge |
| atom | `f_plus` | condensed Fukui electron-acceptor response |
| atom | `f_minus` | condensed Fukui electron-donor response |
| atom | `f_zero` | condensed radical response, if available |
| atom | donor signal | configurable normalized signal derived from atom descriptors |
| atom | acceptor signal | configurable normalized signal derived from atom descriptors |
| pair | distance A | Cartesian distance between two atoms |
| pair | direct WBO | current WBO between the atom pair |
| pair | graph distance | shortest path length in the current WBO graph |
| pair | pair class | configurable class such as organic, metal-substrate, metal-ligand, ligand-framework |
| pair | stretch ratio | distance divided by a chosen size/reference scale |
| pair | charge polarization | `abs(q_i - q_j)` |
| pair | frontier coupling | endpoint Fukui cross-response |
| path | mapped WBO delta | WBO change for an AAM-aligned pair between two structures |

Keep raw component columns in CSV outputs. A final score without components is
not auditable enough for pathway discovery.

## Formation Measurement

Formation measurement asks which absent or weak direct contacts are geometrically
and electronically compatible with bond/contact formation.

Use an adjustable score of the form:

```text
formation_score(i -> j)
  = combine(
      donor_signal(i),
      acceptor_signal(j),
      distance_factor(i,j),
      direct_wbo_absence_factor(i,j),
      optional_terms(i,j)
    )
```

A common multiplicative baseline is:

```text
formation_score(i -> j)
  = donor_signal(i)
  * acceptor_signal(j)
  * distance_factor(i,j)
  * direct_wbo_absence_factor(i,j)
  * optional_terms(i,j)
```

Suggested measurable components:

```text
donor_signal(i)
  = normalize(f_minus(i), q(i), local_softness(i), optional_atom_terms(i))

acceptor_signal(j)
  = normalize(f_plus(j), q(j), local_softness(j), optional_atom_terms(j))

distance_factor(i,j)
  = configurable decreasing function of distance from a chosen contact scale

direct_wbo_absence_factor(i,j)
  = configurable penalty based on current direct WBO(i,j)
```

Record graph distance, but keep its treatment configurable. It is a topology
measurement, not a descriptor of direct WBO.

## Cleavage Measurement

Cleavage measurement asks which existing bonds or contacts are compatible with
weakening or breaking.

Use an adjustable class-aware score of the form:

```text
cleavage_score(i-j)
  = class_normalize(
      combine(
        weak_wbo_component(i,j),
        stretch_component(i,j),
        charge_polarization_component(i,j),
        frontier_component(i,j),
        local_softness_component(i,j),
        optional_terms(i,j)
      )
    )
```

Suggested measurable components:

```text
weak_wbo_component(i,j)
  = WBO weakness measured relative to comparable pair classes

stretch_component(i,j)
  = distance(i,j) / chosen_size_scale(i,j)

charge_polarization_component(i,j)
  = abs(q_i - q_j)

frontier_component(i,j)
  = f_plus(i) * f_minus(j) + f_plus(j) * f_minus(i)

local_softness_component(i,j)
  = f_plus(i) + f_minus(i) + f_plus(j) + f_minus(j)
```

Do not compare cleavage values across heterogeneous pair classes unless they
have been normalized onto a comparable scale.

## Coupled Bond-Change Measurement

Coupled measurement ranks a proposed formation with a proposed cleavage or
weakening event.

Use an adjustable score of the form:

```text
coupled_score(form i->j, cleave k-l)
  = combine(
      formation_score(i->j),
      cleavage_score(k-l),
      coupling_factor(i,j,k,l),
      optional_terms(i,j,k,l)
    )
```

Suggested measurable coupling terms:

```text
shared_endpoint
  = true if {i,j} intersects {k,l}

prox_A
  = minimum Cartesian distance between endpoints of the formation pair
    and endpoints of the cleavage pair

coupling_factor
  = configurable function of shared_endpoint and prox_A
```

Report formation rank, cleavage rank, and coupled rank separately.

## Atom Alignment

Always state the atom-index frame used in a table or sentence.

- Native frame: atom labels from the current XYZ file.
- Reference frame: atom labels after AAM mapping into a chosen structure.
- Local adjacent frame: atom labels from the two structures being directly
  compared.

When deriving pathway order or progress, choose one reference atom frame and use
it consistently. Usually use the head/reactant as the forward reference. For
product-side auditing, use the tail/product as the backward reference.

For multi-step paths:

1. Run local AAM between adjacent structures.
2. Compose adjacent mappings into the chosen reference frame.
3. Re-express WBO deltas and event counts in that reference frame.
4. Use direct head/tail AAM only as an additional diagnostic.

Never compare native indices from different XYZ files as if they are the same
atom identity.

## Verification Measurements

After a proposal is made, verify with mapped WBO changes when comparison
structures are available.

For each mapped pair, report:

| field | meaning |
|---|---|
| pair in reference frame | aligned atom pair in the chosen reference |
| pair in comparison frame | corresponding native or mapped pair |
| WBO before | WBO in earlier/current structure |
| WBO after | WBO in comparison structure |
| delta WBO | `WBO_after - WBO_before` |
| event class | formed/strengthened, broken/weakened, or below threshold |
| pair class | same class scheme used for scoring |

Event thresholds are configurable and should be reported with the table.

## Output Requirements

For each investigated intermediate, write auditable artifacts:

- `atom_stats.csv`: native index, element, charge, Fukui values, donor signal,
  acceptor signal.
- `formation_scores.csv`: raw formation candidates with component columns.
- `cleavage_scores.csv` or `cleavage_scores_by_class.csv`: cleavage candidates
  with component columns and pair classes.
- `coupled_move_scores.csv`: paired formation/cleavage candidates with component
  columns.
- Markdown summary with concise tables and no hidden filters.
- Optional heatmap for matrix-style scores.
- Optional 3Dmol viewer with native 1-based atom labels and any requested guide
  pairs.

If a score or rule is experimental, label it as such and keep the raw
measurements visible.
