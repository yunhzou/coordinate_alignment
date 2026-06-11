# Metric Terms And Current Rules

This folder summarizes the product-blind pathway investigations from the xTB/AAM
work. The goal is to preserve what was learned without hiding raw ranks behind
post-filters.

## Descriptor Terms

| term | meaning | use | caveat |
|---|---|---|---|
| Native index | 1-based atom index in the current XYZ file | Use for visual inspection and proposed moves from one structure | Native indices differ across structures |
| Mapped index | Atom identity transferred by AAM between structures | Use for R/P or adjacent-intermediate WBO changes | Always state the source and target frame |
| WBO | xTB Wiberg/Mayer-like bond order | Measures existing or changing bond/contact strength | Do not compare raw WBO globally across organic and metal bonds |
| `q` | xTB partial charge | Context for electron-rich/electron-poor sites | Charge alone is not a donor/acceptor score |
| `f_minus` | Condensed Fukui donor response | Larger values mark stronger donor-like centers | Compare within the same method/charge protocol |
| `f_plus` | Condensed Fukui acceptor response | Larger values mark stronger acceptor-like centers | Compare within the same method/charge protocol |
| Distance A | Cartesian interatomic distance in Angstrom | Tests whether a contact is geometrically accessible | A long distance can still be possible after conformational change |
| Graph distance | Number of existing WBO-graph edges between atoms | Diagnostic tag for topology/proximity | Do not hard-drop graph-distance-2 candidates |

## Score Terms

| score | built from | answers | interpretation |
|---|---|---|---|
| Raw formation score | donor strength, acceptor strength, distance factor, direct-WBO penalty | Which absent/weak contact wants to form? | Keep raw rank visible. It is not a rate constant. |
| Cleavage score | class-normalized WBO weakness, stretch, charge polarization, Fukui/local softness | Which existing bond/contact wants to weaken or break? | Rank within bond class: organic, metal-substrate, metal-ligand. |
| Coupled score | formation score + cleavage score + endpoint/proximity relation | Which formation-cleavage pair is coherent? | Useful for rearrangements, migrations, and collapses. |
| Prox A | minimum distance between endpoints of a formation and cleavage pair | Whether two events are locally connected | `0.00` often means the formation and cleavage share an atom. |

## Current Formation Rule

Use direct WBO as the hard formation gate:

| condition | treatment |
|---|---|
| Direct WBO substantial | Not a new formation; classify as strengthening or cleavage candidate |
| Direct WBO near zero | Eligible formation candidate |
| Graph distance 1 | Already bonded; exclude from formation |
| Graph distance 2 | Retain and tag as proximal/collapse/cyclization-like |
| Graph distance >= 3 | Retain and tag as remote/new-contact-like |

The key correction from post-`12`: `O23->C1` was the correct next formation.
It had direct WBO `0.000`, distance `2.43 A`, raw rank `#2`, and graph distance
`2`. The earlier hard `graph_distance > 2` filter was too aggressive.

## AAM Rule

Use AAM for verification and identity tracking, not for cheating during the
proposal step.

When deriving pathway order or progress, choose one reference atom frame and use
it consistently. Usually use the head/reactant as the forward reference and the
tail/product as the backward reference. Adjacent local AAM events should be
composed into that reference frame before comparing WBO deltas or event counts.

| use case | preferred method |
|---|---|
| Verify one proposed transition | Local AAM between the two structures |
| Track a multi-step path | Adjacent AAM, then compose mappings |
| Compare to head/tail | Use composed maps; direct endpoint AAM is diagnostic only |
| Detect suspicious path steps | Watch for too many backward WBO events; tolerate 1-2 |

All summaries below use product-blind scoring unless an AAM verification table
is explicitly labeled as downstream verification.
