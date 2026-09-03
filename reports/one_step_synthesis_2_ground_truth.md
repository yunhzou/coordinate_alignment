# Ground-truth check: `one step synthesis_2.tif`

This check compares the reactants drawn in the supplied image with the
1,919-compound inventory, the 155,303-structure merged fast-delivery bank, the
saved R-to-P detection records, and the recommendation outputs. No full-bank
AAM scan was rerun. All eight assembly/ranking jobs were rebuilt from those
saved records with the same occupation-union workflow.

“Detected” means that the saved rough R-to-P scan retained at least one
fragment candidate. It does not mean that the reactant was selected by the
assembly ranking.

| Reaction | Depicted compound | 1,919 inventory | Merged bank | Saved detection | Recommendation result |
|---|---|---:|---:|---:|---:|
| 1 | Cyclohexylphosphine | No | No | — | — |
| 1 | Formaldehyde / paraformaldehyde unit | Yes | Yes | Yes | No |
| 1 | Ethanol | Yes | Yes | Yes | No |
| 1 | Bis(hydroxymethyl)cyclohexylphosphine intermediate | No | No | — | — |
| 1 | 4-(Trifluoromethyl)aniline | No | Yes | Yes | No |
| 2 | Iodoquinoline | No | Yes | Yes | No |
| 2 | 4-Chlorobenzaldehyde | Yes | Yes | Yes | No |
| 3 | 1-Bromo-4-tert-butylbenzene | No | Yes | Yes | No |
| 3 | Magnesium | Yes | Yes | Yes | Rank 5 |
| 4 | p-Toluidine | Yes | Yes | Yes | Expected set recovered exactly |
| 5 | Acetylacetone | Yes | Yes | Yes | No |
| 5 | 2,6-Diisopropylaniline | Yes | Yes | Yes | Ranks 1–4 and 9–12 |
| 5 | Enaminone intermediate | No | No | — | — |
| 5 | 2-Amino-6-cyanobenzothiazole | No | Yes | Yes | No |
| 6 | 1,5-Cyclooctadiene | No | No | — | — |
| 6 | Iodobenzene | Yes | Yes | Yes | No |
| 7 | Phenol | Yes | Yes | Yes | No |
| 7 | Triflic anhydride | Yes | Yes | Yes | No |
| 7 | Phenyl triflate | No | Yes | Yes | No |
| 7 | n-Butyl vinyl ether | Yes | Yes | Yes | No |
| 7 | Enol ether intermediate | No | No | — | — |
| 8 | Cyclohexanone | Yes | Yes | Yes | No |
| 8 | Tosylhydrazide | Yes | Yes | Yes | No |
| 8 | Cyclohexanone tosylhydrazone | No | Yes | Yes | No |
| 8 | n-Butyllithium | Yes | Yes | Yes | No |

## Route-level conclusions

- Reaction 1 cannot be recovered exactly because the immediate phosphorus
  intermediate is absent. The detected aniline has two valid target
  occupation regions, consistent with its repeated use in the product.
- Reaction 2 has both depicted carbon-framework reactants in the merged bank,
  and both are detected. Together they cover 30 of 31 explicit product atoms;
  the remaining atom is the alcohol hydrogen supplied by the unspecified
  reduction chemistry. Direct evaluation confirms that the pair is therefore
  not a complete explicit-atom assembly.
- Reaction 3 has both reactants in the merged bank and detects both. The saved
  rough candidate covers 24 of 25 product atoms because product Br remains
  uncovered even though Br is present in the aryl bromide. This exposes a
  separate residual-fragment attachment question; it is not changed by the
  occupation-union assembly fix evaluated here.
- Reaction 4 is the direct regression for repeated precursor use. AAM gives
  p-toluidine two distinct 15-atom occupation regions. With occupation-aware
  union assembly, two copies cover all 30 explicit product atoms with no
  overlap. Direct evaluation recovers this exact set. The blind search also
  finds its construction pattern, although other molecules in that same
  geometric class are the displayed variants under the current retention
  ranking.
- Reaction 5 cannot recover the depicted final step because the enaminone
  intermediate is absent. Its other final reactant is present and detected.
- Reaction 6 cannot recover the depicted pair because 1,5-cyclooctadiene is
  absent. Iodobenzene is present and detected.
- Reaction 7 cannot recover the depicted final hydrolysis step because the
  enol ether intermediate is absent. Every depicted precursor from the two
  upstream steps is present and detected.
- Reaction 8 contains the hydrazone and n-butyllithium in the merged bank and
  detects both. The rough candidates do not form a strict complete atom cover.
  n-Butyllithium is a reagent rather than a product atom donor, and the
  hydrazone candidate should be refined before judging this route.

## Corrected full-bank assembly run

All searches used the 155,303-structure merged bank and explicit hydrogens.
The recommendation budget is target-scaled and every result reports that its
search was truncated; these are recommendations, not exhaustive rankings.

| Reaction | Post-processing time | Blind complete patterns | Displayed assemblies | Known immediate set |
|---|---:|---:|---:|---|
| 1 | 122.7 s | 0 | 0 | Cannot test exactly: immediate intermediate absent |
| 2 | 53.1 s | 20 | 16 | Not full: 30/31 explicit atoms |
| 3 | 66.2 s | 3 | 10 | Not full: saved Br residual remains uncovered |
| 4 | 42.9 s | 20 | 10 | Recovered: p-toluidine x2, 30/30, zero overlap |
| 5 | 89.4 s | 0 | 0 | Cannot test exactly: immediate intermediate absent |
| 6 | 79.7 s | 11 | 9 | Cannot test exactly: 1,5-cyclooctadiene absent |
| 7 | 36.0 s | 20 | 14 | Cannot test exactly: enol ether intermediate absent |
| 8 | 86.6 s | 20 | 14 | Not full: n-butyllithium is not a product atom donor |

## Main finding

Inventory lookup and R-to-P detection are working for every depicted compound
that is actually present. The poor ground-truth ranking comes from several
separate causes: missing intermediates in the bank, strict explicit-hydrogen
ownership when an unspecified reagent supplies a product hydrogen, the
residual-fragment issue exposed by reaction 3, and the assembly bug that
collapsed repeated target occupations of the same R. The assembly bug is
fixed by combining whole AAM occupation regions, including repeated copies of
one precursor. The reaction-3 residual question remains separate. Assembly is
a set-union operation: overlaps are allowed, but less-overlapping complete
recommendations rank ahead of otherwise equivalent alternatives.
