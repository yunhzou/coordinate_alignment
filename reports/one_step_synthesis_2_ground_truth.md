# Ground-truth check: `one step synthesis_2.tif`

This check compares the reactants drawn in the supplied image with the
1,919-compound inventory, the 155,305-structure merged fast-delivery bank, the
saved R-to-P detection records, and the recommendation outputs. The merged
bank contains the original 155,303 structures plus two clearly labeled raw
material supplements. Every result combines the original full-bank scan with
the supplement scan and uses the same occupation-union workflow.

“Detected” means that the saved rough R-to-P scan retained at least one
fragment candidate. It does not mean that the reactant was selected by the
assembly ranking.

| Reaction | Depicted compound | 1,919 inventory | Merged bank | Saved detection | Recommendation result |
|---|---|---:|---:|---:|---:|
| 1 | Cyclohexylphosphine | No | Yes (supplement) | Yes | Exact set evaluates, not blindly ranked |
| 1 | Formaldehyde / paraformaldehyde unit | Yes | Yes | Yes | No |
| 1 | Ethanol | Yes | Yes | Yes | No |
| 1 | Bis(hydroxymethyl)cyclohexylphosphine intermediate | No | No | — | — |
| 1 | 4-(Trifluoromethyl)aniline | No | Yes | Yes | No |
| 2 | Iodoquinoline | No | Yes | Yes | No |
| 2 | 4-Chlorobenzaldehyde | Yes | Yes | Yes | No |
| 3 | 1-Bromo-4-tert-butylbenzene | No | Yes | Yes | No |
| 3 | Magnesium | Yes | Yes | Yes | Rank 5 |
| 4 | p-Toluidine | Yes | Yes | Yes | Blind candidate rank 3 |
| 5 | Acetylacetone | Yes | Yes | Yes | No |
| 5 | 2,6-Diisopropylaniline | Yes | Yes | Yes | Ranks 1–4 and 9–12 |
| 5 | Enaminone intermediate | No | No | — | — |
| 5 | 2-Amino-6-cyanobenzothiazole | No | Yes | Yes | No |
| 6 | 1,5-Cyclooctadiene | No | Yes (supplement) | Yes | Blind candidate rank 150 |
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

- Reaction 1 can be assembled directly from the raw atom sources without its
  intermediate: cyclohexylphosphine x2, formaldehyde x4, and aniline x2. This
  exact eight-copy set is not returned by the blind geometric ranking because
  lower-complexity covers rank ahead of it.
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
  returns it at blind candidate rank 3.
- Reaction 5 has all depicted raw ingredients in the bank and detects them,
  but their saved occupations do not yet form a complete final-target cover.
- Reaction 6 is recovered directly from 1,5-cyclooctadiene and iodobenzene
  after adding the absent raw material to the supplement. It is blind
  candidate rank 150, outside the ordinary top-20 display but explicitly
  appended and labeled in the viewer.
- Reaction 7 has all depicted raw ingredients in the bank and detects them.
  Phenol plus n-butyl vinyl ether do not yet form a complete final-target
  cover; triflic anhydride is an activating reagent, not an atom-source module.
- Reaction 8 has all depicted raw ingredients in the bank and detects them.
  Cyclohexanone is the final-product atom source, but its saved occupation does
  not yet retain every explicit final-target atom.

## Corrected full-bank assembly run

All searches used the 155,305-structure merged bank and explicit hydrogens.
The recommendation budget is target-scaled and every result reports that its
search was truncated; these are recommendations, not exhaustive rankings.

| Reaction | Post-processing time | Blind complete patterns | Displayed assemblies | Known immediate set |
|---|---:|---:|---:|---|
| 1 | 158.5 s | 57 | 20 | Exact set evaluates; not returned by blind ranking |
| 2 | 47.9 s | 31 | 20 | Not full: 30/31 explicit atoms |
| 3 | 61.0 s | 25 | 20 | Not full: saved Br residual remains uncovered |
| 4 | 73.2 s | 26 | 20 | Blind rank 3: p-toluidine x2, 30/30, zero overlap |
| 5 | 65.4 s | 58 | 20 | Raw ingredients present; complete cover not recovered |
| 6 | 73.9 s | 30 | 21 | Blind rank 150; appended as labeled ground truth |
| 7 | 39.2 s | 20 | 20 | Raw ingredients present; complete cover not recovered |
| 8 | 85.2 s | 20 | 20 | Raw ingredients present; complete cover not recovered |

## Main finding

After adding cyclohexylphosphine and 1,5-cyclooctadiene, every depicted raw
ingredient is present in the merged bank and detected. Intermediates are not
required: raw atom sources are matched directly to the final product. Exact
raw-material assemblies exist for reactions 1, 4, and 6; the blind recommender
returns reactions 4 and 6, while reaction 1 is only recovered by independent
evaluation of its known raw set. Reactions 2,
3, 5, 7, and 8 expose mapping or explicit-atom ownership gaps despite their
raw ingredients being present; these are algorithm gaps, not inventory gaps.
Assembly is a set-union operation over whole AAM occupation regions, including
repeated copies of one precursor. Overlaps are allowed, but less-overlapping
complete recommendations rank ahead of otherwise equivalent alternatives.
