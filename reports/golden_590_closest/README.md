# Case 590: reference-nearest saved representative

**Correction:** this older view constrained the seven carbons collectively, not
the reference's specific five-carbon subset. It does not answer the requested
five-carbon-choice comparison. See `../golden_590_correct_five/viewer.html`.
Checking both required sets separately returned zero eligible representatives
in either archive. This is not a negative verdict on every compressed realization.

Open `viewer.html` offline or at `/golden_590_closest/viewer.html` on port 8765.
Ground truth and the selected saved detection appear simultaneously, each with R and P.

This is **post-hoc diagnostic selection**, not a new search or the ranker's prediction.
Among both random and distance archives, require the seven-carbon reactant to occupy
the same seven product carbon positions collectively as the reference. Select the
representative with the fewest literal heavy-atom reference-pair disagreements;
break ties lexicographically. There were 769 eligible random and 645 eligible
distance terminal representatives. This does not optimize within their compressed
families or normalize endpoint automorphisms, so it is not a claim of the globally
closest symmetry-equivalent mapping.

The selected random terminal is 720938. All 31 heavy atoms and 55 explicit atoms
are mapped. It has 17 literal reference-pair differences, including symmetry
choices; the exact pairs and archive provenance are in `mapping.json`.

Using original RDF map labels (not the viewer's canonical r-indices):

- The detected five-carbon set is 25,26,27,28,29; remaining 30,31 enter product
  positions 28,27. Thus the five-plus-two product allocation is present.
- The three-ring reactant and its two oxygen assignments agree literally with
  the reference in this selected representative.
- Within the five-carbon side product, the detected atom ordering produces
  three bond-order changes, whereas the reference produces one.
- The ten atoms of the other, oxygen-containing reactant are assigned differently.
  In particular, reference C19/O20 supply carbon monoxide; this representative
  instead supplies it from C16/O17. This is not merely a flipped drawing.
- Reference heavy events: 6 breaks + 4 formations + 2 order changes = 12.
  Selected detection: 7 + 5 + 6 = 18. Hydrogen identities are not annotated in the
  reference, so these comparative counts exclude H.

This identifies disagreements in one saved representative, not the cause of
absence from every compressed family. Full reference recovery for 590 remains
unknown because the previous family verification timed out.

Validation: existing rendering tests passed (3); offline browser checked four
R/P diagrams, explicit-H toggle, linked atom selection, and absence of JS errors.
No core AAM modifications or AAM recomputation.
