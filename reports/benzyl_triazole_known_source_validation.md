# Known-source full-coverage validation

## Accepted example

The validation target is 1-benzyl-4-phenyl-1H-1,2,3-triazole.  The documented
three-component synthesis uses benzyl bromide, sodium azide, and
phenylacetylene.  Copper is catalytic and does not contribute product atoms.

## Blind result

The detector received only the final product and the three-source bank.  It
used explicit hydrogens, tolerance 0.5, and branch cap 100.

- all 3 sources detected
- 4 alternative fragment placements generated
- 0 cap hits
- 31 of 31 explicit product atoms occupied
- 18 of 18 heavy product atoms occupied
- 0 overlapping target occupations
- 0 chirality conflicts

The selected mappings agree with the documented source roles:

- benzyl bromide supplies the complete benzyl fragment; bromide is unmatched
- sodium azide supplies the three triazole nitrogens; sodium is unmatched
- phenylacetylene supplies the phenyl group and both triazole carbons; the
  alkyne bond changes during cycloaddition

This is a valid full-source occupation example.  Unlike the retired
vancomycin test, no target atom is attributed to an omitted atom-donating
reagent.

## Evidence

- primary reaction example: https://pubs.acs.org/doi/10.1021/acscatal.1c05610
- independent one-pot synthesis and 98% reported yield:
  https://pubs.rsc.org/en/content/articlehtml/2016/ra/c5ra25116h
- product identity, PubChem CID 11310894:
  https://pubchem.ncbi.nlm.nih.gov/compound/11310894

## Saved artifacts

- source bank: `docs/example_runs/benzyl_triazole_source_bank.csv`
- target: `docs/example_runs/benzyl_triazole_target.smi`
- blind records: `data/retro_runs/benzyl_triazole/blind_fragment_records.jsonl.gz`
- full-cover assembly: `data/retro_runs/benzyl_triazole/assembly_coverage.json`
- assembled viewer: `reports/benzyl_triazole_assembled_view.html`
