# Vancomycin literature-bank cross-validation

## Answer

The earlier workflow did not recover the published Boger route because six of
its seven protected starting subunits were absent from the merged MCule bank.
Only subunit 17, protected N-methyl-D-leucine, was present as the correct
stereoisomer (`MCULE-7118794312`).  The opposite enantiomer
(`MCULE-2876957375`) was also present.

The six absent structures were added to a separate literature-validation bank.
They were not inserted into the commercial inventory because a literature
starting material is not automatically a purchasable catalog item.

## Blind detection result

After repairing the augmented-fragment baseline, one unguided explicit-H pass
recovered all seven literature subunits:

- target: vancomycin aglycon, PubChem CID 445835
- tolerance: 0.5
- branch cap: 100
- recovered: 7/7 subunits, 16 candidate mappings
- cap hits: 0
- elapsed time: 72.52 seconds with seven workers

Before the repair, subunits 14 and 15 produced valid initial fragments but no
final candidate.  The augmented matcher found target-heavy placements for the
residual protecting groups and discarded them at the cut-boundary check.  It
failed to retain the always-valid placement in which every residual atom maps
to its copied competitor fragment.  The repair adds that canonical copied
placement to every augmented search.  Subunits 14 and 15 then recover without
a target-region hint.

## Assembly interpretation

The seven starting subunits collectively map 107 of 132 explicit-H target
atoms, including 64 of 80 heavy atoms.  They do not directly cover the complete
aglycon, which is expected: the cited synthesis is a 17-step route, not a
single-step assembly.  Oxidation, nitrile conversion, aromatic substitution,
chlorination, macrocyclization, deprotection, and other reagents change or add
atoms between these starting subunits and the final aglycon.

Therefore the cross-validation establishes that the fragment detector can now
discover every published starting subunit once it is present.  It does not
claim that the seven structures alone form a balanced one-step reaction.

## Saved artifacts

- validation bank: `docs/example_runs/vancomycin_literature_validation_bank.csv`
- target: `docs/example_runs/vancomycin_aglycon_pubchem_445835.smi`
- complete detection records: `data/retro_runs/vancomycin_literature_validation/blind_fixed_fragment_records.jsonl.gz`
- detection summary: `data/retro_runs/vancomycin_literature_validation/blind_fixed_fragment_records.jsonl.gz.summary.json`
- all-subunit viewer: `reports/vancomycin_literature_validation_results.html`

Primary source: M. J. Moore et al., *J. Am. Chem. Soc.* 2020, 142,
16039-16050, DOI: 10.1021/jacs.0c07433.  Structures 11-17 are transcribed from
Figure 2 and checked against the supporting information where full compound
names are supplied.
