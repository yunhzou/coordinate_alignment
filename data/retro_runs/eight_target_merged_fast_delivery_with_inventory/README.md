# Merged fast-delivery and inventory scan

This run searches eight targets against `merged_fast_delivery_with_inventory`,
a stereochemically deduplicated bank containing 155,303 structures from:

- Mcule fast-delivery class 1
- Mcule augmented fast-delivery class 2
- the 1,919-row local inventory

Each target directory contains the ranked `results.json`, a standalone
symmetry-enabled `viewer.html`, and `postprocess_time.tsv`. The 352 MB of raw
fragment-detection shards is retained locally but excluded from Git.

| Target | Viewer | Result |
|---|---|---|
| Fluorinated phosphorus ligand | [viewer](t01_fluorinated_phosphorus_ligand/viewer.html) | [JSON](t01_fluorinated_phosphorus_ligand/results.json) |
| Chloro quinoline alcohol | [viewer](t02_chloro_quinoline_alcohol/viewer.html) | [JSON](t02_chloro_quinoline_alcohol/results.json) |
| Aryl magnesium bromide | [viewer](t03_aryl_magnesium_bromide/viewer.html) | [JSON](t03_aryl_magnesium_bromide/results.json) |
| Dimethyl azobenzene | [viewer](t04_dimethyl_azobenzene/viewer.html) | [JSON](t04_dimethyl_azobenzene/results.json) |
| Cyano thiazole amidine | [viewer](t05_cyano_thiazole_amidine/viewer.html) | [JSON](t05_cyano_thiazole_amidine/results.json) |
| Phenyl cyclooctadiene | [viewer](t06_phenyl_cyclooctadiene/viewer.html) | [JSON](t06_phenyl_cyclooctadiene/results.json) |
| Acetophenone | [viewer](t07_acetophenone/viewer.html) | [JSON](t07_acetophenone/results.json) |
| Cyclohexene | [viewer](t08_cyclohexene/viewer.html) | [JSON](t08_cyclohexene/results.json) |
