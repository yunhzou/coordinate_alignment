# Exact blind inventory scan: bulky BIAN target

Target: `CC(C)C1=CC=CC(C(C)C)=C1/N=C2/C(C3=C4C(C=CC=C42)=CC=C3)=N/C5=C(C(C)C)C=CC=C5C(C)C`

The complete merged inventory of 1,919 deduplicated compounds was scanned without supplying a known precursor ID. Detection used explicit hydrogens, all source seeds, full augmented-fragment matching, WBO tolerance 0.5, branch cap 100, and candidate cap 100.

- Slurm execution: 28 shards, 48 CPUs per shard
- Full blind-scan wall time: 480.84 seconds
- Compounds searched: 1,919
- Parse failures: 0
- Compounds with retained-fragment matches: 1,840
- Fragment candidates: 5,361
- Compounds reporting a cap event: 7
- Normal shards completed in roughly 24–118 seconds
- The final wall time was controlled by one 886-atom explicit-H peptide

The top modular assembly recovered the expected condensation foundation:

- 2 × `INVENTORY-000436`: 2,6-diisopropylaniline
- 1 × `INVENTORY-001463`: acenaphthenequinone
- Full target coverage
- Set atom retention: 0.9286
- Set heavy-atom retention: 0.9500
- Six broken precursor bonds and two formed product bonds under the geometric edit accounting

Complete compressed AAM detections are stored in `data/retro_runs/bulky_bian_inventory_exact_accelerated/parts/`. The modular assembly report is `data/retro_runs/bulky_bian_inventory_exact_accelerated/assembly.json`.
