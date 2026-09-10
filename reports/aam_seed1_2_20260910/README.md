# One- and two-seed bidirectional AAM ablations

Status: preparing fresh, separately archived full-dataset runs.

Use the same frozen native-reuse AAM engine `98b01b1`, native binary, input
structures and search/profiling/evaluation adapters as the completed three-seed
and saved ten-seed benchmark. Only `seed_count` changes to one or two per cut.
These orders are prefixes of the same deterministic ten-order sequence.

Both directions, full sweep, cap 100, matching tolerance 1.0, explicit H,
root seed 42 and eight CPUs per directional AAM call remain unchanged.
Test all 1,851 Golden reactions and all 140 XYZ/WBO holdout cases separately.
The holdout has no annotated ground truth. Golden recovery uses the strict
heavy-atom relation modulo endpoint chemical symmetry, including unmatched atoms,
among saved compressed families; it is not top-1 accuracy.

Each experiment has a 512-CPU search budget, for 1,024 CPUs combined. Dynamic
queues distribute independent tasks. Search and evaluation are separate phases,
with 300-second search and 240-second evaluation process watchdogs. Save every
cut and final AAM archive. Keep incomplete outcomes in accuracy denominators.

Compare measured search CPU excluding persistence/loading on common completed
cases against three seeds, ten seeds and SLAP+sweep. Publish actual recovery
counts and case-set differences, not an assumption that fewer seeds preserve
accuracy. Cluster execution spans are distinct from CPU cost and exclude only
the initial queue, not persistence or dispatch.

Artifact roots:

- `/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_one_seed_bidirectional_20260910`
- `/project/yunhengzou/coordinate_alignment/aam_benchmarks/aam_two_seeds_bidirectional_20260910`

The three-seed run and ten-seed baseline remain untouched. Package defaults and
the AAM core are unchanged.
