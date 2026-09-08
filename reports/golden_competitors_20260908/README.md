# Golden: competitor reproduction

All 1,851 original Golden records were submitted to each configuration. These
are **fresh runs of released implementations with default mapping settings**, not
copied literature numbers. Evaluation uses our common, strict mapping criterion;
this is not a claim to reproduce each paper's historical evaluation protocol.

| Implementation | First returned mapping correct | Reference in any returned mapping | Mapping exceptions |
|---|---:|---:|---:|
| RXNMapper 0.4.3 | 1,547 / 1,851 (83.58%) | 83.58% | 0 |
| LocalMapper 0.1.5 | 1,605 / 1,851 (86.71%) | 86.71% | 0 |
| Chython 2.18 / chython-rxnmap 2.0 | 1,561 / 1,851 (84.33%) | 84.33% | 5 |
| SLAPMapper 1.0.0, binary (default) | 1,399 / 1,851 (75.58%) | 1,591 / 1,851 (85.95%) | 0 |
| SLAPMapper 1.0.0, weighted | 1,401 / 1,851 (75.69%) | 1,564 / 1,851 (84.49%) | 0 |
| Indigo 1.46.0 | 692 / 1,851 (37.39%) | 37.39% | 133 |
| RDT 4.0.0 | 1,030 / 1,851 (55.65%) | 55.65% | 0 after retry |

SLAP's first output is **not asserted to be a separately ranked top-one**. Its
default run returned 2,860 candidates; the weighted run returned 2,379. We did not
rerank, reference-guide, or repair any competitor's output. Chython's current
ONNX model is identified explicitly; it is not presented as the original
GraphormerMapper paper checkpoint.

## Evaluation

- Inputs come from `golden_original_20260906`, upstream commit
  `793475e54d8b2c7f714165a61e4eb439435d7d92`. All reference atom-map labels are
  removed before mapping. No balancing or byproduct completion is applied.
- Correctness is the exact shared **heavy-atom mapping relation**, including
  unmatched atoms, modulo chemical automorphisms of both endpoints. Endpoint
  colors include element, charge, isotope, hydrogen count and stereochemistry;
  bonds include order and stereochemistry. This is not merely agreement on a
  bond-change count. Explicit reference hydrogen assignments are not evaluated.
- We reuse the existing AAM evaluation's graph construction. Because external
  tools can modify molecular structures, chemical color names/counts accompany
  Nauty's certificate. This prevents changed chemistry from appearing equivalent
  just because two color partitions have the same topology.
- Original reference and label-free input endpoints agree on all 1,851 records.
  Agent-field components in returned `R>agents>P` SMILES are retained on the
  input side for scoring: Golden originally supplied them there. No atoms are
  deleted and mappings are not changed by that representation normalization.
- Duplicate heavy-atom labels are invalid, not silently overwritten. Invalid
  predictions and execution failures remain in the 1,851 denominator.
- Seventeen existing/new evaluator tests plus the agent-field test pass (18 in
  total). No core AAM implementation or competitor mapping algorithm was edited.

## Important qualifications

**Indigo needs careful interpretation.** Its outputs contain 503 predictions
with duplicate heavy-atom map labels, seven with changed endpoint chemistry, and
133 mapping exceptions. For case 0, the duplicate nitrogen map label is already
present in Indigo's native `atomMappingNumber` API, not introduced by our SMILES
parser. The same case reproduces the duplicate with Indigo 1.29.0. Therefore this
low score is not a general claim about all historical Indigo benchmarks.

Chython produced 29 endpoint-modifying outputs. Those are not accepted as valid
mappings of the original input under this strict protocol.

RDT's first attempt had one JVM failure followed by broken-pipe failures in the
same worker. All 232 affected attempts are preserved. We restarted workers after
an error and retried those records; all 1,851 final records have mapped output.
Three final RDT candidates fail output validation. Final successful-attempt
timings exclude the failed cascade. RDT's own published/repository
chemical-change scoring is a different criterion from exact mapping recovery.

SAMMNet is **not benchmarked**: the linked public repository contains training
and model code but no released pretrained checkpoint or release assets. We did
not train a substitute and label it a reproduction.

## Timing and resources

| Configuration | Sum of successful per-reaction mapping seconds | Actual mapping CPU-hours |
|---|---:|---:|
| RXNMapper | 104.55 | 0.02888 |
| LocalMapper | 544.91 | 0.15098 |
| Chython | 148.52 | 0.31371 |
| SLAP binary | 138.87 | 0.03856 |
| SLAP weighted | 107.22 | 0.02975 |
| Indigo | 315.86 | 0.08772 |
| RDT | 545.27 | 0.52477 |

These sums are **not parallel campaign wall times**. Queue/node setup, imports,
model loading, external serialization, artifact saving, and evaluation are
excluded. Tool-native output generation remains part of the mapping call.
Failure costs remain available separately in raw logs and accounting. RDT CPU
time includes its JVM, not just the Python wrapper. First-call model/JIT work
performed inside mapping is included.

The main campaign used 16 one-compute-thread shards per RXNMapper, LocalMapper,
SLAP variant and Indigo; eight eight-thread shards for Chython's default ONNX
configuration; and eight eight-CPU RDT shards. RDT's JVM uses
`-XX:ActiveProcessorCount=8 -Xmx4g`. Slurm can allocate two logical CPUs for a
one-thread request. Hardware, allocations, startup timings, and job accounting
are saved. Setup-stalled shards were relocated without rerunning completed
records. Per-reaction process watchdog: 300 seconds.

No hardware-normalized speedup or SOTA claim follows from this table. Our larger
multi-seed/bidirectional search budget must be compared separately from these
default competitor runs, and single-answer accuracy must not be equated with
alternative recovery.

## Saved artifacts and reproduction

Complete cluster run:

`/project/yunhengzou/coordinate_alignment/aam_benchmarks/golden_competitors_full_20260908`

This report contains `summary.json` and `saved_outputs.tar.gz`. The archive
includes every returned mapping, SLAP compressed labeled graphs, per-record
evaluations, startup/version records, model/template hashes, driver snapshots,
submission commands, logs, accounting, and preserved failed attempts. Model
weights themselves remain in the pinned cluster environments rather than Git.

Entry points:

```bash
python bench/golden_competitors.py init --run RUN
python bench/golden_competitors.py worker --run RUN --method rxnmapper
PYTHONPATH=src:bench .venv/bin/python bench/golden_competitors.py report --run RUN
```

Use each recorded competitor environment for `worker`, and our existing AAM
environment for `report`. Run `environment --run RUN` in each competitor
environment to record package/model provenance. The archive's exact submission
commands specify all shards and environments; completed record files are reused.

SLAP source pin: `ea248fd9494f52f4865193e87a98cc92c62b5f9e`.
RDT source: tag `v4.0.0`, commit `0d6632251108e7a8f625d15c60607a136093c4d8`.
RDT was built with Java 25 and Maven 3.9.11 using `compile
dependency:build-classpath`; its public `RDT.map` API is called by
`bench/RDTGolden.java`. The packaged-jar profile references an unavailable
assembly-plugin version, so no fat-jar packaging was used or needed.

Primary implementation sources:
[RXNMapper](https://github.com/rxn4chemistry/rxnmapper),
[LocalMapper](https://github.com/snu-micc/LocalMapper),
[Chython](https://github.com/chython/chython),
[SLAPMapper](https://github.com/shin1koda/slap-mapper),
[Indigo](https://github.com/epam/Indigo),
[RDT](https://github.com/asad/ReactionDecoder),
[SAMMNet](https://github.com/maryamastero/SAMMNet).
