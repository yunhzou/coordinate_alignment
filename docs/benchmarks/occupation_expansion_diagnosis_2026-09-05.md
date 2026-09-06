# Why inventory sources 538 and 594 exceed the watchdog

The retro consumer expands AAM's compressed symmetry into many explicit
competitor placements, then repeatedly traverses already-tested symmetry
transitions. The cost accumulates over many AAM paths. This investigation did
not change production matching, hydrogen treatment, seeds, or limits.

## Evidence

Two 100-second captures used the same explicit-H BIAN target, branch cap 100,
and no sweep as the full-bank scan. Each source had a 48-CPU allocation.
Instrumentation saved the exact native occupation-call inputs before executing
the unchanged function. These are bounded diagnostic captures, not completed
precursor results. Jobs: 432154_538 and 432154_594.

| Source | Calls captured | Calls returned before timeout | States emitted by returned calls |
| --- | ---: | ---: | ---: |
| INVENTORY-000538 | 1,355 | 1,337 | 706,476 |
| INVENTORY-000594 | 117 | 116 | 34,212 |

None of the captured native input tuples are exactly duplicated within a
source. Caching entire calls by their inputs would not fix these captures.
Those counts do not imply the full precursor contains only that many calls.

### 1. Competitor-only placements are materialized

`match_augmented_residuals` sends complete fragment images in **P plus the
copied residuals** into `materialize_target_coverage_orbit`. Its equivalence key
includes competitor images and bonds, not just the target. Only afterward does
`_augment_initial_family` project the placements onto the real P target.

| Isolated call | Complete occupation states | Target fragment-region relations |
| --- | ---: | ---: |
| 538, worker 2009879 call 32 | 5,460 | 12 |
| 594, worker 2148408 call 6 | 5,616 | 144 |
| 538, worker 2009869 call 39 | 2,640 | 12 |
| 594, worker 2148408 call 16 | 2,592 | 144 |

The last two calls were still in progress when capture stopped. Replaying the
saved inputs completed them; they are not individually ten-minute calls.
Their full output states were saved. Even keeping the **exact source-to-target
atom pairs**, not just target coverage, there are respectively only 12 and 144
different mappings. Each case has one retained-source-atom set. The remaining
distinctions concern competitors.

In the first 538 call, a transposition-connected pool contains 15 equivalent
competitor hydrogens, of which 12 are occupied. There are 455 subsets of size
12; the 5,460 complete states have 12 target-region relations. In the other
538 call, the analogous pool has 12 hydrogens with nine occupied: 220 subsets,
and 2,640 states for 12 exact target mappings. Explicit hydrogens are correct;
**eagerly listing their competitor-side choices is the representation cost**.

The native input archives for the representative calls are only a few KiB.
The consumer expands those compact actions into thousands of states and
subsequently builds per-state Python hierarchy/candidate objects.

### 2. The same state/generator transitions are revisited

The native walk closes successive stages over every state accumulated so far.
Many stages contain the same generators. It does not remember that an earlier
stage already evaluated a particular state/generator pair.

| Call | Generator entries / distinct generators | Operation attempts | Newly added states |
| --- | ---: | ---: | ---: |
| 538 call 32 | 548 / 64 | 1,372,706 | 5,459 |
| 594 call 6 | 727 / 78 | 1,324,817 | 5,615 |

Over 99.5% of attempts add nothing. A changed witness still requires rebuilding,
sorting, and comparing the full occupation key even when the occupation is
unchanged. For these calls those keys describe 171 and 162 source atoms,
respectively, including their fragment image sets and preserved bonds.

### 3. Branch cap 100 does not bound this work

The cap applies to AAM search branching. This later routine receives an
unlimited occupation limit (`-1`), and is called for multiple returned paths
and initial matching alternatives. It does not count those states against the
AAM branch cap. One CPU-expensive worker can also delay a whole augmented
family result while other workers become idle.

## Diagnostic cache experiment — not a production change

`bench/occupation_probe.cpp` replicates the existing native occupation walk,
with counters and an optional cache. For each generator it remembers how many
states were already evaluated in earlier stages. It skips only those old
state/generator pairs; newly discovered states and new generators still run.
No occupation is removed or sampled.

Isolated one-process replays (one CPU per call):

| Call | Instrumented baseline | Cached traversal | Output |
| --- | ---: | ---: | --- |
| 538 call 32 | 8.68 s | 2.81 s | Identical |
| 594 call 6 | 11.64 s | 4.05 s | Identical |
| 538 call 39 | 3.23 s | 1.21 s | Identical |
| 594 call 16 | 5.29 s | 1.84 s | Identical |

The production native routine, profiling baseline, and caching prototype have
exactly equal states, witnesses, permutation actions, and output order on all
four captures. A separate 280-case randomized comparison also agrees exactly.
The existing four native-group-operation tests pass.

This is roughly a threefold improvement of the sampled native occupation
calls, **not a measured full-precursor or full-bank speedup**. It does not remove
the number of output states or subsequent Python materialization costs.

## Fix direction

1. Reuse already evaluated state/generator transitions in the native walk.
   This is the smaller, output-preserving change demonstrated above.
2. Keep competitor-only choices symbolic when passing results toward target
   projection and assembly. Preserve their correlations and feasibility
   constraints without eagerly constructing every competitor image choice.

Do not simply drop competitor positions from the current key before applying
later generators: that projection need not commute with the recorded actions.
A correct compressed projection needs that property established explicitly.
Neither hydrogen removal, an arbitrary occupation cap, nor raising the AAM
branch cap addresses the representation problem.

## Reproduction and saved artifacts

Root: `/h/399/yunhengzou/coordinate_alignment/data/retro_runs/occupation_diagnosis_20260905/`

- `INVENTORY-000538/`, `INVENTORY-000594/`: exact input tuples, context, and
  per-call completion timings. Unfinished calls intentionally have no done file.
- `profile538.json`, `profile594.json`, `pending538.json`, `pending594.json`:
  per-stage counters and isolated timings (jobs 432157, 432158, 432160, 432161).
- `pending538.states.pkl.gz`, `pending594.states.pkl.gz`: complete native
  output states from the interrupted-call replays.
- `logs/`: capture and replay logs.

Build `bench/occupation_probe.cpp` as a pybind11 extension named
`_occupation_probe` into the diagnostic directory, using C++17 and `-O3`.
`bench/profile_occupation_call.py` loads one saved input, checks the production
native output against both probe variants, and writes measurements.
Compact evidence is in
[`occupation_expansion_diagnosis_2026-09-05.json`](occupation_expansion_diagnosis_2026-09-05.json).
