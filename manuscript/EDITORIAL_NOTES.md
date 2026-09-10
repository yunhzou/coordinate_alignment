# Author review for the first draft

The paper is a substantive first draft with generated figures and recorded-event
animations. Resolve the following before treating it as a submission manuscript.

1. Author order supplied by Yunheng Zou: Yunheng Zou; Olalla Nieto Faza; Shifa Hussain;
   Varinia Bernales; Alán Aspuru-Guzik. Bernales and Aspuru-Guzik are the principal
   investigators. The provided institutional reference assigns Zou to 1 and 6,
   Bernales to 2 and 7, and Aspuru-Guzik to 1–9; those affiliations are now in the
   paper. Nieto Faza's University of Vigo affiliation is confirmed and listed
   under number 10 with her Department of Organic Chemistry. Supply Hussain's
   affiliation. Confirm whether both PIs are corresponding authors and supply contact
   details if desired. Equal-contribution symbols from the reference article
   have not been transferred to this author list. Confirm the method name, title,
   individual contributions, funding and competing interests. The proposed method
   name is **GRAFT**, with working title **GRAFT: Symmetry-Aware Atom Mapping through
   Continuous Fragment Growth**. Source reports and software APIs retain their
   original AAM identifiers.
2. Document the provenance and selection of the 140 elementary steps, the precise
   calculation method/version used for their cached WBO matrices, and which
   structures/data may be redistributed publicly. No WBO calculation level or
   dataset citation has been invented. This collection was inspected during
   development and lacks reference mappings; its results are feasibility,
   score and alternative-coverage evidence, not held-out chemical accuracy.
3. Decide which frozen engine is the main public release. The newer seed ablation
   uses `98b01b1` (1,836 verified at ten orders); the earlier publication/collection
   uses `bb0a7d7` (1,839 verified). The draft keeps them separate. Later adaptive
   or shared-policy optimizations are not silently substituted.
4. The located full Golden no-sweep experiment is the adaptive agenda at
   `0ce0f17`, with 1,723 verified, 33 absent and 95 unknown. It also changes
   scheduler and cap policy. It cannot serve as a pure “sweep on/off” ablation
   of the mature scheduler. The draft reports it under its actual scope. If a
   separate matched-policy no-sweep report is intended, identify that artifact
   before substituting numbers.
5. Review the method description against the intended core contribution:
   ongoing weighted fragment extension, conditional saturated-placement
   branching, local boundary-aware canonicalization, exact continuation-state
   joins, and downstream correlated-family equality/containment. The current
   evidence does not isolate their individual performance contributions. The
   algorithm is the central contribution. SLAP sweep results serve as a supporting
   ablation of single-edge constraint relaxation; the sweep unions individual
   cuts and does not establish that one fixed cut is sufficient.
6. Revise related work using the [10 September novelty audit](evidence/novelty-audit-20260910.md).
   The expanded [review-led assessment](evidence/review-led-novelty-20260910.pdf)
   and [coverage ledger](evidence/review-coverage-20260910.csv) add historical
   fragment-enumeration leads, direct Indigo/RDT/SynKit approximate-growth
   comparisons, and AMLGAM's weighted objective/post-processing distinction.
   Use the historical reviews as discovery sources; unresolved originals are
   not evidence that their complete mechanisms have been excluded.
   For the specific algorithmic argument, use the follow-up
   [mechanism comparison](evidence/mechanism-comparison-20260910.md): aggregate
   accept/defer growth, saturation-level whole-mapping branches, whole-island
   absorption, and contextual continuation equality. Compare OpenChemLib's
   immediate neighbor-pair commitment and FMCS's fragment-shape enumeration.
   Incremental embedding reuse alone can reproduce a fresh subgraph search under
   ideal completeness assumptions; distinguish representation efficiency from
   the search policy's substantive decisions.
   The initial five scientific references are insufficient for a novelty argument.
   Add the closest comparisons: SMSD/McSplit, Jaworski's neighborhood expansion
   and cut search, Mann's constraint-based mapping, Ali's symmetry-distinct
   enumeration, Laffitte's partial-map extensions, and the inspected SLAP/SynKit
   implementations. SLAP already branches and suppresses repeated states.
   Emphasize the specific shared-fragment transition, compressed placements,
   contextual continuation joins and correlated families; fragment growth,
   branching and symmetry reduction alone are established ideas. The targeted
   audit found no exact match for the complete design but does not establish
   priority. The related-work revision and matched component ablations remain
   outstanding; no global optimality, physical mechanism or TS-accuracy claim
   is established by this assessment.
7. Supply a public code URL, archive DOI, licensing and reproducible environment
   for the frozen benchmark releases. The included compact evidence supports
   figure regeneration; it is not the entire multi-gigabyte mapping archive.

Numerical interpretation to preserve during editing:

- Golden recovery concerns a verified reference witness among retained families,
  not top-1 prediction. All 1,851 records remain in denominators; unknowns are
  not counted as successes or certified absences.
- One, two, three and ten “seeds” mean seed orderings per cut, not the number of
  atoms attempted or a single atom-pair initialization.
- SLAP sweep recovery is 1,795; union with an older expanded run is 1,797 and
  has a different budget. The 1,661 uncut control combines both directions and
  bond modes; it is not the default single-call result.
  The [complete saved-output recount](evidence/slap-sweep-count-audit-20260910/README.md)
  reproduces these totals. Standalone nonrecoveries are 52 completed misses plus
  four unresolved cases; the expanded union has 50 plus four. Case 1665 timed
  out in both directions. The four unresolved cases contain nine unfinished
  variants, not four jobs, and completing them is an extended-budget follow-up.
  The 59,466.8/176,323.7 three/ten-seed CPU totals use an AAM-only common set of
  1,837 cases; the manuscript's SLAP-inclusive comparison uses 1,807 cases.
- Paired CPU uses 1,807 mutually completed Golden cases. Instrumentation differs
  and interrupted work is omitted. It does not establish latency or total cost.
- The primary 140-case comparison keeps cap 100. Case 123's cap-200 result is
  a separate follow-up. An extra low-event mapping is not automatically a
  chemically correct alternative.
- Figure 1 is schematic. Movies use actual event traces and terminal mappings,
  sample at most five witnesses per event, and illustrate an event-enabled Python
  path. Their timing is chosen for presentation and is unrelated to runtime.

The user previously excluded the synthetic six-vertex exact-oracle diagnostic
from paper comparisons; it remains excluded.
