# Why aggregate reranking failed, and what other mappers add

No AAM search, molecular-energy calculation or new ranking rule was run here. This investigation analyzes the cached ablation and checks primary descriptions of other methods. The core algorithm and production ranker are unchanged.

## Diagnosis from our saved results

Among **1,832 completed search archives**, 1,435 top-ranked representatives match the reference and **397 do not**. The 19 incomplete original searches are outside this failure decomposition, but remain in the full benchmark accuracy denominator.

Of those 397 wrong top-1 representatives:

- **44** have no reference-equivalent *saved representative*. This is not proof that their compressed families lack the reference.
- **225** have a reference representative at the same event count.
- **128** have a reference representative requiring more events: 57 need +1, 55 need +2, 11 need +3, two need +4 and three need +5.
- **80** of the tied cases have identical cached event count, center-component count, edited-atom count, weighted bond cost and missing-energy metadata. Reweighting only these aggregate features cannot separate that reference/top-candidate pair. Atom-index tie-breaking remains arbitrary with respect to chemistry.

Compactness has little discriminative range: **1,785 / 1,832 top candidates already have one connected reaction center** under our definition. Among the 353 wrong top candidates with an available reference representative, 338 have the same center-component count as the reference.

For the 314 wrong/reference pairs with energy features available on both representatives, the reference energy-weighted edit cost is **higher in 163**, equal in 98 and lower in only 53. These pair-level support counts differ from the previous ablation's whole-case support criterion. Simply increasing the weight of this energy feature is not a generally corrective signal.

These findings concern the tested representations and scores. They are not a proof that every graph-based ranker must fail, nor that every reference annotation is chemically correct.

## Concrete assignments to inspect

### Case 9: oxygen origin, not generic bond strength

In the saved [case 9 viewer](../golden_alternatives_20260907/case9/viewer.html):

- Product ring oxygen **p4** comes from ketone oxygen **r8** in top class 1, but from alcohol oxygen **r5** in reference-equivalent class 2.
- Product water oxygen **p14** correspondingly comes from r5 in class 1 and r8 in class 2.
- Class 1 has **8 events** and normalized weighted cost **8.34**. Reference-equivalent class 2 has **9 events** and weighted cost **10.99**.

Thus both generic objectives favor the non-reference assignment. A transformation-aware oxygen-origin preference can distinguish these candidates; a generic preference for fewer or cheaper edits cannot. This is a statement about the supplied reference and saved assignments, not an independently reconstructed experimental mechanism.

### Case 21: water incorporation

In the saved [case 21 viewer](../golden_alternatives_20260907/case21/viewer.html), reference product carbonyl oxygen **p5** comes from water oxygen **r13**. Top class 1 instead uses an oxygen already in the organic precursor, **r2**. It has seven events versus eight for reference-equivalent class 4. Favoring structural retention alone does not encode the role of water in the transformation.

## What other methods use

| Method | Information beyond generic edit totals |
|---|---|
| Jaworski et al. mapper (2019) | Twenty expert reaction heuristics generate intermediate assignments and give recognized transformations lower scores than unrelated individual bond cuts. Categories include pericyclic reactions, carbonyl chemistry, metathesis and rearrangements. This explicitly allows a chemically preferred mapping to require more raw edits. [Paper](https://www.nature.com/articles/s41467-019-09440-2) |
| NameRXN | A mechanism-based rule system. The authors describe applying SMIRKS transformations to the reactants and identifying a mechanism/mapping when the product is obtained. This is not merely MCS plus bond energies. [Author presentation](https://www.nextmovesoftware.com/talks/Sayle_Regioselectivity_ACS_201808.pdf), [product documentation](https://www.nextmovesoftware.com/namerxn.html) |
| AAMFixer | Explicit predefined remapping rules derived from corresponding wrong/correct mapped reactions. This adds correction knowledge downstream of a mapper. [Author repository](https://github.com/cimm-kzn/AAMFixer) |
| RXNMapper | Reaction-language-model attention learned from unmapped reaction data, with element matching and a neighbor-attention multiplier during decoding. Its chemical preferences are learned rather than a library of named-reaction rules. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8026122/) |
| LocalMapper | Graph/message-passing and attention model trained on chemist-labeled mappings through human-in-the-loop curation. Reaction templates are also used in its confidence assessment; it is not simply a hand-written template mapper. Its examples specifically discuss correcting oxygen origins in hydrolysis. [Paper](https://www.nature.com/articles/s41467-024-46364-y) |
| ChemAxon AutoMapper | Its public documentation describes MCS and minimal chemical distance. That description alone does not establish what all internal chemical preferences are. [Documentation](https://docs.chemaxon.com/latest/automapper_user-guide.html) |

## Interpretation and next experiment

The strongest methods are not necessarily solving a more accurate generic energy formula. They add **transformation context**: which atom is the nucleophile, which group leaves, where water contributes an atom, or which correlated edits constitute a recognized rearrangement. Aggregate bond costs discard much of that context. Generic average homolytic bond energies also do not model a catalytic, ionic or multistep pathway's activation barrier.

The clean next experiment is a **separate local reaction-pattern scorer over existing AAM candidates**, not additional constraints inside fragment matching:

1. Build a local before/after graph around each candidate's edited bonds, preserving elements, bond orders, charges, H counts, input-molecule membership and the correlated atom assignment.
2. Score that change pattern against either independently defined chemical transformation rules or statistics learned from a separate curated dataset.
3. Keep all original alternatives; report the preference and its reason rather than silently replacing mappings or deleting unsupported candidates.
4. Evaluate on a frozen split, including gains, losses and unknown patterns. Golden-informed rule development is exploratory tuning, not independent SOTA validation.

This is a proposed experiment, not a newly implemented chemistry rule engine or a demonstrated accuracy gain. It can use the existing shortlist, so no quantum calculations or full AAM reruns are inherently needed. Full model inference would be a separate computational choice.

## Reproducibility

`bench/diagnose_golden_ranking.py` generates `summary.json` and per-failure `cases.json` from `reports/golden_fast_rerank_20260907/features.json.gz`. Case 9 and 21 atom correspondences were checked against their saved mappings and original endpoint bond matrices. All counts refer to the unchanged fixed-search archives. Primary sources above were consulted on 2026-09-07.
