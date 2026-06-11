# Backward Search Note

Question tested: if a product/intermediate `P` is given, can we infer the likely
previous intermediate by reversing the forward logic?

## Result

Backward inference is harder than forward proposal. The simple reverse metric
ranked parent `2` only around rank 3 in the three-child test, which showed that
formation likelihood and cleavage likelihood are not symmetric.

## Why It Is Different

| forward question | backward question |
|---|---|
| Which absent/weak contact wants to form from current geometry? | Which existing bond/contact may have just formed and should be undone? |
| Uses local donor/acceptor readiness plus distance | Needs to decide which strong present bond is historically new |
| Often one-step local | More combinatorial; many possible removals relax to similar structures |

## Brutal Search Observation

A brute-force backward perturbation search was tried by weakening/removing
candidate WBO contacts and relaxing. It was useful diagnostically but not clean
enough to replace forward scoring.

| observation | meaning |
|---|---|
| Many perturbations remain close to the child | Removing one bond/contact can be absorbed by relaxation |
| Energetic ordering alone was not decisive | Relaxed energies reflect many effects, not just pathway history |
| AAM verification is still needed | Mapping identifies whether the recovered event matches the pathway |

## Current Rule

Use backward logic as a verification or diagnostic mode, not as the primary
intermediate generator.

For a known pathway, prefer:

1. Local adjacent AAM.
2. Composed identity mapping.
3. Head/tail progress counts as sanity checks.
4. Backward brute force only for targeted ambiguity.

Raw artifact:

- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/backward_bruteforce_wbo030_120/summary.csv`
