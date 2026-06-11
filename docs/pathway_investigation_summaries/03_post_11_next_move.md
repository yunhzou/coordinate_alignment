# Post-11 Next-Move Investigation

Current structure: intermediate `11`, charge `+1`.

Question: from native `11`, what move should be proposed next?

## Raw Formation Candidates

| raw rank | move in native 11 | score | distance A | WBO now | class | read |
|---:|---|---:|---:|---:|---|---|
| 1 | `Au27->O24` | 0.1609 | 3.15 | 0.000 | metal-substrate | coordination signal |
| 2 | `Au27->C25` | 0.1060 | 4.43 | 0.000 | metal-substrate | long metal contact |
| 3 | `C3->C25` | 0.0934 | 2.25 | 0.000 | organic | raw high, later treated as likely artifact for this step |
| 4 | `O23->C1` | 0.0850 | 2.52 | 0.000 | organic | secondary organic collapse candidate |
| 5 | `C1->O24` | 0.0755 | 2.39 | 0.000 | organic | substrate candidate |
| 6 | `C1->C25` | 0.0596 | 3.46 | 0.000 | organic | lower organic candidate |
| 7 | `C3->O26` | 0.0516 | 2.35 | 0.000 | organic | local substrate candidate |

## Cleavage Candidates By Class

| class rank | bond in native 11 | class score | WBO | distance A | read |
|---:|---|---:|---:|---:|---|
| organic 1 | `C4-O26` | 0.7062 | 0.795 | 1.54 | strongest organic cleavage |
| organic 2 | `C3-O24` | 0.6698 | 0.862 | 1.46 | just-formed 2->11 bond; backward risk |
| organic 3 | `C17-O23` | 0.6323 | 0.914 | 1.45 | O23 framework contact |
| organic 4 | `C2-O23` | 0.6285 | 0.898 | 1.45 | O23 framework contact |
| metal-substrate 1 | `C1-Au27` | 0.3500 | 0.700 | 2.07 | metal contact |
| metal-ligand 1 | `Au27-P28` | 0.3500 | 0.639 | 2.38 | ligand contact |

## Coupled Candidates

| rank | formation | cleavage | coupled score | f score | cleavage score | shared endpoint? |
|---:|---|---|---:|---:|---:|---|
| 1 | `Au27->O24` | `C3-O24` | 0.2353 | 0.1609 | 0.6698 | yes |
| 2 | `Au27->O24` | `C1-Au27` | 0.1710 | 0.1609 | 0.3500 | yes |
| 6 | `C3->C25` | `C3-O24` | 0.1366 | 0.0934 | 0.6698 | yes |
| 10 | `O23->C1` | `C17-O23` | 0.1203 | 0.0850 | 0.6323 | yes |
| 12 | `O23->C1` | `C2-O23` | 0.1199 | 0.0850 | 0.6285 | yes |

## Corrected Product-Blind Call

The raw formation table suggested several possibilities, including `C3->C25`.
Manual/native-index inspection showed `C3-C25` was not the useful next-step
description for this structure. The more robust call was cleavage-centered:

| proposal | support |
|---|---|
| `C4-O26` weakens/breaks | top organic cleavage score, 0.7062 |
| `Au27` moves toward `O24` | top raw formation/coordination score, 0.1609 |
| avoid immediate backward move | `C3-O24` is the just-formed `2->11` organic bond |

## Downstream AAM Verification: 11 -> 12

| event type | pair in 11 | pair in 12 | WBO 11 | WBO 12 | delta |
|---|---|---|---:|---:|---:|
| weakened | `C1-C3` | `C1-C3` | 1.819 | 1.292 | -0.527 |
| broken | `C4-O26` | `C4-O26` | 0.795 | 0.000 | -0.795 |
| strengthened | `C3-C4` | `C3-C4` | 0.964 | 1.518 | +0.553 |

Read: the post-`11` prediction was validated by the observed `11 -> 12`
transition: native `11` `C4-O26` breaks.

## Raw Artifacts

- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current11_intrinsic_forward.summary.md`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current11_intrinsic_forward.formation_scores.csv`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current11_intrinsic_forward.cleavage_scores_by_class.csv`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/aam_11_to_12/summary.md`
