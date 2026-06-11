# Post-3 Next-Move Investigation

Current structure: intermediate `3`, charge `+1`.

Question: from `3` alone, what move should be proposed next?

## One-Body Reactivity Centers

| rank | donor center | donor strength | q | f- | acceptor center | acceptor strength | q | f+ |
|---:|---|---:|---:|---:|---|---:|---:|---:|
| 1 | `Au20` | 1.000 | +0.181 | 0.147 | `P21` | 1.000 | +0.017 | 0.186 |
| 2 | `P21` | 0.721 | +0.017 | 0.106 | `C27` | 0.778 | +0.362 | 0.108 |
| 3 | `O5` | 0.455 | -0.336 | 0.050 | `O26` | 0.508 | -0.212 | 0.096 |
| 4 | `C2` | 0.436 | -0.144 | 0.056 | `O28` | 0.444 | -0.198 | 0.084 |

## Raw Formation Candidates

| raw rank | move | score | distance A | WBO now | class/read |
|---:|---|---:|---:|---:|---|
| 1 | `C2->C27` | 0.1108 | 2.81 | 0.000 | substrate/acyl-fragment engagement |
| 2 | `C2->O26` | 0.0849 | 2.52 | 0.000 | same local acyl/O cluster |
| 3 | `C2->O28` | 0.0803 | 2.43 | 0.000 | same local acyl/O cluster |
| 4 | `Au20->C27` | 0.0706 | 4.88 | 0.000 | metal-substrate long contact |
| 5 | `Au20->O26` | 0.0599 | 4.48 | 0.000 | metal-substrate long contact |
| 6 | `Au20->O28` | 0.0576 | 4.37 | 0.000 | metal-substrate long contact |
| 7 | `C1->C27` | 0.0546 | 2.39 | 0.000 | substrate |
| 8 | `O5->Au20` | 0.0539 | 2.91 | 0.000 | metal coordination |

## Cleavage Candidates

| rank | bond | cleavage score | WBO | distance A | note |
|---:|---|---:|---:|---:|---|
| 1 | `Au20-P21` | 0.7922 | 0.659 | 2.36 | metal-ligand; track separately |
| 2 | `C2-Au20` | 0.6628 | 0.663 | 2.07 | metal-substrate |
| 3 | `C7-O26` | 0.6092 | 0.733 | 1.58 | strongest organic cleavage |
| 4 | `C1-O28` | 0.5098 | 0.844 | 1.46 | organic |
| 5 | `C4-O5` | 0.4987 | 0.909 | 1.45 | organic |
| 6 | `C3-O5` | 0.4930 | 0.921 | 1.43 | organic |

## Product-Blind Call

The signal after `3` is not just a pure formation event. It is a coupled channel:

| component | best metric signal |
|---|---|
| formation cluster | `C2->C27/O26/O28` ranks 1-3 |
| organic cleavage | `C7-O26` is the strongest organic cleavage |
| metal contact | `O5->Au20` exists but is raw formation rank 8 |

Read: propose weakening/cleavage of `C7-O26` while `C2` engages the
`C27/O26/O28` fragment. This was later supported by adjacent AAM in path1:
native `3` `C7-O26` maps to a threshold cleavage in the next step.

## Raw Artifacts

- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current3_intrinsic_forward.summary.md`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current3_intrinsic_forward.formation_scores.csv`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current3_intrinsic_forward.cleavage_scores.csv`
