# Post-12 Next-Move Investigation

Current structure: intermediate `12`, charge `+1`.

Question: from native `12`, what move should be proposed next?

## Previous-Step Guardrail

AAM `11 -> 12` showed:

| pair | status in 11 -> 12 | use as guardrail |
|---|---|---|
| `C4-O26` | broken | re-forming it is likely backward |
| `C1-C3` | weakened | re-strengthening may be backward |
| `C3-C4` | strengthened | breaking it may be backward |

These are guardrails, not hard product knowledge for scoring.

## Raw Formation Candidates

| raw rank | move in native 12 | score | distance A | direct WBO | graph distance | class | current treatment |
|---:|---|---:|---:|---:|---:|---|---|
| 1 | `Au27->O23` | 0.2079 | 2.99 | 0.000 | 3 | metal-substrate | valid coordination candidate |
| 2 | `O23->C1` | 0.1843 | 2.43 | 0.000 | 2 | organic | retain; correct high-rank candidate |
| 3 | `O23->Au27` | 0.1550 | 2.99 | 0.000 | 3 | metal-substrate | valid coordination candidate |
| 4 | `O26->C1` | 0.1474 | 2.61 | 0.000 | 4 | organic | valid but not top raw organic truth |
| 18 | `C1->O23` | 0.0590 | 2.43 | 0.000 | 2 | organic | reverse-direction lower score |

Correction learned here: the old `graph_distance > 2` hard filter was wrong.
`O23->C1` is graph-distance 2 through `C1-C2-O23`, but the direct WBO is zero
and the 3D distance is good. Graph-distance 2 should be retained and annotated
as proximal/collapse-like.

## Cleavage Candidates

| rank | bond in native 12 | class score | WBO | distance A | read |
|---:|---|---:|---:|---:|---|
| 1 | `C2-O23` | 0.6973 | 0.878 | 1.45 | strongest organic cleavage; relevant to O23 collapse |
| 2 | `C17-O23` | 0.6629 | 0.958 | 1.41 | O23 framework contact |
| 3 | `C3-O24` | 0.5909 | 0.964 | 1.41 | organic |
| 4 | `C2-C17` | 0.5518 | 0.820 | 1.53 | organic |
| 5 | `O24-C25` | 0.5002 | 1.066 | 1.38 | organic |
| 6 | `C25-C32` | 0.4975 | 1.005 | 1.50 | organic |
| 7 | `C1-C2` | 0.4744 | 1.140 | 1.46 | local to O23-C1 formation |

## Coupled Candidates

| rank | formation | cleavage | coupled score | f score | cleavage score | prox A | read |
|---:|---|---|---:|---:|---:|---:|---|
| 1 | `Au27->O23` | `C2-O23` | 0.3112 | 0.2079 | 0.6973 | 0.00 | coordination + O23 bond lability |
| 2 | `Au27->O23` | `C17-O23` | 0.3022 | 0.2079 | 0.6629 | 0.00 | coordination + O23 bond lability |
| 3 | `O23->Au27` | `C2-O23` | 0.2320 | 0.1550 | 0.6973 | 0.00 | reverse direction notation for same contact |
| 10 | `O26->C1` | `C1-C2` | 0.1795 | 0.1474 | 0.4744 | 0.00 | old filtered organic candidate |
| 11 | `O26->C1` | `C2-O23` | 0.1765 | 0.1474 | 0.6973 | 1.46 | old filtered organic candidate |

The stored coupled table was generated after the old filter, so it does not rank
`O23->C1` coupled variants correctly. The raw formation table is the more
important artifact for this stage.

## Final Read

| conclusion | support |
|---|---|
| correct next formation is `O23->C1` | raw rank #2, score 0.1843, distance 2.43 A, direct WBO 0.000 |
| graph-distance hard filter failed | it removed graph-distance-2 `O23->C1` even though the move is valid |
| likely local cleavage partners | `C2-O23` and/or `C1-C2`, based on cleavage scores and local topology |
| metal channel remains real | `Au27->O23` is raw rank #1 and may assist or compete |

Read: the detector worked; the post-filter was the weak part. Keep raw ranks in
all future reports.

## Raw Artifacts

- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current12_intrinsic_forward.formation_scores.csv`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current12_intrinsic_forward.cleavage_scores_by_class.csv`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/current12_intrinsic_forward.coupled_move_scores.csv`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/12_native_index_viewer.html`
