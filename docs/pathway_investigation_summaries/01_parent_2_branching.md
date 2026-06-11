# Parent 2 Branching: 2 -> 2p, 3, 11

Current structure: parent `2`, charge `+1`.

Question: can descriptor-only scoring explain why `2p`, `3`, and `11` are
probable children?

## Raw Formation Matrix

| raw rank | proposed move in native 2 | score | distance A | parent WBO | observed child |
|---:|---|---:|---:|---:|---|
| 1 | `O18->Au31` | 0.2383 | 2.73 | 0.118 | latent/pre-existing contact |
| 2 | `O5->Au31` | 0.1163 | 3.21 | 0.000 | `2p` |
| 3 | `O18->C2` | 0.0569 | 2.95 | 0.000 | `11` |
| 4 | `O18->C1` | 0.0414 | 3.79 | 0.000 | `3` |

Interpretation: the three known child channels all appear in the top four raw
formation candidates from parent `2`, without using child structures.

## Observed Child Events

| transition | dominant event | WBO before | WBO after | note |
|---|---|---:|---:|---|
| `2 -> 2p` | `O5-Au31` forms | 0.000 | 0.414 | metal coordination branch |
| `2 -> 3` | `O18-C1` forms | 0.000 | 0.844 | organic formation branch |
| `2 -> 3` | `C1-C2` weakens | 2.516 | 1.848 | coupled bond-order redistribution |
| `2 -> 3` | `C2-Au31` strengthens | 0.232 | 0.663 | metal support/contact |
| `2 -> 11` | `O18-C2` forms | 0.000 | 0.862 | native `11` pair is `C3-O24` |
| `2 -> 11` | `C1-C2` weakens | 2.516 | 1.819 | coupled bond-order redistribution |
| `2 -> 11` | `C1-Au31` strengthens | 0.226 | 0.700 | metal contact strengthens |

## AAM Verification For 2 -> 11

| event type | pair in 2 frame | pair in native 11 | WBO 2 | WBO 11 | delta |
|---|---|---|---:|---:|---:|
| weakened | `C1-C2` | `C1-C3` | 2.516 | 1.819 | -0.696 |
| formed/strengthened | `C1-Au31` | `C1-Au27` | 0.226 | 0.700 | +0.474 |
| formed | `C2-O18` | `C3-O24` | 0.000 | 0.862 | +0.862 |

Read: the parent-2 scorer captured the dominant organic alternatives and the
metal coordination alternative. The remaining WBO changes are coupled response,
not separate hand-written mechanism rules.

## Raw Artifacts

- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/parent2_donor_acceptor_score.md`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/parent2_donor_acceptor_pair_scores.csv`
- `tmp_xtb_single_points/run_chrg_plus1_correct2p_AO9ezU/aam_2_to_11/summary.md`
