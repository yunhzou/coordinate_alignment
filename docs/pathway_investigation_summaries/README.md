# Pathway Investigation Summaries

These files summarize the mechanism-discovery work so far. They prioritize
metric tables, native atom indices, and concise reasoning.

## Files

| file | contents |
|---|---|
| `00_metric_terms_and_rules.md` | annotated definitions and current scoring rules |
| `01_parent_2_branching.md` | parent `2` explanation for children `2p`, `3`, and `11` |
| `02_post_3_next_move.md` | product-blind next-move analysis from `3` |
| `03_post_11_next_move.md` | product-blind next-move analysis from `11`, plus `11->12` AAM verification |
| `04_post_12_next_move.md` | corrected post-`12` read; raw `O23->C1` retained |
| `05_post_13_status.md` | explicit status note: post-`13` metrics not run yet |
| `06_path1_aam_progression.md` | direct vs path-composed AAM pathway-progress metrics |
| `07_backward_search_note.md` | short summary of backward-search lessons |

## Main Lessons

| lesson | current rule |
|---|---|
| Raw formation scores are useful | Always show raw rank and score before filters |
| Graph distance can mislead | Use graph distance as annotation, not as a hard drop except for already-bonded pairs |
| Cleavage differs from formation | Maintain separate cleavage ranking by bond class |
| Metal contacts need their own class | Do not compare Au/P or Au/substrate WBO directly to organic WBO |
| AAM is for verification and identity tracking | Do not use downstream structures to rank a product-blind proposal |
| Adjacent AAM is preferred for paths | Compose local mappings instead of relying only on direct endpoint mapping |

## Storage

This summary bundle is stored in:

```text
docs/pathway_investigation_summaries/
```
