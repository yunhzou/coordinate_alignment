# Path1 AAM Progression Investigation

Path tested: `2 -> 3 -> 4 -> 5 -> 6 -> 7in -> 8in -> 9`, charge `+1`.

Purpose: test whether mapped WBO event counts can quantify pathway progress, and
whether adjacent/composed AAM is more reliable than direct endpoint remapping.

## Direct Endpoint AAM Counts

Direct comparison maps each intermediate independently to head `2` and tail `9`.

| step | events from 2 | organic from 2 | remaining events to 9 | remaining organic |
|---:|---:|---:|---:|---:|
| 3 | 3 | 2 | 11 | 9 |
| 4 | 3 | 3 | 9 | 7 |
| 5 | 6 | 5 | 5 | 5 |
| 6 | 8 | 7 | 7 | 7 |
| 7 | 8 | 7 | 4 | 2 |
| 8 | 8 | 7 | 4 | 2 |

Read: mostly sensible, but endpoint minimum-event mappings can rematch similar
atoms and hide path identity.

## Path-Composed AAM Counts

Adjacent mappings were composed along the path to preserve atom identity.

| step | composed events from 2 | organic from 2 | composed remaining events to 9 | remaining organic |
|---:|---:|---:|---:|---:|
| 3 | 3 | 2 | 12 | 10 |
| 4 | 5 | 5 | 10 | 8 |
| 5 | 8 | 7 | 8 | 6 |
| 6 | 10 | 9 | 8 | 6 |
| 7 | 10 | 9 | 4 | 2 |
| 8 | 10 | 9 | 4 | 2 |

Read: composed AAM is preferred for pathway diagnostics because atom identities
are inherited through local steps.

## Adjacent Break/Form Counts

| transition | forward broken | forward formed | backward broken | backward formed |
|---|---:|---:|---:|---:|
| `2->3` | 1 | 2 | 2 | 1 |
| `3->4` | 1 | 0 | 0 | 1 |
| `4->5` | 2 | 1 | 1 | 2 |
| `5->6` | 0 | 1 | 1 | 0 |
| `6->7` | 2 | 1 | 1 | 2 |
| `7->8` | 0 | 0 | 0 | 0 |
| `8->9` | 2 | 2 | 2 | 2 |

## Adjacent WBO Changes In Original-2 Atom Frame

| transition | events | original-2-frame changes |
|---|---:|---|
| `2->3` | 3 | `C1-C2` 2.516->1.848; `C1-O18` 0.000->0.844; `C2-Au31` 0.232->0.663 |
| `3->4` | 1 | `C7-O16` 0.733->0.000 |
| `4->5` | 3 | `C3-O5` 0.907->0.000; `C2-Au31` 0.385->0.000; `O5-Au31` 0.142->0.628 |
| `5->6` | 1 | `C2-O5` 0.000->0.923 |
| `6->7` | 3 | `C1-C3` 1.723->1.208; `C4-O5` 0.891->0.000; `C3-C4` 0.995->1.612 |
| `7->8` | 0 | no threshold WBO event |
| `8->9` | 4 | `O5-Au31` 0.545->0.000; `C17-O18` 0.809->0.000; `O18-Au31` 0.000->0.445; `O5-C17` 0.000->1.080 |

## Read

The path is not a simple monotonic count of plus/minus events. Local steps can
redistribute WBO while cumulative distance to head/tail changes in jumps. Still,
forward composed distance generally increases from `2`, and remaining distance
to `9` generally shrinks.

This supports using forward progress as a useful intermediate plausibility
metric, with backward tolerance around 1-2 events. Larger backward movement
should trigger closer inspection of AAM quality or branch choice.

## Raw Artifacts

- `tmp_xtb_single_points/path1_direct_aam_chrg_plus1/direct_head_tail_aam/head_tail_aam_distance_report.md`
- `tmp_xtb_single_points/path1_direct_aam_chrg_plus1/path_composed_aam/path_composed_aam_report.md`
- `tmp_xtb_single_points/path1_direct_aam_chrg_plus1/path_composed_aam/adjacent_progression_head2_frame_report.md`
