# Signed-event decoding of saved AAM families

Keep search coverage separate from event equivalence. `iso_tolerance=1` may be useful for finding mappings in full reactions. The decoder scores their original raw WBOs with inclusive event thresholds of 0.5 for ordinary pairs and 0.3 for metal pairs, including hydrogen.

Postprocessing symmetry uses exact event-response labels. Each pair's label records whether it gives a decrease, no event, or increase against the relevant opposite-endpoint WBO values. An automorphism must preserve those labels. An approximate raw-weight tolerance such as 0.2 is insufficient near event thresholds: with reactant weights 1.00 and 0.90, swapping product weights 0.51 and 0.49 changes the event count from zero to one although those product weights differ by only 0.02.

The loose search family remains intact. Shuffles that fail the stricter symmetry check may produce different valid event alternatives, which must still be decoded.

## Running a complete decode

From the repository environment:

```sh
python bench/decode_saved_events.py \
  --archive /path/to/saved/aam.pkl.gz \
  --output /path/to/new/output-directory \
  --max-events 5
```

The command performs no AAM search. It uses one worker, a five-minute process watchdog, a monitored 3 GiB worker RSS limit, and a 6 GiB available-memory reserve. The watchdog is configurable with `--seconds`. RSS checks are polling safeguards. A completed result is recorded in `comparison.json`; a stopped or unresolved result remains explicitly incomplete and exits nonzero.

`--max-events 5` requests every distinct event class with at most five events in the retained full-mapping families. Omit `--max-events` to request all event counts. Completeness of a retained archive does not establish global AAM search completeness, particularly when its search was capped.

New event witnesses and completed-family certificates are saved incrementally. Continue from earlier directories without rerunning their completed families:

```sh
python bench/decode_saved_events.py \
  --archive /path/to/saved/aam.pkl.gz \
  --output /path/to/new/continuation-directory \
  --max-events 5 \
  --resume /path/to/previous/output-directory
```

Repeat `--resume` for every earlier certificate directory in the chain. Resume validates the archive hash, event policy, window and search configuration, and recanonicalizes stored mapping witnesses. A partially written final certificate line is ignored; malformed complete lines are rejected.

## Why the faster checks are safe

1. Reachable-image lower bounds relax correlated actions into sets of possible images. A pair contributes an unavoidable event only when all image pairs produce an event. Exceeding the event window proves the family cannot contribute. Passing the bound does not create a witness or discard correlations in the fallback model.
2. Source/target symmetry certificates prove that every mapping factors as `target symmetry · representative · source symmetry`. This preserves the event class even when literal event edges shuffle. The action order is preserved; arbitrary interleaving is not assumed safe.
3. The full correlated symbolic decoder handles the remaining families. Its event tables factor identical rows into shared selectors and accumulate fixed event counts as integers. Each solver exclusion removes an entire output event class.

`SignedEventIndex` caches exact response labels, source groups and bounded immutable edge-action metadata. `extract_path_events(..., max_patterns=None, on_pattern=callback)` removes the class-count cap and emits each newly discovered class for checkpointing. Library time budgets are soft; the command supplies the process watchdog.

The independent legacy dense canonicalizer and frozen event-table encoder are retained as test oracles, not used in the production decode path. Historical floor-based event scorers implement a different policy and remain separate.
