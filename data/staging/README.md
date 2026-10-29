# data/staging — provenance scripts for `reference_inventories.csv`

Scripts that built or edited `data/reference_inventories.csv`, kept so every row can be
traced to the code that wrote it. They were run once, in date order; **do not re-run
them** on the current CSV (they append or overwrite).

| date | script | what it did |
|---|---|---|
| 2026-09-14 | (hand-written first version of the CSV, commit `6f3b9b1`) | 3 legacy anchors, 12 Balakshin 2016 hardwoods (Table 4 of `manuscript/hf-2014-0328.pdf`), Capanema 2005 *E. grandis* (Table 1 of `manuscript/jf0515330.pdf`) |
| 2026-09-14 | `2026-09-14_stage_candidates.py` | appended 52 rows from `docs/13` with `use_fit = 0`, `vetted = pending` |
| 2026-09-14 | `2026-09-14_defaults_and_vetting.py` | authors' unit defaults (Capanema 5-5 ÷ 2, Peng β-β ÷ 2, anchors as three sources) and the provisional vetting (wood → fit; grasses/bark/cork/phloem/shell → validation-only) |

Later edits should be made by a new dated script here, or by hand with the change
recorded in `docs/15_expanded_calibration.md` §5.
