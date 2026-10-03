# Protocol log (held-back cases)

Events that bear on the claim that held-back outcomes were not seen before predictions were
frozen. Each entry: date, what was exposed, consequence.

- 2026-10-03 — While searching agent working files for an input sentence (Farjami 2009, CoPt
  nucleation temperature), a broad `grep` printed lines from `.tmp/l10_work/part_pd_outcomes.yaml`,
  exposing part of one sealed FePd outcome (a favoured-variant fraction after two-step field
  ordering at 673 K + 773 K) and the authors' notes on its temperatures. Consequence: all FePd
  field-ordering outcome records (route SC+FA) are excluded from the held-out accuracy and
  reported as "exposed" development cases. Searches are now restricted to source-text folders.
- 2026-10-03 — The value-free conditions export (`data/sealed/*.conditions.yaml`) redacted text but
  not numeric fields inside `conditions`. Collectors had placed outcome-related numbers there
  (`Mr_Am2_kg`, `M_30kOe_Am2_kg`, `M_at_max_field_Am2_kg`, `BHmax_kJ_m3`, `Br_before_600C_T`,
  `J_at_10kOe_T`, `Js_same_sample_T`, `Ku_*`, `Hc_*`). While reading the MnAl case list (first 113
  lines) these were seen for: mnal_hd_madugundo2016_perp/_par, mnal_am_tiden2024_ref_mrms,
  mnal_sc_sa_kubo1976_ex5, mnal_sc_sa_kubo1976_s3_*, the (BH)max values of mnal_hd_* (Kubo,
  Sakamoto 1980, Kojima 1979) and of mnal_sc_sa_kubo1976_s5..s33 ((BH)max bounds Br from below),
  and Ku/Hc of mnal_ep_takata2018_*. Consequence: every MnAl case whose record appeared in those
  lines is classed "exposed" (development, not held out). The export now keeps only an allow-list
  of process parameters; all condition files were regenerated before any further reading.
- 2026-10-03 — The process text of mnbi_other_tang2022_fieldsinter250 contains "(270 C: 0.86)", a
  remanence ratio of a sibling sample of the same series. Consequence: both Tang 2022
  field-sintering cases are classed "exposed".

## Freeze before opening batch1 and batch2 (2026-10-03)

Predictions for every sealed case (batch1: 10 files; batch2: ndfeb_routes, borides_feta) and
for the known cases, together with the scoring rules (src/align/scoring.py) and the method
changes (docs/METHOD.md section 8), are committed in the commit that adds this entry, before
`tools/seal.py open` is run for either batch.

    sha256 data/curated/cases.yaml        eb6ee47de1f561568ff545347e7336fb226fee9f205a0a82637db39892252eca
    sha256 data/curated/cases_known.yaml  b10093d6d74e49d65376add9c557c61f464f621244a5748166b8ca9fcdb31499

Development data seen while writing the rules: the value-free condition files only, and the
known cases of PROBLEM.md (data/raw/known_cases_*.yaml), on which the scoring was dry-run.
That dry run exposed a grade-to-Js mapping error (columnar alnico 5DG/5-7 read as alnico 8),
fixed from sourced Js values before this freeze. Exposed sealed cases remain as listed above.
The Mn4N collection (Codex task A3) is still running; its outcomes will be sealed as a later
batch and frozen separately.
