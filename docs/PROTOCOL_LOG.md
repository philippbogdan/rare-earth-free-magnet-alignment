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
