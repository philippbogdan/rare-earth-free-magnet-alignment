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

## Opening and the held-out result (v1, frozen at 15c80aa)

Batches 1 and 2 were opened at 4060e30 (hashes verified) and scored with the frozen rules:
results/scores_v1.json. This is the held-out test; it is never overwritten.

## Post-hoc revision v2 (after opening; not a held-out test)

Rules changed after seeing the outcomes. Each item names the trigger. Scores of v2 are in
results/scores_v2.json and are reported only as post-hoc.

1. Measuring direction 'film normal / perpendicularly magnetized' was read as 'across the
   alignment axis' (encoding error; trigger uni_out_rt_fept_ep_shima2002).
2. Tetrataenite in meteorites orders below the Curie point of its ferromagnetic parent, whose own
   domain magnetisation is a bias (rule 2 of METHOD section 8 already says so; the case had been
   encoded as unbiased; trigger feni_meteorite_poirier2015).
3. Fine polytwinned L1_0 crystals: equal variants can act as one exchange-coupled multiaxial unit,
   as the spinodal networks do (the source itself says so); both readings now given (trigger
   copt_sc_kandaurova1992a).
4. Isotropic cluster-assembled films are flagged as possibly exchange-enhanced, as other
   fine-grained isotropic samples were (trigger hfco7_iso_cl_bal2014).
5. Co-rod assemblies made without an alignment step: the rod arrangement in the measured sample
   is unverified (3-D random, flat-lying, or parallel bundles) -> open; dried agglomerates of
   parallel rods have a net axis and turn in a die field (triggers bousnina2020, maurer2007,
   ener2018).
6. Fe-Cr-Co with < 3 wt% Mo: the nearest-<100> path is sourced at 3 wt% Mo only -> both paths
   (trigger fecrco_FA_sintered_fe28cr7co2mo_vompe2019).
7. Binary MnAl 'DS + hot deformation' with temperature and sequence not stated -> both readings
   (massive / displacive) (trigger mnal_other_jia2022a).
8. MnAl-C single crystals compressed during the transformation: the tetragonal-strain rule left
   2 or 4 variants tied (A <= 0.73), but the crystals reach 0.77-0.89: the transformation shear
   breaks the tie. The case ceiling reverts to the account's SC+SA cell (one variant, A = 1)
   (exposed cases mnal_sc_sa_kubo1976_s9-s15, s27-s29).
9. AlNiCo-5 (Ti-free): the open rod-axis question is settled by the opened data (equiaxed
   field-cooled Alnico 5 reaches Br/Js = 0.895 with the same-sample Js, above the 0.831 that
   rods held on the nearest <100> allow; in columnar crystals A follows the cosine of the field
   angle): rods follow the field. The phase entry cites these outcome records.
