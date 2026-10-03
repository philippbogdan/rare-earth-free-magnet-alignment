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

## Freeze before opening batch3 (Mn4N)

Mn4N curated (cubic, four <111> easy axes: a random polycrystal keeps the nearest one, A = 0.866)
and its six sealed cases predicted with the v2 rules before opening batch3.

    sha256 data/curated/cases_v2.yaml  fed2017ea7a6c8dc82c6ddf87cdef7f4300bca417e3ae55113b0eac10bc391af
10. (account text, no prediction changed) Spinodal phases: every 'cannot align' cell states the
    coupled-network bound that held in all held-out unbiased cases (0.866 for random grains, 0.934
    along the best axis of a <100> fibre); the cell value stays the uniaxial reference 1/2, so
    route comparisons are between alignment ceilings.

## Freeze v3 before opening batch 4 (heldout2_b: 139 records from sources not used before)

Rules v3 = v2 plus the account changes logged above (EP geometric ceiling, physical TGG and
direct-formation reasons, unbiased bound stated) and, for this batch, phase-label mappings,
Br units mT/kGs, same-sample Js links of the form _br_/_js_same_sample_, mass-unit Js not
converted, exchange-enhanced bonded/composite samples flagged, MnBi powder milled from
ribbon open, epsilon-Fe2O3 variant-family fractions not A. No prediction for batches 1-3
changed (checked against cases_v2.yaml).

    sha256 data/curated/cases_v3.yaml  b203b12b3f0b75f6935c9cde3aa4ee783f2c0fb79f406b5557cf7bef55716588
freeze v3: 3d50da0

## After opening batch 4

Frozen v3 result: results/scores_v3.json (held-out = batch 4). A conversion bug found after
opening: for phases other than AlNiCo the specs did not carry the link to a same-sample Js, so
Br was converted with the phase Js (understating A). The frozen result is kept; a rescoring with
the same-sample Js is reported next to it (scores_v3.json, sensitivity_js), and the link is
fixed for later versions. In that rescoring one isotropic bonded Nd-Fe-B print exceeds 1/2 by
exchange enhancement of its rapidly solidified powder (excluded from A by definition); the
v3 rules had not flagged 'Nd-Fe-B/PPS bonded' as such (fixed for later versions).

## Freeze v4 before opening batch 5 (heldout2_a: 60 records from sources not used before)

Rules v4 = v3 plus: same-sample Js linked for every phase; bonded isotropic Nd-Fe-B and
hard/soft mixtures (phase content with alpha-Fe, cubic Co3Pt, soft phases) flagged as
exchange-coupled (not tests); thin films measured in plane not tests (shape anisotropy),
unbiased (111)-fibre L1_0 films capped at 0.577 out of plane; epitaxial L1_0 films measured
along the normal tested; alnico abnormal grain growth to near-[111] grains as a tie; stress-
biased grain growth with unmeasured texture not a test; Fe-Cr-Co cold rolling between ageing
steps = deformation ageing; '8h' in an id no longer read as grade 8H.

    sha256 data/curated/cases_v4.yaml  b8ce592f93b14a5bc774115a5cfe590a3fd1c647e7332174071f77fe0123c91b
freeze v4: 92b950f

## After opening batch 5

Frozen v4 result: results/scores_v4.json (held-out = batch 5). Two 'cannot align' predictions
failed: loose L1_0 FePt nanoparticle powders (sintered on NaCl at 700-750 C, 6-7 nm) measured
as loose powder reach Mr/Ms above the random value. Diagnosis: loose particles can turn in the
measuring field, and particles sintered together can be exchange-coupled; the 'no alignment
step' rule had treated a loose powder as a fixed random sample. Post-hoc rule v5: loose powders
measured without fixing are open (random if held fixed; up to 1 if they turn), decided by
whether the particles were immobilised during the measurement.

## Freeze v5 before opening batch 6 (heldout3: 81 records from sources not used before)

Rules v5 = v4 plus: loose powders and powders whose fixing is not reported are open (random if
fixed, up to 1 if the particles turn in the measuring field); composites with a soft or second
magnetic phase (MnO2, CoFe2O4, Fe) are not tests; substituted M-type ferrite labels map to the
hexaferrite phase; pressing/SPS of platelets is shape orientation (MO); textured hexaferrite
films are tested along the normal only. Account: HDDR from hydride thermodynamics, physical
TGG, textured-parent variant angles (no case prediction of batches 1-5 depends on them).

    sha256 data/curated/cases_v5.yaml  0c7a5fd10bc6252afbdd1f9a24101db312ed32ec76170e6650a0e308afe15583
