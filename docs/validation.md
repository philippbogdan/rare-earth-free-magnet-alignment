# What the frozen predictions survived

The frozen v1 ceilings held for 394 of 441 convertible held-out measurements, or 89.3%. Another 34 straddled the scoring boundary; 13 exceeded it even under the most favourable conversion. These are tests of an upper bound. A ceiling near perfect alignment is easy to satisfy and does not establish an accurate prediction of a sample's remanence.

The present account uses the current physical rules. Historical scores use the saved predictions for each version. No parameter was fitted to the outcome values, but failures did inform later rules; those revisions receive a held-out score only on a later batch.

## Measurements, boundaries and exclusions

Each prediction contains `A_hi`, the sample's ceiling, and `A_base`, its alignment without the route's mechanism. The measured quantity becomes an interval $`[m_{lo},m_{hi}]`$. Remanence ratios convert directly; remanence in tesla requires saturation polarisation and phase fraction. The saturation value comes from the same sample where linked, a sourced grade interval, or the phase value with a 5% range. Supported texture measurements are converted geometrically.

```math
t=\max(0.03,2u),\qquad b=A_{hi}+t.
```

Here $`u`$ is the stated uncertainty expressed in alignment units. `within` means $`m_{hi}\leq b`$; `straddling` means $`m_{lo}\leq b<m_{hi}`$; `violated` means $`m_{lo}>b`$. In the score files these labels are `ok`, `possible violation` and `violated`. Uncertainty enters the conversion interval as well as the tolerance, so this is a deliberately permissive ceiling test.

Of the 772 records classed as held out in the v1 batches, 228 were not alignment tests, 88 could not be converted and 15 were not saturated. The 441 scored measurements are not a count of independent papers or materials; several measurements can come from one study.

Of those scored measurements, 221 had ceilings at or above 0.97. Among the 220 below 0.97, there were 193 within, 14 straddling and all 13 violations. This more restrictive subset is the better check that the mechanism imposes a useful limit.

The separate mechanism tests refuted 12 of 130 `none` predictions and showed alignment beyond the baseline in 254 of 290 `acts` predictions. Failure to demonstrate an acting mechanism is not a refutation of its ceiling. All 10 scored paired comparisons had the expected relation. Open cases retain or reject each stated hypothesis separately.
A post hoc check using the same-sample saturation in tesla that v1 had omitted when converting eight remanence measurements to alignment left all ceiling results unchanged and raised supported `acts` predictions from 254 to 260 of 290; 254 of 290 remains the frozen held-out result.

## Batch order and frozen files

Outcomes were encrypted before prediction. The available records contained process conditions with outcome values removed. Predictions were fixed before the corresponding batch was opened. [MANIFEST.json](../data/heldout/MANIFEST.json) preserves batch membership and the original outcome hashes; the condition files preserve the inputs to the sample rules.

| Rules | Newly opened batch | Condition file stem | Held-out scored | Within | Straddling | Violated |
|---|---|---|---:|---:|---:|---:|
| v1 | 1 and 2 | Initial topic collections | 441 | 394 | 34 | 13 |
| v2 | 3 | `mn4n_outcomes` | 4 | 4 | 0 | 0 |
| v3 | 4 | `heldout2_b_outcomes` | 89 | 89 | 0 | 0 |
| v4 | 5 | `heldout2_a_outcomes` | 40 | 34 | 4 | 2 |
| v5 | 6 | `heldout3_outcomes` | 60 | 60 | 0 | 0 |

The v3, v4 and v5 batches used sources absent from all earlier batches. Scores on an already opened batch are `post hoc`, even when a newer version improves them. These rows must not be pooled as a test of a single unchanged model.

| Version and file under `data/curated/` | SHA-256 |
|---|---|
| v1, `cases.yaml` | `eb6ee47de1f561568ff545347e7336fb226fee9f205a0a82637db39892252eca` |
| v2, `cases_v2.yaml` | `fed2017ea7a6c8dc82c6ddf87cdef7f4300bca417e3ae55113b0eac10bc391af` |
| v3, `cases_v3.yaml` | `b203b12b3f0b75f6935c9cde3aa4ee783f2c0fb79f406b5557cf7bef55716588` |
| v4, `cases_v4.yaml` | `b8ce592f93b14a5bc774115a5cfe590a3fd1c647e7332174071f77fe0123c91b` |
| v5, `cases_v5.yaml` | `0c7a5fd10bc6252afbdd1f9a24101db312ed32ec76170e6650a0e308afe15583` |

These prediction files remain byte-identical. Their original header references to an earlier protocol description and problem statement are retained as part of the hashed record. Those documents are not included; this page explains the record. Hashes identify files; the timing of each freeze and opening rests on the recorded batch sequence. Published outcome quotations were shortened, so current outcome-file hashes cannot authenticate the originally sealed bytes. The manifest is a record of those earlier hashes, not a successful integrity check of the current outcome files.

## Exposure and known cases

Three exposures are excluded from held-out scoring:

- A FePd field-ordering result and its treatment temperatures were visible before freezing. All related FePd field-ordering cases were marked `exposed`.
- Numeric outcome fields survived in an early MnAl condition export. The affected MnAl records, including remanence, energy-product and anisotropy information, were marked `exposed`. The condition files were then restricted to process fields.
- Process text for a MnBi field-sintering sample contained a sibling sample's remanence ratio. Both samples in that series were marked `exposed`.

The exact exclusions are the IDs and prefixes in `cases.EXPOSED`. In v1, 75 exposed measurements were scored separately: 63 within, 5 straddling and 7 violated. They contribute nothing to the held-out result.

The four known examples are field-aged columnar AlNiCo, printed MnAl(C), hot-extruded MnAl-C and field-pressed sintered Nd-Fe-B. They were visible when rules were formulated. Their 107 scored measurements give 106 within, one straddling and no violations; this is a fidelity check. A grade-to-saturation mapping error in the AlNiCo examples was corrected before v1 was frozen.

## v1 failures and the v2 revision

All 13 v1 violations were diagnosed after opening. These explanations are interpretations of failures, not additional validation.

| Failure group | Observation against ceiling | Revision |
|---|---|---|
| FePt epitaxial film | 0.90 against 0.50 | Interpret the film normal as the alignment axis |
| Meteorite FeNi | 0.62 against 0.33 | Include bias from the ferromagnetic parent's own magnetisation |
| Polytwinned CoPt crystal | 0.82 against 0.58 | Include a coupled magnetic unit as an alternative to independent variants |
| Isotropic HfCo7 cluster film | 0.55 against 0.50 | Flag possible exchange enhancement |
| Co nanorod assemblies | Seven cases, 0.54 to 0.70 against 0.50 to 0.53 | Leave unmeasured rod arrangement open; allow aligned bundles to turn as bodies |
| Fe-Cr-Co with Mo | 0.88 to 1.06 against 0.83 | Do not impose the path sourced at 3 wt% Mo on a sample with 2 wt% Mo |
| Binary MnAl | 0.67 against 0.61 | Retain both transformation readings when treatment order and temperature are unspecified |

Converted intervals can extend above unity; they are not clipped to hide inconsistent conversion or measurement assumptions. The case IDs and full intervals are in `scores_v1.json`.

Further v2 changes concerned exposed or already opened cases. Transformation shear resolves ties left by a strain-only rule for stress-biased MnAl-C single crystals. AlNiCo-5's rod direction was settled using opened observations: the rods can follow the field rather than remain on the nearest parent axis. These are learned rule choices and are not independently supported by rescoring the same cases.

Mn4N added the nearest easy-direction rule for a cubic phase. The account also states the coupled-network alternative for unbiased spinodal alloys, while retaining the uniaxial reference for alignment comparisons. The only fresh v2 test was the Mn4N batch.

## v3: a fresh batch and a conversion failure

The next rules added phase aliases, remanence-unit conversions and additional same-sample saturation links. They recognised possible exchange enhancement in bonded or composite samples, left the internal structure of MnBi particles milled from ribbon open, and kept epsilon-Fe2O3 variant-family fractions separate from alignment. Account changes treated epitaxy geometrically and made template-growth and direct-formation prerequisites explicit.

The frozen result was 89 within out of 89. After opening, a conversion defect was found: outside AlNiCo, the frozen specifications could use phase saturation instead of an available same-sample value, understating alignment. The sensitivity calculation covers 24 records. With same-sample saturation, 19 ceilings hold, one is exceeded and four have no scored ceiling result. The exceedance is an isotropic bonded Nd-Fe-B print with exchange-enhanced remanence, a contribution excluded from alignment but missed by the frozen classification. `scores_v3.json` retains the frozen result alongside `sensitivity_js`; the correction does not become a new held-out test.

## v4: corrected conversions, then loose-powder failures

The next freeze linked same-sample saturation for every phase and excluded recognised isotropic bonded Nd-Fe-B and hard/soft mixtures from alignment tests. It separated film-normal measurements from in-plane measurements affected by shape anisotropy and applied the geometric cap to unbiased inclined-fibre films.

Further sample rules covered near-diagonal AlNiCo grain-growth textures, stress-biased growth with unmeasured texture, Fe-Cr-Co rolling between ageing steps, and ambiguous grade text in sample IDs. The fresh batch gave 34 within, four straddling and two violations.

Both failures were loose L1_0 FePt nanoparticle powders measured without fixing: 0.58 and 0.61 against 0.50. Rotation in the measuring field or coupling within sintered clusters can raise remanence. An absence of a deliberate alignment step does not establish that the measured powder is a fixed random assembly.

## v5: measurement mobility becomes a condition

Loose powders, and powders whose fixing is unspecified, became open: immobilised random particles retain the random baseline, while mobile particles can turn towards the field. The deciding measurement is whether the particles were fixed during magnetometry.

The rules also excluded additional composite samples, mapped substituted M-type ferrites, recognised pressing or spark-plasma sintering of platelets as shape orientation, and restricted textured hexaferrite film tests to the normal. Account changes refined hydride thermodynamics, template-growth limits and textured-parent variant angles.

All 60 scored measurements in the fresh batch were within their ceilings. This supports that frozen version on that batch; it does not remove the earlier failures, conversion defect or untested manufacturing prerequisites.
