# Data schema and evidence rules

Two kinds of records, kept apart:

* **Inputs** (`data/raw/<topic>_inputs.yaml`, visible): properties of phases and
  transformations that the account uses to make predictions.
* **Outcomes** (`.tmp/sealed/<topic>_outcomes.yaml`, sealed until predictions are frozen):
  measured remanence, Mr/Ms, or texture of magnets/samples made by a specific route.
  Exception: the four known cases of PROBLEM.md (field-aged columnar AlNiCo, printed
  MnAl(C), hot-extruded MnAl-C, field-pressed sintered NdFeB) go to the visible file
  `data/raw/known_cases.yaml`.

## Evidence rules

1. Never invent a value, a quote or a reference. If not found, write nothing.
2. Every record carries `quote`: the verbatim sentence(s) from the source that contain the
   value (copy exactly, including units). Paraphrase goes in `note`, never in `quote`.
3. `evidence` is one of: `full_text` (you read the paper), `abstract` (only the abstract was
   accessible), `datasheet` (manufacturer/standard), `review` (value quoted by a review or
   handbook, not the original), `book`.
4. `url` is the page or file where you read the quote. `doi` must resolve on
   api.crossref.org; copy `title`, `year`, first author from Crossref.
5. Record what the measurement is: one sample or many, which direction, at what field
   (is it saturating?), at what temperature. No claim stronger than its data.
6. Values read off a figure: `locator: "Fig. 3b (read from plot)"` and an uncertainty.
7. Input quotes must not carry outcome values (remanence/texture of a route) of other
   samples; put those in the sealed outcome file instead.

## Input record (`<topic>_inputs.yaml`)

```yaml
- id: mnbi_Tc_xyz2014          # unique, lowercase
  phase: LTP-MnBi              # phase name as used in the topic file header
  quantity: Tc                 # see list below
  value: 630                   # number, or text for qualitative quantities
  unit: K
  uncertainty: 10              # if stated or estimated (say which in note)
  conditions: "bulk, Mn50Bi50" # composition/temperature/etc.
  source: {doi: 10.xxxx/yyy, title: "...", first_author: "...", year: 2014,
           url: "https://...", evidence: abstract, locator: "abstract"}
  quote: "verbatim text containing the value"
  note: ""
```

Quantities (use these names): `Js` (μ0Ms at 300 K, T), `Ms_mass` (Am²/kg), `K1` (MJ/m³),
`Tc` (K), `density` (g/cm³), `lattice` (a, c in Å), `structure`, `easy_axis`,
`T_stab` (highest temperature the phase survives, with time scale), `formation_path`
(parent phase, mechanism), `T_formation` (K, range), `OR` (orientation relationship,
verbatim), `n_variants`, `tetragonality` (c/a of product, and strain relative to parent),
`twins` (type, misorientation, effect on c-axis), `growth_direction` (preferred
solidification growth direction of the phase or its parent), `primary_crystallisation`
(does the phase or its parent crystallise from the melt as discrete crystals; liquidus
range), `single_variant_particles` (can single-crystal/single-variant particles be made,
grain/twin size), `deformation` (slip/twinning/DRX mechanisms, hot-workability window),
`delta_chi` (paramagnetic susceptibility anisotropy), `field_effect_on_formation`,
`stress_effect_on_formation`, `phase_fraction_limit`, `other`.

## Outcome record (`<topic>_outcomes.yaml`, sealed)

```yaml
- id: mnbi_fp_xyz2016a
  phase: LTP-MnBi
  route: FP                    # FP FS DS AM FA DS+FA AM+FA RX+FA SC+FA SA DS+SA HD EP TGG
                               # or 'isotropic' (no alignment step) or 'other: ...'
  process: "jet-milled powder aligned in 2 T, warm compacted 300 C, 1 h"
  conditions: {field_T: 2, T_process_K: 573, parent_texture: "none",
               magnet_axis: "field direction", density_rel: 0.93}
  quantity: Br                 # Br (T) | Mr_over_Ms | Js_same_sample | texture_metric
  value: 0.0                   # (placeholder)
  unit: T
  uncertainty: 0.0
  measured_along: "alignment direction"
  measurement: "VSM, max field 9 T, 300 K, demagnetisation corrected?"
  n_samples: 1                 # or "grade/datasheet range"
  source: {doi: ..., title: ..., first_author: ..., year: ..., url: ..., evidence: abstract,
           locator: "Table 2"}
  quote: "verbatim"
  note: ""
```

Texture metrics (`quantity: texture_metric`): give the kind in `unit` —
`pole_figure_max_mrd`, `orientation_factor`, `lotgering_f00l`, `rocking_fwhm_deg`,
`variant_fraction`, `A_reported` — and the reflection/pole.
