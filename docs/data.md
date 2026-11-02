# Evidence behind a cell

The account uses 476 sourced records. Supporting quotations are confirmed at the source for 468: 352 online, 68 in saved source text and 48 by inspection of source text or page images. Eight quotations remain unchecked. The verification register records resolving DOIs for 443 of these records. A further 116 values are labelled estimates; none sets a ceiling. These counts are in [summary.json](../results/summary.json), under `evidence`.

These are recorded checks, not a guarantee that every source remains accessible online. In [verification.json](../data/verification.json), a local `quote_where` entry names the saved source text; it is a provenance label, not a repository link. The texts are not distributed. The online checker can produce different results as source access changes, and the published quotation excerpts are shorter than some passages used for the recorded checks.

## Unchecked quotations

Seventy earlier quote checks matched a printout of the records themselves rather than source text. Those matches prove nothing and have been withdrawn with `quote_ok: false`. Eight of the affected records supply the account, all in [alnico_fecrco_inputs.yaml](../data/raw/alnico_fecrco_inputs.yaml):

| Record ID | Input |
|---|---|
| `alpha1_Tc_kaneko1977` | Curie point of the Fe-Co-rich precipitate |
| `fecrco_Mo_anisotropic_sugimoto1991` | Field selection in Mo-containing Fe-Cr-Co |
| `fecrco_Mo_needs_100_homma1980` | Parent-axis requirement in Mo-containing Fe-Cr-Co |
| `fecrco_SC_orientation_indep_okadathesis` | Field response of differently oriented single crystals |
| `fecrco_elastic_beats_field_Mo_ikuta1983` | Competition between elastic and field bias |
| `fecrco_elastic_small_in_ridge_ikuta1983` | Elastic bias near the spinodal ridge |
| `fecrco_shape_anisotropy_sugimoto1987` | Shape anisotropy of the precipitates |
| `fecrco_vs_alnico_selection_okada1978` | Variant selection in Fe-Cr-Co compared with AlNiCo |

Their values and citations stand, but their quotations are unchecked. A retained `quote_score` or `numbers_found` from a withdrawn match is not confirmation; read `quote_ok` and any entry in [verification_manual.yaml](../data/verification_manual.yaml). Confirming these passages requires the cited source text or page image. When rechecking, populate the optional source-text cache only with material taken from the cited publications, never exports of these records.

## Follow a result to its source

Start with a row of [account.csv](../results/account.csv). Its `refs` are record IDs in [raw](../data/raw/README.md) or [curated records](../data/curated/records.yaml). The [phase descriptions](../data/curated/phases.yaml) show which records supply each structure, easy axis, transformation and processing prerequisite.

For a validation point, use its `case` ID in `results/scores_v*.json` to find the measured sample in [outcomes](../data/outcomes/README.md). The corresponding frozen file in [curated](../data/curated/README.md) supplies the ceiling and conversion rule. `topic` identifies the collection; `status` determines whether the comparison is held out, exposed, known or post hoc.

## Record fields

| Field | What a reader can check |
|---|---|
| `id`, `phase`, `quantity` | Stable record identity, material and measured or described property |
| `value`, `unit`, `uncertainty` | The value and its stated or estimated uncertainty |
| `conditions` | Composition, temperature and other conditions needed to interpret it |
| `source` | DOI or URL, title, first author, year, evidence level and passage or figure locator |
| `quote` | Source wording around the supporting value |
| `note` | Conversions, interpretation, qualifications and transcription details |

Outcome records also describe the route, process, measurement direction, measurement method and sample count where available. A property of a phase is kept separate from the remanence or texture of a particular processed sample. The known-case files use the outcome schema even though they sit under `raw/`.

`full_text`, `abstract`, `datasheet`, `review` and `book` describe the evidence actually available. A review citation is secondary evidence. Values read from a plot require a locator and uncertainty. A missing value remains missing; an estimate must be identified as such and explain its basis. Source wording belongs in `quote`; interpretation belongs in `note`.

Physical magnitudes in a phase description resolve through `{ref: record_id}`. `{est: value, why: explanation}` marks an estimate. Plain crystallographic indices and symmetry names are conventions. This separation lets a reader replace a physical input without changing the orientation operators.

The coverage list includes all 24 candidate phases and all 16 route classes named in the selected literature. DS+SA and RX+SA extend that route list. This is coverage of the cited selection, not proof that every possible material or route is represented. The account includes Nd2Fe14B as a rare-earth reference; some supporting outcome records also concern substituted ferrites containing rare earths. Inspect composition when selecting a strictly rare-earth-free measurement set.

## Quotations and frozen records

Published excerpts contain at most 300 source characters; the `[...]` omission markers are additional to that count. Shortening the excerpts changes the outcome-file bytes, so those files do not match the original sealed hashes in [MANIFEST.json](../data/heldout/MANIFEST.json). The manifest preserves the batch membership, case IDs and hashes of the files as sealed. Those original sealed files are not included in this snapshot.

The frozen `cases*.yaml` files are unchanged, including their header comments. Some headers refer to an earlier protocol description and problem statement that are not included here. [Validation](validation.md) gives the complete reader-facing account of the frozen batches and exclusions. The preserved headers are part of the hashed record, not instructions to follow missing documents.

## Rights and dependencies

The code, curation, original text and figures use [CC0 1.0](../LICENSE). Numerical facts are accompanied by source citations; the short quotations remain attributed to their publications. This dedication does not relicense the source papers.

Dependencies are installed separately: NumPy uses BSD-3-Clause, PyYAML uses MIT, and Matplotlib uses its PSF-based Matplotlib licence. Their code is not bundled. Source PDFs, publisher pages and local text caches are not included.
