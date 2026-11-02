# Reproduce the account

From the repository root, in a Python environment with [requirements.txt](requirements.txt) installed, run:

```sh
python reproduce.py
```

Allow about 40 seconds on a laptop CPU, an estimate that depends on the machine. The calculation is offline and uses the included data. Install the pinned dependencies first if needed:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

## What is checked

The command runs every `tests/test_*.py` script, rebuilds `results/account.csv`, scores the frozen predictions into `results/scores_v1.json` through `scores_v5.json`, and writes `results/summary.json`. Geometry, magnetic-energy calculations and representative route rules are checked before regeneration.

It then compares every regenerated result file with the SHA-256 in [tests/expected.json](tests/expected.json). A differing file produces `MISMATCH` and a failing exit status. This verifies numerical and textual reproducibility of the included results; it does not revalidate the papers online or recreate the historical sealing of outcomes. Figures and figure CSVs are generated separately.

Expected output:

```text
tests passed
wrote results/account.csv (558 cells, 31 phases including Nd2Fe14B reference), results/scores_v1.json to scores_v5.json and results/summary.json
account: 558 cells (31 phases x 18 routes, including Nd2Fe14B reference): full 103, capped 30, none 196, open 229
MnAl(C) from a perfectly textured parent: easy axes at 52.2 deg, A_max 0.6133 against 0.5 for random grains
v1 held-out: 441 measurements, within 394, straddling 34, violated 13
v2 held-out: 4 measurements, within 4, straddling 0, violated 0
v3 held-out: 89 measurements, within 89, straddling 0, violated 0
v4 held-out: 40 measurements, within 34, straddling 4, violated 2
v5 held-out: 60 measurements, within 60, straddling 0, violated 0
MATCH: every regenerated file is identical to the expected one
```

The checksums describe the supplied result files. `python reproduce.py --write-expected` replaces those checksums with hashes of the files currently on disk; it does not run tests or regenerate anything. Use that option only after independently checking an intentional change. A mismatch should be investigated, not hidden by replacing the baseline.

## Re-check the public sources

The numerical inputs are the cited facts in the YAML records. The source checker fetches DOI metadata, abstracts and the linked pages or PDFs to compare their wording and numbers. Install Poppler's `pdftotext` for PDF extraction, then run:

```sh
python tools/verify_sources.py data/raw/*.yaml data/outcomes/*.yaml data/curated/records.yaml
```

This online step is separate from the short offline run and can take much longer. It updates `data/verification.json` and caches fetched material under `.cache/fetched/`. Optional source-text files can be placed under `.cache/texts/`; use only text extracted from the cited publications, never printouts or exports of the records being checked. No source cache is bundled. Local `quote_where` filenames in the supplied register identify past source texts, not bundled files. Current access failures, shortened excerpts and differences between PDF text and page images can require manual checking. The [evidence guide](docs/data.md#unchecked-quotations) lists the account's eight unchecked quotations.

To inspect the closest cached passage for a record, pass its ID to `python tools/quote_context.py RECORD_ID`. Record a checked passage and the reason for any automatic mismatch in `data/verification_manual.yaml`. To retry only failed quotations, use `--recheck-failed` before the same file list in the verification command. Rechecking sources may change the evidence fields in the summary, so a subsequent checksum mismatch can be legitimate and still needs inspection.

## Conditions to predictions

The value-free `data/heldout/*.conditions.yaml` records carry the sample's process, composition and measuring direction. `align.cases.build()` maps them to phase and route rules using current curated phase data. It restores case IDs from the manifest, links same-sample saturation records by ID, and attaches the exposure status and paired comparisons. Outcome values are used later by the scorer.

The existing `data/curated/cases*.yaml` files are the frozen predictions used for historical scores. Keep them unchanged. Current phase data and rules include revisions informed by earlier outcomes, so rebuilding now cannot recreate every earlier freeze. To inspect predictions made by the current rules without overwriting a frozen file:

```sh
PYTHONPATH=src python - <<'PY'
from pathlib import Path
import yaml
from align import cases

predictions = cases.build()
Path('results/cases_current.yaml').write_text(
    yaml.safe_dump(predictions, sort_keys=False, allow_unicode=True, width=110)
)
PY
```

`cases.build_known()` applies the current rules to the known examples. Avoid running `python -m align.cases` on this snapshot: its default writer targets an existing versioned filename. The [validation record](docs/validation.md) supplies batch order, frozen hashes, exposure rules and every revision's limitations.

## Figure data and rendering

Regenerate the orientation populations, texture-spread curve and measurement comparisons with:

```sh
PYTHONPATH=src python -m align.figdata
```

Allow about 2 seconds on a laptop CPU, an estimate. The CSVs appear in [figures/data](figures/data/README.md); the atlas uses `results/account.csv` directly. Their columns and interpretation are documented there. Regenerate result files before exporting measurement comparisons, because `figdata` reads the saved scores.

Render the figures from the repository root:

```sh
python tools/render_figures.py
```

Rendering requires only the pinned Python requirements; Pillow, used for RGB export, is installed with Matplotlib. The font preference is Helvetica Neue, Helvetica, Arial, then DejaVu Sans; another machine may reproduce the content without identical pixels. On a laptop, the static command completed in 1.8 s wall time (1.5 s CPU time), and two runs produced byte-identical PNGs.

The command writes `mnal-inheritance.png`, `alignment-atlas.png`, `heldout-ceilings.png` and `mnal-texture-spread.png` under [figures](figures/README.md). Successful rendering ends with `Wrote four static alignment figures to figures/.` The renderer checks the hero geometry against the summary, includes every account cell and matches each plotted held-out interval to its saved score.

To trace a row in the held-out figure, `tools/render_figures.py::heldout_rows()` returns the plotted case IDs and intervals in order: ceiling first, interval midpoint relative to the boundary second, then case ID. The gap in the image separates ceilings below 0.97 from those at or above it; it does not omit measurements.
