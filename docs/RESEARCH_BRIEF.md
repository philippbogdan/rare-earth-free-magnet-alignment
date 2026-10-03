# Research brief for literature collection

Project (see PROBLEM.md, GOAL.md, docs/METHOD.md): for each rare-earth-free hard magnetic
phase, which processing routes can align its easy axes, the remanence ceiling of each
route, and why. The account is a physics engine: it takes *inputs* (crystallography,
magnetism, transformation paths of each phase) and predicts the alignment factor
A = <|cos θ|> and Br = Js·φ·A for every route. *Outcomes* (measured remanence/texture of
real samples) are only used to test it and are sealed until predictions are frozen.

## Rules (all mandatory)

* Work only inside the working directory `/Users/phil/magnet-align-frontier`; do not
  read or search local files elsewhere. The internet is open.
* Contact nobody, post nothing, submit no forms, no logins.
* Never invent values, quotes or references. Every value carries a verbatim `quote`
  from a page you actually fetched. If you cannot find something, say so (in the
  `gaps` list) rather than filling it.
* Follow `data/SCHEMA.md` exactly (YAML lists of records). Keep files small (text only;
  no PDFs in the repository: if you download PDFs, put them in `.tmp/` and delete them
  when done).
* Inputs go to `data/raw/<topic>_inputs.yaml`. Outcomes go to
  `.tmp/sealed/<topic>_outcomes.yaml`, except the four known cases (columnar field-aged
  AlNiCo, printed MnAl(C), hot-extruded MnAl-C, field-pressed sintered NdFeB) which go to
  `data/raw/known_cases_<topic>.yaml`.
* **Do not put outcome values (remanence, Mr/Ms, texture of samples) in your final reply.**
  Your final reply lists: files written, number of records by quantity/route, the
  sources used (citation + DOI), what could not be found (gaps), and any physics
  insight about mechanisms (orientation relationships, variant counts, formation
  temperature vs Curie temperature) — mechanisms are inputs, so they may be reported.
* Do not edit PROBLEM.md, GOAL.md or files of other topics. Do not commit (the lead
  commits).

## Access tips

* Metadata/abstracts: `https://api.crossref.org/works/<doi>` (title, authors, year; some
  abstracts), `https://api.openalex.org/works/doi:<doi>` (abstract as
  `abstract_inverted_index`; open-access locations), Semantic Scholar
  `https://api.semanticscholar.org/graph/v1/paper/DOI:<doi>?fields=title,abstract,openAccessPdf`.
* Full text: OSTI.gov (DOE-lab papers: Ames, ORNL, PNNL, NREL, LANL...), arXiv, PMC,
  institutional repositories, open-access journals. MDPI pages block scripts but PDFs
  load from `https://mdpi-res.com/d_attachment/<journal>/<journal>-<vv>-<nnnnn>/article_deploy/<journal>-<vv>-<nnnnn>.pdf`.
  `pdftotext` (poppler) and the Python venv `.venv/bin/python` (numpy, scipy, pyyaml,
  pypdf) are available.
* Publisher pages (ScienceDirect, Springer, Wiley) often refuse scripts; then use the
  abstract from Crossref/OpenAlex/Semantic Scholar and mark `evidence: abstract`.
* Datasheets (Arnold Magnetic Technologies, Electron Energy, MMPA standard 0100-00,
  TDK, Proterial/Hitachi Metals, Shin-Etsu, ...) are acceptable with `evidence: datasheet`.
* Use the YAML loader (`.venv/bin/python -c "import yaml; yaml.safe_load(open(...))"`) to
  check your files parse before finishing.
