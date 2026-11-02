# Source checks and figures

`verify_sources.py` compares DOI metadata and supporting excerpts with the cited online material, storing results in the verification register. `quote_context.py` shows the closest passage in a cached source when an automatic match needs inspection. Neither tool changes the measured values.

The online checks need source access and Poppler's `pdftotext` for PDFs. [REPRODUCE.md](../REPRODUCE.md) gives the full file list, cache locations and instructions for recording manual checks.

`render_figures.py` reads the included account and figure CSVs to draw the inheritance geometry, verdict atlas, held-out comparison and parent-spread curve as static PNGs. Rendering uses Matplotlib; commands and font fallbacks are in [REPRODUCE.md](../REPRODUCE.md#figure-data-and-rendering).
