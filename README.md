# Rare-earth-free magnet alignment

![MnAl(C) geometry and ideal alignment: inherited parent texture puts the product easy axes 52.2 degrees off the build direction, giving A = 0.61 versus 0.50 for random grains and 1 for perfect alignment.](figures/mnal-texture.png)

**A textured parent need not make an aligned magnet.** For the printed MnAl(C) texture, even perfect orientation inheritance gives an ideal remanence-to-saturation ratio of about 0.61. This assumes equal variant populations and no field or stress bias. [Calculation and assumptions](results/account.csv).

This project models how processing sets magnetic alignment across 31 phase entries and 18 routes, using published crystal structures, transformations and magnetic properties. It calculates conditional ceilings and identifies measurements that could refute them. It does not predict coercivity, cost or a working manufacturing process.

![All 558 phase-route combinations. Squares mark 103 high alignment ceilings, circles 30 geometric caps, crosses 196 routes without an alignment gain, and question marks 229 unresolved cases. Nd2Fe14B is included as a reference.](figures/route-map.png)

Each cell has its mechanism, assumptions and proposed test in the [data table](results/account.csv). A square means the model permits an alignment factor of at least 0.95 under ideal conditions. Some routes concern films or single crystals. The [findings](FINDINGS.md) carry the details.

![Five separate held-out batches: v1 has 394 within-ceiling measurements, 34 conversion-straddling cases and 13 clear violations; v2 has 4 within; v3 has 89 within; v4 has 34 within, 4 straddling and 2 violations; v5 has 60 within.](figures/validation.png)

These checks ask whether published measurements exceed a frozen ceiling. Passing does not establish an accurate point prediction. The [report](REPORT.md) includes failures, conversion corrections and results after rule changes. [Scoring data](results/) and [method](docs/METHOD.md).

## Reproduce

```sh
python3 -m venv .venv
.venv/bin/pip install numpy==2.4.6 PyYAML==6.0.3 matplotlib==3.11.1
.venv/bin/python tools/make_figures.py
```

The figures read the committed result files. Run `./regenerate.sh` to rebuild the results and reports from the source records. Figure inputs, hashes and transformations are in [provenance.json](figures/provenance.json). [SVG figures](figures/) are also included.
