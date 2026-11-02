# Ceilings and their tests

`account.csv` contains 558 phase and route combinations, including Nd2Fe14B as a rare-earth reference. Every row names its mechanism and a measurement that tests it; open rows give the deciding quantity and alternative ceilings. Read `axis`, `chain` and `prereq` alongside `A_max`.

`scores_v1.json` through `scores_v5.json` retain individual measurement intervals, scoring tolerances and statuses. The fresh held-out subset differs by version. `scores_v3.json` also preserves the same-sample saturation sensitivity check. `summary.json` collects the headline counts and frozen prediction hashes. [Validation](../docs/validation.md) explains what the scores establish and where they failed.
