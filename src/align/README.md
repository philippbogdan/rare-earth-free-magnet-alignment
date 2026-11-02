# One set of orientation rules

`orient.py` and `variants.py` define lattices, rotations, point groups and transformation variants. `engine.py` applies texture, inheritance, selection and body rotation. `physics.py` supplies magnetic-energy bounds; `routes.py` and `account.py` combine them into the route account.

`data.py` and `curate.py` resolve evidence records and phase descriptions. `cases.py` converts sample conditions into current predictions; `scoring.py` converts measurements and tests ceilings; `evaluate.py` scores the frozen files. `accuracy.py` provides a lower-level conditional-interval interface and is not the frozen-file scorer. `validate.py` checks reference orientation geometry.

`cli.py` writes the account, historical scores and `summary.py` headline fields. `figdata.py` exports figure CSVs. To add a material, supply sourced phase data and use the same operators. Start with the [method](../../docs/method.md); [REPRODUCE.md](../../REPRODUCE.md) shows how to run without replacing frozen predictions.
