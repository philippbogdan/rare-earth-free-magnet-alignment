# Checks before regeneration

The plain Python scripts check analytic orientation geometry, magnetic-energy limits, variant selection and representative route rules. The account fixture uses deliberately simple values to isolate rule behaviour; it is not additional experimental evidence.

`expected.json` records SHA-256 values for the supplied results. The root reproduction command runs the scripts before regenerating and comparing those files. See [REPRODUCE.md](../REPRODUCE.md) for expected output and the distinction between a checksum match and scientific validation.
