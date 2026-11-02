# Physics and frozen predictions

`phases.yaml` links each physical input to a source record or a labelled estimate. `records.yaml` holds additional cited inputs. `universe.yaml` groups the candidate phases and routes named in the selected literature.

The `cases*.yaml` files preserve the frozen predictions for held-out and known examples. Their bytes, including old header references to documents absent from this snapshot, are retained for hash checks. Use the [validation record](../../docs/validation.md) to interpret versions and exclusions; use the safe example in [REPRODUCE.md](../../REPRODUCE.md) to build current predictions into a separate file.
