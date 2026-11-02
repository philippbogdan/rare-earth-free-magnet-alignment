# Data behind the figures

| CSV | Columns and interpretation |
|---|---|
| `mnal_c_axes_fibre_parent.csv` | `parent`, `variant`, `x`, `y`, `z`, `weight`: product easy axes from a perfect parent fibre |
| `mnal_c_axes_random_parent.csv` | The same fields for random parent orientations |
| `mnal_c_A_vs_parent_spread.csv` | `parent_fibre_spread_deg`, `A_along_fibre`: alignment along the parent fibre at each angular spread; `random` is a separate reference |
| `measurements_vs_ceiling.csv` | `version`, `status`, `case`, `phase`, `route`, `predicted_class`, `A_ceiling`, `A_base`, `A_measured_lo`, `A_measured_hi`, `tol`, `ceiling`: saved sample comparisons |

Axes are folded onto the upper hemisphere because an easy axis is unchanged by reversing its sign. Each parent contributes equally; `weight` is the variant fraction within that parent. The two populations have different sample counts, so normalise them separately before comparing densities.

In the spread CSV, the numeric coordinate is the standard deviation of each Cartesian component of the rotation vector applied to a parent orientation. It is a simulation parameter, not a measured rocking-curve width. The `random` row comes from a separate random-parent population and has no angular coordinate.

For the held-out figure, select `version == v1` and `status == heldout`. Known cases recur across versions and must not be counted as fresh tests. Measurement endpoints relative to the tested boundary are `A_measured_lo - A_ceiling - tol` and `A_measured_hi - A_ceiling - tol`. Preserve intervals and out-of-range values. The atlas comes from [account.csv](../../results/account.csv), whose Nd2Fe14B row is a rare-earth reference.
