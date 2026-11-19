# Finding archived results

Paths below are relative to the repository root. Molecule codes are MeOH
(methanol), EG (ethylene glycol), Gly (glycerol), DMF (N,N-dimethylformamide),
DMSO (dimethyl sulfoxide), and DEC (diethyl carbonate).

## Final force-field validation

For each molecule `MOL`, look under
`Opt_ES/analysis/at_00/MOL/ExpVal/opt_res/`:

| Path | Contents |
| --- | --- |
| `ms_val_no_opt/ms_data.csv` | Base force-field simulation results, with replicate rows |
| `ms_val_opt/ms_data.csv` | GP-Optimized force-field simulation results, with replicate rows |
| `ms_val_no_opt/` and `ms_val_opt/` | Property tables and validation plots |
| `best_per_run.csv` | Best parameter vector from each optimization run |
| `pareto_info.csv` | Non-dominated parameter candidates |
| `unique_best_set.csv` | Selected optimized parameters |
| `prop_pred/` | GP predictions and simulation/reference comparisons |

[Opt_ES/analysis/AT-0/](../Opt_ES/analysis/AT-0/) contains aggregate
property-comparison outputs. The names `no_opt` and `opt` refer to Base and
GP-Optimized, respectively; `at_00` identifies the atom-typing scheme used
in these results.

Use the [validation tables](../reproducibility/validation/README.md) for usable
replicate counts, per-state-point means and variation across replicate means.
Counts vary by property and temperature. These tables distinguish missing
values from zero physical uncertainty and document the methanol exclusion.

## Calibration data and GP models

[Build_GPs/analysis/](../Build_GPs/analysis/) groups results by molecule, with
`ld_iters/` for liquid-density iterations and `vle_iters/` for the surface-tension
iteration. These directories include candidate parameter tables, simulation
outputs, train/test partitions, saved models and diagnostics. The `vle_iters`
directory name is retained from the workflow; this calibration stage fits
surface tension, rather than performing the final Cassandra VLE validation.

The [GP analysis guide](../reproducibility/analysis/README.md) describes saved
model settings, grouped validation, sensitivity analyses and regeneration
commands. Its [reference directory](../reproducibility/analysis/reference/)
contains the numerical comparison results. For the five-candidate DMF figure,
use the [dedicated plotting guide](../reproducibility/dmf/README.md).

These archives provide processed numerical results and selected simulation
inputs, not complete trajectories or every generated restart file. For a new
simulation, start with the [force-field files](../Opt_ES/FF_params/README.md)
and [simulation input guide](../reproducibility/README.md).
