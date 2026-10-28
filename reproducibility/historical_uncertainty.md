# Interpreting historical surface-tension uncertainties

The original ST-iteration processing code read the GROMACS `Err.Est.` column
into a variable called `vars_est` in
[`calc_blk_avg_energ`](../Build_GPs/vle_iters/project.py), then calculated
`max(sqrt(vars_est))`. GROMACS reports block error estimates in that column,
not variances. The extra square root therefore changes the meaning and units
of the stored `surf_tens_unc` field. The corresponding simulation averages
were calculated separately and are unaffected by this transformation.

Original CSV values are preserved so historical analyses remain inspectable.
Do not interpret those stored within-run ST uncertainty fields as standard
deviations or confidence intervals without checking their parsing lineage.
Do not apply the correction to density uncertainties or replicate standard
deviations.

## Recoverable DMF quantity

For the five-candidate DMF data used in the manuscript, the code converted the
GROMACS pressure-times-length output from bar nm to mN/m and divided by two
for the two interfaces. The numerical conversion factor is 0.05. If `E` is the
largest original GROMACS block error estimate and `u` is the stored value,

```text
u = 0.05 * sqrt(E)
corrected_max_error_mN_per_m = 0.05 * E = 20 * u**2
plotted_halfwidth_mN_per_m = 1.96 * corrected_max_error_mN_per_m
```

The original estimates used three through ten blocks. Taking the maximum and
square root commutes for these nonnegative values. The original DMF per-block
output files were not recovered, so this inverse recovers only the maximum
error estimate, not the eight individual estimates or a complete block series.
The 1.96 factor preserves the manuscript plotting convention; it does not
demonstrate calibrated 95% confidence coverage.

The [corrected DMF workflow](dmf/README.md) checks all 25 plotted points against
the approved numerical table. The inverse is specific to this established
field/conversion pipeline; check units and code history before applying it to
another dataset.

## GP noise and final-validation bars are separate

The ST GP noise calculation in
[`get_Eotvos_scale`](../Build_GPs/utils/models.py) fits an Eotvos-style
temperature relation to the simulated surface-tension means, using molar
volumes from a density GP. It aggregates the mean absolute residuals (mean
across parameter vectors for glycerol, maximum for the other molecules),
squares the aggregate and scales the variance to GP output units. It does not
use the stored within-run `surf_tens_unc` fields. Thus the parser issue does
not change the simulation averages or this residual-derived GP noise.

The final-validation property comparison bars instead use variation across
independent replicate means, as documented in
[validation/README.md](validation/README.md). They must not be corrected using
the DMF inverse transformation. The corrected five-candidate DMF figure uses
within-run block error estimates and is a different plot.

White-kernel variance, Gaussian likelihood variance, and numerical Cholesky
jitter are separate quantities. The historical helper's comments call the
likelihood initialization `jitter`; that label does not make it the numerical
jitter used in a covariance factorization. Saved-model checks are needed to
establish the actual trained/fixed parameters of each reported GP.
