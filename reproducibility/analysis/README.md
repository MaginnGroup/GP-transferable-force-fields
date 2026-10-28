# Additional GP validation and sensitivity analyses

These commands use committed public data and saved GPs. They do not require a
manuscript checkout, the private development repository, REFPROP, signac jobs,
GROMACS or Cassandra. Write regenerated results to a separate directory.
Use the [lightweight or GP environment](../environments/README.md) as indicated.

## Verify saved numerical results (lightweight environment)

```sh
python reproducibility/analysis/verify_saved_results.py --output-dir /tmp/es-ffo-saved-check
```

This independently verifies all 745 out-of-fold predictions from 30 grouped
folds, their MAE/RMSE/nominal-95% predictive coverage, complete-vector fold
exclusion, the seed-42 assignments, and three-start selection receipts. It also
recalculates ranks for all 48 saved local Jacobian blocks and checks 24 public
sampled-domain SVD decompositions for orthonormal right vectors and consistent
singular-value-weighted rankings. No fitted model is loaded by this command.

The fixed [reference data](reference/) come from the retained v13 audit package
in `dowlinglab/ESS-FO-Audit`, checkpoint `2256590`. They include the final
**three-start** grouped results, not the earlier exploratory one-start run.
`grouped_validation_3starts/protocol.json` gives the fit protocol. The
manuscript-parameter reference contains the rounded approved values; the exact
validation vectors remain in the original public CSVs. `model_specs.json`
records constants and their source hash. Property bounds and spans preserve
the original float32 calculations; subtracting JSON endpoints in float64 can
change the physical-unit conversion slightly.

## Regenerate local diagnostics from saved models (GP environment)

```sh
python reproducibility/analysis/inspect_saved_models.py --output-dir /tmp/es-ffo-model-settings
python reproducibility/analysis/local_estimability.py inventory --output-dir /tmp/es-ffo-inventory
python reproducibility/analysis/local_estimability.py local --output-dir /tmp/es-ffo-local
```

The inventory checks full-precision vectors against the rounded manuscript
parameters, split overlap and multistart results. The local command loads the
final density RQ and ST Matern32 models, evaluates predictions and derivatives
at Base/GP-Optimized vectors and the five fitting temperatures, then writes
property-specific and combined Jacobians, singular values, ranks at three
thresholds, predictions, and finite-difference derivative diagnostics. Fixed
hydroxyl-hydrogen LJ coordinates and the temperature column are excluded from
parameter ranks. The combined design has ten response rows.

Objective-diagonal weighting uses the residual scale with its 2% floor. The
alternative measurement-plus-GP diagonal weighting adds the saved GP
prediction variance to that scale; it is a separate diagnostic, not a new fit.
The matrices are surrogate derivatives, not evidence of structural
identifiability or parameter confidence intervals. Rank and conditioning depend
on scaling, thresholds, parameterization and GP fidelity.

## Refit grouped validation (GP environment)

```sh
python reproducibility/analysis/grouped_GP_validation.py --restarts 3 --maxiter 2000 --output-dir /tmp/es-ffo-grouped
```

Use a fresh output directory. All five temperatures of a complete scaled LJ
vector stay in one test fold; the 25 vectors for each of the first five
molecules and 24 for DEC give 125 and 120 eligible rows, respectively. The
Matérn-3/2 ST kernel and saved mean-function form are retained. Each fold uses
three deterministic starts and chooses the smallest training loss among
converged fits, or among all fits if none converges. Held-out ST responses do
not initialize or select hyperparameters. A seed-42 permutation assigns folds;
restart draws use their recorded integer seeds.

The fixed white variance is recomputed from **training-fold ST vectors** using
the Eötvös residual rule. The original density surrogate supplies molar
volumes; it is not refitted within each ST fold. Thus this is a test of ST
prediction conditional on the available density model, not a joint holdout of
all density and ST information. The Gaussian likelihood variance remains
1e-5 in scaled ST-output units. Coverage uses `predict_y` variance and is
nominal predictive coverage, not replicate-SD coverage or a parameter interval.

Compare the newly generated results independently:

```sh
python reproducibility/analysis/verify_gp_reruns.py --local-results /tmp/es-ffo-local --grouped-results /tmp/es-ffo-grouped --output-dir /tmp/es-ffo-gp-check
```

The comparator uses relative tolerance 1e-8 and absolute tolerance 1e-10,
checks row identities, and excludes runtime seconds and optimizer prose.
The tested existing macOS/GPflow environment reproduced all 82 compared files
exactly, including 30 sets of selected fitted parameters. Other numerical
backends can converge differently; a failed comparison must be investigated
and reported rather than silently replacing the manuscript reference results.

## Saved-model semantics and historical helpers

- All six saved final density models have an RQ kernel, no White term (`c=0`),
  a linear mean, and trained Gaussian likelihood variance.
- All six saved ST models have Matern32 plus fixed White variance and fixed
  Gaussian likelihood variance 1e-5. Five use an epsilon/temperature-masked
  linear mean; DMF uses the full linear mean including sigma inputs.
- White-kernel variance, likelihood variance and numerical factorization jitter
  are distinct. See the [historical uncertainty explanation](../historical_uncertainty.md).
- The generic `Build_GPs/utils/models.py` helper currently contains a different
  glycerol-density configuration, and its generic custom ST mean does not
  encode the saved DMF exception. Running that helper unchanged does not
  reproduce every saved model. The commands above explicitly load the final
  saved models or preserve their forms during the added grouped refits.
- Existing `make_svd_sensitivity.py` describes **sampled-domain** gradients at
  different parameter vectors. Those rankings differ from **local final-FF**
  Jacobian ranks. Neither establishes unique physical LJ parameters.
- Per-molecule LHS plots summarize atom-count-weighted parameter sums; their
  correlations do not establish independence of all individual coordinates.

Only load the trusted pickles supplied by the repository. Pickles execute
Python code and are not an interchange format for untrusted downloaded models.
The public audit distinguishes model/analysis reproduction from molecular
simulation reproduction, which was not attempted on this machine.
