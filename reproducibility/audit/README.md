# Reproducibility audit for the public cleanup PR

The audit checks retained public inputs, simulation results, saved GP models
and added analyses. It does **not** reproduce the molecular simulations.
GROMACS/Cassandra preprocessing and execution, trajectory reanalysis, and
recovery of missing simulations were not run on this machine.

## Recorded run

[audit.json](audit.json) records the commands, environment, full tested source
revision, UTC time, exit codes and 586 original CSV/model hashes. The run used
a fresh `git archive` of `c21f0fac23af6039b9133ce7b35e5f7ebc4109f4`, with
`PYTHONPATH` unset. Scripts resolved their inputs within that export; no sibling
private development or audit repository was required. The original input
hashes were unchanged after all commands.

The lightweight Python 3.12 environment was installed fresh from the pinned
[requirements](../environments/analysis-py312.txt). GP calculations used the
separate existing, tested Python 3.10/GPflow environment; a fresh GP installation
was not validated. The short revision argument recorded by the driver is
resolved to the full archived commit in `audit.json`.

| Check | Result | Evidence |
|---|---|---|
| Twelve variant exports, 132 LJ values, charges, units, scaling and two input examples | Pass | [inputs.txt](inputs.txt) |
| All 380 property/state-point replicate counts and raw/processed row comparisons | Pass; zero discrepancies | [validation.txt](validation.txt), [counts and summaries](../validation/README.md) |
| Corrected DMF figure | Pass; all 25 data rows match | [dmf.txt](dmf.txt), [figure provenance](../dmf/provenance.json) |
| Saved grouped predictions and numerical SVD results | Pass; 30 folds, 48 local blocks, 24 global decompositions | [saved-result report](saved_result_verification.json) |
| Saved kernel/mean/noise settings | Pass for all 12 final GPs | [model settings](saved_model_settings.json) |
| Full-precision vectors, split overlap and multistart inventory | Pass | [inventory.txt](inventory.txt) |
| Saved-model local predictions/Jacobians | Pass; exactly match references | [local.txt](local.txt), [comparison report](gp_rerun_comparison.json) |
| Sampled-domain GP gradients, singular values and rankings | Pass; numerical agreement within 4.55e-13 | [global.txt](global.txt), [comparison report](gp_rerun_comparison.json) |
| Five-fold grouped refits, three starts per fold | Pass; all 30 selected folds converged | [grouped_refit.txt](grouped_refit.txt) |
| All 130 regenerated result-file comparisons | Pass; maximum absolute difference 4.547473508864641e-13 | [comparison report](gp_rerun_comparison.json) |

The numerical comparator uses relative tolerance 1e-8 and absolute tolerance
1e-10. Saved-result arithmetic uses relative and absolute tolerances of 1e-10.
The DMF reference-data comparison uses 1e-12. These tolerances are explicit in
the scripts and reports. Runtime seconds and optimizer prose are excluded from
numerical equivalence checks. Local results and grouped-refit predictions,
selected parameters, noise values and metrics matched exactly in this run;
the small difference arose in global SVD arithmetic.

Of 90 individual grouped-fit starts, 83 converged. Each fold selected a
converged candidate using the documented training-loss rule. TensorFlow emitted
Cholesky warnings at rejected optimization trials; these are retained in the
logs. No fit failure or warning was silently counted as a successful individual
start. The final selected fits and their comparisons passed.

## Checks exposed by using a clean export

The first export lacked fold-2 reference files because the repository's old
`*old2*` ignore pattern matched their names. That failed the saved-result check.
Narrow exceptions now include the reference CSV/JSON files, and the complete
export passed on repetition. The failed exploratory run is retained in the
separate audit repository rather than presented as a passing result here.

[input_counterchecks.json](input_counterchecks.json) records that the checker
rejects a missing manifest variant and refuses optimized Python (`-O`), which
would disable its assertions. A disposable parameter-package copy regenerated
byte-identically. Fixed-width PDB padding and copied simulator formatting are
preserved; `.gitattributes` excludes only those archived formats from trailing
space diagnostics. Derived template whitespace was cleaned without changing
any parsed settings.

To repeat the audit, use [run_audit.py](../run_audit.py) with the
[documented environments and command](../environments/README.md#run-the-audit-as-one-command).
Subsequent report/documentation commits do not change the audited numerical
code or original results. Four historical plotting launchers now resolve their
working directories from their own file locations instead of a collaborator's
cluster directory. Only directory selection changed; their full workflow plots
were not rerun and are not covered by the numerical pass claims here. The PR diff preserves all original simulation CSVs
and saved GP pickles. The corrected DMF data are derived outputs with explicit
provenance, not replacements for historical CSV fields.

## Scope and remaining publication steps

This PR targets public `main` and retains its prior author changes, including
the merged SVD PR and later force-field exports. It adds the input package,
usable-count documentation, corrected DMF plotting command, historical parser
explanation, final grouped/local analysis package, tested-environment records,
and independent numerical checks. It fixes README navigation, stale figure
numbers, the unavailable branch reference, and the DMF/saved-model distinctions.

The usable replicate counts stay in the public repository, matching the
approved manuscript reference. Generated Cassandra execution inputs, complete
restart histories, original DMF block outputs, and the excluded methanol
trajectory/diagnostic were not recovered. Their absence is documented; files
were not fabricated as substitutes.

Merge, release/tag creation and verification of final stable URLs are separate
publication steps. R1-15/R1-16 must be reconciled with the final integrated
revision. R1-16 should describe variant-specific parameter exports, stage
settings, selected examples and generators, disclose reconstructed Base
exports and unavailable execution/restart assets, and avoid claiming that all
executed simulation inputs were recovered. No manuscript or reviewer-response
source is changed by this PR.
