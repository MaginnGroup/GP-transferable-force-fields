# Tested analysis environments

Use separate environments for retained-data plotting/checks and historical
TensorFlow/GPflow models. No simulation engine is required for these commands.

## Lightweight checks and plots

Python 3.12 with [analysis-py312.txt](analysis-py312.txt) was installed in a
fresh isolated environment. The file pins the resolved plotting dependencies,
including NumPy 2.3.5, pandas 2.2.3 and Matplotlib 3.11.2.

```sh
python3.12 -m venv .venv-analysis
.venv-analysis/bin/python -m pip install -r reproducibility/environments/analysis-py312.txt
.venv-analysis/bin/python reproducibility/check_inputs.py
.venv-analysis/bin/python reproducibility/validation_counts.py --output-dir /tmp/es-ffo-validation
.venv-analysis/bin/python scripts/plot_dmf_corrected.py --output-dir /tmp/es-ffo-dmf
.venv-analysis/bin/python reproducibility/analysis/verify_saved_results.py --output-dir /tmp/es-ffo-saved-check
```

The first two commands require only Python's standard library. The figure and
saved-analysis checks use the third-party environment. If Matplotlib cannot
write its default cache, set `MPLCONFIGDIR` to a writable temporary directory.

## Historical GP models and refits

The existing tested GP environment uses Python 3.10.21, NumPy 1.26.4,
TensorFlow-macos 2.15.0, TensorFlow Probability 0.23.0 and GPflow 2.9.1 on macOS
arm64. [gp-py310.txt](gp-py310.txt) pins its direct analysis dependencies;
[gp-installed-packages.json](gp-installed-packages.json) records installed
distribution versions. The model loading, local derivatives and full grouped
refits were actually executed in that environment. A fresh GP installation
was not tested. The non-macOS TensorFlow requirement is provided as a setup
starting point, not evidence of Linux numerical equivalence.

With a compatible Python 3.10 environment, run the commands in the
[analysis README](../analysis/README.md). GPflow 2.9.1 emits a deprecation warning
for `pkg_resources`; the recorded setuptools version remains below 81. The
lightweight and GP environments deliberately use different NumPy versions.
Do not install the lightweight NumPy 2 pin over the TensorFlow 2.15 environment.

## Simulation prerequisites and external reference software

The simulation campaign used GROMACS 2022.1 with PLUMED 2.8, Cassandra 1.3.1,
and the signac-flow/mBuild/Foyer/MoSDeF-Cassandra workflow. Those engines were
not installed, compiled or executed during this audit. The historical
`hfcs-fffit.yml`/`.yaml` files are preserved; they do not establish one verified
environment for every reported software version. Original timing logs and some
generated execution/restart inputs are unavailable.

Retained numerical reference values are already included. The commands above
do not call licensed REFPROP or regenerate experimental reference data.
Regenerating those sources is a distinct task with the applicable source
software and permissions. No claim of molecular-simulation or trajectory
reproduction follows from successful analysis checks.

## Run the audit as one command

```sh
.venv-analysis/bin/python reproducibility/run_audit.py --output-dir /tmp/es-ffo-audit
```

To include saved-model settings, inventory, local/global derivatives and all
30 three-start grouped refits, append `--gp-python /path/to/gp-env/bin/python
--refit`. The audit keeps simulation/result/model inputs unchanged and hashes
them before and after the run. It records the tested revision, environment,
commands, exit codes, logs and explicit not-run simulation checks. For a source
export without `.git`, supply `--revision FULL_COMMIT_SHA`. Use a fresh output
directory for each audit.
