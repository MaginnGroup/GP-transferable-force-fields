# Running the calibration and validation workflow

Start with the [simulation input guide](../reproducibility/README.md) if you
want to use the supplied force fields in a new simulation. This page covers
running the calibration campaign: simulate candidate parameters, fit Gaussian
processes, optimize parameters, and validate the final force fields.

## Software and job setup

The simulation campaign used GROMACS 2022.1 with PLUMED 2.8, Cassandra 1.3.1,
signac-flow, mBuild/Foyer, and MoSDeF-Cassandra. The historical environment
exports are [hfcs-fffit.yml](../hfcs-fffit.yml) and
[hfcs-fffit.yaml](../hfcs-fffit.yaml). They are starting points, not a tested
single installation for all software versions. See the
[environment guide](../reproducibility/environments/README.md) for separate
analysis environments. Full molecular-simulation reruns have not been tested
with this documentation update.

The scripts use [signac](https://signac.io/) job workspaces and signac-flow
operations. An initializer creates jobs from molecule, temperature and parameter
settings; `project.py status` reports their state; `project.py run` executes an
operation locally; `project.py submit` submits it to a configured scheduler.
The `qsub` examples below use the original PBS cluster setup. Adapt scheduler
and resource settings for your machine before submitting.

Use a separate working copy for a new campaign so saved results remain
available for comparison. Inspect molecule selections, temperatures, sample
counts, iteration limits and executable paths before initializing jobs. Existing
analysis folders affect iteration selection: decide whether to start a new
campaign or extend an existing one. Do not simply run the distributed campaign
settings as a single-job example.

The saved workspaces lack some coordinates, restarts and Cassandra fragment
inputs needed to resume the original jobs. The simulation input guide identifies
the available templates and missing inputs. Activate a suitable simulation
environment before using the commands below. Each code block starts at the
repository root; parentheses keep directory changes local to that block.

## 1. Liquid-density iterations

[init_LD.py](../Build_GPs/init_LD.py) initializes density jobs using the molecule
data, parameter bounds and candidate samples. The `LD` operation group creates
force fields and systems, minimizes/equilibrates them, runs production and
extracts properties.

```sh
(cd Build_GPs && python init_LD.py)
(cd Build_GPs/ld_iters && python project.py status)
(cd Build_GPs/ld_iters && python project.py submit -o LD --bundle=12 --parallel)
```

After jobs complete, analyze the results and select candidates for subsequent
iterations using [post_ld.py](../Build_GPs/post_ld.py):

```sh
(cd Build_GPs && qsub -N postLD ../submit_job post_ld.py)
```

The campaign also uses [gen_lhs_samples.py](../Build_GPs/gen_lhs_samples.py) to
generate Latin hypercube samples. Check its settings and the initializer's
expected input files before starting a new sampling design.

## 2. Surface-tension iteration and GP fitting

The surface-tension project is named `vle_iters` in the source tree. Review the
selected LD candidates before initializing its jobs. The `IFT` group includes
packing, minimization, equilibration, density estimation, slab preparation,
interface simulations and property extraction.

```sh
(cd Build_GPs && python init_IFT.py)
(cd Build_GPs/vle_iters && python project.py status)
(cd Build_GPs/vle_iters && python project.py submit -o IFT)
```

After jobs complete, [post_vle.py](../Build_GPs/post_vle.py) processes the results
and fits/analyzes the surface-tension models:

```sh
(cd Build_GPs && qsub -N postVLE ../submit_job post_vle.py)
```

The saved DMF surface-tension GP uses a linear mean over all input features;
the historical fitting helper applies a masked linear mean to all six solvents.
Use the [saved-model settings and refit guide](../reproducibility/analysis/README.md#saved-model-semantics-and-historical-helpers)
when reproducing the reported models. Running the historical helper unchanged
does not reproduce DMF's saved mean-function specification.

## 3. Parameter optimization

[init_opt_at.py](../Opt_ES/init_opt_at.py) initializes the optimization jobs.
`OptGenFF` generates non-dominated candidates and optimizes the objective with
the fitted density and surface-tension GPs.

```sh
(cd Opt_ES && python init_opt_at.py)
(cd Opt_ES/opt_at_params_new && python project.py status)
(cd Opt_ES/opt_at_params_new && python project.py submit -o OptGenFF --bundle=12 --parallel)
```

After jobs complete:

```sh
(cd Opt_ES && qsub -N postOptAnaly ../submit_job post_analysis_opt.py)
```

## 4. Final simulation validation

Review the selected force-field vectors and state points in
[init_gemc_val.py](../Opt_ES/init_gemc_val.py) and
[init_ift_val.py](../Opt_ES/init_ift_val.py) before creating validation jobs:

```sh
(cd Opt_ES && python init_gemc_val.py)
(cd Opt_ES && python init_ift_val.py)
```

Validation projects ending in `_no_opt` use Base parameters; projects ending
in `_opt` use GP-Optimized parameters. For example, the GP-Optimized projects
are submitted as follows. Complete and inspect each Cassandra stage before
submitting the next stage; these are separate job submissions, not a command
that waits for all stages to finish.

```sh
(cd Opt_ES/gemc_val_opt && python project.py status)
(cd Opt_ES/gemc_val_opt && python project.py submit -o NPT_NVT --bundle=6 --parallel)
# After initialization finishes:
(cd Opt_ES/gemc_val_opt && python project.py submit -o EQ --bundle=6 --parallel)
# After equilibration finishes:
(cd Opt_ES/gemc_val_opt && python project.py submit -o PROD --bundle=6 --parallel)

# Interface validation is a separate project:
(cd Opt_ES/ift_val_opt && python project.py submit -o IFT)
```

Repeat with the corresponding `_no_opt` directories for Base validation. After
the simulations finish, process their results:

```sh
(cd Opt_ES && qsub -N postValAnaly ../submit_job post_analysis_val.py)
```

For figure regeneration, validation counts, saved-model diagnostics and grouped
validation without new molecular simulations, use the
[analysis guide](../reproducibility/analysis/README.md),
[validation guide](../reproducibility/validation/README.md), and
[DMF figure guide](../reproducibility/dmf/README.md).
