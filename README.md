# GP-transferable-force-fields

Simulation inputs, saved results, and calibration code for **Gaussian Process
Optimization of Transferable Force Fields: A Case Study for Six Common Solvents
from Liquid Density and Surface Tension Data**, by Montana N. Carlozo,
Alexander W. Dowling, and Edward J. Maginn.

The study fits Lennard–Jones parameters using Gaussian-process models of liquid
density and surface tension for six solvents. Start with the parameter files
below if you want to use the force fields; the calibration software is needed
only if you want to fit or analyze models yourself.

## 1. Simulation input files

**[Force-field files and usage notes](Opt_ES/FF_params/README.md)** are the main
entry point. Each solvent includes both a **Base** model (selected from the
sampled candidates) and a **GP-Optimized** model (obtained by optimizing the GP
objective). The study's recommended variant differs by solvent:

| Solvent | Input files | Recommended variant |
| --- | --- | --- |
| Methanol (MeOH) | [Base](Opt_ES/FF_params/MeOH/Base/) · [GP-Optimized](Opt_ES/FF_params/MeOH/GP-Optimized/) | GP-Optimized |
| Ethylene glycol (EG) | [Base](Opt_ES/FF_params/EG/Base/) · [GP-Optimized](Opt_ES/FF_params/EG/GP-Optimized/) | GP-Optimized |
| Glycerol (Gly) | [Base](Opt_ES/FF_params/Gly/Base/) · [GP-Optimized](Opt_ES/FF_params/Gly/GP-Optimized/) | Base |
| N,N-dimethylformamide (DMF) | [Base](Opt_ES/FF_params/DMF/Base/) · [GP-Optimized](Opt_ES/FF_params/DMF/GP-Optimized/) | Base |
| Dimethyl sulfoxide (DMSO) | [Base](Opt_ES/FF_params/DMSO/Base/) · [GP-Optimized](Opt_ES/FF_params/DMSO/GP-Optimized/) | Base |
| Diethyl carbonate (DEC) | [Base](Opt_ES/FF_params/DEC/Base/) · [GP-Optimized](Opt_ES/FF_params/DEC/GP-Optimized/) | GP-Optimized |

Each variant provides:

- `system.top` — GROMACS topology.
- `ff.xml` — Foyer force-field definitions.
- `species1.mcf` — Cassandra connectivity and force-field parameters.
- `species1.pdb` — single-molecule coordinates and atom ordering.

The [simulation input guide](reproducibility/README.md) covers units, combining
rules, nonbonded settings, [GROMACS stage templates](reproducibility/templates/gromacs/),
[Cassandra generator settings](reproducibility/templates/cassandra/settings.json),
and [two saved GROMACS job examples](reproducibility/examples/).

Base parameter files are reconstructed from the saved validation parameter
vectors. The topology molecule count must match your coordinates; the
single-molecule PDB is not a liquid or slab configuration. The examples are
input collections, not complete simulation restarts. Read the input guide
before starting a simulation.

To download the repository and check parameter-file consistency (Python 3 only;
no simulation engine required):

```sh
git clone https://github.com/MaginnGroup/GP-transferable-force-fields.git
cd GP-transferable-force-fields
python3 reproducibility/check_inputs.py
```

## 2. Archived results

| Results | Where to look |
| --- | --- |
| Final Base and GP-Optimized simulation results | [Opt_ES/analysis/](Opt_ES/analysis/); see the [results guide](docs/results.md) for file locations |
| Usable replicate counts, exclusions, means and uncertainty-bar definitions | [Validation results](reproducibility/validation/README.md) |
| Density and surface-tension calibration iterations, training data and saved GP models | [Build_GPs/analysis/](Build_GPs/analysis/) |
| Grouped GP validation and parameter-sensitivity results | [Saved analyses](reproducibility/analysis/README.md) and [numerical reference files](reproducibility/analysis/reference/) |
| Five-candidate DMF figure data and plotting command | [DMF figure](reproducibility/dmf/README.md) |

These are processed results, parameter records and selected job inputs. Complete
trajectories and all simulation restart files are not included. Replicate
availability varies by state point; use the validation tables for counts.

## 3. Code to rerun the workflow

The [workflow guide](docs/workflow.md) explains the simulation software,
signac job management, scheduler setup and commands for:

1. Generating candidate parameters and running liquid-density iterations.
2. Running surface-tension simulations and fitting GP models.
3. Optimizing parameters and validating the final force fields.

[Build_GPs/](Build_GPs/) contains candidate-generation, simulation and GP-fitting
code. [Opt_ES/](Opt_ES/) contains parameter optimization and final validation.
The [analysis guide](reproducibility/analysis/README.md) provides commands to
check saved numerical results or rerun GP diagnostics without new molecular
simulations. Use the [environment guide](reproducibility/environments/README.md)
to choose the required software for each task.

## Manuscript and related papers

This repository accompanies the six-solvent manuscript named above, currently
under revision. Cite the manuscript and the repository commit or release used
in your work. A publication DOI will be added when available.

The workflow builds on these earlier papers:

- Wang, N.; Carlozo, M. N.; Marin-Rimoldi, E.; Befort, B. J.; Dowling, A. W.;
  Maginn, E. J. [Machine Learning-Enabled Development of Accurate Force Fields
  for Refrigerants](https://doi.org/10.1021/acs.jctc.3c00338).
  *Journal of Chemical Theory and Computation* **2023**, *19*, 4546–4558.
- Carlozo, M. N.; Wang, N.; Dowling, A. W.; Maginn, E. J.
  [Development of accurate transferable hydrofluorocarbon refrigerant force
  fields using a machine learning and optimization approach](https://doi.org/10.1039/D5DD00537J).
  *Digital Discovery* **2026**, *5*, 1650–1674.

The force fields supplied here are for the six solvents listed above; the
refrigerant studies are methodological background.

## Credits and contact

Supported by NSF award EEC-2330175 (ERC EARTH). Computing resources were provided
by the University of Notre Dame Center for Research Computing. Montana Carlozo
acknowledges the Department of Education GAANN fellowship, grant P200A210048.

For questions, contact Montana Carlozo (mcarlozo@nd.edu), or
[open an issue](https://github.com/MaginnGroup/GP-transferable-force-fields/issues).
