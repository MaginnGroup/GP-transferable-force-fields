# Gaussian Process Optimization of Transferable Force Fields: A Case Study for Six Common Solvents from Liquid Density and Surface Tension Data
Authors: Montana N. Carlozo, Alexander W. Dowling, and Edward J. Maginn
<!-- Introduction: Provide a brief introduction to the project, including its purpose, goals, and any key features or benefits. -->
## Introduction
**ES-FFO** is a repository used to rapidly calibrate the LJ parameters of electrolyte solvent force fields given experimental data. The key feature of this work is using machine learning tools in the form of Gaussian processes (GPs) which allow us to cheaply estimate the results of a molecular simulation given temperature state points and thermophysical property data. This workflow combines the approaches of Wang et al., 2023, J. Chem. Theory Comput. and Carlozo et al., 2025, Digital Discovery to create one repository to build GP models and optimize parameters for electrolyte solvent FFs. **ES-FFO-Public** is based on the "for_sub" branch of private repository **ES-FFO**.

**NOTE**: We use Signac and signac flow (`<https://signac.io/>`) to manage the setup and execution of each workflow. These instructions assume a working knowledge of that software. Use the files in this public checkout; a private branch is not required for the retained-data checks below.

## Force fields and simulation inputs

Both Base and GP-Optimized parameter files for all six molecules are indexed in [Opt_ES/FF_params](Opt_ES/FF_params/README.md), with explicit provenance and recommended-variant labels. [Reproducibility instructions](reproducibility/README.md) provide GROMACS stage inputs, two archived examples, Cassandra generator settings, units, 1-4 scaling and the limits of the locally archived inputs. Run `python3 reproducibility/check_inputs.py` to check their consistency without running simulations.

## Start with a reproducibility task

| Task | Entry point | Software needed |
|---|---|---|
| Select Base or GP-Optimized force-field inputs | [Variant inventory](Opt_ES/FF_params/README.md) | None to inspect; Python for consistency checks |
| Inspect usable replicate counts and exclusions | [Validation results](reproducibility/validation/README.md) | Python standard library |
| Regenerate the corrected DMF figure | [DMF workflow](reproducibility/dmf/README.md) | NumPy, pandas, Matplotlib |
| Verify saved analyses or rerun GP diagnostics | [Additional analyses](reproducibility/analysis/README.md) | Lightweight or GP environment, depending on task |
| Set up new molecular simulations | [Simulation inputs and limitations](reproducibility/README.md) | Matching simulation engines and workflow dependencies |

See [tested environments](reproducibility/environments/README.md) and the
[reproducibility audit](reproducibility/audit/README.md). Historical simulation
CSVs are preserved. Pooled cross-molecule correlation plots are retained as
historical scripts but were omitted from the revised paper. The corrected
DMF command supersedes the old uncertainty-bar plotting convention.

## Citation
This work is currently in preparation and has not yet been published. Please cite the following when referring to this work:

Carlozo, M. N., Dowling, A. W., & Maginn, E. J. Gaussian Process Optimization of Transferable Force Fields: A Case Study for Six Common Solvents from Liquid Density and Surface Tension Data. In preparation. Cite the exact repository revision used alongside the manuscript; a stable release identifier will be added after review.
   
## Available Data

### Repository Organization
The repository is organized as follows: <br />
ES-FFO/ is the top level directory. It contains: <br />
1. .gitignore prevents large files from the signac workflow and plots from being tracked by git and prevents  tracking of other unimportant files. <br />
2. utils/ contains functions used by all pieces of this workflow including molecule experimental data and reference parameters <br />
3. block_average/ contains code for block averaging via the methods in H. Flyvbjerg and H.G. Petersen. Error estimates on averages of correlated data. J. Chem. Phys. 91:461-466, 1989.<br />
4. fffit/ contains another set of utility functions used by multiple aspects of this workflow. <br />
5. Build_GPs/ contains the workflow that runs and analyzes results from applying the methods modified from those in Wang et al. 2023 to create Base FFs <br />
6. Opt_ES/ contains the workflow that runs and analyzes results from applying the methods in Carlozo et al. 2025. to create GP-optimized FFs <br />
7. Opt_ES/ Also contains the molecular simulation results and analysis of the FF(s) developed in steps 5 and 6.  <br />
8. hfcs-fffit.yml is a historical environment export; tested retained-analysis environments are documented separately. <br />
9. submit_job is a shell script for submitting jobs to the cluster. <br />
10. create_analysis_figs.ipynb is a jupyter notebook for creating the final analysis figures including all figures from running the scripts in 11-17. <br />
11. make_corr_figs_all_molec.py is a python script for creating correlation figures for all molecules. These pooled cross-molecule plots are historical and are omitted from the revised manuscript. <br />
12. make_corr_figs_one_molec.py is a python script for creating correlation figures for individual molecules. Includes LHS weighted-parameter-sum correlations and per-molecule property correlations. <br />
13. make_GP_vs_sim_and_sens.py is a python script for creating figures comparing GP and actual predictions and generating the results for individual molecule sensitivity analyses. For the current complete-set and local analyses, use reproducibility/analysis/ instead. <br />
14. make_ift_val_figs.py is a python script for creating visualizations of the interfacial tension boxes for all ST simulations. <br />
15. make_param_comp_figs.py is a python script for creating parameter comparison figures. Produces parallel-coordinate parameter comparisons. <br />
16. make_pareto_comp_figs.py is a python script for visualizing differences between pareto-optimal LJ parameter sets in ST (and LD) iterations. For the corrected five-candidate DMF plot, use scripts/plot_dmf_corrected.py. <br />
17. make_svd_sensitivity.py is a python script for creating visualizations of the SVD sensitivity analysis. Produces sampled-domain singular values, right vectors and rankings; see reproducibility/analysis/ for interpretation. <br />
18. Claude_Audit/ contains the files generated by Claude to audit the experimental data in esolvs.py and generate reference-provenance and property-comparison data. <br />
Note that all workspace/ subdirectories below include one example signac workspace with some example input data and key simulation validation figures where appropriate.


### utils/
This directory contains:
1. molec_class_files/esolvs.py: A .py file which loads all electrolyte solvent data into a class
2. prep_ms_data.py: Contains scripts used to prepare data for GP training and get error data


### fffit/fffit/
This directory contains functions used by both workflows:
1. models.py: Functions related to building GP models
2. pareto.py: Functions for determining pareto efficient points
3. plot.py: Functions for plotting data
4. signac.py: Functions for extracting signac data
5. utils.py: Functions for scaling and unscaling data and shuffling/splitting GP data

### Build_GPs/
This directory contains all data related to running the GP building workflow, including:
1. analysis/mol/ld_iters: Results for liquid density (LD) iterations for each molecule "mol"
2. analysis/mol/vle_iters: Results for surface tension (ST) iterations
3. utils/*.py: Functions for analyzing the results of LD and ST iterations
4. ld_iters/ : The signac project directory for running the LD iterations including project.py, templates/, and workspace/ 
5. vle_iters/ : The signac project directory for running the ST iterations including project.py, templates/, and workspace/
6. gen_lhs_samples.py: A script for generating the LHS data necessary for LD and ST iterations
7. init_LD.py: A script to generate the signac project for LD iterations
8. init_IFT.py: A script to generate the signac project for ST iterations
9. post_ld.py: A script to analyze LD iteration data
10. post_vle.py: A script to analyze ST iteration data
11. get_CPU_hrs.py : A script to calculate the total CPU hours used by the workflow

### Opt_ES/
This directory contains all data related to running the optimization workflow of Carlozo et al., 2025, Digital Discovery including:
1. analysis/: Data from analyzing the results of the optimization workflow
2. utilsOpt/: Functions for running the optimization workflow
3. gemc_val_opt/: The signac project directory for running the gemc validation for GP-optimized FF including project.py, templates/, and workspace/ 
4. gemc_val_no_opt/: The signac project directory for running the gemc validation for the Base FF including project.py, templates/, and workspace/
5. ift_val_opt/: The signac project directory for running the interfacial tension validation for GP-optimized FF including project.py, templates/, and workspace/
6. ift_val_no_opt/: The signac project directory for running the interfacial tension validation for the Base FF including project.py, templates/, and workspace/
7. opt_at_params_new/: The signac project directory for running the optimization of FF parameters with updated GP models including project.py, templates/, and workspace/
8. FF_params/: The directory for compiling FF parameters for each molecule for the Base and GP-optimized FFs including ff.xml, species1.mcf, species1.pdb, and system.top files; use the explicit Base/ and GP-Optimized/ subdirectories
9. compress_gemc.py: A script to compress the gemc validation data
10. compress_ift.py: A script to compress the interfacial tension validation data
11. get_CPU_hrs.py: A script to calculate the total CPU hours used by the optimization and validation workflows
12. get_GPU_hrs.py: A script to calculate the total GPU hours used by the optimization and validation workflows
13. init_gemc_val.py: A script to initialize the gemc validation workflow
14. init_ift_val.py: A script to initialize the interfacial tension validation workflow
15. init_opt_at.py: A script to initialize the optimization workflow
16. post_analysis_opt.py: A script to analyze the results of the optimization workflow
17. post_analysis_val.py: A script to analyze the results of the validation workflows. Generates property-comparison plots
18. get_RD2_data.py: A script to extract select RD2 property data for each molecule for the Base and GP-optimized FFs. Generates simulation-design reference values.
19. all_properties.csv: A csv file containing select RD2 property data for each molecule for the Base and GP-optimized FFs. Simulation-design reference values.
20. compile_ff_params.py: Historical parameter-export helper. Use reproducibility/export_force_fields.py for the documented variant package.

#### Opt_ES/analysis/
This directory contains data from analyzing the results of the optimization workflow for each molecule "mol".:
1. at_00/mol/ExpVal/opt_res/jac_approx:  Includes jacobian approximations
2. at_00/mol/ExpVal/opt_res/hess_approx:  Includes hessian approximations
3. at_00/mol/ExpVal/opt_res/MAPD:  Includes Mean Absolute Percentage Deviation values
4. at_00/mol/ExpVal/opt_res/ms_val_opt:  Includes the error data and validation results for the GP-Optimized FF
5. at_00/mol/ExpVal/opt_res/ms_val_no_opt:  Includes the error data and validation results for the Base FF
6. at_00/mol/ExpVal/opt_res/prop_pred/:  Includes comparisons of GP, simulated, and experimental values for the GP-optimized and Base FFs
7. at_00/mol/ExpVal/opt_res/best_per_run.csv:  The best FF optimization results for each optimization run
8. at_00/mol/ExpVal/opt_res/pareto_info.csv:  The pareto-optimal FF optimization results for each molecule
9. at_00/mol/ExpVal/opt_res/unique_best_set.csv:  The unique LJ parameter sets identified for each molecule during optimization
10. AT-0/ms_val_no_opt/: Base FF property prediction plots for heat of vaporization and vapor pressure (h_p_vap.pdf), surface tension (surf_tens.pdf), and density (vle.pdf)
11. AT-0/ms_val_opt/: GP-optimized FF property prediction plots for heat of vaporization and vapor pressure (h_p_vap.pdf), surface tension (surf_tens.pdf), and density (vle.pdf)
12. AT-0/ms_val_opt_comp/: GP-optimized FF and Base FF property prediction plots. (h_p_vap.pdf, surf_tens.pdf, vle.pdf)
13. comp_err_data.csv: Data comparing the errors of different FFs
14. lit_err_data.csv: Data comparing the errors of the literature FFs
15. lit_ff_data.csv: Data comparing the properties of the literature FFs
16. lit_ff_data_w_no_opt.csv: Data comparing the properties of the literature FFs with the Base FFs
17. lit_ff_data_w_opt.csv: Data comparing the properties of the literature FFs with the GP-optimized FFs
18. lit_Hvap_est_w_no_opt.csv: Data comparing the heat of vaporization of the literature FFs and Base FFs
19. lit_Hvap_est_w_opt.csv: Data comparing the heat of vaporization of the literature FFs and GP-optimized FFs

#### Opt_ES/utilsOpt/
This directory contains all functions related to running the optimization workflow of Carlozo et al., 2025, Digital Discovery including:
1. utilsOpt/atom_type.py: Class definitions for atom typing schemes
2. utilsOpt/opt_atom_types.py: Other functions necessary for optimizing FF parameters and analyzing workflow results
3. utilsOpt/plot.py: Plotting functions for analyzing optimization results
4. utilsOpt/signac.py: Functions for extracting signac data

#### Claude_Audit/
This directory contains all functions related to the Claude Code Audit of this repository including:
1. scripts/: Includes relevant scripts generated by Claude Code
2. refprop_output/: Includes ethyleneglycol_refprop_montana.txt and methanol_refprop_montana.txt, the REFPROP property data for EG and MeOH
3. figures/: Figures (baseline_vs_exp_*.png/.pdf) comparing experimental data in esolvs.py to other data sources for each molecule
4.audit_data/: All data necessary for the files in scripts/ and all data generated from the files in scripts/
5. refs.bib: References used in the main text

#### Claude_Audit/audit_data
This directory contains all functions related to the Claude Code Audit of this repository including:
1. baseline_agreement_stats.py: Comparison of experimental data sources
2. digitized_litff.csv: The digitized values for the Vahid & Maginn, Jahn et al., and Gonzalez-Salgado & Vega sources
3. litexp_*.csv: The literature experimental data from various FF development sources
4. litff_*.csv: The literature FF data from various FF development sources
5. master_ff_predictions.csv: Predictions of property data for all investigated FFs
6. montana_reference_values.csv: The values from esolvs.py laid out in csv format
7. provenance_map.csv: Data from the sources of experimental data used in esolvs.py
8. provenance_table.tex: Historical reference-provenance table
9. recomputed_metrics.csv: RD1 and RD2 property-comparison metrics
10. reference_values_all.csv: Compiled reference values used for rulers R1 and R2
11. refprop_validation.csv: Comparison of esolvs.py and REFPROP values for EG and MeOH
12. rename_map.csv: Information related to data sources
13. student_ff_errors_repo.csv: Compiled error data from ES-FFO/Opt_ES/analysis/
14. student_ff_recompute: Error data generated by the Claude Audit. Retained property-comparison data
15. tables_6_7_reconstructed.tex: Historical comparison-table LaTeX
16. audit_table3_refdata.csv: Claude audit for the historical reference-data table

#### Claude_Audit/scripts
This directory contains all functions related to the Claude Code Audit of this repository including:
1. build_master_ff_predictions.py: Creates a consolidated list of all FF predictions
2. build_reference_values_all.py: Generates RD1 and RD2 comparison-reference values
3. extract_student_ff_values.py: Extracts FF ms data from ES-FFO/Opt_ES/analysis/
4. generate_provenance_table.py: Generates the reference-provenance table
5. generate_tables.py: Generates comparison-table LaTeX
6. plot_baseline_vs_experiment.py: Generates experimental data comparison figures
7. recompute_metrics_all.py: Recomputes error metrics for predictions. Generates reference-based comparison metrics
8. recompute_student_metrics.py:  Recomputes error metrics for predictions from the Base and GP-Optimized FFs with data from ES-FFO/Opt_ES/analysis/
9. audit_table3_refdata.py: Generates audit for the historical reference-data table
10. unit_value_audit.py: Cross-checks esolvs.py's expt_* dictionaries (the manuscript's reference values) against ruler R1
11. validate_refprop.py: Generates comparison of esolvs.py and REFPROP values for EG and MeOH


## Installation

For retained-data checks and figure/GP analysis, use the separate
[tested analysis environments](reproducibility/environments/README.md).
The old `hfcs-fffit.yml` and `hfcs-fffit.yaml` files document historical workflow
requirements; neither is a verified single environment for all reported
simulation and analysis versions.

The historical simulation campaign used GROMACS 2022.1 with PLUMED 2.8 and
Cassandra 1.3.1, plus signac-flow, mBuild/Foyer and MoSDeF-Cassandra. Those engines
were not installed or executed during the retained-results audit. The cluster
commands below are historical workflow guidance and require a suitable
simulation host, scheduler and complete inputs; the selected example collections
are not complete standalone restarts.

<!-- Usage: Provide instructions on how to use the project, including any configuration or customization options. Examples of usage scenarios can also be added. -->
## Usage

### Liquid Density Optimization (LD Iterations)
To run liquid density iterations, follow the following steps:
1. Initialize Signac workflow
   ```
     conda activate hfcs-fffit
     cd Build_GPs/
     python init_LD.py
   ```    
2. Run LD iterations
   - **Note: rm -r workspace/ signac_project_document.json signac.rc will remove everything and allow you to start fresh if you mess up**
  ```
     cd Build_GPs/ld_iters/
     python project.py submit -o LD --bundle=12 --parallel
   ```  
3. Check status a few times throughout the process
   ```  
     python project.py status
   ```   
Note: Step 2 operation LD runs multiple operations in series. Alternatively, the following can be run one at a time:
4. Alternative submission
   ```
     python project.py run -o create_forcefield
     python project.py submit -o create_system
     python project.py run -o fix_topology
     python project.py submit -o em_sim --bundle=12 --parallel
     python project.py submit -o nvt_eq_sim --bundle=12 --parallel
     python project.py submit -o npt_eq_sim --bundle=12 --parallel
     python project.py submit -o npt_prod_sim --bundle=12 --parallel
     python project.py submit -o calculate_props --bundle=12 --parallel
   ```
     
5. Extract densities, run GP optimization and get samples for the next iteration in the Build_GPs/ directory
   - **Note: vle_iters/ will be populated with parameter sets automatically once the termination criteria for LD iters is satisfied**
   ```
     qsub -N postLD submit_job post_ld.py
   ```       

### Surface Tension Optimization (ST Iterations)

To run surface tension iterations, follow the following steps:
1. Initialize Signac workflow
   ```
     conda activate hfcs-fffit
     cd Build_GPs/
     python init_IFT.py
   ```    
2. Run ST iterations
   - **Note: rm -r workspace/ signac_project_document.json signac.rc will remove everything and allow you to start fresh if you mess up**
  ```
     cd Build_GPs/vle_iters/
     python project.py submit -o IFT #--bundle=2 --parallel
   ```  
3. Check status a few times throughout the process
   ```  
     python project.py status
   ```   
Note: Step 2 operation IFT runs multiple operations in series. Alternatively, the following can be run one at a time:
4. Alternative submission
   ```
     python project.py run -o create_forcefield
     python project.py submit -o create_system
     python project.py run -o fix_topology
     python project.py submit -o em_sim #--bundle=2 --parallel
     python project.py submit -o nvt_eq_sim #--bundle=2 --parallel
     python project.py submit -o npzzat_eq_sim #--bundle=2 --parallel
     python project.py submit -o npzzat_prod_sim #--bundle=2 --parallel
     python project.py submit -o init_inter_eq_sim #--bundle=2 --parallel
     python project.py submit -o inter_eq_sim #--bundle=2 --parallel
     python project.py submit -o inter_prod_sim #--bundle=2 --parallel
     python project.py submit -o calculate_props #--bundle=2 --parallel
   ```
     
5. Extract properties, run GP optimization and get samples for the next iteration in the Build_GPs/ directory
   ```
     qsub -N postVLE submit_job post_vle.py
   ```  

> **Note — regenerating the surface-tension GPs will not reproduce the archived DMF model.**
> Five of the six archived surface-tension GPs use a masked linear mean function, but DMF's
> uses a plain linear mean function; the current code masks all six. The archived DMF model
> is the one behind the published results. See the "Mean functions, and a known exception for
> DMF" section of [SENSITIVITY_ANALYSIS.md](SENSITIVITY_ANALYSIS.md) for details.

### Distinct Atom Type Optimization
To run generalized FF calibration, follow the following steps:
1. Initialize Signac workflow
   ```
     conda activate hfcs-fffit
     cd Opt_ES/
     python init_opt_at.py
   ```    
2. Run LD iterations
   - **Note: rm -r workspace/ signac_project_document.json signac.rc will remove everything and allow you to start fresh if you mess up**
  ```
     cd Opt_ES/opt_ff/
     python project.py submit -o OptGenFF --bundle=12 --parallel
   ```  
3. Check status a few times throughout the process
   ```  
     python project.py status
   ```   
Note: Step 2 operation OptGenFF runs multiple operations in series. Alternatively, the following can be run one at a time:
4. Alternative submission
   ```
     python project.py submit -o gen_pareto_sets --bundle=12 --parallel
     python project.py submit -o run_obj_alg --bundle=12 --parallel
   ```
     
5. Extract and analyze data:
   ```
     qsub -N postOptAnaly ../submit_job post_analysis_opt.py
   ``` 

### GP-optimized and Base FF Validation
To validate the generalized FF from above, follow the following steps:
1. Initialize Signac workflow
   ```
     conda activate hfcs-fffit
     cd Opt_ES/
   ```    
2. Run GEMC validation simulations (Pvap, rho_v, Hvap) - Simultaneously with Steps 5-7
   - **Note: rm -r workspace/ signac_project_document.json signac.rc will remove everything and allow you to start fresh if you mess up**
  ```
     cd Opt_ES/gemc_val_xxx/
     python project.py submit -o NPT_NVT --bundle=6 --parallel
     python project.py submit -o EQ --bundle=6 --parallel
     python project.py submit -o PROD --bundle=6 --parallel
     
   ``` 
Note: In step 2, "xxx" represents which FF is being validated. use "no_opt" for the Base FF and  "opt" for the GP-optimized FF.
3. Check status a few times throughout the process
   ```  
     python project.py status
   ```   
Note: Step 2 operation is multiple operations in series. Alternatively, the following can be run one at a time:
4. Alternative submission
   ```
     python project.py run -o create_forcefield
     python project.py run -o calc_boxes
     python project.py submit -o NVT_liqbox --bundle=4 --parallel
     python project.py run -o extract_final_NVT_config
     python project.py submit -o NPT_liqbox --bundle=4 --parallel
     python project.py run -o extract_final_NPT_config
     python project.py submit -o run_gemc_eq --bundle=4 --parallel
     python project.py submit -o gemc_eq_restart --bundle=4 --parallel (optional)
     python project.py submit -o run_gemc_prod --bundle=4 --parallel
     python project.py submit -o check_prod_data --bundle=4 --parallel
     python project.py submit -o calculate_props --bundle=4 --parallel
     python project.py run -o del_job (optional)
   ```
     
5. Run validation MD interface simulations (rho_l and surface tension) - Simultaneously with Steps 2-4
   - **Note: rm -r workspace/ signac_project_document.json signac.rc will remove everything and allow you to start fresh if you mess up**
  ```
     cd Opt_ES/ift_val_xxx/
     python project.py submit -o IFT #--bundle=2 --parallel
   ```  
Note: In step 5, "xxx" represents which FF is being validated. use "no_opt" for the Base FF and  "opt" for the GP-optimized FF.
6. Check status a few times throughout the process
   ```  
     python project.py status
   ```   
Note: Step 2 operation IFT runs multiple operations in series. Alternatively, the following can be run one at a time:
7. Alternative submission
   ```
     python project.py run -o create_forcefield
     python project.py submit -o create_system
     python project.py run -o fix_topology
     python project.py submit -o em_sim #--bundle=2 --parallel
     python project.py submit -o nvt_eq_sim #--bundle=2 --parallel
     python project.py submit -o npzzat_eq_sim #--bundle=2 --parallel
     python project.py submit -o npzzat_prod_sim #--bundle=2 --parallel
     python project.py submit -o init_inter_eq_sim #--bundle=2 --parallel
     python project.py submit -o inter_eq_sim #--bundle=2 --parallel
     python project.py submit -o inter_prod_sim #--bundle=2 --parallel
     python project.py submit -o calculate_props #--bundle=2 --parallel
   ```
     
8. Extract and analyze data:
   ```
     cd Opt_ES/
     qsub -N postOptAnaly ../submit_job post_analysis_val.py

     cd ../
     qsub -N makeCorrAll submit_job make_corr_figs_all_molec.py
     qsub -N makeCorrInd submit_job make_corr_figs_one_molec.py
     qsub -N makeGPSimSense submit_job make_GP_vs_sim_and_sens.py
     qsub -N makeIFTValfigs submit_job make_ift_val_figs.py
     qsub -N makeParamComp submit_job make_param_comp_figs.py
     qsub -N makeParetoComp submit_job make_pareto_comp_figs.py
     qsub -N makeSVDSens submit_job make_svd_sensitivity.py

     cd Claude_Audit/scripts
     python extract_student_ff_values.py
     python recompute_student_metrics.py
     python plot_baseline_vs_experiment.py
     python generate_provenance_table.py
     python build_master_ff_predictions.py
     python build_reference_values_all.py
     python recompute_metrics_all.py
     python generate_tables.py
   ``` 
Note: The second set of jobs create final figures for the SI and main text. Alternatively, these can be created by running the jupyter notebook `create_analysis_figs.ipynb`.

### Known Issues
The instructions outlined above seem to be system-dependent. In some cases, users have the following error:
```
ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.29' not found
```
If you observe this, please try the following in the terminal
```
export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
```
which should fix the problem. This is not an optimal solution and is something we would like to address. We found that related projects [1](https://github.com/openmm/openmm/issues/3943), [2](https://github.com/conda/conda/issues/12410) have similar issues.
If you are aware of a robust solution to this issue, please let us know by raising an issue or sending an email!

## Credits
This research is based upon work supported by the National Science Foundation under award number EEC-2330175 for the Engineering Research Center EARTH. Computing resources were provided by the Center for Research Computing (CRC) at the University of Notre Dame. MC acknowledges support from the Graduate Assistance in Areas of National Need fellowship from the Department of Education, grant number P200A210048.

## Contact
Please contact Montana Carlozo (mcarlozo@nd.edu) with any questions, suggestions, or issues.
