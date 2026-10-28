# Corrected five-candidate DMF figure

This package reproduces the approved manuscript figure at full text width,
with larger lettering and corrected within-run surface-tension error bars.
The molecule/parameter identities, simulation means, temperatures, colors and
RD2 reference values are preserved. Density bars use 1.96 times the stored
within-run block-average estimates. Surface-tension bars use 1.96 times the
recovered maximum GROMACS block error estimate over three through ten blocks.

From the repository root, with NumPy, pandas and Matplotlib installed:

```sh
python scripts/plot_dmf_corrected.py --output-dir /tmp/es-ffo-dmf
```

The script locates this public repository relative to its own path. It also
accepts `--repo-root` explicitly. It reads only committed public files and
writes PNG, SVG, numerical plotted data, and source hashes to the output
directory. No molecular simulation software, GP model loading, REFPROP
installation, sibling private repository or manuscript checkout is required.

The generator checks candidate IDs 11, 13, 15, 17 and 22, their five
temperatures, and all 25 numerical rows against
[reference_figure_data.csv](reference_figure_data.csv) at relative and absolute
tolerances of 1e-12. That table and the lettering/layout originate from audit
commit `0e2766c` in `dowlinglab/ESS-FO-Audit`. The figure was approved for the
manuscript in commit `a78d6b5` of its separate Overleaf repository.

[provenance.json](provenance.json) records the public input hashes and tested
plotting package versions. The regenerated PNG was also compared pixel by
pixel with the approved manuscript image in that environment and matched.
Other Matplotlib versions can render differently while preserving numerical
content. Compare data and inspect the plot rather than relying solely on
image-byte equality.

See [historical_uncertainty.md](../historical_uncertainty.md) for the parser
history and inverse conversion. Original CSV uncertainty fields are unchanged.
The original DMF block outputs were unavailable: only their maximum error
estimate is recoverable. These bars are not asserted to be calibrated 95%
confidence intervals. Use this command for the current manuscript figure;
historical plotting scripts can reproduce the old stored-field convention.
