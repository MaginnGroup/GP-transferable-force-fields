# Usable validation replicates

The final validation design specified three independent replicates per force
field and state point. The retained usable count varies by property and
temperature. [replicate_counts.csv](replicate_counts.csv) lists all 380
property/state-point combinations for the six molecules and both force-field
variants, including zero counts. Each row links to its repository-relative
source CSV and lists the retained `restart` identifiers. These identifiers are
replicate labels in the results tables, not complete simulator job histories.

Regenerate the counts and independently calculated summaries using Python 3.10
or later, without third-party packages or simulation software:

```sh
python3 reproducibility/validation_counts.py --output-dir /tmp/es-ffo-validation
```

The command does not modify input files. It checks unique temperature/replicate
keys and compares each retained raw property with the corresponding `sim_*`
value in `ms_data.csv` using relative and absolute tolerances of 1e-12. The
committed [report.json](report.json) records the 24 input-file hashes and the
comparison outcome. `ms_data.csv` retains replicate rows; it is not a table of
replicate means. The summaries are recalculated directly from those rows.

## Means and plotted error bars

A usable value is a finite entry for the particular property. Missing one
property does not discard the other properties from that replicate. The
arithmetic mean and population standard deviation (`ddof=0`) are calculated
across available replicate values at each temperature. The property comparison
plots use 1.96 times that standard deviation, as implemented in
[Opt_ES/utilsOpt/plot.py](../../Opt_ES/utilsOpt/plot.py). These bars describe
variation across replicate means; they are not confidence intervals for the
mean. A single usable replicate has population SD zero, which does not establish
zero physical uncertainty. No usable replicates produces blank summary values.

[replicate_summaries.csv](replicate_summaries.csv) records the independently
recomputed means, population SDs and plotted halfwidths. These calculations do
not use the historical within-run `surf_tens_unc` fields. Critical-point fits and
their uncertainties are separate calculations and are not covered by this
counting script.

## Incomplete results and the methanol exclusion

- GP-Optimized methanol at 500 K retains density and surface-tension replicate
  labels 2 and 3. The manuscript reports discarding one interface run with a
  2.35 nm liquid region, compared with the stated 4.5 nm criterion. Label 1 is
  absent from these retained property values. The corresponding trajectory and
  density-profile diagnostic were not recovered here; the CSV alone cannot
  verify the physical exclusion criterion or identify a signac job directory.
- Neither methanol variant retains GEMC property values at 500 K; neither
  ethylene-glycol variant retains them at 700 K.
- Glycerol and DMF have no retained GEMC property values in these tables.
- GP-Optimized DMSO has two usable GEMC replicates at 318.15 K.
- Base DEC has two at 353.19 K. GP-Optimized DEC has three at 333.18 K and two
  at 353.19 K; both variants have no retained GEMC results at lower temperatures.

The CSV gives the property-specific counts and identifiers. Empty entries show
unavailable retained results; they do not by themselves identify why a
simulation failed. No simulations, trajectory reanalysis, or recovery of
missing results were performed by this audit.
