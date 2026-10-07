# metaSPARSim simulations

This directory contains the 38-setting simulation and complete taxon-level
inputs for the external simulation results in the manuscript. Each replication
has 152 controls, 91 cases, and a fixed testing family of 100 taxa. There are
two global-null settings with 500 replications each and 36 alternative settings
with 100 replications each. Alternatives affect 20 taxa.

## Recompute the reported results

From the repository root, with R and `data.table` installed:

```sh
Rscript Analysis/MetaSPARSimSimulation/summarize_results.R
```

The script reads every saved p-value family, recalculates raw-p BH and
null-calibrated BH, and writes `tables/`. It does not fit models or generate
simulated data. `tables/main_plot_data.csv` contains the 180 points for panel B
of the main simulation figure.

| Table | Contents |
| --- | --- |
| `type1_error.csv` | Global-null rejection at raw p < 0.05 and Monte Carlo standard errors |
| `power.csv` | Power after raw-p BH and calibrated-p BH |
| `fdr.csv` | Mean false discovery proportions after each BH procedure |
| `availability.csv` | Overall, signal-taxon, and null-taxon test availability |
| `setting_summary.csv` | Complete setting-level metrics |
| `replication_metrics.csv.gz` | Metrics for each replication and test |
| `null_calibration.csv` | Empirical-null calibration checks |
| `design.csv`, `software_versions.csv` | Design and software used for the saved results |

All rates are proportions. Power and FDR are averaged over replications;
comparisons across settings weight settings equally. Unavailable tests receive
decision p = 1 and remain in the 100-taxon family. BH rejections use q < 0.05.

Calibration is separate for each method, test, and covariate condition. For an
evaluation replication, all 100 null tests with the same replication index are
held out. With the remaining 49,900 null p-values, the calibrated p-value is
`(1 + number of reference p-values <= evaluation p) / 49901`. BH is then applied
within each replication, setting, and test.

## Saved inputs

The four `data/raw_*.csv.gz` files contain all 6,420,000 p-value records:
1,380,000 SACAT results, 1,840,000 comparison-component results, 2,760,000
general-abundance results, and 440,000 native joint results. No unavailable
taxa or null families are removed. Each row records the setting, replication,
taxon, test, p-value, availability, component truth, and any-effect truth.
`tests.csv` maps test identifiers to methods and components; `design.csv`
maps settings to effects and covariates. `taxon_template.csv`,
`configuration.csv`, and `replication_index.csv` record the generator parameters
and seeds.

The two `native_joint_*.rds` files preserve the original component p-values and
taxon zero rates without loss of numeric precision and are read with `readRDS`.
To reconstruct the native joint p-values, install ZINQ 2.0 and maaslin3 1.5.3,
then run:

```sh
Rscript Analysis/MetaSPARSimSimulation/reconstruct_joint_tests.R
Rscript Analysis/MetaSPARSimSimulation/summarize_results.R
```

ZINQ uses its native Cauchy combination of the original prevalence and three
quantile p-values, including the taxon zero rate. MaAsLin 3 uses its native
joint calculation and component-error handling. Their native adjusted p-values
are not used: the common calibration and BH procedures follow the joint raw
p-values. `tables/native_joint_checks.csv` records availability and fallback
counts.

## Generate a new simulation run

The saved results use the versions listed in `data/software_versions.csv`,
including SACAT 0.7.0 and metaSPARSim 1.1.2. Full simulation also requires the
comparison packages listed there. The supplied metaSPARSim source package uses
R-seeded sampling and can be installed into the active R library with:

```r
install.packages("Analysis/MetaSPARSimSimulation/dependencies/metaSPARSim_1.1.2.tar.gz",
                 repos = NULL, type = "source")
```

The runner retains the six SA effects (.10, .15, .20, .25, .30, .40), six RA log
effects (.20, .40, .50, .60, .80, 1.20), and their paired joint effects. It writes
to `output/` by default; `METASPAR_ROOT` can select another output directory.

```sh
cd Analysis/MetaSPARSimSimulation
Rscript run_simulation.R preflight
mkdir -p logs
sbatch submit.slurm
```

The Slurm array has 200 tasks and uses the current R environment. Adjust Slurm
resource or partition options for the cluster. After all tasks finish, replace
the analysis inputs and regenerate tables with:

```sh
Rscript run_simulation.R summarize
Rscript prepare_inputs.R
Rscript reconstruct_joint_tests.R
Rscript summarize_results.R
```

`prepare_inputs.R` accepts an optional completed simulation output directory.
It reads `summary/` and the 2,200 original null/joint replication objects in
`results/`; the complete simulation archive is not required to reanalyse the
bundled inputs.
