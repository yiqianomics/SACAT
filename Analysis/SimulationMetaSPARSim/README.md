# metaSPARSim simulation

This study uses metaSPARSim's Gamma-based abundance variation and sampling
without replacement to evaluate structural-absence (SA), relative-abundance
(RA), and joint effects under a generator external to the SACAT model.

## Design

Each replication has 152 controls, 91 cases, and 100 focal taxa plus a
background category. Target median depths are 5,855.5 and 3,678 reads.
Alternatives affect 20 focal taxa. The 38 settings comprise two global nulls
with 500 replications each and 36 alternatives with 100 each.

SA effects are 0.10, 0.15, 0.20, 0.25, 0.30, and 0.40. RA log effects are
0.20, 0.40, 0.50, 0.60, 0.80, and 1.20. Joint effects pair these values,
with ten reinforcing and ten opposing taxa. Settings include either no
covariate or an independent standard-normal covariate included in the
generator and fitted models. Exact settings and generator parameters are
in `data/design.csv` and `data/configuration.csv`.

## Figures

| PDF in `figures/` | Contents |
| --- | --- |
| [power.pdf](figures/power.pdf) | Structural-only, abundance-only, and joint-effect power |
| [fdr.pdf](figures/fdr.pdf) | FDR for the same comparisons |
| [type1_error.pdf](figures/type1_error.pdf) | Raw global-null Type I error for all 15 tests |

Power and FDR each have two pages: raw p-values followed by BH, then
null-calibrated p-values followed by BH. Columns distinguish the three
effect directions; rows distinguish settings without and with a covariate.
Effect levels are equally spaced and labeled. Joint tick labels give SA
above RA. The FDR reference line is 0.05; vertical ranges differ between
pages.

SA panels compare SACAT SA with the prevalence tests of ZINQ and MaAsLin 3.
RA panels use the corresponding abundance tests and six general-abundance
methods. Joint panels use SACAT's Bonferroni omnibus, the package-native
ZINQ and MaAsLin 3 joint tests, and the six general-abundance methods.

Type I error uses original raw p-values under the global null. Bars show
the estimate plus or minus 1.96 Monte Carlo standard errors across 500
replications. The vertical reference line is 0.05; the values on the right
are the point estimates.

## Calculations

Power is the proportion of affected taxa rejected. FDR is the mean
replication-level false discovery proportion, with zero for no discoveries.
Raw Type I error averages the proportion of null taxa with raw p < 0.05
across replications. Monte Carlo standard errors use variation across
replications.

Unavailable tests receive p = 1 and remain in the 100-taxon family. BH
rejects at adjusted p < 0.05. Calibration is separate for each method, test,
and covariate condition. The null replication sharing the evaluation
replication's index is held out. Using the remaining 49,900 null p-values,
the calibrated p-value is
`(1 + number of reference p-values <= evaluation p) / 49901`.
BH is then applied within each replication and test.

## Files and code

- `code/`: simulation, native joint-test reconstruction, input preparation, summaries, and plotting.
- `data/`: complete raw p-value families, simulated truth, designs, generator parameters, and software versions. The two `native_joint_*.rds` files contain the inputs for reconstructing the native joint tests.
- `tables/`: power, FDR, raw Type I error, availability, and complete setting- and replication-level summaries.
- `figures/`: the three PDFs above.

`tables/main_plot_data.csv` supplies the metaSPARSim panel of the
[main simulation figure](../MainFigures/README.md). The study PDFs also
include general-abundance methods in the abundance-only panels.

With `data.table`, `R.utils`, `ggplot2`, and `systemfonts` installed, regenerate
tables and figures from this study directory:

```sh
Rscript code/summarize_results.R
Rscript code/plot_results.R
```

Scripts locate study files from their own paths. For new simulations,
install the comparison packages in `data/software_versions.csv`.
`code/install_metasparsim.R` installs metaSPARSim 1.1.2 with R-seeded
sampling into the active R library. Submit from this study directory:

```sh
Rscript code/run_simulation.R preflight
mkdir -p logs
sbatch code/submit.slurm
```

New runs write to `output/`, or `METASPAR_ROOT` when specified. After all
tasks finish, run `code/run_simulation.R summarize`, then
`code/prepare_inputs.R`, `code/reconstruct_joint_tests.R`,
`code/summarize_results.R`, and `code/plot_results.R`.
Native joint reconstruction uses ZINQ 2.0 and maaslin3 1.5.3.
