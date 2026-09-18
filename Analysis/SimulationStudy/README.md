# SACAT Simulation Study

This directory contains the simulation program, cluster submission template, figures, numerical summaries, and tabular results for the SACAT simulation study.

## Study design

Each simulation contains 50 focal taxa in two equally sized groups, with 10 taxa designated for perturbation. The design crosses three sample sizes (60, 80, or 120 observations per group) and two covariate settings. Without a covariate, the data-generating model and fitted models contain only the group effect. With a covariate, a standard-normal variable is generated independently of group, standardized across samples, and included in both data generation and analysis. Library sizes have a target median of 8,000 reads in both groups.

The experiments cover:

- observed-prevalence perturbations with fixed structural-absence probabilities;
- structural-absence perturbations calibrated to preserve expected observed prevalence;
- present-conditional abundance perturbations;
- structural-only perturbations for component specificity.

| Experiment | Effect levels |
|---|---|
| Observed-prevalence difference with fixed structural absence | 0, 0.08, 0.16, 0.24, 0.32 |
| Structural-absence difference with matched observed prevalence | 0, 0.10, 0.20, 0.35, 0.50 |
| Absolute present-conditional mean log-abundance difference | 0, 0.20, 0.40, 0.75, 1.40 |
| Structural-only difference | 0, 0.10, 0.20, 0.35, 0.50 |

The 20 experiment-effect combinations and six design strata give 120 settings, each evaluated in 100 independent replications. Abundance effects are positive for five target taxa and negative for five. Each experiment includes its zero-effect endpoint. Baseline profiles retain a fixed 20-taxon template pool, of which 10 taxa are perturbed.

## Programs

- `run_simulation.R` defines the data-generating mechanisms, runs one replication, applies the analysis methods, and aggregates the simulation results.
- `submit.slurm` is the SLURM array template used for the 100 replications.
- `summarize_results.R` validates the completed grid and creates the figures and CSV tables.

The simulation program uses SACAT 0.7.1 and the packages listed in its `replicate_packages` and `summary_packages` objects. The supplied SLURM file provides a portable 100-task array template. Load the required R environment before submission and add any cluster-specific partition or resource options locally. All methods include the covariate only when it is present; MaAsLin 3 and ZINQ also adjust for standardized log sequencing depth.

ZINQ 2.0 is available from the [ZINQ-v2 repository](https://github.com/wdl2459/ZINQ-v2), and MaAsLin3 1.5.3 is available from the [MaAsLin3 repository](https://github.com/biobakery/maaslin3).

The figure and table script requires `data.table`, `ggplot2`, `patchwork`, `scales`, and an R build with Cairo graphics.

From `Analysis/SimulationStudy`, set the repository and output directories and run:

```bash
export SACAT_PROJECT_ROOT=/path/to/SACAT
export SACAT_SIMULATION_ROOT=/path/to/sacat_internal_validation_output
Rscript run_simulation.R preflight
Rscript run_simulation.R design
sbatch submit.slurm
```

After the array has completed, aggregate its outputs with:

```bash
Rscript run_simulation.R summarize
```

## Simulation results

The full simulation result archive contains taxon-level results, simulation truth, replication-level metrics, setting summaries, completion records, method-status summaries, and software information. Counts and latent states are retained in the dataset files.

With `SACAT_SIMULATION_ROOT` set to the run directory, create figures and tables from its `summary/` directory with:

```bash
Rscript summarize_results.R
```

The script checks the setting grid against the saved design and verifies all 100 replications. New figures and tables are written to `report/` within the run directory.

The figures and tables already in this directory accompany the SACAT 0.6.0 result archive. Their 540-setting design, covariate structures, and software versions are recorded in `tables/`. To reproduce those outputs, extract the archive into `results_data/source/` and run the figure script without `SACAT_SIMULATION_ROOT` set; the reproduced files are written to `simulation_report/`. The figure script reads the design and software information supplied with each set of results. `SACAT_SIMULATION_SOURCE` and `SACAT_SIMULATION_REPORT_ROOT` can specify other input and output directories.

## Figures

Core figures:

| File | Content |
|---|---|
| [`main_mechanism_separation.pdf`](figures/main_mechanism_separation.pdf) | Structural-absence and observed-prevalence rejection curves in the two estimand-separation experiments |
| [`main_abundance_power.pdf`](figures/main_abundance_power.pdf) | Present-conditional abundance power across all sample sizes, signal fractions, covariate structures, and nonzero effect sizes |

Additional figures:

| File | Content |
|---|---|
| [`estimand_calibration.pdf`](figures/estimand_calibration.pdf) | Target and achieved probability contrasts in the estimand-separation experiments |
| [`component_specificity.pdf`](figures/component_specificity.pdf) | Marginal Type I error for the non-target SACAT component |
| [`abundance_fdr_signal20.pdf`](figures/abundance_fdr_signal20.pdf) | Present-conditional abundance false discovery rates with 10 perturbed taxa |
| [`abundance_fdr_signal40.pdf`](figures/abundance_fdr_signal40.pdf) | Present-conditional abundance false discovery rates with 20 perturbed taxa |
| [`global_null_family_rejection_structural.pdf`](figures/global_null_family_rejection_structural.pdf) | Structural and observed-prevalence familywise Type I error under three global-null generators |
| [`global_null_family_rejection_abundance.pdf`](figures/global_null_family_rejection_abundance.pdf) | Abundance familywise Type I error under three global-null generators |
| [`correlated_community_abundance_power.pdf`](figures/correlated_community_abundance_power.pdf) | Abundance power under correlated-community perturbations |
| [`correlated_community_abundance_fdr.pdf`](figures/correlated_community_abundance_fdr.pdf) | Abundance false discovery rates under correlated-community perturbations |
| [`correlated_community_structural_performance.pdf`](figures/correlated_community_structural_performance.pdf) | Structural and observed-prevalence power and false discovery rates under correlated-community structural-absence perturbations |

Each curve is evaluated within its displayed sample-size, signal-fraction, covariate, and effect stratum. Points and intervals summarize 100 independent replications. The vector PDFs use a consistent publication-scale layout.

## Tables

Design and implementation records:

- `simulation_constants.csv` and `simulation_design.csv`;
- `method_targets.csv` and `software_versions.csv`;
- `confounding_design.csv` (covariate-generation summaries);
- `test_availability.csv` and `test_unavailability_reasons.csv`.

Numerical values underlying the figures:

- `estimand_calibration_summary.csv`;
- `structural_estimand_summary.csv`;
- `component_specificity_summary.csv`;
- `abundance_performance_summary.csv`;
- `global_null_family_rejection_summary.csv`;
- `correlated_community_summary.csv`.

## Performance measures

Within each replication, setting, method, and component, Benjamini-Hochberg adjustment is applied to the fixed family of 50 focal taxa. Each unavailable taxon test enters the multiplicity family with a p-value of one and contributes a non-rejection.

- Power is the proportion of directly perturbed taxa rejected after adjustment.
- The false discovery proportion is the number of false discoveries divided by the total number of discoveries, with zero assigned when no discovery occurs. The empirical false discovery rate is its mean over 100 replications.
- Under a global null, familywise Type I error is the probability of at least one rejection. It equals the empirical false discovery rate because all 50 focal taxa are null.
- `type1_error` and `type1_error_bh` report the mean rejection proportion among null taxa before and after BH adjustment. `familywise_type1_error_raw` and `familywise_type1_error_bh` record whether any null taxon is rejected. Each is computed within replication before averaging.
- Marginal Type I error in the component-specificity experiment is the mean rejection probability for the non-target SACAT component among taxa designated for the other component.
- Monte Carlo intervals equal the replication-level mean plus or minus 1.96 Monte Carlo standard errors, truncated to the parameter range.

The global-null figures use the zero-effect abundance experiment. For result archives containing multiple signal-template strata, those strata are averaged within replication before the Monte Carlo mean and interval are computed.
