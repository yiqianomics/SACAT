# DASRA Simulation Study

This directory contains the simulation program, cluster submission template, manuscript figures, numerical summaries, and tabular results for the DASRA simulation study.

## Study design

Each simulation contains 50 focal taxa in two equally sized groups. The design crosses three sample sizes (60, 80, or 120 observations per group), two signal fractions (20% or 40%, corresponding to 10 or 20 perturbed taxa), and unconfounded or confounded covariate structures. The 45 base settings and 12 design strata define 540 settings, each evaluated in 100 independent Monte Carlo replications.

The experiments cover:

- observed-prevalence perturbations with fixed structural-absence probabilities;
- structural-absence perturbations calibrated to preserve expected observed prevalence;
- present-conditional abundance perturbations;
- structural-only perturbations for component specificity;
- a global null with unequal library-size distributions;
- correlated-community global-null, abundance, and structural experiments.

The effect grids and common design constants are recorded in [`tables/simulation_design.csv`](tables/simulation_design.csv) and [`tables/simulation_constants.csv`](tables/simulation_constants.csv).

## Programs

- `dasra_formal_hpc_simulation.R` defines the data-generating mechanisms, runs one replication, applies the analysis methods, and aggregates the simulation results.
- `submit.slurm` is the SLURM array template used for the 100 replications.
- `summarize_results.R` validates the completed grid and creates the manuscript figures and CSV tables.

The simulation program uses DASRA 0.6.0 and the packages listed in its `replicate_packages` and `summary_packages` objects. The supplied SLURM file records the cluster configuration used for the reported study. Adapt its partition, resources, notification settings, environment location, working directory, log directory, and output root to the target cluster.

The figure and table script requires `data.table`, `ggplot2`, `patchwork`, `scales`, and an R build with Cairo graphics.

The simulation program reads its output root from `DASRA_FORMAL_ROOT`:

```bash
export DASRA_FORMAL_ROOT=/path/to/dasra_formal_simulation
mkdir -p "$DASRA_FORMAL_ROOT/logs"
Rscript dasra_formal_hpc_simulation.R preflight
Rscript dasra_formal_hpc_simulation.R design
sbatch submit.slurm
```

After the array has completed, aggregate its outputs with:

```bash
Rscript dasra_formal_hpc_simulation.R summarize
```

## Simulation results

The full simulation result archive is stored outside version control because of its size and is available from the authors upon request. The archive contains the complete taxon-level results, simulation truth, replication-level metrics, setting summaries, completion records, method-status summaries, and software information.

Extract the supplied archive into `results_data/`, or place the output of the `summarize` step under `results_data/source/`, and run:

```bash
Rscript summarize_results.R
```

The script verifies all 54,000 setting-replication combinations before creating the figures and tables below.

## Figures

The two main-text figures are:

| File | Content |
|---|---|
| [`main_mechanism_separation.pdf`](figures/main_mechanism_separation.pdf) | Structural-absence and observed-prevalence rejection curves in the two estimand-separation experiments |
| [`main_abundance_power.pdf`](figures/main_abundance_power.pdf) | Present-conditional abundance power across all sample sizes, signal fractions, covariate structures, and nonzero effect sizes |

The supplementary figures are:

| File | Content |
|---|---|
| [`estimand_calibration.pdf`](figures/estimand_calibration.pdf) | Target and achieved probability contrasts in the estimand-separation experiments |
| [`component_specificity.pdf`](figures/component_specificity.pdf) | Rejection proportions for the non-target DASRA component |
| [`abundance_fdr_signal20.pdf`](figures/abundance_fdr_signal20.pdf) | Present-conditional abundance false discovery rates with 10 perturbed taxa |
| [`abundance_fdr_signal40.pdf`](figures/abundance_fdr_signal40.pdf) | Present-conditional abundance false discovery rates with 20 perturbed taxa |
| [`global_null_family_rejection_structural.pdf`](figures/global_null_family_rejection_structural.pdf) | Structural and observed-prevalence family-rejection probabilities under three global-null generators |
| [`global_null_family_rejection_abundance.pdf`](figures/global_null_family_rejection_abundance.pdf) | Abundance family-rejection probabilities under three global-null generators |
| [`correlated_community_abundance_power.pdf`](figures/correlated_community_abundance_power.pdf) | Abundance power under correlated-community perturbations |
| [`correlated_community_abundance_fdr.pdf`](figures/correlated_community_abundance_fdr.pdf) | Abundance false discovery rates under correlated-community perturbations |

Each curve is evaluated within its displayed sample-size, signal-fraction, covariate, and effect stratum. Points and intervals summarize 100 independent replications. The vector PDFs are formatted for placement with manuscript captions.

## Tables

Design and implementation records:

- `simulation_constants.csv` and `simulation_design.csv`;
- `method_targets.csv` and `software_versions.csv`;
- `confounding_design.csv`;
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
- Under a global null, the family-rejection probability is the probability of at least one rejection. It equals the empirical false discovery rate because all 50 focal taxa are null.
- Structural-estimand and component-specificity curves report the proportion of designated taxa rejected after adjustment.
- Monte Carlo intervals equal the replication-level mean plus or minus 1.96 Monte Carlo standard errors, truncated to the parameter range.

For the global-null summaries, the two signal-template strata are averaged within replication before the Monte Carlo mean and interval are computed.
