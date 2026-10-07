# Internal simulations

This directory contains the model-based simulations used in the manuscript,
including complete saved p-value families, design and truth records, simulation
programs, and reproducible numerical summaries. The manuscript figure is drawn
by [MainFigures](../MainFigures/README.md) from `tables/main_plot_data.csv`.

## Design

Every dataset has 50 focal taxa. Each setting has 100 independent replications.
The main experiments use 60, 80, or 120 observations per group and either no
covariate or an independently generated standard-normal covariate included in
the data-generating and fitted models. The target median library size is 8,000
reads in both groups. Nonzero settings designate 10 taxa for perturbation.

| Experiment | Settings | Effect levels |
|---|---:|---|
| Observed prevalence, fixed structural absence | 30 | 0, 0.08, 0.16, 0.24, 0.32 |
| Structural absence, matched observed prevalence | 30 | 0, 0.10, 0.20, 0.35, 0.50 |
| Present-conditional abundance | 30 | 0, 0.20, 0.40, 0.75, 1.40 |
| Structural-only component specificity | 30 | 0, 0.10, 0.20, 0.35, 0.50 |
| Joint structural and abundance effects | 42 | SA/RA: 0/0, .10/.20, .15/.40, .20/.50, .25/.60, .35/.80, .50/1.40 |
| Additional abundance effects at n = 80 | 6 | 0.30, 0.50, 0.60, with and without a covariate |

The joint experiment includes five reinforcing and five opposing effect taxa.
SACAT uses its Bonferroni omnibus test; ZINQ and MaAsLin 3 use their native joint
tests. The other six methods use their abundance tests. The abundance supplement
also repeats ZINQ on the ten original n = 80 abundance settings, including the
two null settings, so its original and added effect levels use the same analysis
environment and null references.

## Saved inputs

`data/` contains five compressed taxon-result files. Each contains all 50 taxa
within every included method-result-replication family, including unavailable
tests. The fields retain raw p values, availability, component-specific truth,
designated-taxon labels, and status. Joint results also identify the primary
comparisons. P-value precision is preserved; no model fitting is required to
recompute the tables.

| File | Taxon results |
|---|---:|
| `base_n060.csv.gz`, `base_n080.csv.gz`, `base_n120.csv.gz` | 1,700,000 each |
| `joint.csv.gz` | 3,150,000 |
| `abundance.csv.gz` | 160,000 |

The supplement contains 90,000 tests at the added effects, 50,000 locally rerun
ZINQ tests from the original grid, and 20,000 saved SACAT/MaAsLin 3 null tests.
The base files retain the complete original grid; the manuscript table uses
the supplement's ZINQ abundance results at n = 80. Its displayed abundance
effects are 0.20, 0.30, 0.40, 0.50, 0.60, and 0.75.

The three `*_design.csv` files record the exact settings. `base_truth.csv.gz`
and `joint_truth.csv.gz` retain the full simulation truth records.

Recorded R and package versions are in `data/software_versions.csv`.
The supplied base and joint runners require SACAT 0.7.1 or later.

## Recompute the tables

With R, `data.table`, and `R.utils` installed, run from any working directory:

```sh
Rscript Analysis/SimulationStudy/summarize_results.R
```

Set `SACAT_SIMULATION_DATA` or `SACAT_SIMULATION_TABLES` to use other input or
output directories. All reported rates are proportions between zero and one.

- `setting_summary.csv` contains raw-p + BH and calibrated-p + BH power, FDR,
  designated-taxon rejection, availability, and Monte Carlo standard errors.
- `replication_metrics.csv.gz` retains all 168,200 result families used in
  those summaries.
- `type1_error.csv` reports global-null Type I error from raw p < 0.05 before
  multiplicity adjustment, averaged over taxa and replications.
- `test_availability.csv` gives complete test counts by result family and status.
- `main_plot_data.csv` contains the 204 manuscript Panel A points.
- `manuscript_summary.csv` gives the equally weighted setting means for SACAT
  structural FDR (48 nonzero structural-only/matched-prevalence settings) and
  abundance FDR (24 nonzero abundance settings).

For calibration, null p values are matched by scenario, sample size, covariate,
method, and component. The entire null replication sharing the evaluation
replication index is excluded. The calibrated p value is
`(1 + number of reference p values <= p) / 4951`, using 4,950 reference tests.
Unavailable tests enter with p = 1. BH adjustment is then applied within each
fixed family of 50 taxa at 0.05. FDR is the mean replication-level false discovery
proportion, with zero when there are no discoveries. The first two manuscript
scenarios report designated-taxon rejection; abundance and joint scenarios
report power. Raw and calibrated results have separate columns throughout.

## Run the simulations

The runners list their required R packages. ZINQ 2.0 is available from
[ZINQ-v2](https://github.com/wdl2459/ZINQ-v2) and MaAsLin 3 1.5.3 from
[MaAsLin3](https://github.com/biobakery/maaslin3). Use the recorded package
versions when comparing numerical fits across environments.

From this directory, set `SACAT_PROJECT_ROOT` to the repository and
`SACAT_SIMULATION_ROOT` and `SACAT_JOINT_ROOT` to the desired run directories.
Load the R environment and add cluster-specific options to the SLURM templates.

```sh
mkdir -p logs
Rscript run_simulation.R preflight
Rscript run_simulation.R design
sbatch submit.slurm
# After all 100 array tasks finish:
Rscript run_simulation.R summarize

Rscript run_joint.R check
Rscript run_joint.R design
sbatch submit_joint.slurm
# After all 100 array tasks finish:
Rscript run_joint.R summarize
```

The abundance supplement uses SACAT 0.7.1 and the saved base design and null
references. Set `SACAT_ABUNDANCE_ROOT` to its output directory, then run:

```sh
Rscript run_abundance_supplement.R prepare
Rscript run_abundance_supplement.R preflight
Rscript run_abundance_supplement.R worker 1 1
Rscript run_abundance_supplement.R zinq 1 1
```

The two final arguments are worker index and total workers. Different worker
indices can be run concurrently. To prepare another completed run for the same
table workflow, collect its inputs into a separate directory:

```sh
Rscript prepare_results.R base "$SACAT_SIMULATION_ROOT/summary" /path/to/new_inputs
Rscript prepare_results.R joint "$SACAT_JOINT_ROOT/summary" /path/to/new_inputs
Rscript prepare_results.R abundance "$SACAT_ABUNDANCE_ROOT" /path/to/new_inputs
SACAT_SIMULATION_DATA=/path/to/new_inputs Rscript summarize_results.R
```

Preparation and summarization only read saved results. Counts and fitted
objects remain in the simulation run directories.
