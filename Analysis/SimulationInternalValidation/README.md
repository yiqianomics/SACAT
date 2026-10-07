# Internal model validation

Model-based simulations evaluate structural absence (SA), observed prevalence,
present-conditional relative abundance (RA), and joint effects.

## Design

Each setting has 100 replications and 50 focal taxa, with ten designated
taxa under nonzero effects. Sample sizes are 60, 80, and 120 per group.
Settings include either no covariate or an independent standard-normal
covariate included in the generator and fitted models. Target median
sequencing depth is 8,000 reads in both groups.

The base experiments vary prevalence with SA fixed, SA with matched expected
prevalence, SA alone, or RA alone. Joint effects include five reinforcing
and five opposing taxa. The abundance supplement adds effects of 0.30,
0.50, and 0.60 at n = 80. Exact settings are in `data/*_design.csv`.

## Figures

| PDF in `figures/` | Contents |
| --- | --- |
| [mechanism_separation.pdf](figures/mechanism_separation.pdf) | Designated-taxon rejection with SA fixed or expected prevalence matched |
| [power.pdf](figures/power.pdf) | Structural-only, abundance-only, and joint-effect power |
| [fdr.pdf](figures/fdr.pdf) | FDR for the same comparisons |
| [type1_error.pdf](figures/type1_error.pdf) | Raw global-null Type I error and Monte Carlo intervals |

Mechanism, power, and FDR figures show n = 80. Their first page uses raw
p-values followed by BH; the second uses null-calibrated p-values followed
by BH. Rows distinguish settings without and with a covariate. Power and
FDR effect levels are equally spaced and labeled. Joint labels give SA/RA
effect pairs. The FDR reference line is 0.05; vertical ranges differ between
pages.

SA panels compare SACAT SA with ZINQ and MaAsLin 3 prevalence tests.
RA panels use the corresponding abundance tests. Joint panels use SACAT's
Bonferroni omnibus and the package-native ZINQ and MaAsLin 3 joint tests.
The six general-abundance methods appear in RA and joint power/FDR panels.
Their RA curves use the original effect grid; the three component methods
also use the supplementary effects.

Type I error includes all three sample sizes. Page 1 shows structural and
abundance tests; page 2 shows joint-effect comparisons. Structural tests
use the structural-only null, RA tests use the abundance-only null, and
joint comparisons use the joint null. ZINQ's n = 80 RA results and null
reference come from the supplementary run. Bars are the estimate plus or
minus 1.96 Monte Carlo standard errors; the reference line is 0.05.

## Calculations

Power is the proportion of affected taxa rejected, using component-specific
truth. FDR is the mean replication-level false discovery proportion, with
zero for no discoveries. Raw Type I error averages the proportion of null
taxa with raw p < 0.05 across replications. Monte Carlo standard errors use
variation across replications.

Unavailable tests receive p = 1 and remain in the 50-taxon family. BH rejects
at adjusted p < 0.05. Calibration uses null p-values matched by scenario,
sample size, covariate, method, and component, excluding the null replication
with the same index as the evaluation replication. The calibrated p-value is
`(1 + number of reference p-values <= evaluation p) / 4951`.

## Files and code

- `code/`: simulation runners, input preparation, summaries, and plotting.
- `data/`: complete raw p-value families, simulated truth, designs, and software versions. `abundance_reference_truth.csv` supplies the supplementary runner's preflight reference.
- `tables/`: all sample sizes, setting and replication summaries, raw Type I error, availability, and manuscript summaries.
- `figures/`: the four PDFs above.

`tables/main_plot_data.csv` supplies the internal-validation panel of the
[main simulation figure](../MainFigures/README.md).

With `data.table`, `R.utils`, and `ggplot2` installed, regenerate tables
and figures from this study directory:

```sh
Rscript code/summarize_results.R
Rscript code/plot_results.R
```

Scripts locate study files from their own paths. For new simulations, install
the packages listed by the runners. Set `SACAT_PROJECT_ROOT` to the repository,
and `SACAT_SIMULATION_ROOT` and `SACAT_JOINT_ROOT` to separate output
directories. From this study directory:

```sh
Rscript code/run_simulation.R preflight
Rscript code/run_joint.R check
mkdir -p logs
sbatch code/submit.slurm
sbatch code/submit_joint.slurm
```

Each runner's `summarize` mode collects completed replications.
`code/prepare_results.R base|joint|abundance INPUT_DIRECTORY [DATA_DIRECTORY]`
prepares inputs for the summary script. The supplementary runner supports
`prepare`, `preflight`, `worker <index> <workers>`, and
`zinq <index> <workers>` modes.
