# Internal model validation

This study examines how SACAT separates structural absence (SA) from relative
abundance (RA) conditional on presence. Counts are generated under the
latent-state model, allowing each source of variation to be controlled.

## Design

Each setting has 100 replications, 50 focal taxa, and ten affected taxa under
the alternatives. Sample sizes are 60, 80, and 120 per group, with target
median sequencing depth of 8,000 reads in both groups. Settings either omit
the covariate or include a standard-normal covariate generated independently
of group and standardized within each replication.

The experiments examine:

- abundance-driven detection changes while SA remains fixed;
- SA changes while expected detection rates remain matched between groups;
- SA-only and RA-only effects; and
- joint effects, with five reinforcing and five opposing taxa.

Additional RA effects of 0.30, 0.50, and 0.60 are evaluated at 80 samples
per group. Exact grids are in `data/*_design.csv`.

## Figures and results

| PDF in `figures/` | Contents |
| --- | --- |
| [mechanism_separation.pdf](figures/mechanism_separation.pdf) | Distinguishing abundance-driven detection changes from SA changes |
| [power.pdf](figures/power.pdf) | SA-only, RA-only, and joint-effect power |
| [fdr.pdf](figures/fdr.pdf) | False discovery rate (FDR) for the same comparisons |
| [type1_error.pdf](figures/type1_error.pdf) | Raw global-null Type I error and Monte Carlo intervals |

The mechanism, power, and FDR figures use 80 samples per group. Their first
page uses original p-values followed by Benjamini–Hochberg (BH) adjustment at
0.05; the second uses matched-null p-value calibration followed by BH.
Rows distinguish settings without and with a covariate.

SA panels compare SACAT SA with ZINQ and MaAsLin 3 prevalence tests. RA panels
use abundance tests. Joint panels use the package-native omnibus or joint
tests of SACAT, ZINQ, and MaAsLin 3, together with six abundance methods.
The six abundance methods use the original RA effect grid; SACAT, ZINQ, and
MaAsLin 3 also use the supplementary effects.

Power is the proportion of affected taxa detected; FDR is the mean false
discovery proportion across replications. Type I error is the mean proportion
of null taxa with raw p < 0.05, shown for all three sample sizes with intervals
of ±1.96 Monte Carlo standard errors.

## Files and code

- `code/`: simulation, summary, and plotting scripts.
- `data/`: raw p-values, component-specific truth, designs, and software versions.
- `tables/`: setting-level and replication-level results.
- `figures/`: the four PDFs above.

`tables/main_plot_data.csv` supplies the internal-validation panel of the
[main simulation figure](../MainFigures/README.md).

With `data.table`, `R.utils`, and `ggplot2` installed, regenerate the tables
and figures from this directory:

```sh
Rscript code/summarize_results.R
Rscript code/plot_results.R
```
