# Main manuscript figures

Each R script generates one PDF from saved numerical inputs.
The scripts resolve paths relative to their own location and can be run from
any working directory. Model fitting is not repeated.

## Generate the figures

From the repository root:

```sh
Rscript Analysis/MainFigures/make_simulation_figure.R
Rscript Analysis/MainFigures/make_real_data_figure.R
```

Use R 4.1 or later, `ggplot2` 3.5 or later, `patchwork` 1.3 or later, and
`systemfonts`. The scripts use Arial when available, then Liberation Sans or
the system sans-serif font. PDFs are written with Quartz on macOS and Cairo
on other platforms.

| Script | Output in `figures/` |
| --- | --- |
| `make_simulation_figure.R` | `main_simulation_calibrated_bh_power_combined.pdf` |
| `make_real_data_figure.R` | `main_real_data_combined.pdf` |

## Simulation figure

Panel A shows model-based simulations with 80 samples per group. Panel B
shows metaSPARSim simulations. Rows distinguish settings without and with a
covariate. Component comparisons use SACAT, ZINQ, and MaAsLin 3; joint-effect
comparisons also include the six abundance methods.

The figure reads `tables/main_plot_data.csv` from
[`SimulationStudy`](../SimulationStudy/README.md) and
[`MetaSPARSimSimulation`](../MetaSPARSimSimulation/README.md). Together, these
tables contain the 384 plotted values and their design labels.
`plotted_value` uses matched-null p-value calibration followed by BH
adjustment at 0.05 within each replication and method-result family.
`original_bh_value` records the corresponding result without null calibration
and is not plotted. Panel A's first two scenarios report rejection among the
designated taxa; its abundance and joint scenarios, and all of Panel B,
report power.

Both simulation directories include the complete raw p-value inputs and scripts
for recomputing their summaries.

## Real-data figure

Panels A and B show mock-community signal recovery. Panels C–F show method
comparisons and component profiles for *Coprococcus*, *Enterococcus*, and
*Faecalibacterium* in the Schubert CDI cohort. Inputs are in `data/real_data/`:

- `rauer_mock_recovery.csv` contains discovery totals over 52 allocations for
  eight structural targets or three abundance targets at each sequencing depth.
- `schubert_cdi_figure_method_results.csv` contains all 88 taxa in each
  method-result family, including p-values, BH-adjusted p-values, and signed
  statistics.
- `schubert_cdi_figure_sample_values.csv` contains counts, library sizes, and
  posterior SA probabilities (`gamma`) for the three displayed taxa. These
  probabilities come from the SA score-test fit with its SA group coefficient
  fixed at zero. Panel E averages `gamma` and `I(count = 0) - gamma` over all
  samples in each group to obtain the estimated SA and ND shares.
- `schubert_cdi_figure_fitted_abundance.csv` contains the fitted abundance
  markers in Panel F. Each is 100 times the exponential of the fitted mean log
  relative abundance conditional on presence, standardized over all 243
  covariate profiles with group set to Healthy or CDI. The boxplots use positive
  observed relative abundances; the axis is logarithmic.

Data preparation and analysis are documented under
[`cdi_schubert`](../RealDataAnalysis/cdi_schubert/README.md) and
[`rauer_mock`](../RealDataAnalysis/rauer_mock/README.md).
