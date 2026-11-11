# Main manuscript figures

These figures bring together the simulation studies and the real-data
applications. The scripts create PDFs from the saved numerical results.

## Simulation figures

The figures show internal model validation at 80 samples per group and
metaSPARSim simulations, available as combined or separate PDFs.

The comparisons cover structural-absence (SA), relative-abundance (RA), and
joint results. Component panels use the corresponding SACAT, ZINQ, and
MaAsLin 3 tests; joint panels also include six abundance methods. Rows
distinguish settings without and with a covariate. Plotted values use
matched-null p-value calibration followed by Benjamini–Hochberg (BH)
adjustment at 0.05.

Inputs are the `tables/main_plot_data.csv` files in
[Internal model validation](../SimulationInternalValidation/README.md) and
[metaSPARSim simulation](../SimulationMetaSPARSim/README.md).

## Real-data figures

In the combined figure, panels A and B show signal recovery in the spike-in
experiment as sequencing depth decreases. Panels C–F examine *Coprococcus*, *Enterococcus*,
and *Faecalibacterium* in the Schubert *Clostridioides difficile* infection
(CDI) cohort:

- C compares BH-adjusted discoveries across methods.
- D shows SACAT's signed SA and RA statistics.
- E places observed detection alongside model-estimated sources of zeros.
- F shows observed positive abundances and fitted abundance conditional on
  presence, standardized over the cohort's covariates.

The estimated SA and nondetection (ND) shares average posterior zero-source
probabilities from the fitted model allowing group and covariate effects over
all samples in each group.

Numerical inputs are in `data/real_data/`. The analyses are described under
[Schubert CDI](../RealDataAnalysis/cdi_schubert/README.md) and
[Rauer spike-in benchmark](../RealDataAnalysis/rauer_mock/README.md).
The separate spike-in figure also reads
`../RealDataAnalysis/rauer_mock/table/rauer_mock_null_results.csv`.
After rerunning those analyses, refresh these inputs from the completed
results and saved CDI fit:

```sh
Rscript Analysis/MainFigures/refresh_real_data_inputs.R
```

The exported CSVs contain the plotted zero-source shares, positive abundances,
fitted abundances, full method-result families, and spike-in recovery counts,
so generating the manuscript PDFs does not require a saved fit object.

## Generate the figures

Use R 4.1 or later, `ggplot2` 3.5 or later, `patchwork` 1.3 or later, and
`systemfonts`. From the repository root:

```sh
Rscript Analysis/MainFigures/make_simulation_figure.R
Rscript Analysis/MainFigures/make_real_data_figure.R
```

| Script | Output in `figures/` |
| --- | --- |
| `make_simulation_figure.R` | `main_simulation_calibrated_bh_power_combined.pdf` |
| `make_real_data_figure.R` | `main_real_data_combined.pdf` |

### Separate simulation figures

```sh
Rscript Analysis/MainFigures/make_simulation_figure.R --split
```

The outputs are `figures/main_simulation_internal.pdf` and
`figures/main_simulation_metasparsim.pdf`.

### Separate application figures

To save the spike-in experiment and CDI cohort as separate PDFs:

```sh
Rscript Analysis/MainFigures/make_real_data_figure.R --split
```

`figures/main_spike_in.pdf` shows target recovery in A–B and negative-control
rejection rates in C–D. Recovery uses BH-adjusted discoveries. Negative-control
rates use raw p-values at 0.05, averaging the rejection fraction over all 11
taxa within each of the 111 partitions; unavailable tests count as
non-rejections. The partitions enumerate balanced splits of the same 32 libraries.

`figures/main_cdi.pdf` contains the method comparison and component displays
for the three CDI taxa, labeled A–D.
