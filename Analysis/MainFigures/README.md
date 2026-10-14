# Main manuscript figures

These two figures bring together the simulation studies and the real-data
applications. Each script creates one PDF from the saved numerical results.

## Simulation figure

Panel A shows internal model validation at 80 samples per group; Panel B shows
the metaSPARSim simulations. Together they display component separation and
signal recovery under model-based and externally generated data.

The figure compares structural-absence (SA), relative-abundance (RA), and
joint results. Component panels use the corresponding SACAT, ZINQ, and
MaAsLin 3 tests; joint panels also include six abundance methods. Rows
distinguish settings without and with a covariate. Plotted values use
matched-null p-value calibration followed by Benjamini–Hochberg (BH)
adjustment at 0.05.

Inputs are the `tables/main_plot_data.csv` files in
[Internal model validation](../SimulationInternalValidation/README.md) and
[metaSPARSim simulation](../SimulationMetaSPARSim/README.md).

## Real-data figure

Panels A and B show signal recovery in the defined spike-in experiment as
sequencing depth decreases. Panels C–F examine *Coprococcus*, *Enterococcus*,
and *Faecalibacterium* in the Schubert *Clostridioides difficile* infection
(CDI) cohort:

- C compares BH-adjusted discoveries across methods.
- D shows SACAT's signed SA and RA statistics.
- E places observed detection alongside model-estimated sources of zeros.
- F shows observed positive abundances and fitted abundance conditional on
  presence, standardized over the cohort's covariates.

Numerical inputs are in `data/real_data/`. The analyses are described under
[Schubert CDI](../RealDataAnalysis/cdi_schubert/README.md) and
[Rauer spike-in benchmark](../RealDataAnalysis/rauer_mock/README.md).

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
