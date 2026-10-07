# Schubert *Clostridioides difficile* infection stool cohort

This analysis uses the cohort reported by Schubert et al. in “Microbiome data
distinguish patients with *Clostridioides difficile* infection and
non-*C. difficile*-associated diarrhea from healthy controls”
([doi:10.1128/mBio.01021-14](https://doi.org/10.1128/mBio.01021-14)). The public
count table and metadata are distributed in MicrobiomeHD version 3,
[Zenodo record 1146764](https://doi.org/10.5281/zenodo.1146764), as
`cdi_schubert_results.tar.gz`.

## Final prepared analysis

- **Comparison:** 152 nondiarrheal healthy controls (reference) versus 91
  participants with *C. difficile* infection (comparison), 243 participants
  in total. The contrast contains the healthy-control and *C. difficile*
  infection groups.
- **Tested panel:** 88 genus or deepest-resolved higher-rank taxon bins.
- **Adjustment model:** `~ age_z + sex + antibiotics_3mo + group`, with age
  standardized within the final cohort.
- **Independent unit:** one stool profile per retained participant in the
  public cohort.
- **Original library size:** the column sum of the complete MicrobiomeHD count
  table before feature filtering, taxonomic aggregation, prevalence
  filtering, or the two-group support filter.

Download `cdi_schubert_results.tar.gz` from the Zenodo record and extract it
as `raw/cdi_schubert_results/`. The preparation uses the RDP-assigned 100% de
novo feature table and the accompanying metadata. It retains the publication's
minimum depth of 1,450 reads and the unrarefied counts. One healthy sample with
missing sex is excluded from the complete-case model.

Features are aggregated to genus or the deepest resolved higher rank. Taxa
must be observed in at least 5% of the final cohort and have positive counts in
both groups.

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/cdi_schubert/prepare_data.R
```

The SACAT input and cohort, count, taxonomy, and preprocessing summaries are
written to `processed/`.

## CDI and mock-community figure

`figs/schubert_cdi_rauer_mock.pdf` combines recovery curves from the Rauer
mock-community benchmark with method comparisons and component profiles for
*Coprococcus*, *Enterococcus*, and *Faecalibacterium* in the CDI cohort.

Regenerate the PDF from the included numerical inputs, without rerunning model
fits:

```sh
Rscript Analysis/RealDataAnalysis/cdi_schubert/make_schubert_cdi_rauer_mock_figure.R
```

The script uses `ggplot2` and `patchwork`, with Arial for figure text, and writes
only the PDF. Its CDI inputs are in `table/`:

- `schubert_cdi_figure_method_results.csv` contains all 88 taxa in each
  method-result family, including the p-values, BH-adjusted p-values, and signed
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

The mock-community panels read
`../rauer_mock/table/rauer_mock_recovery.csv`, which contains discovery totals
over 52 allocations for eight structural targets or three abundance targets at
each depth.
