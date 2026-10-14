# Schubert *Clostridioides difficile* infection stool cohort

This analysis uses the stool cohort reported by
[Schubert et al. (2014)](https://doi.org/10.1128/mBio.01021-14). The count table
and metadata are available in MicrobiomeHD version 3,
[Zenodo record 1146764](https://doi.org/10.5281/zenodo.1146764).

## Analysis

- Comparison: 152 nondiarrheal healthy controls (reference) versus 91
  participants with *C. difficile* infection; 243 participants in total.
- Tested taxa: 88 genus or deepest-resolved higher-rank bins.
- Adjustment: `~ age_z + sex + antibiotics_3mo + group`, with age
  standardized within the final cohort.
- Selection: one stool profile per participant, at least 1,450 reads, and
  complete adjustment variables. One healthy sample with missing sex is excluded.

## Data and code

Download `cdi_schubert_results.tar.gz` from Zenodo and extract it into this
directory's `raw/cdi_schubert_results/`. Preparation uses the RDP-assigned 100%
de novo feature table and accompanying metadata. Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/cdi_schubert/prepare_data.R
```

The SACAT input and cohort, count, taxonomy, and preprocessing summaries are
written to `processed/`. After the [shared analysis](../README.md), generate
selected-taxon component profiles and discovery overlap across ten result sets:

```sh
Rscript Analysis/RealDataAnalysis/cdi_schubert/make_detailed_figures.R
```

Outputs are `figs/schubert_cdi_component_examples.pdf` and
`figs/schubert_cdi_method_intersections.pdf`, with numerical summaries and
intersection membership in `table/`.

The combined CDI and spike-in figure is documented in
[MainFigures](../../MainFigures/README.md).
