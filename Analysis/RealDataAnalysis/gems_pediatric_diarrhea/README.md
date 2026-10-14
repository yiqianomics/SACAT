# GEMS pediatric diarrhea

This analysis uses genus-count data from the Global Enteric Multicenter Study
reported by [Pop et al. (2014)](https://doi.org/10.1186/gb-2014-15-6-r76).
The cohort includes children from Bangladesh, Kenya, Mali, and The Gambia.
Data are available in the authors'
[MSD1000 repository](https://github.com/jnpaulson/MSD1000).

## Analysis

- Comparison: 484 matched controls (reference) versus 508 children with
  moderate-to-severe diarrhea (MSD); 992 children in total.
- Tested taxa: 78 published genera.
- Adjustment: `~ age_z + country + group`, where age in months is
  standardized within the final cohort and country is categorical.
- Sampling: one published observation per enrolled child.
- Library size: the published `totalCounts` value, before genus filtering.

## Data and code

Download `forserveroptim.rdata` from MSD1000 and place it in this directory's
`raw/`. It contains the published genus counts, metadata, and library sizes.
Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/gems_pediatric_diarrhea/prepare_data.R
```

The SACAT input and count, metadata, taxonomy, and preprocessing tables are
written to `processed/`. Use the [shared analysis](../README.md) to fit the
methods and produce results in `table/` and `figs/`.
