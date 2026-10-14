# Baxter colorectal cancer stool cohort

This analysis uses the colorectal-cancer stool cohort reported by
[Baxter et al. (2016)](https://doi.org/10.1186/s13073-016-0290-3). The count
table and metadata are available in MicrobiomeHD version 3,
[Zenodo record 1146764](https://doi.org/10.5281/zenodo.1146764).

## Analysis

- Comparison: 172 participants without a colorectal lesion (reference)
  versus 120 with colorectal cancer; 292 participants in total.
- Tested taxa: 126 genus or deepest-resolved higher-rank bins.
- Adjustment: `~ age_z + sex + group`, with age standardized within
  the final cohort.
- Sampling: one stool profile per participant in the cross-sectional cohort.

## Data and code

Download `crc_baxter_results.tar.gz` from Zenodo and extract it into this
directory's `raw/crc_baxter_results/`. Preparation uses the RDP-assigned 100%
de novo feature table and accompanying metadata. Unlike the publication's
97% OTU table rarefied to 10,000 reads, this analysis uses unrarefied counts.

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/crc_baxter/prepare_data.R
```

The SACAT input and cohort, count, taxonomy, and preprocessing summaries are
written to `processed/`. Use the [shared analysis](../README.md) to fit the
methods and produce results in `table/` and `figs/`.
