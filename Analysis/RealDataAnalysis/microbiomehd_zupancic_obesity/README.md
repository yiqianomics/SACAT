# Zupancic Old Order Amish obesity

This analysis uses the Old Order Amish cohort reported by
[Zupancic et al. (2012)](https://doi.org/10.1371/journal.pone.0043052).
The count table and metadata are available in MicrobiomeHD version 3,
[Zenodo record 1146764](https://zenodo.org/records/1146764). The underlying
sequence study is NCBI SRA project SRP002465.

## Analysis

- Comparison: 88 participants in the healthy category (H; reference)
  versus 95 in the obesity category (OB); 183 participants in total.
- Tested taxa: 41 genus or deepest-resolved higher-rank bins.
- Adjustment: `~ sex + group`, with sex treated as categorical.
- Selection: one visit-1 sample per participant. Among repeated eligible
  samples, the largest original library is retained; sample identifier breaks
  ties.

## Data and code

Download `ob_zupancic_results.tar.gz` from Zenodo and extract it into this
directory's `raw/ob_zupancic_results/`. Keep the archive's structure, including
`RDP/ob_zupancic.otu_table.100.denovo.rdp_assigned` and
`ob_zupancic.metadata.txt`. MicrobiomeHD count filtering is applied before
selecting the analysis cohort.

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/microbiomehd_zupancic_obesity/prepare_data.R
```

The SACAT input and count, metadata, taxonomy, and preprocessing tables are
written to `processed/`. Use the [shared analysis](../README.md) to fit the
methods and produce results in `table/` and `figs/`.
