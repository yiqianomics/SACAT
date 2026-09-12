# Zupancic Old Order Amish obesity

This analysis uses the Old Order Amish cohort from Zupancic et al., *PLOS ONE*
2012
([doi:10.1371/journal.pone.0043052](https://doi.org/10.1371/journal.pone.0043052)).
The processed count table and metadata are distributed in MicrobiomeHD v3,
[Zenodo record 1146764](https://zenodo.org/records/1146764), as
`ob_zupancic_results.tar.gz`. The underlying sequence study is NCBI SRA project
SRP002465.

## Final prepared analysis

- **Comparison:** 88 participants in the healthy category (H; reference)
  versus 95 participants in the obesity category (OB; comparison), 183
  participants in total.
- **Tested panel:** 41 genus or deepest-resolved higher-rank taxon bins.
- **Adjustment model:** `~ sex + group`, with sex treated as categorical.
- **Independent unit:** one visit-1 sample per participant. When more than one
  eligible sample was available, the sample with the largest original library
  size was selected, with sample identifier used to break ties.
- **Original library size:** the column sum of the complete RDP-assigned de
  novo count table, calculated before the MicrobiomeHD feature filter,
  taxonomic aggregation, prevalence filtering, or the two-group support
  filter.

Download `ob_zupancic_results.tar.gz` from the Zenodo record and extract it as
`raw/ob_zupancic_results/`. The preparation requires
`raw/ob_zupancic_results/RDP/ob_zupancic.otu_table.100.denovo.rdp_assigned`
and `raw/ob_zupancic_results/ob_zupancic.metadata.txt`. The primary
comparison is obesity versus the healthy category at visit 1. The model
adjusts for sex. Visit-1 records are reduced to one sample per public
participant using the rule above.

The standard MicrobiomeHD count filter is applied to the complete cohort. The
analysis uses unrarefied counts. Features are aggregated to genus or the
deepest resolved higher rank, retained when present in at least 5% of the final
cohort, and required to occur in both groups.

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/microbiomehd_zupancic_obesity/prepare_data.R
```

The SACAT input and readable count, metadata, taxonomy, and preprocessing
tables are written to `processed/`.
