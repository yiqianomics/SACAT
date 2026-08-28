# Baxter colorectal cancer stool cohort

This analysis uses the colorectal-cancer cohort reported by Baxter et al. in
“Microbiota-based model improves the sensitivity of fecal immunochemical test
for detecting colonic lesions” ([doi:10.1186/s13073-016-0290-3](https://doi.org/10.1186/s13073-016-0290-3)).
The public count table and metadata are distributed in MicrobiomeHD version 3,
[Zenodo record 1146764](https://doi.org/10.5281/zenodo.1146764), as
`crc_baxter_results.tar.gz`.

## Final prepared analysis

- **Comparison:** 172 participants without a colorectal lesion (reference)
  versus 120 participants with colorectal cancer (comparison), 292
  participants in total. Adenoma groups are not part of this contrast.
- **Tested panel:** 126 genus or deepest-resolved higher-rank taxon bins.
- **Adjustment model:** `~ age_z + sex + group`, with age standardized within
  the final cohort.
- **Independent unit:** one stool profile per retained participant in the
  public cross-sectional cohort.
- **Original library size:** the column sum of the complete MicrobiomeHD count
  table before feature filtering, taxonomic aggregation, prevalence
  filtering, or the two-group support filter.

Download `crc_baxter_results.tar.gz` from the Zenodo record and extract it as
`raw/crc_baxter_results/`. The preparation uses the RDP-assigned 100% de novo
feature table and the accompanying metadata. The published analysis used a
different 97% OTU table and rarefied samples to 10,000 reads; this analysis
does not rarefy because the original library size is modeled explicitly.

The preparation applies the documented MicrobiomeHD count filters, retains
the no-lesion and colorectal-cancer groups, aggregates features to genus or
the deepest resolved higher rank, and retains taxa observed in at least 5% of
the final cohort and in both groups.

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/crc_baxter/prepare_data.R
```

The DASRA input and readable cohort, count, taxonomy, and preprocessing audit
tables are written to `processed/`.
