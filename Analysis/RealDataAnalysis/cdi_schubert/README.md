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
  in total. The non-*C. difficile* diarrhea group is not part of this contrast.
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
minimum depth of 1,450 reads but does not rarefy the counts. One healthy sample
with missing sex is excluded from the complete-case model.

Features are aggregated to genus or the deepest resolved higher rank. Taxa
must be observed in at least 5% of the final cohort and have positive counts in
both groups.

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/cdi_schubert/prepare_data.R
```

The DASRA input and readable cohort, count, taxonomy, and preprocessing audit
tables are written to `processed/`.
