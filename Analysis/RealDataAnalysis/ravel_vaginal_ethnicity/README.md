# Ravel cross-sectional vaginal microbiome cohort

This analysis uses the cross-sectional cohort reported by Ravel et al. in
“Vaginal microbiome of reproductive-age women”
([doi:10.1073/pnas.1002611107](https://doi.org/10.1073/pnas.1002611107)).
The original study enrolled 396 women, and the public count table contains 394
microbiome profiles. The prespecified binary comparison retains 97 White
participants as the reference group and 104 Black participants as the
comparison group.

## Final prepared analysis

- **Comparison:** 97 White participants (reference) versus 104 Black
  participants (comparison), 201 participants in total.
- **Tested panel:** 70 published taxon labels after the 5% prevalence and
  two-group support filters.
- **Adjustment model:** `~ group`; no covariates are included.
- **Independent unit:** one vaginal sample per participant, as supplied by the
  cross-sectional source dataset.
- **Original library size:** the exact column sum of the published unrarefied
  count table, before prevalence and two-group support filtering. This agrees
  with the source `Depth` field.

The preparation requires the published taxon counts at `raw/counts.tsv` and
participant metadata at `raw/metadata.tsv`. The accompanying
`raw/taxonomy.tsv` is retained for source documentation; the preparation does
not need it because the count-table rows already carry the published taxon
labels. These files are available in the paper's
supplementary material and are also mirrored in the
[ETH Zurich Microbiome Data Analysis workshop](https://www.gdc-docs.ethz.ch/MDA/reports/ravel.html).
The count table is unrarefied, and its exact column total is retained as each
participant's library size.

The primary model contains the ethnicity-group indicator. Vaginal pH, Nugent
score, and community state type are not adjustment variables because they are
microbiome-related outcomes rather than pre-exposure confounders. Taxa are
retained when they are observed in at least 5% of the selected cohort and in
both groups.

From the repository root, run:

```sh
Rscript Analysis/RealDataAnalysis/ravel_vaginal_ethnicity/prepare_data.R
```

The resulting RDS and readable audit tables are written to `processed/`.
