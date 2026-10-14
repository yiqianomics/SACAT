# Ravel cross-sectional vaginal microbiome cohort

This analysis uses the cross-sectional cohort reported by
[Ravel et al. (2011)](https://doi.org/10.1073/pnas.1002611107). The original
study enrolled 396 women; the public count table contains 394 microbiome
profiles.

## Analysis

- Comparison: 97 White participants (reference) versus 104 Black participants;
  201 participants in total.
- Tested taxa: 70 published taxon labels.
- Adjustment: `~ group`; no covariates are included. Vaginal pH, Nugent score,
  and community state type are treated as microbiome-related outcomes and
  omitted from the adjustment set.
- Sampling: one vaginal profile per participant.
- Library size: the published `Depth` field, which equals the full count-table
  column sum.

## Data and code

Place the published taxon counts at `raw/counts.tsv` and participant metadata
at `raw/metadata.tsv` within this directory. These files are available in the
paper's supplementary material and mirrored in the
[ETH Zurich Microbiome Data Analysis workshop](https://www.gdc-docs.ethz.ch/MDA/reports/ravel.html).

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/ravel_vaginal_ethnicity/prepare_data.R
```

The SACAT input and summary tables are written to `processed/`. Use the
[shared analysis](../README.md) to fit the methods and produce results in
`table/` and `figs/`.
