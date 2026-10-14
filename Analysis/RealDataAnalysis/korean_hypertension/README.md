# Korean hypertension gut microbiome cohort

This analysis uses the Korean adult cohort reported by
[Song et al. (2023)](https://doi.org/10.3389/frmbi.2023.1072059). The article's
public supplementary workbook contains the ASV counts, taxonomy, and clinical
metadata used here.

## Analysis

- Comparison: 503 adults with normotension (reference) versus 120 with
  hypertension; 623 participants in total.
- Tested taxa: 73 genus or deepest-resolved higher-rank bins.
- Adjustment: `~ age_z + sex + bmi_z + group`; age and BMI are
  standardized within the final cohort.
- Sampling: one microbiome profile per participant in the cross-sectional
  workbook.
- Systolic and diastolic blood pressure define the study groups and are not
  adjustment variables.

## Data and code

Download the supplementary `Table_1.xlsx` workbook and place it in this
directory's `raw/`. Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/korean_hypertension/prepare_data.R
```

The SACAT input and summary tables are written to `processed/`. Use the
[shared analysis](../README.md) to produce results in `table/` and `figs/`.

The SACAT profile labels are “Normal BP” and “High BP”; exported tables retain
`Normotension` and `Hypertension`.
