# Korean hypertension gut microbiome cohort

This analysis uses the 120 adults with hypertension and 503 adults with
normotension reported in “The association between gut microbiome and
hypertension varies according to enterotypes: a Korean study”
([doi:10.3389/frmbi.2023.1072059](https://doi.org/10.3389/frmbi.2023.1072059)).
The public supplementary workbook contains the ASV count table, taxonomy, and
clinical metadata used here.

## Final prepared analysis

- **Comparison:** 503 adults with normotension (reference) versus 120 adults
  with hypertension (comparison), 623 participants in total.
- **Tested panel:** 73 genus or deepest-resolved higher-rank taxon bins.
- **Adjustment model:** `~ age_z + sex + bmi_z + group`; age and BMI are
  standardized within the final cohort.
- **Independent unit:** the public cross-sectional workbook contributes one
  microbiome sample per participant; no repeated observations are retained.
- **Original library size:** the sum of all ASV counts for the sample in the
  published unrarefied table, before taxonomic aggregation, prevalence
  filtering, or the two-group support filter.

Download the article's supplementary `Table_1.xlsx` workbook and place it at
`raw/Table_1.xlsx`, then run
from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/korean_hypertension/prepare_data.R
```

Normotension is the reference group and hypertension is the comparison group.
Age, sex, and BMI are included as adjustment variables. Systolic and diastolic
blood pressure are not included because they define the study group. The
script keeps the original, unrarefied ASV library totals, aggregates ASVs to
genus or the deepest resolved higher rank, and retains taxa present in at
least 5% of the cohort and with positive counts in both groups.

The reproducible analysis RDS and readable audit tables are written to
`processed/`.

For readability, the DASRA profile figure uses the display labels “Normal BP”
and “High BP”; the model and exported tables retain the formal group labels
`Normotension` and `Hypertension`.
