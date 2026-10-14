# RISK pediatric Crohn's disease

This analysis uses the RISK cohort reported by
[Gevers et al. (2014)](https://doi.org/10.1016/j.chom.2014.02.005).
Public data are available as Qiita study 1939 and European Nucleotide Archive
project PRJEB13679.

## Analysis

- Comparison: 101 non-IBD controls (reference) versus 143 children with
  Crohn's disease; 244 children in total.
- Tested taxa: 212 genus or deepest-resolved higher-rank bins.
- Adjustment: `~ age_z + sex + group`, with age standardized within
  the final cohort.
- Selection: RISK rectal mucosal samples with `diagnosis = no` or `CD` and
  `antibiotics = FALSE`. One sample per child is retained, choosing the largest
  original library and breaking ties by sample identifier.

## Data and code

Fetch the taxonomy-annotated BIOM table and matching metadata using redbiom
0.3.9 with the context below. Run from the repository root:

```sh
redbiom fetch qiita-study --study-id 1939 --context Pick_closed-reference_OTUs-Greengenes-Illumina-16S-V4-90nt-44feac --resolve-ambiguities most-reads --fetch-taxonomy --remove-blanks --output-basename Analysis/RealDataAnalysis/qiita_1939_pediatric_crohn/raw/qiita_1939
```

This writes `raw/qiita_1939.biom` and `raw/qiita_1939.tsv` within this dataset
directory. Prepare the inputs with:

```sh
Rscript Analysis/RealDataAnalysis/qiita_1939_pediatric_crohn/prepare_data.R
```

The SACAT input and summary tables are written to `processed/`. Use the
[shared analysis](../README.md) to fit the methods and produce results in
`table/` and `figs/`.
