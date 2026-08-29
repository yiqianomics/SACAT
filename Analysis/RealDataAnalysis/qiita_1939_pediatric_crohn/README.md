# RISK pediatric Crohn's disease

This dataset comes from the RISK cohort reported by Gevers et al. in “The
treatment-naive microbiome in new-onset Crohn's disease”
([doi:10.1016/j.chom.2014.02.005](https://doi.org/10.1016/j.chom.2014.02.005)).
The public data are available as Qiita study 1939 and from the European
Nucleotide Archive under PRJEB13679.

Place the taxonomy-annotated public BIOM table at `raw/qiita_1939.biom` and the
matching Qiita metadata at `raw/qiita_1939.tsv`. These are the files written by
a redbiom Qiita-study fetch with output basename `qiita_1939` and taxonomy
enabled.

## Final prepared analysis

- **Comparison:** 101 non-IBD controls (reference) versus 143 children with
  Crohn's disease (comparison), 244 children in total.
- **Tested panel:** 212 genus or deepest-resolved higher-rank taxon bins.
- **Adjustment model:** `~ age_z + sex + group`, with age standardized within
  the final cohort.
- **Independent unit:** one rectal mucosal sample per child, retaining the
  largest original library size when repeated eligible samples are available;
  sample identifier breaks an exact tie. The tie-breaking rule is independent
  of diagnosis and association results.
- **Original library size:** the column sum of the complete downloaded BIOM
  feature table before taxonomic aggregation or any taxon filter.

The prespecified cohort contains RISK rectal mucosal samples from children
recorded as having Crohn's disease (`diagnosis = CD`) or no IBD
(`diagnosis = no`). Samples recorded as collected during antibiotic exposure
are excluded. When a child has more than one eligible sample, the sample with
the largest original library size is retained. The model compares Crohn's
disease with the non-IBD control group and adjusts for standardized age and
sex. Age is standardized after cohort construction and subject deduplication.

Run the preparation from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/qiita_1939_pediatric_crohn/prepare_data.R
```

The script does not rarefy counts. It preserves each sample's full downloaded
BIOM total as the original library size, aggregates features to genus or the
deepest resolved higher rank, retains taxa present in at least 5% of the final
cohort and in both groups, and writes the analysis RDS and readable audit CSV
files to `processed/`.
