# Real-data analyses

These analyses examine group differences in seven public human microbiome
cohorts and a spike-in experiment with known differences between communities.

SACAT's structural-absence (SA), conditional relative-abundance (RA), and
omnibus results are compared with other microbiome association methods.

## Human cohorts

The first group in each contrast is the reference group. Numeric variables
ending in `_z` are standardized within the final cohort.

| Dataset | Reference / comparison | Samples | Taxa |
|---|---|---:|---:|
| [Baxter colorectal cancer](crc_baxter/README.md) | No-lesion control / colorectal cancer | 172 / 120 | 126 |
| [Schubert CDI](cdi_schubert/README.md) | Healthy control / *C. difficile* infection | 152 / 91 | 88 |
| [GEMS pediatric diarrhea](gems_pediatric_diarrhea/README.md) | Matched control / moderate-to-severe diarrhea | 484 / 508 | 78 |
| [Korean hypertension](korean_hypertension/README.md) | Normotension / hypertension | 503 / 120 | 73 |
| [Amish obesity](microbiomehd_zupancic_obesity/README.md) | Healthy category / obesity category | 88 / 95 | 41 |
| [RISK pediatric Crohn's disease](qiita_1939_pediatric_crohn/README.md) | Non-IBD control / Crohn's disease | 101 / 143 | 212 |
| [Ravel vaginal microbiome](ravel_vaginal_ethnicity/README.md) | White / Black participants | 97 / 104 | 70 |

Each dataset README provides the source, cohort selection, adjustment model,
and preparation command. The human-cohort analyses use unrarefied counts and original library
sizes recorded before taxon filtering. Tested taxa occur in at least 5% of the
final cohort and have positive counts in both groups.

## Spike-in benchmark

The [Rauer spike-in benchmark](rauer_mock/README.md) uses 32 libraries from two
defined communities to assess SA/RA signal recovery, false positives, and test
availability as sequencing depth decreases.

## Comparisons

The human-cohort outputs contain 15 result sets: SA, RA, and combined tests from
SACAT; prevalence, abundance, and combined tests from MaAsLin3 and ZINQ; and
abundance tests from ANCOM-BC2, LinDA, corncob, edgeR, DESeq2, and metagenomeSeq.
Benjamini–Hochberg (BH) adjustment is applied separately to each result set over
all tested taxa. Discoveries use an adjusted p-value at or below 0.05. SACAT's
combined result uses the Bonferroni minimum-p test.

## Results

Each human cohort has numerical results in `table/`:

- `*_method_results_all_taxa.csv`: taxon-level results for all 15 result sets.
- `*_discovery_summary.csv`: discovery counts for each result set.
- `*_upset_intersections.csv` and `*_upset_intersection_membership.csv`:
  discovery overlap and the taxa in each intersection.
- `*_analysis_input_summary.csv` and `*_method_status.csv`: cohort summaries
  and method-level information, including software versions.

The `figs/` directories contain SACAT component profiles and method-overlap
figures. The [CDI analysis](cdi_schubert/README.md) also provides detailed
component examples and a ten-result-set intersection figure. The
[main manuscript figures](../MainFigures/README.md) combine CDI profiles with
spike-in recovery curves.

The `summary/` directory contains cross-cohort tables and PDFs of test
availability, SA/RA discovery categories, and method overlap.

## Run the analyses

Run commands from the repository root. First download and prepare the inputs
using each dataset's README. Then fit the seven human cohorts and summarize
their results:

```sh
Rscript Analysis/RealDataAnalysis/analysis.R
Rscript Analysis/RealDataAnalysis/summarize_datasets.R
```

To fit one cohort, add its directory name, for example:

```sh
Rscript Analysis/RealDataAnalysis/analysis.R --dataset=cdi_schubert
```

Use `--plot-only` to regenerate intersection tables and the six non-CDI UpSet
figures from saved results. Detailed CDI figures and the spike-in analysis have
their own commands in the linked dataset READMEs.

Source files belong in `raw/`; prepared inputs and intermediate files are
written to `processed/` and `work/`. Fitting requires SACAT 0.7.1 or later.
Versions used for the saved results are recorded in
[`software_versions.csv`](software_versions.csv), the cohort method-status
tables, and the [spike-in version table](rauer_mock/software_versions.csv).
