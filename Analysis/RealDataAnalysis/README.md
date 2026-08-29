# Real-data analyses

This directory contains seven completed analyses of public human microbiome
cohorts spanning gastrointestinal and vaginal communities and several disease
or exposure contrasts. Every dataset follows the same preparation, modeling,
multiplicity, and reporting workflow.

The purpose is to determine whether an association is better described as a
change in structural absence, a change in abundance among samples in which the
taxon is present, or evidence from both components. The comparison methods
place those findings in the context of established differential abundance
workflows. DASRA's contribution in these analyses is the explicit
decomposition of a combined association into two biologically distinct signal
types.

## Completed datasets

The first group in each contrast is the reference group. Numeric variables
ending in `_z` are standardized within the final cohort.

| Dataset | Analysis contrast | Samples, reference / comparison | Tested taxa | Primary DASRA adjustment variables | Publication |
|---|---|---:|---:|---|---|
| [Baxter colorectal cancer](crc_baxter/README.md) | No-lesion control / colorectal cancer | 172 / 120 | 126 | `age_z`, sex | Baxter et al. 2016, [doi:10.1186/s13073-016-0290-3](https://doi.org/10.1186/s13073-016-0290-3) |
| [Schubert *C. difficile* infection](cdi_schubert/README.md) | Nondiarrheal healthy control / CDI | 152 / 91 | 88 | `age_z`, sex, recent antibiotics | Schubert et al. 2014, [doi:10.1128/mBio.01021-14](https://doi.org/10.1128/mBio.01021-14) |
| [GEMS pediatric diarrhea](gems_pediatric_diarrhea/README.md) | Matched control / moderate-to-severe diarrhea | 484 / 508 | 78 | `age_z`, country | Pop et al. 2014, [doi:10.1186/gb-2014-15-6-r76](https://doi.org/10.1186/gb-2014-15-6-r76) |
| [Korean hypertension](korean_hypertension/README.md) | Normotension / hypertension | 503 / 120 | 73 | `age_z`, sex, `bmi_z` | Song et al. 2023, [doi:10.3389/frmbi.2023.1072059](https://doi.org/10.3389/frmbi.2023.1072059) |
| [Zupancic Old Order Amish obesity](microbiomehd_zupancic_obesity/README.md) | Healthy category / obesity category | 88 / 95 | 41 | sex | Zupancic et al. 2012, [doi:10.1371/journal.pone.0043052](https://doi.org/10.1371/journal.pone.0043052) |
| [RISK pediatric Crohn's disease](qiita_1939_pediatric_crohn/README.md) | Non-IBD control / Crohn's disease | 101 / 143 | 212 | `age_z`, sex | Gevers et al. 2014, [doi:10.1016/j.chom.2014.02.005](https://doi.org/10.1016/j.chom.2014.02.005) |
| [Ravel vaginal microbiome](ravel_vaginal_ethnicity/README.md) | White / Black participants | 97 / 104 | 70 | group only | Ravel et al. 2011, [doi:10.1073/pnas.1002611107](https://doi.org/10.1073/pnas.1002611107) |

Each linked dataset README records the public source, cohort construction,
independent sampling unit, filtering, original library-size definition,
adjustment variables, and preparation command.

## What the completed analyses show

Across the seven datasets, DASRA identified 157 taxa with a significant
combined BH-adjusted result. Both component results were available for 144 of
these discoveries; their mechanism categories were:

| Signal category | Taxa | Percentage of classified discoveries |
|---|---:|---:|
| Structural absence only | 71 | 49.31% |
| Present-conditional abundance only | 57 | 39.58% |
| Both components | 15 | 10.42% |
| Omnibus only | 1 | 0.69% |

Thus, 128 of 144 classified discoveries (88.89%) were significant in exactly
one component analysis. Among the remaining 13 combined discoveries, only the
structural result formed for 2 taxa and only the present-conditional abundance
result formed for 11 taxa.

The number of combined discoveries was 41 for CDI, 5 for colorectal cancer,
32 for GEMS diarrhea, 55 for pediatric Crohn's disease, and 24 for the Ravel
cohort. The Korean hypertension and Zupancic obesity analyses had no combined
BH discoveries.

Among the 157 combined discoveries, 156 were also significant under at least
one comparison method. Full dataset-level counts and percentages are in
[`summary/`](summary/).

## Input construction

All analyses begin with public nonnegative integer count tables. Counts are not
rarefied and are not converted to relative abundance before modeling. Each
sample's `library_size` is the source-recorded read depth or the column total
of the complete source count table before analysis-level taxonomic aggregation
and filtering. The sum of the retained analysis taxa can therefore be smaller
than the original library size.

After the cohort and adjustment set are fixed, taxa are retained when they are
observed in at least 5% of the pooled analysis cohort. A separate support rule
requires a positive count in both comparison groups. These rules do not use
association results. Dataset-specific preparation scripts document source
quality filters, taxonomic aggregation, complete-case restrictions, and any
handling of repeated observations.

Some comparison methods require a complete composition for normalization. For
those methods only, [`analysis.R`](analysis.R) adds an `Other_unmodeled` row
equal to the part of the original library not represented by the tested taxa.
DASRA receives only the retained tested taxa together with the original
library size; `Other_unmodeled` is neither tested nor included in the DASRA
target-excluded abundance reference.

## Reproducing the workflow

Run commands from the repository root. Public source files are not committed
to version control because of their size and source-specific redistribution
terms. Download and place them exactly as described in each dataset README.

Prepare all seven inputs:

```sh
for dataset in \
  crc_baxter \
  cdi_schubert \
  gems_pediatric_diarrhea \
  korean_hypertension \
  microbiomehd_zupancic_obesity \
  qiita_1939_pediatric_crohn \
  ravel_vaginal_ethnicity
do
  Rscript "Analysis/RealDataAnalysis/${dataset}/prepare_data.R"
done
```

Run all completed-analysis configurations:

```sh
Rscript Analysis/RealDataAnalysis/analysis.R
```

Run one dataset independently:

```sh
Rscript Analysis/RealDataAnalysis/analysis.R \
  --dataset=gems_pediatric_diarrhea
```

Regenerate the cross-dataset descriptive summaries without rerunning any
statistical method:

```sh
Rscript Analysis/RealDataAnalysis/summarize_datasets.R
```

The local `raw/`, `processed/`, and `work/` directories and the optional
`R_lib/` compatibility library are excluded from version control. They are
regenerable analysis materials, not reported results. Formal reviewer-facing
files are in each dataset's `table/` and `figs/` directories and in the
cross-dataset `summary/` directory. Random seeds are fixed by dataset in
`analysis.R`, so changing the order of the dataset list does not change a
method's random-number stream.

## Methods and multiplicity

Every retained taxon remains in each result set's multiplicity family. An
unavailable result is assigned an operational p-value of one for harmonized
BH adjustment, and a discovery requires both an available result and a
BH-adjusted q-value no greater than 0.05.

The DASRA combined result is available when at least one component forms. Its
`components_used` field records whether one or both components entered the
Bonferroni minimum-p omnibus test; an unformed component contributes its
operational p-value of one.

The 15 exported result sets are:

1. DASRA structural absence, present-conditional abundance, and combined;
2. MaAsLin3 prevalence, abundance, and combined;
3. ZINQ prevalence, abundance, and combined; and
4. the primary differential abundance result from ANCOM-BC2, LinDA, corncob,
   edgeR, DESeq2, and metagenomeSeq.

Component definitions differ among methods. Comparisons across methods are
therefore descriptive comparisons of their prespecified outputs, not claims
that every method estimates the same biological quantity. Method-specific
library-size terms, offsets, and normalization conventions are recorded in the
taxon-level outputs and method-status tables.

## Per-dataset outputs

Each completed dataset contains six CSV files in `table/`:

1. `*_analysis_input_summary.csv` — sample counts, tested taxa, retained-read
   fractions, and the multiplicity definition;
2. `*_method_results_all_taxa.csv` — all taxon-level results for the 15 result
   sets, including a plain-language availability explanation;
3. `*_method_status.csv` — method completion, version, discovery-row count,
   elapsed time, and notes; for a multi-component family, the discovery count
   is the sum of significant result rows across its result sets, while
   `*_discovery_summary.csv` gives the separate count for each result set;
4. `*_discovery_summary.csv` — BH discovery counts by result set;
5. `*_upset_intersections.csv` — every discovery-set intersection; and
6. `*_upset_intersection_membership.csv` — the taxa in every intersection.

Each `figs/` directory contains two one-page PDF figures:

- `*_upset.pdf` compares discoveries across the 15 result sets. It displays
  the 20 largest intersection patterns; the two UpSet CSV files retain the
  complete intersection inventory and taxon membership.
- `*_dasra_profile.pdf` displays the DASRA structural and
  present-conditional components. It displays up to 24 significant combined
  results; when no combined discovery exists, it displays the 10 strongest
  available combined results as a descriptive profile.

## Cross-dataset summaries

[`summarize_datasets.R`](summarize_datasets.R) reads only completed result
tables and writes:

- dataset-level DASRA availability and human-readable unavailability reasons;
- structural-only, present-conditional-only, both-component, and omnibus-only
  counts among discoveries with both components available;
- the number of combined discoveries with one available component;
- descriptive overlap with the comparison methods;
- equal-dataset descriptive percentages;
- a PDF overview of signal categories by dataset.

The pooled percentage weights each classified discovery equally. The
equal-dataset percentage averages within-dataset percentages among datasets
with at least one classified combined DASRA discovery.

## Interpretation boundaries

- These are observational associations. Adjustment reduces measured
  confounding but does not establish causality.
- Taxonomic resolution and available covariates differ across public studies;
  comparisons should therefore emphasize recurring patterns rather than a
  universal effect size.
- The Ravel White/Black comparison uses the categories reported in the source
  study. It can reflect social, environmental, behavioral, and clinical
  differences and must not be interpreted as an intrinsic biological effect
  of race or ethnicity.
- Structural-absence and present-conditional abundance findings answer
  different questions. Their separation is the primary scientific purpose of
  this analysis collection.
