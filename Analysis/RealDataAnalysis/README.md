# Real-data analyses

This directory contains seven completed analyses of public human microbiome
cohorts spanning gastrointestinal and vaginal communities and several disease
or exposure contrasts. Every dataset follows the same preparation, modeling,
multiplicity, and reporting workflow.

## Completed datasets

The first group in each contrast is the reference group. Numeric variables
ending in `_z` are standardized within the final cohort.

| Dataset | Analysis contrast | Samples, reference / comparison | Tested taxa | Primary SACAT adjustment variables | Publication |
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

## Defined mock-community benchmark

The [Rauer mock-community benchmark](rauer_mock/README.md) uses 32 libraries
from two defined source communities to evaluate component-specific recovery,
calibration, and availability across native and controlled sequencing depths.
The cross-cohort summaries cover the seven human cohorts.

## Human-cohort input construction

The seven human-cohort analyses use public unrarefied integer count tables.
Each sample's `library_size` is the source-recorded read depth or the column
total of the complete source count table before analysis-level taxonomic
aggregation and filtering. The sum of the retained analysis taxa can therefore
be smaller than the original library size.

After the cohort and adjustment set are fixed, taxa are retained when they are
observed in at least 5% of the pooled analysis cohort. A separate support rule
requires a positive count in both comparison groups. Both retention rules are
applied before association testing. Dataset-specific preparation scripts
document source quality filters, taxonomic aggregation, complete-case
restrictions, and any handling of repeated observations.

Some comparison methods require a complete composition for normalization.
For these methods, [`analysis.R`](analysis.R) adds an `Other_unmodeled` row
equal to the part of the original library represented outside the tested taxa.
SACAT receives the retained tested taxa and the original library size. Its
target-excluded abundance reference is formed from the retained tested taxa.

## Reproducing the workflow

Run commands from the repository root. Download and place the public source
files exactly as described in each dataset README.

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

Regenerate the seven standard UpSet figures and their intersection tables from
the completed taxon-level result tables:

```sh
Rscript Analysis/RealDataAnalysis/analysis.R --plot-only
```

Regenerate the cross-dataset descriptive summaries without rerunning any
statistical method:

```sh
Rscript Analysis/RealDataAnalysis/summarize_datasets.R
```

After preparing and analyzing the Schubert CDI dataset, regenerate its
detailed component and method-intersection outputs with:

```sh
Rscript Analysis/RealDataAnalysis/cdi_schubert/make_detailed_figures.R
```

Prepare, analyze, and summarize the defined mock-community benchmark with:

```sh
Rscript Analysis/RealDataAnalysis/rauer_mock/prepare_data.R
Rscript Analysis/RealDataAnalysis/rauer_mock/analysis.R
Rscript Analysis/RealDataAnalysis/rauer_mock/make_figures.R
```

Preparation and fitting write source inputs and intermediate files to `raw/`,
`processed/`, and `work/`. Reported tables and figures are in each dataset's
`table/` and `figs/` directories and in the cross-dataset `summary/` directory.
Dataset-specific fixed seeds provide reproducible random-number streams
independent of execution order.

[`software_versions.csv`](software_versions.csv) records the R version and
direct package versions used to generate the reported outputs. Each
`*_method_status.csv` also records the fitted method package version for that
dataset.

## Human-cohort methods and multiplicity

Every retained taxon remains in each result set's multiplicity family. An
unavailable result is assigned an operational p-value of one for harmonized
BH adjustment, and a discovery requires both an available result and a
BH-adjusted p-value no greater than 0.05.

The SACAT combined result is available when at least one component forms. Its
`components_used` field records whether one or both components entered the
Bonferroni minimum-p omnibus test; an unformed component contributes its
operational p-value of one.

The 15 exported result sets are:

1. SACAT structural absence, present-conditional abundance, and combined;
2. MaAsLin3 prevalence, abundance, and combined;
3. ZINQ prevalence, abundance, and combined; and
4. the primary differential abundance result from ANCOM-BC2, LinDA, corncob,
   edgeR, DESeq2, and metagenomeSeq.

Component definitions differ among methods. Cross-method comparisons summarize
these prespecified outputs while retaining their method-specific estimands.
Method-specific library-size terms, offsets, and normalization conventions are
recorded in the taxon-level outputs and method-status tables.

## Human-cohort outputs

Each completed human cohort contains six standard CSV files in `table/`:

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

Each `figs/` directory contains two standard PDF figures:

- `*_upset.pdf` compares discoveries across the 15 result sets. Every
  nonempty intersection pattern is displayed in one continuous matrix.
- `*_sacat_profile.pdf` displays the SACAT structural and
  present-conditional components. It displays up to 24 significant combined
  results; when no combined discovery exists, it displays the 10 strongest
  available combined results as a descriptive profile.

The Schubert CDI directory also contains a detailed component figure, a
10-result-set method-intersection figure, and three corresponding data tables:

- `schubert_cdi_component_examples.pdf` and
  `schubert_cdi_component_examples.csv`;
- `schubert_cdi_method_intersections.pdf` and
  `schubert_cdi_method_intersections.csv`; and
- `schubert_cdi_method_intersection_membership.csv`.

## Cross-dataset summaries

[`summarize_datasets.R`](summarize_datasets.R) reads only completed result
tables and writes:

- dataset-level SACAT availability and human-readable unavailability reasons;
- pooled availability for the primary result from each method, using the common
  taxon-by-dataset result grid;
- SACAT component and omnibus unavailability reasons, with counts, percentages
  of all results, and percentages among unavailable results;
- structural-only, present-conditional-only, both-component, and omnibus-only
  counts among discoveries with both components available;
- the number of combined discoveries with one available component;
- descriptive overlap with the comparison methods;
- dataset-averaged descriptive percentages;
- PDF overviews of method availability and signal categories by dataset.

The availability outputs are `method_availability.csv`,
`sacat_unavailability_reasons.csv`, and `method_availability.pdf`. MaAsLin3 and
ZINQ contribute their package-provided combined results; SACAT contributes its
omnibus result, and the other methods contribute their primary
differential-abundance results.

Availability percentages weight each of the 688 prespecified taxon-by-dataset
results equally. For signal-category summaries, the pooled percentage weights
each classified discovery equally. The dataset-averaged percentage averages
within-dataset percentages among datasets with at least one classified
combined SACAT discovery.
