# Real-data negative-control analysis

This directory contains reproducible negative-control analyses for 20 public microbiome datasets. Each dataset directory stores the standardized input, an executable entry point, and the complete result tables. Shared analysis code is kept at the top level, and `combined/analysis.R` validates the full experiment and produces the cross-dataset comparison.

## Directory structure

Each dataset directory contains:

- `analysis_input.rds`: the complete standardized input used in the experiment, including 200 sample-by-genus compositions, aligned sample identifiers and metadata, source information, and the 30-genus analysis panel;
- `analysis.R`: the dataset-specific entry point;
- `results/`: taxon-level results, randomization summaries, diagnostics, depth summaries, and the run manifest.

`run_negative_control.R` is the shared analysis engine used by every dataset-specific `analysis.R` entry point. It generates the randomized counts, runs DASRA, ZINQ, and MaAsLin3 under the common design, validates each completed unit, and writes the result tables.

The stored `analysis_input.rds` files are sufficient to rerun the analyses. They contain processed genus-level analysis data rather than raw sequencing reads, so no download or preparation step is required for ordinary replication.

## Datasets

Each analysis input represents 200 independent empirical profiles. When a public resource contains repeated observations, one profile per participant is represented. Exact sample identifiers and source metadata are stored in the corresponding `analysis_input.rds` file.

| Dataset | Public source, study overview, and stored analysis input |
|---|---|
| `NogueraJulianHIV` | Stool 16S profiles from the [Noguera-Julian HIV microbiome study](https://doi.org/10.1016/j.ebiom.2016.01.032) and its [public data archive](https://doi.org/10.5281/zenodo.1146764); 200 baseline profiles from 240 Barcelona and Stockholm participants. |
| `BaxterE_2016` | Genus-level 16S stool data from the [Baxter colorectal-cancer study](https://doi.org/10.1186/s13073-016-0290-3); 200 profiles from the 292-sample primary cohort. |
| `HMPV35Throat` | Public Human Microbiome Project V3-5 throat 16S data distributed by [HMP16SData](https://waldronlab.io/HMP16SData/reference/V35.html); 200 profiles sampled from a 202-participant earliest-visit cohort. |
| `Qiita13631ASD` | Public [Qiita study 13631](https://qiita.ucsd.edu/public/?study_id=13631), a gut microbiome study of children with autism spectrum disorder and controls; 200 fecal profiles from 286 participants. |
| `ArtPrize_2015_Forehead` | Central-forehead 16S profiles from the cross-sectional [ArtPrize study](https://doi.org/10.1128/mBio.00839-19) and its [public phyloseq object](https://doi.org/10.6084/m9.figshare.8127287.v1); 200 profiles from 372 participants. |
| `MehtaRS_2018` | Stool relative-abundance profiles from curatedMetagenomicData and the [Mehta study of fecal-microbiome stability in adult men](https://pubmed.ncbi.nlm.nih.gov/29335554/); 200 profiles, one per participant, from 308 participants represented by 928 stool assays. |
| `VatanenT_2016` | Stool relative-abundance profiles from curatedMetagenomicData and the [DIABIMMUNE study of early-life microbiome LPS immunogenicity](https://pubmed.ncbi.nlm.nih.gov/27133167/); 200 profiles, one per participant, from 212 infants represented by 785 stool assays. |
| `LeChatelierE_2013` | Stool relative-abundance profiles from curatedMetagenomicData and the [Le Chatelier gut microbial richness study](https://pubmed.ncbi.nlm.nih.gov/23985870/); 200 profiles from 292 independent stool samples. |
| `ORIGINS_2022_Healthy_Plaque` | Healthy-site subgingival-plaque 16S profiles from the public [ORIGINS study](https://doi.org/10.1038/s41522-022-00289-w); 200 profiles from 668 participants, with one healthy-site profile per participant. |
| `NHANESOral_2011_2012` | CDC [NHANES Oral Microbiome](https://wwwn.cdc.gov/Nchs/Nhanes/omp/) genus-level oral-rinse data from the 2011-2012 survey cycle; 200 profiles from 4,887 participants with unique SEQN identifiers. |
| `NielsenHB_2014` | Stool relative-abundance profiles from curatedMetagenomicData and the [Nielsen metagenome assembly study](https://pubmed.ncbi.nlm.nih.gov/24997787/); 200 profiles from 318 participants. |
| `NHANESOral_2009_2010` | CDC [NHANES Oral Microbiome](https://wwwn.cdc.gov/Nchs/Nhanes/omp/) genus-level oral-rinse data from the 2009-2010 survey cycle; 200 profiles from 4,960 participants with unique SEQN identifiers. |
| `JieZ_2017` | Stool relative-abundance profiles from curatedMetagenomicData and the [Jie atherosclerotic cardiovascular-disease study](https://pubmed.ncbi.nlm.nih.gov/29018189/); 200 profiles from 385 independent stool samples. |
| `Qiita11993Colombia` | Public [Qiita study 11993](https://qiita.ucsd.edu/public/?study_id=11993), examining gut microbiota and cardiometabolic health in a Colombian population undergoing Westernization; 200 fecal profiles from 441 participants. |
| `LiJ_2014` | Stool relative-abundance profiles from curatedMetagenomicData and the [Li integrated human gut microbial gene catalog](https://pubmed.ncbi.nlm.nih.gov/24997786/); 200 profiles from 249 participants. |
| `Atlas1006` | The CC0 [HITChip Atlas 1006 dataset](https://doi.org/10.5061/dryad.pk75d), a population resource of adult stool microbiome profiles; 200 time-zero profiles from 1,006 participants. |
| `SchirmerM_2016` | Stool relative-abundance profiles from curatedMetagenomicData and the [Schirmer study of the gut microbiome and inflammatory cytokine production](https://doi.org/10.1016/j.cell.2016.10.020); 200 profiles from 471 independent stool samples. |
| `ZeeviD_2015` | Stool relative-abundance profiles from curatedMetagenomicData and the [Zeevi personalized-nutrition study of postprandial glycemic responses](https://pubmed.ncbi.nlm.nih.gov/26590418/); 200 profiles from 900 participants. |
| `VilaAV_2018` | Stool relative-abundance profiles from curatedMetagenomicData and the [Vila study of gut microbiome composition and function in IBD and IBS](https://pubmed.ncbi.nlm.nih.gov/30567928/); 200 profiles from 355 participants. |
| `QinJ_2012` | Stool relative-abundance profiles from curatedMetagenomicData and the [Qin metagenome-wide association study of type 2 diabetes](https://pubmed.ncbi.nlm.nih.gov/23023125/); 200 profiles from 363 participants. |

## Analysis design

For each dataset, genera with prevalence of at least 0.10 and mean relative abundance of at least `1e-5` in the 200 stored compositions are ranked by decreasing mean abundance, with taxon name as the deterministic tie-break. The first 30 genera form the analysis panel. Community mass outside the panel is retained as `Other_unmodeled` during multinomial sampling.

Each of 100 randomizations assigns 100 samples to `H` and 100 to `Case`. The same group assignment is used for both library-depth settings:

| Setting | Median depth in H | Median depth in Case |
|---|---:|---:|
| Balanced | 4,000 | 4,000 |
| Fourfold | 1,500 | 6,000 |

Library sizes follow a log-normal distribution with `sdlog = 0.45` and are limited to 300-30,000 reads. Counts are sampled from the 30 target genera together with `Other_unmodeled`; the latter preserves the remaining community mass but is not tested. DASRA, ZINQ, and MaAsLin3 analyze the same 30 target taxa in every randomization.

DASRA uses its public Bonferroni omnibus result. An unavailable component contributes a conservative p-value of 1, and the omnibus is formed when at least one component is formed. MaAsLin3 requires both component fits and a combined p-value. An unavailable method result is represented by an analysis p-value of 1 in the common 30-taxon denominator. Type I error is the proportion of analysis p-values at or below 0.05.

Randomization and count-generation seeds are deterministic functions of the dataset, randomization, and depth setting. Each dataset analysis uses seven independent R workers.

## Software requirements

The analysis requires R, DASRA 0.6.0, ZINQ, maaslin3, and ggplot2.

## Running the analyses

From this directory, run one dataset through its entry point, for example:

```sh
Rscript SchirmerM_2016/analysis.R
```

An interrupted run can be restarted; completed randomization-setting units are reused. Temporary checkpoints and method work files are removed after all 200 units pass validation.

Each completed `results` directory contains:

- `taxon_pvalues.csv`: taxon-level p-values, availability, and status;
- `replicate_metrics.csv`: Type I error and availability by randomization, setting, and method;
- `diagnostics.csv`: generated positive-count and DASRA component diagnostics;
- `depth_diagnostics.csv`: realized group sizes and library-depth summaries;
- `run_manifest.csv`: design constants and software versions.

After all 20 dataset analyses finish, run:

```sh
Rscript combined/analysis.R
```

The combined script validates the complete 20-dataset x 100-randomization x 2-setting x 3-method x 30-taxon result grid. It writes randomization-level and dataset-level Type I error summaries, method-availability and completeness tables, formal fourfold support summaries, and `combined/negative_control_type1_error.pdf`. Figure intervals are 95% Monte Carlo intervals based on variation across the 100 randomization-level Type I error estimates.
