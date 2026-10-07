# Native-count negative-control analysis

This directory contains SACAT negative-control analyses of 20 public human microbiome datasets. Each analysis uses 200 participants, a fixed panel of 30 genera, and 100 random assignments to two groups of 100. The published counts and sample library sizes remain fixed throughout the analysis.

## Data and analysis design

The supplied `analysis_input.rds` files contain the genus counts, sample identifiers, library sizes, selected taxa, and source information needed to run each analysis. Library size is taken from the full published count table before genus and panel selection; GEMS uses the published `totalCounts` vector. The `source` field records the public data location and dataset-specific preparation, and `metadata` identifies the selected profiles.

Before group labels are assigned, genera with prevalence of at least 10% and mean relative abundance of at least `1e-5` are ranked by decreasing mean relative abundance, with genus name breaking ties. The first 30 genera form the fixed panel. Relative abundance is computed using the full library size.

For each randomization, SACAT 0.7.1 fits `~ group` to groups `A` and `B`, with `A` as the reference. The analysis reports the structural-absence (SA) component, relative-abundance (RA) component, and Bonferroni omnibus test. Group assignments use the seed `202609190 + dataset_index * 100000 + replicate * 1009`; the dataset index is stored in each input.

Global-null Type I error is the number of raw p-values below 0.05 divided by 30 in each randomization, averaged over the 100 randomizations. An unavailable test contributes a p-value of 1, keeping the denominator fixed. BH-adjusted p-values are also stored in the taxon-level tables. Overall means give each dataset equal weight.

## Datasets

Each input contains 200 selected profiles. Repeated visits are reduced to one profile per participant; the Goodrich and iMSMS selections also account for family and household membership.

| Directory | Study and selected profiles |
|---|---|
| `BaxterE_2016` | [Baxter et al. (2016)](https://doi.org/10.1186/s13073-016-0290-3): stool 16S profiles from the 292-participant primary colorectal-cancer cohort. |
| `cdi_schubert` | [Schubert et al. (2014)](https://doi.org/10.1128/mBio.01021-14): stool 16S profiles from 243 participants with CDI or healthy-control status and complete analysis covariates. |
| `gems_pediatric_diarrhea` | [Pop et al. (2014)](https://doi.org/10.1186/gb-2014-15-6-r76): stool 16S profiles from 992 children in GEMS. |
| `korean_hypertension` | [Song et al. (2023)](https://doi.org/10.3389/frmbi.2023.1072059): stool 16S profiles from 623 participants. |
| `qiita_1939_pediatric_crohn` | [Gevers et al. (2014)](https://doi.org/10.1016/j.chom.2014.02.005): rectal mucosal profiles from 244 children with Crohn's disease or no IBD, complete analysis covariates, and no recorded antibiotic exposure during collection. |
| `NogueraJulianHIV` | [Noguera-Julian et al. (2016)](https://doi.org/10.1016/j.ebiom.2016.01.032): baseline stool profiles from the Barcelona and Stockholm HIV cohorts. |
| `ArtPrize_2015_Forehead` | [Dimitriu et al. (2019)](https://doi.org/10.1128/mBio.00839-19): forehead-skin profiles from 372 participants. |
| `HMPV35Throat` | [HMP16SData](https://doi.org/10.1093/aje/kwz006): V3-5 throat profiles, using the earliest visit from each of 202 participants. |
| `Qiita11993Colombia` | [de la Cuesta-Zuluaga et al. (2018)](https://doi.org/10.1038/s41598-018-29687-x): fecal profiles from 441 Colombian participants. |
| `Qiita13631ASD` | [Dan et al. (2020)](https://doi.org/10.1080/19490976.2020.1747329): fecal profiles from 286 participants in an autism-spectrum-disorder study. |
| `NHANESOral_2009_2010` | [Vogtmann et al. (2023)](https://doi.org/10.1016/S2666-5247(22)00333-0): oral-rinse profiles from the 2009-2010 survey cycle, identified by unique SEQN values. |
| `NHANESOral_2011_2012` | [Vogtmann et al. (2023)](https://doi.org/10.1016/S2666-5247(22)00333-0): oral-rinse profiles from the 2011-2012 survey cycle, identified by unique SEQN values. |
| `Goodrich_2014` | [Goodrich et al. (2014)](https://doi.org/10.1016/j.cell.2014.09.053): stool profiles selected from 532 twin families, with one profile per family. |
| `Zeevi_2015_16S` | [Zeevi et al. (2015)](https://doi.org/10.1016/j.cell.2015.11.001): stool 16S profiles from 729 participants in the public source metadata. |
| `Singh_2015_EDD` | [Singh et al. (2015)](https://doi.org/10.1186/s40168-015-0109-2): baseline stool profiles selected from 201 patient identifiers in the public source metadata. |
| `AmericanGut_2018` | [McDonald et al. (2018)](https://doi.org/10.1128/mSystems.00031-18): human stool profiles selected from 9,511 participants in the deposited Deblur table. |
| `HillBurns_2017` | [Hill-Burns et al. (2017)](https://doi.org/10.1002/mds.26942): stool profiles from 332 participants in a Parkinson's disease study. |
| `Heinzel_2021_TREND` | [Heinzel et al. (2021)](https://doi.org/10.1002/ana.26128): stool profiles from 663 TREND participants. |
| `Wallen_2020_Birmingham` | [Wallen et al. (2020)](https://doi.org/10.1038/s41531-020-0112-6): stool profiles from 507 participants in the Birmingham cohort. |
| `iMSMS_2022` | [International Multiple Sclerosis Microbiome Study (2022)](https://doi.org/10.1016/j.cell.2022.08.021): healthy controls selected from 575 households, with one participant per household and the authors' participant-level pooled counts. |

The inputs include the authors' published count-table processing. Details such as taxonomic assignments, source objects, and upstream feature filtering are recorded in each input's `source` field.

## Running the analyses

Model fitting requires R and SACAT 0.7.1. Summaries and the figure require `data.table`, `ggplot2`, `scales`, and an R build with Cairo graphics.

From this directory, run one dataset with seven workers:

```sh
Rscript BaxterE_2016/analysis.R 7
```

To run all datasets sequentially, using seven workers for the randomizations within each dataset:

```sh
Rscript run_negative_control.R all 7
```

Completed fits are saved in each dataset's `results/replicates/` directory and reused when a run is resumed. These local fit objects are excluded from Git; the complete result tables are included. Windows runs use one worker.

To reproduce the combined summaries and PDF from the supplied result tables:

```sh
Rscript combined/analysis.R
```

## Files and outputs

Every dataset directory contains `analysis_input.rds`, an `analysis.R` entry point, and five result tables:

- `taxon_pvalues.csv`: raw and BH-adjusted p-values, test availability, and raw-p rejection decisions for every genus, component, and randomization.
- `replicate_metrics.csv`: rejection counts, Type I error, and availability for each randomization.
- `diagnostics.csv`: SACAT's taxon-level fit diagnostics.
- `depth_diagnostics.csv`: realized group sizes, library-depth summaries, and randomization seeds.
- `analysis_settings.csv`: software version and fitted analysis settings.

The shared fitting and export functions are in `run_negative_control.R`. The combined script writes `dataset_type1_error.csv`, `replicate_type1_error.csv`, `overall_type1_error.csv`, `availability_summary.csv`, `dataset_design.csv`, and `negative_control_type1_error.pdf` to `combined/`.

The dataset-level table includes the Monte Carlo standard deviation across randomizations and its standard error, computed as that standard deviation divided by the square root of 100. Figure intervals are the mean plus or minus 1.96 Monte Carlo standard errors, limited to the interval from 0 to 1.
