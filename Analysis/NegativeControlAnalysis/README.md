# Native-count negative-control analysis

This analysis evaluates SACAT's null calibration using 20 public microbiome
datasets. Each contributes 200 profiles and 30 genera. Native counts and
original library sizes remain fixed while labels are randomly assigned to
two groups of 100, repeated 100 times.

Genera are selected before randomization: prevalence at least 10%, mean
relative abundance at least `1e-5`, and the highest 30 mean abundances.
Library sizes come from the complete source tables; GEMS supplies `totalCounts`.
Each fit uses `~ group`, with group `A` as reference, and evaluates structural
absence (SA), conditional relative abundance (RA), and the Bonferroni omnibus.

## Type I error

Type I error is the number of genera with raw p < 0.05 divided by 30, averaged
over randomizations. Overall estimates give each dataset equal weight.

## Datasets

| Directory | Public study and selected profiles |
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

## Running the analysis

Fitting requires R and SACAT 0.7.1. Summaries and figures also use
`data.table`, `ggplot2`, `scales`, and Cairo graphics support.

From this directory, run one dataset or all 20 datasets with seven workers:

```sh
Rscript BaxterE_2016/analysis.R 7
Rscript run_negative_control.R all 7
```

Recreate the combined summaries and figure from the supplied result tables with:

```sh
Rscript combined/analysis.R
```

Each dataset's `results/` contains `taxon_pvalues.csv` and
`replicate_metrics.csv`. The `combined/` directory contains dataset-level and
overall Type I error, availability and design summaries, and
`negative_control_type1_error.pdf`.
