# Real-data analyses

This directory contains the shared analysis implementation and the dataset-specific preparation code for two public microbiome studies:

- `crc_baxter/`: colorectal cancer versus healthy controls.
- `cdi_schubert/`: *Clostridioides difficile* infection versus healthy controls.

`analysis.R` is the common entry point for DASRA and all comparison methods. Dataset-specific models, reference groups, and input names are declared near the beginning of that script. Each dataset directory contains its own `prepare_data.R` and separates source data, prepared inputs, method workspaces, tables, and figures.

## Source data and preparation

Both source archives come from MicrobiomeHD version 3, [Zenodo record 1146764](https://doi.org/10.5281/zenodo.1146764):

- `crc_baxter/raw/crc_baxter_results.tar.gz`
- `cdi_schubert/raw/cdi_schubert_results.tar.gz`

From the repository root, extract the archives and prepare the analysis inputs with:

```sh
tar -xzf Analysis/RealDataAnalysis/crc_baxter/raw/crc_baxter_results.tar.gz \
  -C Analysis/RealDataAnalysis/crc_baxter/raw
tar -xzf Analysis/RealDataAnalysis/cdi_schubert/raw/cdi_schubert_results.tar.gz \
  -C Analysis/RealDataAnalysis/cdi_schubert/raw

Rscript Analysis/RealDataAnalysis/crc_baxter/prepare_data.R
Rscript Analysis/RealDataAnalysis/cdi_schubert/prepare_data.R
```

The preparation scripts document the cohort definitions, sample exclusions, taxonomic aggregation, filtering rules, and use of original library sizes. They write the validated analysis objects and supporting audit tables to each dataset's `processed/` directory.

## Analysis

Run both datasets from the repository root with:

```sh
Rscript Analysis/RealDataAnalysis/analysis.R
```

To run one dataset independently, use:

```sh
Rscript Analysis/RealDataAnalysis/analysis.R --dataset=crc_baxter
Rscript Analysis/RealDataAnalysis/analysis.R --dataset=cdi_schubert
```

The random seed is fixed within `analysis.R`. Package versions, method status, elapsed time, and method-specific conventions are recorded in the exported status tables.

## Outputs

Each dataset produces six reviewer-facing CSV files in `table/`:

1. `*_analysis_input_summary.csv`: sample counts, taxon counts, retained-read summaries, and multiplicity definition.
2. `*_method_results_all_taxa.csv`: taxon-level results for every method and component.
3. `*_method_status.csv`: package versions, completion status, discoveries, elapsed time, and analysis conventions.
4. `*_discovery_summary.csv`: discovery counts by method.
5. `*_upset_intersections.csv`: intersection definitions and sizes used in the comparison figure.
6. `*_upset_intersection_membership.csv`: taxon membership in the plotted intersections.

The corresponding UpSet figure is written to `figs/` as `*_upset.pdf`.

`processed/` contains reproducible inputs derived from the public archives. `work/` contains regenerable method-specific intermediate files. `R_lib/` is an optional machine-local compatibility library and is not part of the analysis output.
