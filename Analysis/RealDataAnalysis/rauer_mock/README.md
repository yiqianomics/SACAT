# Rauer defined mock-community benchmark

This analysis uses the mock-community sequencing experiment reported by Rauer
et al. in “De-biasing microbiome sequencing data: bacterial morphology-based
correction of extraction bias and correlates of chimera formation”
([doi:10.1186/s40168-024-01998-4](https://doi.org/10.1186/s40168-024-01998-4)).
Raw reads are available through ENA accession
[PRJEB67827](https://www.ebi.ac.uk/ena/browser/view/PRJEB67827) and
[OSF](https://osf.io/ykrbp/). The sample metadata and original study code are
available from the authors’
[GitHub repository](https://github.com/LuiseRauer/Extraction-bias-correction).

## Analysis design

- **Libraries:** 32 cell-based mock-community profiles, comprising 16 even
  D6300 libraries and 16 D6321 spike-in libraries. Each source community is
  represented by two input levels and all eight combinations of extraction
  kit, lysis condition, and extraction buffer.
- **Tested panel:** eight D6300 taxa with known structural differences and
  three D6321 taxa with known abundance differences.
- **Alternative design:** 52 covariate-balanced allocations. The two groups
  contain 12 D6300 and 4 D6321 libraries, and 4 D6300 and 12 D6321 libraries,
  respectively.
- **Complete-null design:** all 111 distinct covariate-balanced partitions with
  eight libraries from each source community in both groups.
- **Sequencing depths:** native depth and controlled downsampling to 1,000,
  750, 500, 250, 100, and 50 reads per library.
- **Adjustment variables:** extraction kit, lysis condition, extraction buffer,
  and input stratum. MaAsLin3 and ZINQ also include standardized log library
  size when library sizes vary.
- **Multiplicity:** BH adjustment at level 0.05 is applied separately within
  each 11-taxon method-component family. Unavailable tests remain in the
  family and are counted as non-discoveries.
- **Library size:** the total number of DADA2-denoised forward-read ASV counts
  in a library before reference-taxon selection.

Downsampling provides a controlled sequencing-depth perturbation. At each
fixed depth, every library has the stated total count.

## Source files

Place the 94 forward-read FASTQ files from the OSF archive in `raw/fastq/`.
All 94 files are used to learn the DADA2 error model. The benchmark uses the 32
libraries described above. Copy `Input/Metafile.csv` from the authors' GitHub
repository to `raw/Metafile.csv`.

Download the D6300 references from the
[ZymoBIOMICS reference archive](https://doi.org/10.5281/zenodo.3935737) and the
D6321 references from the
[official Zymo Research archive](https://s3.amazonaws.com/zymo-files/BioPool/D6321.refseq.zip).
Place their 16S FASTA files at:

```text
raw/reference/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/
raw/reference/D6321.refseq/16S/
```

The benchmark uses the eight bacterial D6300 references and the three D6321
references.

## Reproducing the analysis

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/rauer_mock/prepare_data.R
Rscript Analysis/RealDataAnalysis/rauer_mock/analysis.R
Rscript Analysis/RealDataAnalysis/rauer_mock/make_figures.R
```

Set `DASRA_WORKERS` to use more than one local process during preparation or
analysis. The analysis can also be resumed one design at a time with
`--design=nonnull` or `--design=null`.

Prepared inputs and allocation-level checkpoints are written to `processed/`
and `work/`.

## Outputs

- `table/rauer_mock_nonnull_results.csv` and
  `table/rauer_mock_null_results.csv` contain the complete taxon-level results.
- `table/rauer_mock_native_taxon_summary.csv` and
  `table/rauer_mock_depth_summary.csv` contain the reported numerical
  summaries.
- `figs/rauer_mock_main.pdf` is the benchmark figure.
