# Rauer spike-in benchmark

This benchmark uses the defined spike-in communities of
[Rauer et al.](https://doi.org/10.1186/s40168-024-01998-4) to assess structural-
absence (SA) and conditional relative-abundance (RA) signal recovery as
sequencing depth decreases. Balanced null partitions assess false positives.

## Analysis design

- **Libraries:** 16 D6300-plus-spike-in and 16 spike-in-only profiles, covering
  two input levels and eight extraction-kit, lysis, and buffer combinations.
- **Tested panel:** eight D6300 taxa with known structural differences and
  three D6321 taxa with known abundance differences.
- **Allocation rule:** pair the two source communities within each of the 16
  experimental cells, placing one library from each pair in each group.
- **Alternative design:** 52 covariate-balanced allocations. The two groups
  contain 12 D6300 and 4 D6321 libraries, and 4 D6300 and 12 D6321 libraries,
  respectively.
- **Complete-null design:** 111 covariate-balanced partitions under this
  paired-cell rule, with eight libraries from each source community in both
  groups; complementary group assignments are counted once.
- **Sequencing depths:** native depth and controlled downsampling to 1,000,
  750, 500, 250, 100, and 50 reads per library.
- **Adjustment variables:** extraction kit, lysis condition, extraction buffer,
  and input stratum. MaAsLin3 and ZINQ also include standardized log library
  size when library sizes vary.
- **Multiplicity:** Benjamini–Hochberg (BH) adjustment at level 0.05 is applied
  separately within each 11-taxon method-component family.
- **Library size:** the total number of DADA2-denoised forward-read ASV counts
  in a library before reference-taxon selection.

## Source files

Place the 94 forward-read FASTQ files from [OSF](https://osf.io/ykrbp/) in
`raw/fastq/`; all 94 files are processed with DADA2. Reads are also
deposited under [PRJEB67827](https://www.ebi.ac.uk/ena/browser/view/PRJEB67827).
Copy `Input/Metafile.csv` from the authors'
[repository](https://github.com/LuiseRauer/Extraction-bias-correction)
to `raw/Metafile.csv`.

Download the D6300 references from the
[ZymoBIOMICS reference archive](https://doi.org/10.5281/zenodo.3935737) and the
D6321 references from the
[official Zymo Research archive](https://s3.amazonaws.com/zymo-files/BioPool/D6321.refseq.zip).
Place their 16S FASTA files at:

```text
raw/reference/ZymoBIOMICS.STD.refseq.v2/ssrRNAs/
raw/reference/D6321.refseq/16S/
```

## Code and outputs

The analysis requires SACAT 0.7.1 or later. Versions used for the saved results
are listed in [`software_versions.csv`](software_versions.csv).

Run from the repository root:

```sh
Rscript Analysis/RealDataAnalysis/rauer_mock/prepare_data.R
Rscript Analysis/RealDataAnalysis/rauer_mock/analysis.R
Rscript Analysis/RealDataAnalysis/rauer_mock/make_figures.R
```

Set `SACAT_WORKERS` to use more than one local process during preparation or
analysis. To run one design, use
`--design=nonnull` or `--design=null`.

Prepared inputs and intermediate files are written to `processed/` and `work/`.

- `table/rauer_mock_nonnull_results.csv` and
  `table/rauer_mock_null_results.csv` contain the complete taxon-level results.
- `table/rauer_mock_native_taxon_summary.csv` and
  `table/rauer_mock_depth_summary.csv` contain the reported numerical
  summaries.
- `figs/rauer_mock_main.pdf` is the benchmark figure.

The combined recovery and CDI figure is in
[`Analysis/MainFigures`](../../MainFigures/README.md).
