# SACAT

[![R-CMD-check](https://github.com/yiqianomics/SACAT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yiqianomics/SACAT/actions/workflows/R-CMD-check.yaml)

**Structural-Absence and Conditional-Abundance Testing** provides depth-aware taxon-level inference for microbiome count data. It separates a group association into two complementary components:

- a **structural-absence component**, which tests whether the probability that a taxon is absent differs between groups; and
- a **relative-abundance component**, which tests the covariate-adjusted group difference in mean log relative abundance conditional on taxon presence.

SACAT models microbial counts together with their original sequencing depths,
allowing zeros observed at different depths to provide different evidence about
structural absence. For samples in which a taxon is present, the abundance
component estimates the adjusted group difference in mean log relative
abundance and centers that difference against a target-excluded background of
taxa. The reported abundance effect is therefore a present-conditional relative
contrast.

## Installation

Install the runtime dependencies, then install SACAT from GitHub:

```r
install.packages(c("Rcpp", "statmod", "remotes"))
remotes::install_github("yiqianomics/SACAT")
```

From a local source checkout, run `R CMD INSTALL .` at the package root.

## Input

`counts` is a raw integer count matrix. By default, taxa are rows and samples
are columns. `metadata` must have sample identifiers as row names.
`library_size` must contain the original total number of sequencing reads for
each sample.

The formula is one-sided, contains the tested binary `group` variable as an
additive main effect, and may include additive fixed-effect adjustment terms.
Variables used by the formula, the tested group, and a metadata-based library
size must be complete. Prepare a complete analysis set while keeping counts,
metadata, and any named depth vector aligned. Missing values in unused metadata
columns are allowed.

The main analysis controls are:

- `min_positive_samples = 3L`, the minimum number of samples with a positive
  count required by the preliminary retention rule;
- `min_reference_taxa = 4L`, the number of target-excluded taxa required to
  form an abundance reference;
- `structural_conditional_present_starts = "adaptive"`, with `"full"`
  available for an immediate five-start structural fit; the abundance arm
  always uses its full five-start bank;
- `structural_quadrature_points = 1001L`, the structural quadrature order;
- `abundance_quadrature_points = 41L`, the abundance quadrature order;
- `store_plot_data = FALSE`; set this to `TRUE` when the fitted object will be
  used to draw the dual-component association profile; and
- `workers = 1L` and `verbose = FALSE`. Increasing `workers` uses an ordered,
  cross-platform process cluster, while `verbose = TRUE` reports arm-level
  progress.

The retention and reference thresholds define the tested families and should
be selected before inspecting results. Set `full_output = TRUE` to retain fitted
objects and detailed numerical diagnostics.

## Example

```r
set.seed(2026)
n <- 80
taxa <- paste0("Taxon_", seq_len(8))
samples <- paste0("Sample_", seq_len(n))
group <- rep(c(0, 1), each = n / 2)
library_size <- sample(seq(8000L, 12000L, by = 500L), n, replace = TRUE)

baseline <- seq(-6.0, -5.2, length.out = length(taxa))
abundance_shift <- c(0.25, -0.20, rep(0, length(taxa) - 2L))
absence_reference <- c(0.15, 0.20, 0.18, 0.22, 0.16, 0.24, 0.19, 0.21)
absence_comparison <- c(0.30, 0.20, 0.10, 0.22, 0.16, 0.24, 0.19, 0.21)

probability <- matrix(0, nrow = n, ncol = length(taxa))
for (j in seq_along(taxa)) {
    absent_probability <- ifelse(
        group == 0, absence_reference[j], absence_comparison[j]
    )
    present <- runif(n) > absent_probability
    latent_abundance <- baseline[j] + abundance_shift[j] * group +
        rnorm(n, sd = 0.35)
    probability[, j] <- present * plogis(latent_abundance)
}

count_by_sample <- t(vapply(seq_len(n), function(i) {
    draw <- rmultinom(
        1,
        size = library_size[i],
        prob = c(probability[i, ], 1 - sum(probability[i, ]))
    )
    draw[seq_along(taxa), 1]
}, numeric(length(taxa))))
counts <- t(count_by_sample)
rownames(counts) <- taxa
colnames(counts) <- samples

metadata <- data.frame(
    group = factor(group, levels = c(0, 1), labels = c("control", "case")),
    reads = library_size,
    row.names = samples
)

fit <- sacat(
    counts = counts,
    metadata = metadata,
    formula = ~ group,
    group = "group",
    library_size = "reads",
    component = "all",
    store_plot_data = TRUE
)

fit
head(fit$results)
```

## Visualization

`plot(fit)` draws a taxon-aligned summary of the two complementary components.
The structural-absence triangle and present-conditional abundance circle share
one component Z-statistic axis. Their shared sequential color scale encodes the
corresponding component adjusted p-values, while superscript stars beside each
feature encode the primary omnibus adjusted p-value. Group colors remain
visually distinct from this evidence scale, and the reference and comparison
labels come directly from the fitted contrast by default.

For the default BH adjustment, short gates mark the family-wide discovery
boundary in each component lane at the selected `alpha`. The upper triangle
lane represents structural absence, and the lower circle lane represents
present-conditional abundance. Each boundary uses the complete fitted
component family and is invariant to the displayed feature subset. A component
gate appears when that component has at least one BH discovery. A warning
identifies any component whose boundary cannot be displayed. Set
`show_component_guides = FALSE` to hide the gates.

By default, the plot displays up to 24 omnibus-adjusted discoveries at 0.05,
ordered by their adjusted and raw omnibus p-values. Use `selection = "top"` for
a fixed top-results view, `selection = "all"` for every formed result, or
`features` to supply an explicit order.

```r
plot(fit)

plot(fit, selection = "top", max_features = 20L)

plot(
    fit,
    features = c("Taxon_2", "Taxon_1"),
    p_color_limits = c(1e-8, 1)
)

plot(fit, group_labels = c("Control", "Disease"))

plot(fit, show_component_guides = FALSE)
```

The default adaptive color range is shared by both components. Its light end
is fixed at p = 1, and its dark-plum end is the smallest positive component
adjusted p-value among the displayed markers rounded down to a power of ten.
Both components therefore use a directly comparable scale that remains
distinct from the reference/comparison colors.
Supplying numeric `p_color_limits` fixes the range when feature subsets or
analyses must be compared directly.

Feature width adapts to the longest displayed name and its significance
superscript. The side profiles are covariate-standardized model summaries for
the fitted groups; the central symbols carry the fitted structural score
statistic and the signed square root of the adjusted abundance likelihood-ratio
statistic. The abundance side panel shows the fitted present-conditional
geometric mean relative
abundance, standardized over the observed covariate distribution, as a
percentage of total reads. A log scale with paired group points represents
these positive values. Side-panel segments provide descriptive fitted group
summaries, and the central Z-statistics carry component inference. An omitted
central marker or `--` identifies an unavailable component or side summary.
Optional `group_labels` shortens long fitted group names while preserving the
analysis and group ordering.

The same call can save a vector or raster figure. Default PDF,
PNG, and SVG dimensions use a 180-mm manuscript width and an adaptive height:

```r
plot(fit, file = "sacat-profile.pdf")
plot(fit, file = "sacat-profile.png", width = 7.2, height = 6.5, dpi = 600)
```

SVG output requires an R build with Cairo support.

## Results

For the comparison-minus-reference contrast:

- positive `z_structural_absence` indicates greater structural absence in the comparison group;
- positive `estimate_relative_abundance` indicates that the taxon's present-conditional mean log-relative-abundance contrast exceeds the target-excluded compositional background;
- `p_omnibus` is the Bonferroni minimum-p combination of the two component tests; and
- `p_omnibus_cauchy` is an equal-weight Cauchy sensitivity combination.

Adjusted p-value columns use the `p_adj_` prefix. The adjustment method is
recorded in `fit$settings$p_adjust_method`.

Use `fit$diagnostics` to review positive-count support in each group and
component formation. Retained taxa remain in the corresponding testing family.
An unavailable component receives an operational p-value of one, with its cause
recorded in the corresponding `reason_*` field. Eligible abundance results are
reported independently of structural component formation. When the
present-conditional abundance test forms for fewer than 80% of retained taxa,
users should assess null calibration before interpreting its discoveries.

## Interpreting abundance results

1. State the estimand as a present-conditional, reference-centered relative
   contrast, and identify the comparison and reference groups from
   `fit$settings$contrast`.
2. Review group-specific positive-count support and component formation in
   `fit$diagnostics`.
3. With `full_output = TRUE`, inspect
   `fit$fits$relative_abundance$taxon`. This table separates the raw taxon
   estimate, target-excluded background estimate, and corrected estimate, and
   records the reference size and model diagnostics.

SACAT reports the reference-centered test as its abundance result. Detailed
output also records the raw taxon estimate and the target-excluded background.
The abundance test compares the fitted likelihood with a constrained fit whose
raw contrast equals that background. Its likelihood-ratio statistic is scaled
by the raw contrast's inverse-information variance divided by the corrected
sandwich variance, then compared with a chi-squared distribution with one degree
of freedom. `z_relative_abundance` reports the signed square root of this
statistic. With `full_output = TRUE`, the abundance output includes the profile
statistics and constrained fits.
