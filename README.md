# DASRA

**DASRA** provides depth-aware taxon-level inference for microbiome count data. It separates a group association into two complementary components:

- a **structural-absence component**, which tests whether the probability that a taxon is absent differs between groups; and
- a **relative-abundance component**, which tests the covariate-adjusted group difference in mean log relative abundance conditional on taxon presence.

Finite sequencing depth is incorporated through a latent-state binomial count
model. Within the fitted conditional factorization, conditioning on a positive
count removes the structural-presence probability from the abundance
likelihood. Every retained taxon is therefore analyzed with one zero-truncated
present-conditional model: zero counts carry no abundance score, and the
abundance fit does not require a structural-nuisance estimate or a
data-dependent model switch. The fixed observed-design abundance effect is
then centered against a target-excluded cross-taxon reference background to
account for compositional closure.

The two components answer complementary questions. The structural-absence
component concerns the probability of taxon absence, whereas the
relative-abundance component concerns the covariate-standardized mean log
relative abundance among samples in which the taxon is present. The reported
abundance effect is the target taxon's contrast relative to a target-excluded
shared background, rather than an absolute-abundance change. Reference
centering is intended for analyses in which the eligible taxa provide a stable,
clearly separated common-background cluster.

For a finite regular interior structural nuisance solution, the structural test
uses the sum of nuisance-orthogonalized per-sample estimating-function
contributions together with the empirical sandwich variance formed from those
same contributions. At an exact nuisance root, this numerator is algebraically
identical to the restricted target score. A returned finite nuisance fit is not
used for regular inference if the exact zero structural-absence-probability
limit has a strictly lower detection objective; this comparison establishes
that the returned fit is nonoptimal without classifying the global optimum as
a boundary solution.

## Installation

Install the dependencies and then install DASRA from the package source
directory:

```r
install.packages(c("Rcpp", "statmod"))
devtools::install()
```

The released GitHub version can be installed with
`remotes::install_github("yiqianomics/DASRA")`.

## Input

`counts` is a raw integer count matrix. By default, taxa are rows and samples
are columns. `metadata` must have sample identifiers as row names.
`library_size` must contain the original total number of sequencing reads for
each sample.

The formula is one-sided, contains the tested binary `group` variable as an
additive main effect, and may include adjustment terms. Group interactions,
offsets, and random-effect terms are not supported. Variables used by the
formula, the tested group, and a metadata-based library size must be complete.
DASRA does not silently remove samples; handle missing values first while
keeping counts, metadata, and any named depth vector aligned. Missing values in
unused metadata columns are allowed.

The main analysis controls are:

- `min_positive_samples = 3L`, the minimum number of samples with a positive
  count required by the preliminary retention rule;
- `min_reference_taxa = 4L`, the number of target-excluded taxa required to
  form an abundance reference;
- `structural_conditional_present_starts = "adaptive"`, with `"full"`
  available for an immediate five-start structural fit; the abundance arm
  always uses its full five-start bank;
- `structural_quadrature_points = 1001L`, the validated structural quadrature
  default;
- `abundance_quadrature_points = 41L`, the validated abundance quadrature
  default;
- `store_plot_data = FALSE`; set this to `TRUE` when the fitted object will be
  used to draw the dual-component association profile; and
- `workers = 1L` and `verbose = FALSE`. Increasing `workers` uses an ordered,
  cross-platform process cluster, while `verbose = TRUE` reports arm-level
  progress.

The retention and reference thresholds define the tested families and should
be selected before inspecting results. The two arms expose separate numerical
controls because their validated paths use different quadrature orders and
start strategies. With `full_output = TRUE`, DASRA records fixed-fit
higher-order comparisons and detailed numerical diagnostics without refitting
the model or changing primary inference. For abundance fits, the detailed
taxon table records whether the comparison was performed and succeeded, the
two quadrature orders, and the absolute conditional-log-likelihood and effect
discrepancies.

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

fit <- dasra(
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
one signed-evidence axis. Their shared sequential color scale encodes the
corresponding component adjusted p-values, while superscript stars beside each
feature encode the primary omnibus adjusted p-value. Group colors are kept
visually distinct from this evidence scale, and the reference and comparison
labels come directly from the fitted contrast by default.

By default, the plot displays omnibus-adjusted discoveries at 0.05, ordered by
their adjusted and raw omnibus p-values, with at most 24 rows. It does not fill
unused rows with nonsignificant features. Use `selection = "top"` for a fixed
top-results view, `selection = "all"` for every formed result, or `features`
to supply an explicit order.

```r
plot(fit)

plot(fit, selection = "top", max_features = 20L)

plot(
    fit,
    features = c("Taxon_2", "Taxon_1"),
    p_color_limits = c(1e-8, 1)
)

plot(fit, group_labels = c("Control", "Disease"))
```

The default adaptive color range is shared by both components. Its light end
is fixed at p = 1, and its dark-plum end is the smallest positive component
adjusted p-value among the displayed markers rounded down to a power of ten.
This uses the available color range without giving the two arms incomparable
scales or reusing the reference/comparison colors.
Supplying numeric `p_color_limits` fixes the range when feature subsets or
analyses must be compared directly.

Feature space is calculated from the longest displayed name, including its
significance superscript, rather than reserved at a fixed width. The side
profiles are covariate-standardized model summaries for the fitted groups; the
central symbols carry the prespecified structural score statistic and the
reference-corrected abundance Wald statistic. The abundance side panel shows
the fitted present-conditional geometric mean relative abundance, standardized
over the observed covariate distribution, as a percentage of total reads. It
uses a log scale and paired group points rather than bars, because a
logarithmic scale has no meaningful zero baseline. Side-panel segments are
descriptive fitted group summaries, not effect estimates or confidence
intervals; component inference is carried by the central signed statistics.
An omitted central marker or `--` means that component or side summary is
unavailable, not that its effect is zero. Optional `group_labels` can shorten
long fitted group names for display without changing the analysis or ordering.

The same call can save a manuscript-ready vector or raster figure. Without
explicit dimensions, PDF, PNG, and SVG output use a 180-mm manuscript width
and an adaptive height:

```r
plot(fit, file = "dasra-profile.pdf")
plot(fit, file = "dasra-profile.png", width = 7.2, height = 6.5, dpi = 600)
```

## Results

For the comparison-minus-reference contrast:

- positive `z_structural_absence` indicates greater structural absence in the comparison group;
- positive `estimate_relative_abundance` indicates that the taxon's present-conditional mean log-relative-abundance contrast exceeds the target-excluded compositional background;
- `p_omnibus` is the Bonferroni minimum-p combination of the two component tests; and
- `p_omnibus_cauchy` is an equal-weight Cauchy sensitivity combination.

Adjusted p-value columns use the `p_adj_` prefix. The adjustment method is
recorded in `fit$settings$p_adjust_method`.

Begin interpretation with `fit$diagnostics`. The columns
`n_positive_reference` and `n_positive_comparison` summarize the observed
positive-count support in the two groups. The `formed_*`, `reason_*`, and
`warning_*` columns show whether each requested result was produced and record
its numerical status.

Taxa with fewer positive counts than `min_positive_samples` are not retained.
A retained taxon whose requested component cannot be formed remains in that
testing family with an operational p-value of one. Use the formation indicator
and reason to identify this family-bookkeeping value; it is not evidence of no
association.

Two structural statuses are recorded as nonregular. A taxon with no observed
zeros has no variation in its absence indicator (`no_observed_zeros`). When the
intercept-only structural nuisance equation has $C_0 \leq 0$, its exact
solution lies at zero structural-absence probability
(`structural_absence_boundary_at_zero`). Their signed statistic is undefined,
and they are excluded from the Cauchy sensitivity combination. A finite
covariate-adjusted nuisance fit is not used when the exact zero-limit detection
objective is lower (`structural_absence_nonoptimal_nuisance_fit`). An otherwise
supported all-positive taxon can still enter the abundance arm because its
conditional-mark likelihood has no structural-nuisance parameter.

The abundance standard error is a centered, sample-aligned HC0 sandwich
estimate, and its two-sided p-value uses a first-order standard-normal Wald
reference. Component-use columns in `fit$results` distinguish the components
used by the primary and Cauchy omnibus analyses. Requested components return
lightweight warning codes even when `full_output = FALSE`; setting
`full_output = TRUE` additionally retains detailed fitted objects.

## Recommended abundance workflow

For a reportable abundance analysis:

1. State the estimand as a present-conditional, reference-centered relative
   contrast, and identify the comparison and reference groups from
   `fit$settings$contrast`.
2. Report group-specific positive-count support together with component
   formation, reason, and warning fields from `fit$diagnostics`.
3. With `full_output = TRUE`, inspect
   `fit$fits$relative_abundance$taxon`. This table separates the raw taxon
   estimate, target-excluded background estimate, and corrected estimate, and
   records reference size, bandwidth, curvature, convergence, and quadrature
   diagnostics.
4. When coordinated community-wide changes are scientifically plausible,
   present the cross-taxon raw-effect distribution and reference summaries
   alongside the corrected results. Use these quantities for transparent
   interpretation rather than as outcome-dependent filtering thresholds.

The raw and corrected quantities answer different reference questions; the
raw p-value is a diagnostic companion rather than a replacement for the
reported reference-centered test.
