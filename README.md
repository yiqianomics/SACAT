# DASRA

**DASRA** provides depth-aware taxon-level inference for microbiome count data. It separates a group association into two complementary components:

- a **structural-absence component**, which tests whether the probability that a taxon is absent differs between groups; and
- a **relative-abundance component**, which tests the covariate-adjusted group difference in mean log relative abundance conditional on taxon presence.

Finite sequencing depth is incorporated through a latent-state binomial count
model. The abundance arm conditions on a positive count and fits the resulting
zero-truncated mark likelihood to every taxon. Under the model's separable
structural-presence gate, the structural-presence probability cancels from this
likelihood exactly, so zero counts carry zero abundance score and no
structural-nuisance estimate or data-dependent switch is needed. The fixed
observed-design abundance effect is then centered against a target-excluded
cross-taxon reference background to account for compositional closure.

The reference correction has a pointwise fixed-taxon justification when each
target-excluded reference set contains a separated strict majority sharing one
common background and the selected kernel mode is stable and isolated. The
reported effect is therefore a reference-centered, present-conditional relative
contrast. A count majority by itself is not a finite-sample guarantee: many
same-direction changes can overlap and shift the reference mode. This is the
main practical limitation of the abundance arm.

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
- `conditional_present_starts = "adaptive"`, with `"full"` available for an
  immediate five-start structural fit;
- `structural_quadrature_points = 1001L`, the validated structural quadrature
  default; and
- `workers = 1L` and `verbose = FALSE`. Increasing `workers` uses an ordered,
  cross-platform process cluster, while `verbose = TRUE` reports arm-level
  progress.

The retention and reference thresholds change the tested families and should
be selected before inspecting results. Lower quadrature orders trade numerical
accuracy for speed; with `full_output = TRUE`, DASRA reports comparison with a
strictly higher-order rule as a sensitivity diagnostic rather than an error
bound.

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
    component = "all"
)

fit
head(fit$results)
```

## Results

For the comparison-minus-reference contrast:

- positive `z_structural_absence` indicates greater structural absence in the comparison group;
- positive `estimate_relative_abundance` indicates that the taxon's present-conditional mean log-relative-abundance contrast exceeds the target-excluded compositional background;
- `p_omnibus` is the Bonferroni minimum-p combination of the two component tests; and
- `p_omnibus_cauchy` is an equal-weight Cauchy sensitivity combination.

Adjusted p-value columns use the `p_adj_` prefix. The adjustment method is
recorded in `fit$settings$p_adjust_method`.

Taxa with fewer positive counts than `min_positive_samples` are not retained. A
requested component that cannot be formed for a retained taxon enters its
multiplicity family with p-value one, and the reason is recorded in
`fit$diagnostics`.

Two structural outcomes are nonregular. A taxon with no observed zeros has no
variation in its absence indicator (`no_observed_zeros`). When the
intercept-only structural nuisance equation has $C_0 \leq 0$, its exact
solution lies at zero structural-absence probability
(`structural_absence_boundary_at_zero`). Both outcomes use one as an
operational conservative value in the primary Bonferroni family and are
excluded from the Cauchy sensitivity combination; this value is not a
calibrated boundary p-value. Their signed statistic is undefined. A returned
finite covariate-adjusted structural nuisance fit is unavailable with reason
`structural_absence_nonoptimal_nuisance_fit` if the exact zero
structural-absence-probability limit has a strictly lower detection objective.

An otherwise supported all-positive taxon can enter the abundance arm because
the conditional-mark likelihood has no structural-nuisance parameter. The
abundance standard error is a centered, sample-aligned HC0 sandwich estimate,
and its p-value uses a first-order normal Wald approximation rather than a
claimed exact small-sample pivot. Component-use columns in `fit$results`
distinguish the formed components used by the primary omnibus from the regular
components used by the Cauchy sensitivity omnibus.
Requested components return lightweight warning codes even when
`full_output = FALSE`. Structural analyses include
`fit$diagnostics$warning_structural_absence` and
`fit$diagnostics$nonregular_structural_absence`; abundance analyses include
`fit$diagnostics$warning_relative_abundance`. Setting `full_output = TRUE`
additionally retains detailed fitted objects.
