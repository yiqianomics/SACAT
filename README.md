# DASRA

**DASRA** provides depth-aware taxon-level inference for microbiome count data. It separates a group association into two complementary components:

- a **structural-absence component**, which tests whether the probability that a taxon is absent differs between groups; and
- a **relative-abundance component**, which tests the covariate-adjusted group difference in mean log relative abundance conditional on taxon presence.

Finite sequencing depth is incorporated through a latent-state binomial count model. A zero count can therefore contribute to the relative-abundance component when the fitted model supports presence followed by nondetection. The abundance effect is subsequently centered against a target-excluded robust cross-taxon reference background to account for compositional closure.
The correction is designed for analyses in which a strict majority of eligible taxa share a common compositional background, and the reported abundance effect is interpreted as a reference-centered relative contrast.

## Installation

Install the dependencies and then install DASRA from the package source
directory:

```r
install.packages(c("Rcpp", "nleqslv", "statmod"))
devtools::install()
```

The released GitHub version can be installed with
`remotes::install_github("yiqianomics/DASRA")`.

## Input

`counts` is a raw integer count matrix. By default, taxa are rows and samples are columns. `metadata` must have sample identifiers as row names. `library_size` must contain the original total number of sequencing reads for each sample.

The formula is one-sided, contains the tested binary `group` variable as an additive main effect, and may include adjustment terms. Group interactions, offsets, and random-effect terms are not supported.

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

Taxa with fewer than three positive counts are not retained. A requested component that cannot be formed for a retained taxon enters its multiplicity family with p-value one, and the reason is recorded in `fit$diagnostics`.

Two structural outcomes are nonregular. A taxon with no observed zeros has no
variation in its absence indicator (`no_observed_zeros`). When the
intercept-only structural nuisance equation has $C_0 \leq 0$, its exact
solution lies at zero structural-absence probability
(`structural_absence_boundary_at_zero`). Both outcomes return a conservative
p-value of one, remain in the primary Bonferroni family, and are excluded from
the Cauchy sensitivity combination. Their signed score statistic is undefined.

For an otherwise supported all-positive taxon, the full relative-abundance
estimating system has a structural-nuisance boundary. The relative-abundance
component is unavailable with reason
`no_observed_zeros_structural_nuisance_boundary`. Component-use columns in
`fit$results` distinguish the formed components used by the primary omnibus
from the regular components used by the Cauchy sensitivity omnibus.
Requested components return lightweight warning codes even when
`full_output = FALSE`. Structural analyses include
`fit$diagnostics$warning_structural_absence` and
`fit$diagnostics$nonregular_structural_absence`; abundance analyses include
`fit$diagnostics$warning_relative_abundance`. Setting `full_output = TRUE`
additionally retains detailed fitted objects.
