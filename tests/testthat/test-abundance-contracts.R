.make_abundance_correction_fits <- function() {
    phi <- cbind(
        c(-0.30, -0.12, -0.03, 0.07, 0.14, 0.24),
        c(-0.12, -0.04, -0.01, 0.03, 0.05, 0.09),
        c(0.08, -0.10, 0.02, 0.05, -0.07, 0.02),
        c(-0.05, 0.06, -0.02, 0.03, 0.01, -0.03),
        c(0.03, -0.02, 0.04, -0.06, 0.02, -0.01)
    )
    raw_delta <- c(0.18, 0.10, 0.11, 0.09, 0.105)
    raw_se <- sqrt(colSums(phi^2))

    lapply(seq_along(raw_delta), function(index) {
        fit <- DASRA:::.dasra_abundance_empty_fit("ok", nrow(phi))
        fit$available <- TRUE
        fit$status <- "ok"
        fit$raw_delta <- raw_delta[index]
        fit$raw_se <- raw_se[index]
        fit$raw_p <- 2 * pnorm(-abs(raw_delta[index] / raw_se[index]))
        fit$phi <- phi[, index]
        fit
    })
}

.make_abundance_reference_fits <- function(raw_delta, raw_se, n = 120L) {
    p <- length(raw_delta)
    phi <- matrix(0, nrow = n, ncol = p)
    scale <- raw_se / sqrt(2)
    phi[cbind(seq_len(p), seq_len(p))] <- scale
    phi[p + 1L, ] <- -scale

    lapply(seq_len(p), function(index) {
        fit <- DASRA:::.dasra_abundance_empty_fit("ok", n)
        fit$available <- TRUE
        fit$status <- "ok"
        fit$raw_delta <- raw_delta[index]
        fit$raw_se <- raw_se[index]
        fit$raw_p <- 2 * pnorm(-abs(raw_delta[index] / raw_se[index]))
        fit$phi <- phi[, index]
        fit
    })
}

test_that("all-positive abundance taxa stop at the nuisance boundary", {
    y <- c(3, 5, 4, 6, 7, 4, 8, 5)
    group <- rep(c(0, 1), each = 4L)
    depth <- rep(100L, length(y))
    gh <- DASRA:::make_count_gh_rule(3L)
    solver_entered <- FALSE

    fit <- with_mocked_bindings(
        DASRA:::.dasra_abundance_fit_taxon(
            y = y,
            N = depth,
            group = group,
            z = NULL,
            gh_fit = gh,
            gh_effect = gh,
            control = DASRA:::.dasra_abundance_control()
        ),
        .dasra_abundance_initial_mark = function(...) {
            solver_entered <<- TRUE
            stop("The solver should not be entered.")
        },
        .package = "DASRA"
    )

    expect_false(solver_entered)
    expect_false(fit$available)
    expect_identical(
        fit$status,
        "no_observed_zeros_structural_nuisance_boundary"
    )
    expect_equal(fit$phi, rep(NA_real_, length(y)))

    corrected <- DASRA:::.dasra_abundance_correct(
        fits = list(fit),
        taxa = "Taxon_1",
        n_samples = length(y),
        keep_diagnostics = FALSE
    )
    expect_false(corrected$formed[1L])
    expect_equal(corrected$p[1L], 1)
    expect_identical(
        corrected$reason[1L],
        "no_observed_zeros_structural_nuisance_boundary"
    )
})

test_that("public all-positive results preserve both omnibus contracts", {
    n <- 8L
    samples <- paste0("Sample_", seq_len(n))
    counts <- matrix(
        seq_len(6L * n), nrow = 6L,
        dimnames = list(paste0("Taxon_", 1:6), samples)
    )
    metadata <- data.frame(
        group = factor(rep(c("reference", "comparison"), each = n / 2L)),
        reads = rep(1000L, n),
        row.names = samples
    )

    fit <- dasra(
        counts, metadata, ~ group, "group", "reads", component = "all"
    )

    expect_true(all(fit$diagnostics$formed_structural_absence))
    expect_true(all(fit$diagnostics$nonregular_structural_absence))
    expect_false(any(fit$diagnostics$formed_relative_abundance))
    expect_true(all(
        fit$diagnostics$reason_relative_abundance ==
            "no_observed_zeros_structural_nuisance_boundary"
    ))
    expect_equal(fit$results$p_omnibus, rep(1, 6L))
    expect_identical(
        fit$results$components_used,
        rep("structural_absence", 6L)
    )
    expect_equal(fit$results$p_omnibus_cauchy, rep(1, 6L))
    expect_identical(
        fit$results$components_used_cauchy,
        rep("none", 6L)
    )
})

test_that("one observed zero continues to the abundance solver", {
    y <- c(3, 5, 4, 6, 7, 4, 8, 0)
    group <- rep(c(0, 1), each = 4L)
    gh <- DASRA:::make_count_gh_rule(3L)
    initialization_entered <- FALSE

    fit <- with_mocked_bindings(
        DASRA:::.dasra_abundance_fit_taxon(
            y = y,
            N = rep(100L, length(y)),
            group = group,
            z = NULL,
            gh_fit = gh,
            gh_effect = gh,
            control = DASRA:::.dasra_abundance_control()
        ),
        .dasra_abundance_initial_mark = function(...) {
            initialization_entered <<- TRUE
            NULL
        },
        .package = "DASRA"
    )

    expect_true(initialization_entered)
    expect_identical(fit$status, "positive_mark_initialization_failed")
})

test_that("abundance support failures take priority over the boundary", {
    gh <- DASRA:::make_count_gh_rule(3L)
    solver_entered <- FALSE
    run_case <- function(y, group, z = NULL) {
        with_mocked_bindings(
            DASRA:::.dasra_abundance_fit_taxon(
                y = y,
                N = rep(100L, length(y)),
                group = group,
                z = z,
                gh_fit = gh,
                gh_effect = gh,
                control = DASRA:::.dasra_abundance_control()
            ),
            .dasra_abundance_initial_mark = function(...) {
                solver_entered <<- TRUE
                stop("The solver should not be entered.")
            },
            .package = "DASRA"
        )
    }

    insufficient <- run_case(c(2, 3), c(0, 1))
    insufficient_mark_support <- run_case(c(2, 3, 4), c(0, 0, 1))
    one_group <- run_case(rep(2, 6L), rep(0, 6L))
    group <- rep(c(0, 1), each = 4L)
    rank_deficient <- run_case(rep(2, 8L), group, z = matrix(group))

    expect_identical(
        insufficient$status, "fewer_than_three_positive_counts"
    )
    expect_identical(
        insufficient_mark_support$status,
        "positive_mark_initialization_failed"
    )
    expect_identical(one_group$status, "positive_counts_in_one_group_only")
    expect_identical(rank_deficient$status, "rank_deficient_design")
    expect_false(solver_entered)
})

test_that("cross-taxon correction uses sample-aligned covariance", {
    fits <- .make_abundance_correction_fits()
    taxa <- paste0("Taxon_", seq_along(fits))
    corrected <- DASRA:::.dasra_abundance_correct(
        fits = fits,
        taxa = taxa,
        n_samples = length(fits[[1L]]$phi),
        keep_diagnostics = TRUE
    )

    reference <- 2:5
    background <- DASRA:::.dasra_abundance_background(
        values = vapply(
            fits[reference], `[[`, numeric(1), "raw_delta"
        ),
        standard_errors = vapply(
            fits[reference], `[[`, numeric(1), "raw_se"
        ),
        n_samples = length(fits[[1L]]$phi),
        order_key = taxa[reference]
    )
    phi_reference <- vapply(
        fits[reference], `[[`, numeric(6L), "phi"
    )
    phi_background <- as.numeric(
        phi_reference %*% background$influence_weight
    )
    kernel_background <- as.numeric(
        phi_reference %*%
            (background$kernel_weight / sum(background$kernel_weight))
    )
    expected_estimate <- fits[[1L]]$raw_delta - background$estimate
    expected_se <- sqrt(sum((fits[[1L]]$phi - phi_background)^2))
    kernel_weight_se <- sqrt(sum(
        (fits[[1L]]$phi - kernel_background)^2
    ))
    expected_z <- expected_estimate / expected_se
    expected_p <- 2 * pnorm(-abs(expected_z))
    covariance_free_se <- sqrt(
        sum(fits[[1L]]$phi^2) + sum(phi_background^2)
    )

    expect_true(corrected$formed[1L])
    expect_identical(corrected$reason[1L], "ok")
    expect_equal(corrected$estimate[1L], expected_estimate, tolerance = 1e-14)
    expect_equal(corrected$se[1L], expected_se, tolerance = 1e-14)
    expect_equal(corrected$z[1L], expected_z, tolerance = 1e-14)
    expect_equal(corrected$p[1L], expected_p, tolerance = 1e-14)
    expect_gt(abs(expected_se - covariance_free_se), 1e-4)
    expect_gt(abs(expected_se - kernel_weight_se), 1e-6)
    expect_setequal(
        corrected$diagnostics$taxon$reference_taxa[[1L]],
        taxa[reference]
    )
    expect_identical(names(corrected$diagnostics$taxon), c(
        "taxon", "raw_estimate", "raw_standard_error", "raw_p_value",
        "corrected_estimate", "corrected_standard_error",
        "corrected_p_value", "formed", "reason", "background_size",
        "background_pilot", "background_estimate", "score_residue",
        "scaled_score_residue", "root_step", "root_count",
        "bound_expansions", "numerical_warning", "jacobian_condition",
        "equilibrated_jacobian_condition", "jacobian_backward_error",
        "jacobian_rank", "mean_presence_weight",
        "mean_zero_presence_weight", "reference_taxa"
    ))
})

test_that("kernel-mode influence weights match finite differences", {
    values <- c(-0.050, -0.021, 0.002, 0.018, 0.041, 0.067, 0.530, 0.610)
    initial <- 0.01
    bandwidth <- 0.09
    mode <- DASRA:::.dasra_abundance_kernel_mode(
        values, initial, bandwidth
    )
    step <- 1e-6
    derivative <- vapply(seq_along(values), function(index) {
        plus <- minus <- values
        plus[index] <- plus[index] + step
        minus[index] <- minus[index] - step
        plus_mode <- DASRA:::.dasra_abundance_kernel_mode(
            plus, initial, bandwidth
        )
        minus_mode <- DASRA:::.dasra_abundance_kernel_mode(
            minus, initial, bandwidth
        )
        (plus_mode$estimate - minus_mode$estimate) / (2 * step)
    }, numeric(1))

    expect_true(mode$formed)
    expect_gt(mode$relative_curvature, 0)
    residual <- (values - mode$estimate) / bandwidth
    expected_curvature <- sum(
        mode$kernel_weight * (1 - residual^2)
    ) / sum(mode$kernel_weight)
    expect_equal(mode$relative_curvature, expected_curvature)
    expect_equal(
        sum(mode$kernel_weight * (values - mode$estimate)),
        0,
        tolerance = 1e-10
    )
    expect_equal(sum(mode$influence_weight), 1, tolerance = 1e-12)
    expect_equal(derivative, mode$influence_weight, tolerance = 5e-6)
})

test_that("kernel background is equivariant and order invariant", {
    values <- c(-0.10, 0, 0, 0, 0.03, 0.50, 0.60)
    standard_errors <- c(0.08, 0.04, 0.05, 0.06, 0.04, 0.30, 0.35)
    taxa <- paste0("Taxon_", seq_along(values))
    permutation <- c(6L, 2L, 7L, 4L, 1L, 5L, 3L)
    fit <- DASRA:::.dasra_abundance_background(
        values, standard_errors, 120L, taxa
    )
    shifted <- DASRA:::.dasra_abundance_background(
        values + 2.4, standard_errors, 120L, taxa
    )
    reflected <- DASRA:::.dasra_abundance_background(
        -values, standard_errors, 120L, taxa
    )
    reordered <- DASRA:::.dasra_abundance_background(
        values[permutation], standard_errors[permutation], 120L,
        taxa[permutation]
    )
    restored <- match(taxa, taxa[permutation])
    lts <- DASRA:::.dasra_abundance_lts_reference(values, taxa)
    expected_bandwidth <- median(standard_errors[lts$index]) *
        sqrt(2 * log(log(120)))

    expect_true(all(c(
        fit$formed, shifted$formed, reflected$formed, reordered$formed
    )))
    expect_equal(fit$bandwidth, expected_bandwidth)
    expect_equal(shifted$estimate, fit$estimate + 2.4)
    expect_equal(reflected$estimate, -fit$estimate)
    expect_equal(shifted$influence_weight, fit$influence_weight)
    expect_equal(reflected$influence_weight, fit$influence_weight)
    expect_equal(reordered$estimate, fit$estimate)
    expect_equal(
        reordered$influence_weight[restored], fit$influence_weight
    )
})

test_that("kernel background is invariant to labels at exact ties", {
    values <- c(rep(0, 5L), 0.4, 0.5)
    standard_errors <- c(0.02, 0.03, 0.04, 0.30, 0.40, 0.30, 0.30)
    labels <- paste0("Taxon_", seq_along(values))
    renamed <- labels
    renamed[c(1L, 5L)] <- renamed[c(5L, 1L)]

    original <- DASRA:::.dasra_abundance_background(
        values, standard_errors, 120L, labels
    )
    relabeled <- DASRA:::.dasra_abundance_background(
        values, standard_errors, 120L, renamed
    )

    expect_true(original$formed)
    expect_true(relabeled$formed)
    expect_equal(relabeled$pilot, original$pilot)
    expect_equal(relabeled$bandwidth, original$bandwidth)
    expect_equal(relabeled$estimate, original$estimate)
    expect_equal(relabeled$kernel_weight, original$kernel_weight)
    expect_equal(relabeled$influence_weight, original$influence_weight)
})

test_that("kernel background resists same-direction high-SE signals", {
    raw_delta <- c(-0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03,
        0.55, 0.60, 0.65)
    raw_se <- c(rep(0.04, 7L), rep(0.35, 3L))
    fits <- .make_abundance_reference_fits(raw_delta, raw_se)
    taxa <- paste0("Taxon_", seq_along(fits))
    corrected <- DASRA:::.dasra_abundance_correct(
        fits, taxa, 120L, keep_diagnostics = TRUE
    )

    expect_true(all(corrected$formed))
    expect_lte(max(abs(
        corrected$diagnostics$taxon$background_estimate[1:7]
    )), 0.006)
    expect_true(all(corrected$estimate[8:10] > 0.54))
})

test_that("kernel mode reports genuine mode failures", {
    not_converged <- DASRA:::.dasra_abundance_kernel_mode(
        c(-0.1, 0, 0.1), initial = 2, bandwidth = 0.1,
        max_iterations = 1L
    )
    nonpositive_curvature <- DASRA:::.dasra_abundance_kernel_mode(
        c(-1, -1, 1, 1), initial = 0, bandwidth = 0.5
    )

    expect_false(not_converged$formed)
    expect_identical(
        not_converged$reason, "background_mode_not_converged"
    )
    expect_false(nonpositive_curvature$formed)
    expect_identical(
        nonpositive_curvature$reason,
        "background_mode_nonpositive_curvature"
    )
})

test_that("target exclusion leaves its background unchanged", {
    fits <- .make_abundance_correction_fits()
    changed <- fits
    changed[[1L]]$raw_delta <- changed[[1L]]$raw_delta + 0.5
    taxa <- paste0("Taxon_", seq_along(fits))
    original <- DASRA:::.dasra_abundance_correct(
        fits, taxa, length(fits[[1L]]$phi), keep_diagnostics = TRUE
    )
    updated <- DASRA:::.dasra_abundance_correct(
        changed, taxa, length(fits[[1L]]$phi), keep_diagnostics = TRUE
    )

    expect_equal(
        updated$diagnostics$taxon$background_estimate[1L],
        original$diagnostics$taxon$background_estimate[1L]
    )
    expect_equal(updated$se[1L], original$se[1L])
    expect_equal(updated$estimate[1L] - original$estimate[1L], 0.5)
    expect_false(taxa[1L] %in%
        original$diagnostics$taxon$reference_taxa[[1L]])
    expect_null(DASRA:::.dasra_abundance_correct(
        fits, taxa, length(fits[[1L]]$phi), keep_diagnostics = FALSE
    )$diagnostics)
})

test_that("full correction is shift and sign equivariant", {
    fits <- .make_abundance_correction_fits()
    taxa <- paste0("Taxon_", seq_along(fits))
    shifted <- lapply(fits, function(fit) {
        fit$raw_delta <- fit$raw_delta + 0.4
        fit$raw_p <- 2 * pnorm(-abs(fit$raw_delta / fit$raw_se))
        fit
    })
    reflected <- lapply(fits, function(fit) {
        fit$raw_delta <- -fit$raw_delta
        fit$phi <- -fit$phi
        fit
    })
    correct <- function(input) {
        DASRA:::.dasra_abundance_correct(
            input, taxa, length(input[[1L]]$phi), keep_diagnostics = TRUE
        )
    }
    original <- correct(fits)
    translated <- correct(shifted)
    reversed <- correct(reflected)

    expect_equal(translated$estimate, original$estimate)
    expect_equal(translated$se, original$se)
    expect_equal(translated$z, original$z)
    expect_equal(translated$p, original$p)
    expect_equal(
        translated$diagnostics$taxon$background_estimate,
        original$diagnostics$taxon$background_estimate + 0.4
    )
    expect_equal(reversed$estimate, -original$estimate)
    expect_equal(reversed$se, original$se)
    expect_equal(reversed$z, -original$z)
    expect_equal(reversed$p, original$p)
})

test_that("cross-taxon correction is invariant to taxon order", {
    fits <- .make_abundance_correction_fits()
    taxa <- paste0("Taxon_", seq_along(fits))
    permutation <- c(4L, 1L, 5L, 2L, 3L)

    original <- DASRA:::.dasra_abundance_correct(
        fits, taxa, length(fits[[1L]]$phi), keep_diagnostics = TRUE
    )
    reordered <- DASRA:::.dasra_abundance_correct(
        fits[permutation], taxa[permutation],
        length(fits[[1L]]$phi), keep_diagnostics = TRUE
    )
    restored <- match(taxa, taxa[permutation])

    expect_identical(original$formed, reordered$formed[restored])
    expect_identical(original$reason, reordered$reason[restored])
    expect_equal(original$estimate, reordered$estimate[restored])
    expect_equal(original$se, reordered$se[restored])
    expect_equal(original$z, reordered$z[restored])
    expect_equal(original$p, reordered$p[restored])
    for (index in seq_along(taxa)) {
        expect_setequal(
            original$diagnostics$taxon$reference_taxa[[index]],
            reordered$diagnostics$taxon$reference_taxa[[restored[index]]]
        )
    }
})

test_that("count kernels return the same marginal log likelihood", {
    gh <- DASRA:::make_count_gh_rule(41L)
    y <- c(0, 1, 3, 12, 0, 2, 50, 900, 999)
    depth <- c(100, 100, 250, 500, rep(10000, 4L), 1000)
    eta <- c(-8, -5, -3, -1, -12, -9, -5, 2, 20)
    sigma <- 1.8

    log_only <- DASRA:::dasra_count_log_hy_adaptive_cpp(
        y, depth, eta, sigma, gh$node, gh$log_raw_weight
    )
    without_moments <- DASRA:::dasra_count_moments_adaptive_cpp(
        y, depth, eta, sigma, gh$node, gh$log_raw_weight,
        need_moments = FALSE
    )
    with_moments <- DASRA:::dasra_count_moments_adaptive_cpp(
        y, depth, eta, sigma, gh$node, gh$log_raw_weight,
        need_moments = TRUE
    )

    expect_named(without_moments, "log_h")
    expect_equal(without_moments$log_h, log_only, tolerance = 1e-12)
    expect_equal(with_moments$log_h, log_only, tolerance = 1e-12)
})

test_that("no-moment marginal evaluations use the log-only kernel", {
    gh <- DASRA:::make_count_gh_rule(3L)
    log_only_called <- FALSE
    result <- with_mocked_bindings(
        DASRA:::.dasra_abundance_marginal(
            y = c(0, 1), N = c(100, 100), eta = c(-5, -4),
            sigma = 0.8, gh = gh, need_moments = FALSE
        ),
        dasra_count_log_hy_adaptive_cpp = function(...) {
            log_only_called <<- TRUE
            c(-0.2, -1.1)
        },
        dasra_count_moments_adaptive_cpp = function(...) {
            stop("The moments kernel should not be called.")
        },
        .package = "DASRA"
    )

    expect_true(log_only_called)
    expect_identical(result, list(log_h = c(-0.2, -1.1)))
})
