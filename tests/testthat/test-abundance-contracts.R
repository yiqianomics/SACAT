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
    phi_background <- rowMeans(vapply(
        fits[reference], `[[`, numeric(6L), "phi"
    ))
    expected_estimate <- fits[[1L]]$raw_delta - mean(vapply(
        fits[reference], `[[`, numeric(1), "raw_delta"
    ))
    expected_se <- sqrt(sum((fits[[1L]]$phi - phi_background)^2))
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
    expect_setequal(
        corrected$diagnostics$taxon$reference_taxa[[1L]],
        taxa[reference]
    )
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
