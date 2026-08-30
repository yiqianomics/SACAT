test_that("the default structural engine forms a regular result", {
    n <- 80L
    group <- rep(c(0, 1), each = n / 2L)
    depth <- rep(c(4000, 5500, 7000, 8500), length.out = n)
    index <- seq_len(n)
    y <- as.integer(7 + ((index * 7L) %% 13L) + 3L * group)
    y[index %% 5L == 0L | index %% 11L == 0L] <- 0L

    fit <- DASRA:::zt_count_structural_test(
        y = y,
        N = depth,
        g = group,
        keep_fit = TRUE
    )

    expect_true(fit$tested)
    expect_true(fit$regular)
    expect_identical(fit$reason, "ok")
    expect_true(is.finite(fit$p) && fit$p >= 0 && fit$p <= 1)
    expect_true(is.finite(fit$diagnostics$score_z))
    expect_gt(fit$diagnostics$V, 0)
    expect_true(all(is.finite(fit$gamma)))
    expect_true(all(fit$gamma >= 0 & fit$gamma <= 1))
    expect_true(all(is.finite(fit$rho)))
    expect_true(all(fit$rho > 0 & fit$rho < 1))
    expect_identical(fit$diagnostics$fit$quadrature$Q, 1001L)
    expect_equal(
        fit$diagnostics$U,
        sum(fit$diagnostics$fit$adjusted_score),
        tolerance = 1e-12
    )
    expect_equal(
        fit$diagnostics$V,
        sum(fit$diagnostics$fit$adjusted_score^2),
        tolerance = 1e-12
    )
    expect_equal(
        fit$diagnostics$U_raw,
        sum(fit$diagnostics$fit$target_score),
        tolerance = 1e-12
    )
    expect_lt(
        abs(
            fit$diagnostics$nuisance_score_correction -
                sum(fit$diagnostics$fit$estimating_functions %*%
                    fit$diagnostics$fit$nuisance_adjustment)
        ),
        1e-12
    )
    expect_equal(
        fit$diagnostics$score_z,
        fit$diagnostics$U / sqrt(fit$diagnostics$V),
        tolerance = 1e-14
    )
})

test_that("structural inference uses one orthogonalized contribution for U and V", {
    run_candidate <- function(psi, target, adjustment) {
        with_mocked_bindings(
            DASRA:::zt_structural_inference_candidate(
                beta = 0, alpha = 0, y = 0, N = 1, g = 0,
                X_rho = matrix(1), X_eta = matrix(1), gh = list(),
                detection_component = list(),
                step_info = list(
                    base = 1, step = 1, scheme = "central"
                ),
                multiplier = 1
            ),
            zt_structural_linearization = function(...) {
                list(
                    ok = TRUE,
                    psi = psi,
                    A = diag(ncol(psi)),
                    target = target,
                    M = adjustment
                )
            },
            zt_strict_projection = function(...) {
                list(
                    ok = TRUE,
                    adjustment = adjustment,
                    backward_error = 0,
                    rank = ncol(psi),
                    condition = 1,
                    equilibrated_condition = 1
                )
            },
            .package = "DASRA"
        )
    }

    target <- c(0.4, -0.2, 0.3, 0.1)
    adjustment <- c(0.25, -0.5)
    exact_psi <- rbind(
        c(1, 0), c(-1, 1), c(2, -1), c(-2, 0)
    )
    exact <- run_candidate(exact_psi, target, adjustment)
    exact_influence <- target - as.numeric(exact_psi %*% adjustment)

    expect_equal(colSums(exact_psi), c(0, 0))
    expect_equal(exact$U_raw, sum(target))
    expect_equal(exact$U, exact$U_raw, tolerance = 1e-15)
    expect_equal(exact$U, sum(exact_influence), tolerance = 1e-15)
    expect_equal(exact$V, sum(exact_influence^2), tolerance = 1e-15)

    approximate_psi <- rbind(
        c(1, 0), c(-0.5, 1), c(0.25, -0.25), c(0.5, 0.5)
    )
    approximate <- run_candidate(approximate_psi, target, adjustment)
    approximate_influence <- target -
        as.numeric(approximate_psi %*% adjustment)

    expect_equal(approximate$U_raw, 0.6, tolerance = 1e-15)
    expect_equal(approximate$U, 0.9125, tolerance = 1e-15)
    expect_equal(approximate$V, 0.26640625, tolerance = 1e-15)
    expect_equal(approximate$score_z, 1.7679121, tolerance = 1e-7)
    expect_equal(approximate$p, 0.0770756, tolerance = 1e-7)
    expect_equal(
        approximate$nuisance_score_correction,
        sum(approximate_psi %*% adjustment),
        tolerance = 1e-15
    )
    expect_equal(approximate$U, sum(approximate_influence))
    expect_equal(
        approximate$p,
        pchisq(approximate$U^2 / approximate$V, 1,
               lower.tail = FALSE),
        tolerance = 1e-15
    )
    expect_false(isTRUE(all.equal(
        approximate$p,
        pchisq(approximate$U_raw^2 / approximate$V, 1,
               lower.tail = FALSE)
    )))

    summary <- DASRA:::zt_structural_candidate_summary(
        list(approximate)
    )
    expect_equal(summary$U_raw, approximate$U_raw)
    expect_equal(
        summary$nuisance_score_correction,
        approximate$nuisance_score_correction
    )
    expect_equal(summary$score_z_raw, approximate$score_z_raw)
})

test_that("the exact all-zero structural competitor uses the same detection objective", {
    log_r <- log(c(0.15, 0.35, 0.65, 0.85))
    y <- c(0, 2, 0, 1)
    expected <- -sum(c(
        DASRA:::zt_log1mexp(log_r[1]),
        log_r[2],
        DASRA:::zt_log1mexp(log_r[3]),
        log_r[4]
    ))
    zero_nll <- DASRA:::zt_structural_zero_limit_nll(log_r, y)
    relative_tolerance <- sqrt(.Machine$double.eps)
    comparison_tolerance <- relative_tolerance * max(1, abs(zero_nll))

    expect_equal(zero_nll, expected, tolerance = 1e-15)

    dominated <- DASRA:::zt_structural_zero_limit_comparison(
        zero_nll + 1e-3, log_r, y,
        relative_tolerance = relative_tolerance
    )
    within_tolerance <- DASRA:::zt_structural_zero_limit_comparison(
        zero_nll + comparison_tolerance / 2, log_r, y,
        relative_tolerance = relative_tolerance
    )
    finite_better <- DASRA:::zt_structural_zero_limit_comparison(
        zero_nll - 1e-3, log_r, y,
        relative_tolerance = relative_tolerance
    )

    expect_true(dominated$dominated)
    expect_false(within_tolerance$dominated)
    expect_false(finite_better$dominated)
    expect_lt(abs(dominated$improvement - 1e-3), 1e-14)
    expect_equal(
        dominated$tolerance,
        relative_tolerance * max(
            1, abs(dominated$finite_nll), abs(dominated$zero_limit_nll)
        )
    )
})

test_that("a finite structural fit dominated by the exact zero limit is unavailable", {
    y <- c(1, 0, 1, 0, 1, 0, 1, 0)
    N <- rep(100L, length(y))
    g <- rep(c(0, 1), each = length(y) / 2)
    z <- matrix(rep(c(-1, 0, 1, 2), 2), ncol = 1L)
    adaptive_called <- FALSE

    result <- with_mocked_bindings(
        DASRA:::zt_count_structural_test(
            y = y, N = N, g = g, z = z, Q = 3L,
            min_positive_samples = 3L, keep_fit = TRUE
        ),
        make_structural_gh_rule = function(Q) list(Q = Q),
        zt_fit_beta = function(...) {
            list(
                ok = TRUE, reason = "ok", par = rep(0, 4),
                numerical_warnings = character()
            )
        },
        zt_beta_detection_components = function(...) {
            list(log_r = rep(log(0.5), length(y)))
        },
        zt_fit_alpha = function(...) {
            list(
                ok = TRUE, reason = "ok", par = c(0, 0),
                numerical_warnings = character()
            )
        },
        zt_detection_nll_from_components = function(...) 10,
        zt_adaptive_structural_inference = function(...) {
            adaptive_called <<- TRUE
            stop("adaptive inference must not be reached")
        },
        .package = "DASRA"
    )

    expect_false(adaptive_called)
    expect_false(result$tested)
    expect_false(result$regular)
    expect_equal(result$p, 1)
    expect_identical(
        result$reason, "structural_absence_nonoptimal_nuisance_fit"
    )
    expect_gt(
        result$diagnostics$structural_absence_zero_limit_improvement,
        result$diagnostics$structural_absence_objective_tolerance
    )
    expect_identical(
        result$diagnostics$fit$structural_absence$reason, "ok"
    )
})

test_that("the Q=1001 structural quadrature rule preserves its invariants", {
    gh <- DASRA:::make_structural_gh_rule(1001L)
    midpoint <- (gh$Q + 1L) / 2L
    largest_log_weight <- max(gh$log_weight)
    log_weight_sum <- largest_log_weight + log(sum(exp(
        gh$log_weight - largest_log_weight
    )))

    expect_identical(gh$Q, 1001L)
    expect_length(gh$node, 1001L)
    expect_length(gh$log_weight, 1001L)
    expect_true(all(is.finite(gh$node)))
    expect_true(all(is.finite(gh$log_weight)))
    expect_true(all(diff(gh$node) > 0))
    expect_equal(gh$node, -rev(gh$node), tolerance = 1e-10)
    expect_equal(
        gh$log_weight, rev(gh$log_weight), tolerance = 1e-8
    )
    expect_equal(gh$node[midpoint], 0, tolerance = 1e-12)
    expect_equal(log_weight_sum, 0, tolerance = 1e-12)
    expect_equal(sum(gh$weight), 1, tolerance = 1e-12)
    expect_equal(gh$weight, rev(gh$weight), tolerance = 1e-12)
    expect_equal(sum(gh$weight * gh$node), 0, tolerance = 1e-13)
    expect_equal(sum(gh$weight * gh$node^2), 0.5, tolerance = 1e-12)
})
