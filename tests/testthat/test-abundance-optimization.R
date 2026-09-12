.make_mark_fit_fixture <- function(all_positive = FALSE, with_covariate = TRUE) {
    set.seed(801)
    n <- 80L
    group <- rep(0:1, each = n / 2L)
    depth <- rep(c(3000L, 8000L), length.out = n)
    covariate <- as.numeric(scale(seq_len(n)))
    latent <- -4.8 + 0.30 * group + 0.15 * covariate +
        rnorm(n, sd = 0.45)
    y <- rbinom(n, depth, plogis(latent))
    if (all_positive) {
        y <- pmax(y, 1L)
    } else {
        y[c(3L, 9L, 44L, 58L)] <- 0L
    }
    list(
        y = y,
        N = depth,
        group = group,
        z = if (with_covariate) matrix(covariate, ncol = 1L) else NULL
    )
}

.fit_mark_fixture <- function(all_positive = FALSE, with_covariate = TRUE) {
    input <- .make_mark_fit_fixture(all_positive, with_covariate)
    control <- SACAT:::.sacat_abundance_control()
    gh_fit <- SACAT:::.sacat_make_abundance_gh_rule(
        control$quadrature_Q
    )
    gh_effect <- gh_fit
    fit <- SACAT:::.sacat_abundance_fit_taxon(
        y = input$y,
        N = input$N,
        group = input$group,
        z = input$z,
        gh_fit = gh_fit,
        gh_effect = gh_effect,
        control = control
    )
    list(input = input, fit = fit, gh_fit = gh_fit, gh_effect = gh_effect)
}

test_that("the conditional-mark likelihood cancels structural presence", {
    y <- c(0, 2, 5, 0, 7, 3)
    N <- c(100, 120, 150, 200, 250, 300)
    group <- rep(0:1, each = 3L)
    X_eta <- cbind(Intercept = 1, Group = group)
    beta <- c(-4.5, 0.25, log(0.7))
    gh <- SACAT:::make_count_gh_rule(31L)

    component <- SACAT:::zt_beta_components(beta, y, N, X_eta, gh)
    conditional <- SACAT:::zt_beta_loglik_by_sample_inference(
        beta, y, N, X_eta, gh
    )
    positive <- y > 0
    expected <- numeric(length(y))
    expected[positive] <- component$log_hy[positive] -
        component$log_r[positive]

    expect_equal(conditional, expected, tolerance = 1e-13)
    expect_equal(conditional[!positive], rep(0, sum(!positive)))

    rho <- c(0, 0.20, 0.95, 0.35)
    joint_positive <- log1p(-rho) + component$log_hy[positive]
    detection_positive <- log1p(-rho) + component$log_r[positive]
    expect_equal(
        joint_positive - detection_positive,
        conditional[positive],
        tolerance = 1e-13
    )
})

test_that("the mark fit closes one likelihood and one influence calculation", {
    result <- .fit_mark_fixture()
    input <- result$input
    fit <- result$fit

    expect_true(fit$available)
    expect_identical(fit$status, "ok")
    expect_identical(
        fit$solver_diagnostics$estimator,
        "zero_truncated_conditional_mark"
    )
    expect_equal(
        fit$solver_diagnostics$positive_count,
        sum(input$y > 0)
    )
    expect_equal(
        fit$solver_diagnostics$score[input$y == 0, , drop = FALSE],
        matrix(
            0,
            nrow = sum(input$y == 0),
            ncol = length(fit$theta)
        )
    )
    expect_lte(fit$scaled_score_residue, 1e-8)
    expect_lte(fit$root_step, 1e-7)
    expect_equal(mean(fit$phi), 0, tolerance = 1e-15)
    expect_equal(fit$raw_se^2, sum(fit$phi^2), tolerance = 1e-14)
    expect_equal(
        fit$raw_p,
        2 * pnorm(-abs(fit$raw_delta / fit$raw_se)),
        tolerance = 1e-15
    )
    expect_true(all(
        fit$solver_diagnostics$information_eigenvalues > 0
    ))
    expect_lte(fit$jacobian_backward_error, 1e-8)
    expect_identical(fit$quadrature_Q, 41L)
    expect_false(fit$quadrature_checked)
    expect_identical(fit$quadrature_check_succeeded, NA)
})

test_that("the mark influence matches a case-weight perturbation", {
    result <- .fit_mark_fixture()
    input <- result$input
    fit <- result$fit
    index <- which(input$y > 0)[17L]
    step <- fit$solver_diagnostics$derivative_step
    X_eta <- fit$solver_diagnostics$X_eta
    X_b <- SACAT:::.sacat_abundance_designs(
        input$group, input$z
    )$X_b

    refitted_effect <- function(case_weight) {
        weight <- rep(1, length(input$y))
        weight[index] <- weight[index] + case_weight
        beta <- fit$theta
        for (iteration in seq_len(8L)) {
            weighted_loglik <- function(value) {
                weight * SACAT:::zt_beta_loglik_by_sample_inference(
                    value, input$y, input$N, X_eta, result$gh_fit
                )
            }
            linearization <- SACAT:::zt_beta_linearization(
                beta, weighted_loglik, step$step, step$scheme
            )
            score_sum <- colSums(linearization$score)
            if (max(abs(score_sum)) < 1e-10) break
            beta <- beta + as.numeric(solve(
                linearization$information, score_sum
            ))
        }
        SACAT:::.sacat_abundance_mark_effect(
            beta, X_b, result$gh_effect
        )
    }

    epsilon <- 1e-3
    derivative <- (
        refitted_effect(epsilon) - refitted_effect(-epsilon)
    ) / (2 * epsilon)
    expect_lt(abs(derivative - fit$phi[index]), 1e-6)
})

test_that("all-positive taxa use the same conditional-mark estimator", {
    result <- .fit_mark_fixture(all_positive = TRUE)
    fit <- result$fit

    expect_true(all(result$input$y > 0))
    expect_true(fit$available)
    expect_identical(fit$status, "ok")
    expect_identical(
        fit$solver_diagnostics$estimator,
        "zero_truncated_conditional_mark"
    )
    expect_equal(
        fit$solver_diagnostics$positive_count,
        length(result$input$y)
    )
    expect_true(all(is.finite(c(
        fit$raw_delta, fit$raw_se, fit$raw_p
    ))))
})

test_that("mark formation reports the mathematical support failures", {
    control <- SACAT:::.sacat_abundance_control()
    gh <- SACAT:::make_count_gh_rule(control$quadrature_Q)
    run_case <- function(y, group, z = NULL, minimum = 3L) {
        SACAT:::.sacat_abundance_fit_taxon(
            y = y,
            N = rep(1000L, length(y)),
            group = group,
            z = z,
            gh_fit = gh,
            gh_effect = gh,
            control = control,
            min_positive_samples = minimum
        )
    }

    insufficient <- run_case(c(2, 0, 0, 3), c(0, 0, 1, 1))
    one_group <- run_case(
        c(2, 3, 4, 5, 0, 0, 0, 0),
        rep(0:1, each = 4L)
    )
    information <- run_case(c(2, 3, 4, 0), c(0, 0, 1, 1))
    group <- rep(0:1, each = 4L)
    rank_deficient <- run_case(
        rep(2, 8L), group, z = matrix(group, ncol = 1L)
    )

    expect_identical(insufficient$status, "insufficient_positive_support")
    expect_identical(one_group$status, "positive_counts_in_one_group_only")
    expect_identical(
        information$status,
        "insufficient_positive_mark_information"
    )
    expect_identical(
        rank_deficient$status,
        "rank_deficient_positive_mark_design"
    )
})

test_that("structural start modes have explicit public meanings", {
    expect_identical(
        SACAT:::.sacat_resolve_structural_conditional_present_starts(
            "adaptive"
        ),
        list(mode = "adaptive", count = 1L)
    )
    expect_identical(
        SACAT:::.sacat_resolve_structural_conditional_present_starts(
            "full"
        ),
        list(mode = "full", count = 5L)
    )
    expect_identical(
        SACAT:::.sacat_resolve_structural_conditional_present_starts(1L),
        list(mode = "adaptive", count = 1L)
    )
    expect_identical(
        SACAT:::.sacat_resolve_structural_conditional_present_starts(5L),
        list(mode = "full", count = 5L)
    )
    expect_error(
        SACAT:::.sacat_resolve_structural_conditional_present_starts(
            "primary"
        ),
        "adaptive.*full"
    )
})

test_that("abundance quadrature preserves its default numerical rule", {
    expect_identical(
        SACAT:::.sacat_make_abundance_gh_rule(41L),
        SACAT:::make_count_gh_rule(41L)
    )

    stable <- SACAT:::.sacat_make_abundance_gh_rule(81L)
    expect_identical(stable$Q, 81L)
    expect_true(all(is.finite(stable$node)))
    expect_true(all(is.finite(stable$log_raw_weight)))
    expect_false(is.unsorted(stable$node, strictly = TRUE))
    expect_equal(sum(exp(stable$log_weight)), 1, tolerance = 1e-12)

    expect_identical(
        vapply(
            c(41L, 1001L, 1002L, 2001L),
            SACAT:::.sacat_higher_order_quadrature_points,
            integer(1)
        ),
        c(2001L, 2001L, 2003L, 4001L)
    )
})

test_that("abundance quadrature sensitivity is a fixed-fit diagnostic", {
    y <- c(0, 1, 2, 1, 3, 0, 2, 4)
    N <- c(40, 60, 80, 100, 40, 60, 80, 100)
    group <- rep(0:1, each = 4L)
    X_eta <- cbind(Intercept = 1, Group = group)
    X_b <- matrix(1, nrow = length(y), ncol = 1L)
    beta <- c(-8, 0.25, log(8))
    gh <- SACAT:::.sacat_make_abundance_gh_rule(41L)
    effect <- SACAT:::.sacat_abundance_mark_effect(beta, X_b, gh)

    diagnostic <- SACAT:::.sacat_abundance_quadrature_diagnostic(
        beta, y, N, X_eta, X_b, gh, effect,
        check_quadrature = TRUE
    )
    expect_true(diagnostic$quadrature_checked)
    expect_true(diagnostic$quadrature_check_succeeded)
    expect_identical(diagnostic$quadrature_comparison_Q, 2001L)
    expect_true(all(is.finite(c(
        diagnostic$quadrature_conditional_max_abs,
        diagnostic$quadrature_effect_abs
    ))))
    expect_gte(diagnostic$quadrature_conditional_max_abs, 0)
    expect_gte(diagnostic$quadrature_effect_abs, 0)

    failed <- with_mocked_bindings(
        .sacat_abundance_quadrature_diagnostic(
            beta, y, N, X_eta, X_b, gh, effect,
            check_quadrature = TRUE
        ),
        make_structural_gh_rule = function(...) stop("controlled failure"),
        .package = "SACAT"
    )
    expect_true(failed$quadrature_checked)
    expect_false(failed$quadrature_check_succeeded)
    expect_match(failed$quadrature_check_error, "controlled failure")

    boundary_beta <- c(-38.3, 0, log(3))
    boundary_y <- c(0, 1)
    boundary_N <- c(5, 5)
    boundary_group <- 0:1
    boundary_X_eta <- cbind(Intercept = 1, Group = boundary_group)
    boundary_X_b <- matrix(1, nrow = 2L, ncol = 1L)
    boundary_effect <- SACAT:::.sacat_abundance_mark_effect(
        boundary_beta, boundary_X_b, gh
    )
    boundary_diagnostic <-
        SACAT:::.sacat_abundance_quadrature_diagnostic(
            boundary_beta, boundary_y, boundary_N,
            boundary_X_eta, boundary_X_b, gh, boundary_effect,
            check_quadrature = TRUE
        )
    expect_true(boundary_diagnostic$quadrature_checked)
    expect_false(boundary_diagnostic$quadrature_check_succeeded)
    expect_match(
        boundary_diagnostic$quadrature_check_error,
        "non-finite"
    )
})
