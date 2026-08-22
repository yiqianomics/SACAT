.abundance_optimization_inputs <- function() {
    n <- 8L
    list(
        y = c(1, 2, 0, 1, 0, 2, 0, 1),
        N = rep(100, n),
        group = rep(c(0, 1), each = n / 2L)
    )
}

.abundance_expansion_input <- function() {
    n <- 72L
    group <- rep(c(0, 1), each = n / 2L)
    depth <- rep(10000L, n)
    index <- seq_len(n)
    y <- as.integer(round(exp(
        log(35) + 0.5 * sin(index * 1.7) + 0.1 * group
    )))
    y[index %% 9L == 0L] <- 0L
    list(y = y, N = depth, group = group)
}

.run_abundance_classification_case <- function(mode) {
    stopifnot(mode %in% c(
        "strict_first", "plateau_only", "reject",
        "effect_disagreement"
    ))
    input <- .abundance_optimization_inputs()
    control <- DASRA:::.dasra_abundance_control()
    gh <- DASRA:::make_count_gh_rule(3L)
    start_first <- switch(
        mode,
        strict_first = 0,
        effect_disagreement = 0,
        plateau_only = 1,
        reject = 2
    )

    mock_mark <- function(data, gh) {
        list(par = c(start_first, 0, 0), value = 0)
    }
    mock_detection <- function(data, b, omega, zeta, gh) {
        list(par = c(0, 0), value = 0)
    }
    mock_state <- function(theta, data, gh, return_psi = TRUE) {
        centered <- rep(c(-1, 1), length.out = length(data$y))
        centered <- centered - mean(centered)
        equation_mean <- if (abs(theta[1L]) < 0.25) {
            5e-8
        } else if (abs(theta[1L] - 1) < 0.25) {
            1e-12
        } else {
            1e-10
        }
        psi <- matrix(
            0,
            nrow = length(data$y),
            ncol = data$layout$dimension
        )
        psi[, 1L] <- centered + equation_mean
        psi[, data$layout$zeta] <- centered
        list(
            psi = psi,
            presence_weight = rep(0.5, length(data$y))
        )
    }
    mock_jacobian <- function(theta, fn, ...) {
        jacobian <- diag(length(theta))
        jacobian[1L, 1L] <- if (abs(theta[1L]) < 0.25) {
            1
        } else if (abs(theta[1L] - 1) < 0.25) {
            1e-8
        } else {
            5e-10
        }
        jacobian
    }
    mock_gradient <- function(theta, fn, ...) {
        gradient <- numeric(length(theta))
        gradient[length(theta)] <- 1
        gradient
    }
    mock_effect <- if (identical(mode, "effect_disagreement")) {
        function(theta, data, gh_effect) 0.25 + theta[1L]
    } else {
        function(theta, data, gh_effect) {
            0.25 + theta[data$layout$zeta]
        }
    }
    mock_nleqslv <- function(x, fn, ...) {
        if (mode %in% c("strict_first", "effect_disagreement") &&
            abs(x[1L]) < 0.25) {
            x[1L] <- 1
        }
        list(x = x, termcd = 1L)
    }

    with_mocked_bindings(
        with_mocked_bindings(
            DASRA:::.dasra_abundance_fit_taxon(
                y = input$y,
                N = input$N,
                group = input$group,
                z = NULL,
                gh_fit = gh,
                gh_effect = gh,
                control = control
            ),
            nleqslv = mock_nleqslv,
            .package = "nleqslv"
        ),
        .dasra_abundance_initial_mark = mock_mark,
        .dasra_abundance_initial_detection = mock_detection,
        .dasra_abundance_state = mock_state,
        .dasra_abundance_numeric_jacobian = mock_jacobian,
        .dasra_abundance_numeric_gradient = mock_gradient,
        .dasra_abundance_effect = mock_effect,
        .package = "DASRA"
    )
}

.run_persistent_boundary_case <- function() {
    input <- .abundance_optimization_inputs()
    control <- DASRA:::.dasra_abundance_control()
    control$max_bound_expansions <- 1L
    gh <- DASRA:::make_count_gh_rule(3L)
    design <- DASRA:::.dasra_abundance_designs(input$group, NULL)
    layout <- DASRA:::.dasra_abundance_layout(
        design$X_b, design$X_rho
    )
    bounds <- DASRA:::.dasra_abundance_bounds(layout)
    target <- numeric(layout$dimension)
    target[layout$b[1L]] <- bounds$upper[layout$b[1L]] +
        bounds$upper[layout$b[1L]] - bounds$lower[layout$b[1L]]

    mock_mark <- function(data, gh) list(par = c(0, 0, 0), value = 0)
    mock_detection <- function(data, b, omega, zeta, gh) {
        list(par = c(0, 0), value = 0)
    }
    mock_state <- function(theta, data, gh, return_psi = TRUE) {
        psi <- matrix(
            rep(theta - target, each = length(data$y)),
            nrow = length(data$y)
        )
        list(
            psi = psi,
            presence_weight = rep(0.5, length(data$y))
        )
    }
    mock_jacobian <- function(theta, fn, ...) diag(length(theta))
    mock_nleqslv <- function(x, fn, ...) list(x = x, termcd = 1L)
    mock_optim <- function(par, fn, method, lower, upper, control, ...) {
        point <- pmin(pmax(target, lower), upper)
        list(par = point, value = fn(point), convergence = 0L)
    }

    with_mocked_bindings(
        with_mocked_bindings(
            DASRA:::.dasra_abundance_fit_taxon(
                y = input$y,
                N = input$N,
                group = input$group,
                z = NULL,
                gh_fit = gh,
                gh_effect = gh,
                control = control
            ),
            nleqslv = mock_nleqslv,
            .package = "nleqslv"
        ),
        .dasra_abundance_initial_mark = mock_mark,
        .dasra_abundance_initial_detection = mock_detection,
        .dasra_abundance_state = mock_state,
        .dasra_abundance_numeric_jacobian = mock_jacobian,
        optim = mock_optim,
        .package = "DASRA"
    )
}

test_that("a finite root outside initial bounds is recovered by expansion", {
    input <- .abundance_expansion_input()
    control <- DASRA:::.dasra_abundance_control()
    gh <- DASRA:::make_count_gh_rule(control$quadrature_Q)
    original_bounds <- DASRA:::.dasra_abundance_bounds
    narrow_bounds <- function(layout) {
        bounds <- original_bounds(layout)
        bounds$lower[layout$b[1L]] <- -5.4
        bounds$upper[layout$b[1L]] <- -4.0
        bounds
    }

    fit <- with_mocked_bindings(
        DASRA:::.dasra_abundance_fit_taxon(
            y = input$y,
            N = input$N,
            group = input$group,
            z = NULL,
            gh_fit = gh,
            gh_effect = gh,
            control = control
        ),
        .dasra_abundance_bounds = narrow_bounds,
        .package = "DASRA"
    )

    expect_true(fit$available)
    expect_identical(fit$status, "ok")
    expect_gte(fit$bound_expansions, 1L)
    expect_true("expanded_numerical_bounds" %in% fit$numerical_warning)
    expect_true(fit$theta[1L] < -5.4 || fit$theta[1L] > -4.0)
    expect_lt(fit$solver_diagnostics$working_lower[1L], -5.4)
    expect_true(all(
        fit$theta >= fit$solver_diagnostics$working_lower &
            fit$theta <= fit$solver_diagnostics$working_upper
    ))
    expect_false(isTRUE(fit$solver_diagnostics$boundary_following))
})

test_that("a candidate that follows the final boundary is rejected", {
    fit <- .run_persistent_boundary_case()

    expect_false(fit$available)
    expect_identical(fit$status, "persistent_numerical_boundary")
    expect_true(isTRUE(fit$solver_diagnostics$boundary_following))
    expect_identical(fit$bound_expansions, 1L)
    expect_gte(fit$root_count, 1L)
    distance <- pmin(
        fit$solver_diagnostics$selected_theta -
            fit$solver_diagnostics$working_lower,
        fit$solver_diagnostics$working_upper -
            fit$solver_diagnostics$selected_theta
    )
    expect_true(any(distance <= 1e-12))
})

test_that("strict roots take priority and plateau-only roots are bounded", {
    strict_first <- .run_abundance_classification_case("strict_first")
    strict_repeat <- .run_abundance_classification_case("strict_first")
    plateau_only <- .run_abundance_classification_case("plateau_only")
    plateau_repeat <- .run_abundance_classification_case("plateau_only")
    rejected <- .run_abundance_classification_case("reject")

    expect_true(strict_first$available)
    expect_true(any(strict_first$solver_diagnostics$root_strict))
    expect_true(any(
        !strict_first$solver_diagnostics$root_strict &
            strict_first$solver_diagnostics$root_plateau
    ))
    expect_gt(
        min(strict_first$solver_diagnostics$root_raw_residue[
            strict_first$solver_diagnostics$root_strict
        ]),
        min(strict_first$solver_diagnostics$root_raw_residue[
            !strict_first$solver_diagnostics$root_strict &
                strict_first$solver_diagnostics$root_plateau
        ])
    )
    expect_equal(strict_first$theta[1L], 0)
    expect_false(
        "weakly_identified_root_plateau" %in%
            strict_first$numerical_warning
    )
    expect_identical(
        strict_first$solver_diagnostics$root_strict,
        strict_repeat$solver_diagnostics$root_strict
    )
    expect_identical(
        strict_first$solver_diagnostics$root_plateau,
        strict_repeat$solver_diagnostics$root_plateau
    )

    expect_true(plateau_only$available)
    expect_false(any(plateau_only$solver_diagnostics$root_strict))
    expect_true(any(plateau_only$solver_diagnostics$root_plateau))
    expect_true(
        "weakly_identified_root_plateau" %in%
            plateau_only$numerical_warning
    )
    expect_identical(
        plateau_only$solver_diagnostics$root_strict,
        plateau_repeat$solver_diagnostics$root_strict
    )
    expect_identical(
        plateau_only$solver_diagnostics$root_plateau,
        plateau_repeat$solver_diagnostics$root_plateau
    )

    expect_false(rejected$available)
    expect_identical(rejected$status, "estimating_equation_root_step")
    root_step_limit <- DASRA:::.dasra_abundance_control()$root_step_limit
    expect_gt(rejected$root_step, root_step_limit)
})

test_that("disagreeing multiple-root effects leave the test unformed", {
    fit <- .run_abundance_classification_case("effect_disagreement")

    expect_false(fit$available)
    expect_identical(fit$status, "multiple_root_effect_disagreement")
    expect_gte(fit$root_count, 2L)
    effects <- fit$solver_diagnostics$root_effects
    expect_gt(
        diff(range(effects)),
        1e-4 * (1 + max(abs(effects)))
    )

    corrected <- DASRA:::.dasra_abundance_correct(
        fits = list(fit),
        taxa = "Taxon_1",
        n_samples = length(fit$phi),
        keep_diagnostics = FALSE
    )
    expect_false(corrected$formed[1L])
    expect_equal(corrected$p[1L], 1)
    expect_identical(
        corrected$reason[1L],
        "multiple_root_effect_disagreement"
    )
})
