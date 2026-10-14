.shared_mark_fixture <- function() {
    set.seed(904)
    group <- rep(0:1, each = 40L)
    z <- matrix(rep(seq(-1, 1, length.out = 40L), 2L), ncol = 1L)
    N <- rep(c(4000L, 8000L), length.out = length(group))
    y <- rbinom(length(group), N,
                plogis(-5 + 0.3 * group + 0.2 * z[, 1L] +
                           rnorm(length(group), sd = 0.5)))
    y[seq(1L, length(y), by = 5L)] <- 0L
    list(y = y, N = N, group = group, z = z,
         X_eta = SACAT:::.sacat_abundance_designs(group, z)$X_eta,
         gh = SACAT:::.sacat_make_abundance_gh_rule(41L))
}

test_that("both components reuse one five-start search", {
    x <- .shared_mark_fixture()
    Y <- matrix(x$y, ncol = 1L, dimnames = list(NULL, "Taxon"))
    optimize <- SACAT:::zt_fit_beta
    calls <- 0L
    local_mocked_bindings(
        zt_fit_beta = function(...) {
            calls <<- calls + 1L
            optimize(...)
        },
        .package = "SACAT"
    )
    shared <- SACAT:::.sacat_conditional_present_fits(
        Y, x$N, x$group, x$z, quadrature_points = 41L
    )
    expect_identical(calls, 1L)
    expect_true(shared[[1L]]$ok)
    expect_identical(shared[[1L]]$optimization$n_starts, 5L)
    expect_false(any(vapply(shared[[1L]], is.function, logical(1))))

    control <- SACAT:::.sacat_abundance_control(41L)
    ra <- SACAT:::.sacat_abundance_fit_taxon(
        x$y, x$N, x$group, x$z, x$gh, x$gh, control,
        conditional_present_fit = shared[[1L]]
    )
    sa <- SACAT:::zt_count_structural_test(
        x$y, x$N, x$group, x$z, keep_fit = TRUE,
        check_quadrature = FALSE,
        conditional_present_fit = shared[[1L]]
    )
    expect_true(ra$available)
    expect_true(sa$tested)
    expect_identical(calls, 1L)

    independent <- SACAT:::.sacat_abundance_fit_taxon(
        x$y, x$N, x$group, x$z, x$gh, x$gh, control
    )
    expect_identical(ra$theta, independent$theta)
    expect_identical(ra$raw_delta, independent$raw_delta)
    expect_identical(ra$raw_se, independent$raw_se)
    expect_identical(ra$phi, independent$phi)
})

test_that("SA refinement and diagnostics use the structural quadrature", {
    x <- .shared_mark_fixture()
    gh <- SACAT:::make_structural_gh_rule(1001L)
    control <- SACAT:::.sacat_abundance_control(1001L)
    search <- SACAT:::zt_fit_beta(x$y, x$N, x$X_eta, x$gh)
    refined <- SACAT:::.sacat_refine_structural_beta(
        search, x$y, x$N, x$X_eta, gh, control
    )
    direct <- SACAT:::zt_fit_beta(x$y, x$N, x$X_eta, gh)
    direct_refined <- SACAT:::.sacat_refine_structural_beta(
        direct, x$y, x$N, x$X_eta, gh, control
    )
    expect_true(refined$ok)
    expect_identical(refined$quadrature_Q, 1001L)
    expect_equal(refined$par, direct_refined$par, tolerance = 1e-7)
    expect_lt(max(abs(refined$gradient)), control$score_tolerance)
    expect_true(all(refined$eigenvalues > 0))
    expected <- SACAT:::zt_beta_nll(refined$par, x$y, x$N, x$X_eta, gh)
    expect_equal(refined$nll(refined$par), expected, tolerance = 1e-12)
    expect_equal(refined$optimization_objective(refined$par),
                 expected / sum(x$y > 0), tolerance = 1e-12)
    expect_equal(refined$hessian,
                 refined$root_polish$linearization$information /
                     sum(x$y > 0), tolerance = 1e-12)
    expect_false(refined$bound_tracking)
    expect_null(refined[["optimization"]])
})

test_that("SA can fit independently when the search or root refinement fails", {
    x <- .shared_mark_fixture()
    gh <- SACAT:::make_structural_gh_rule(81L)
    control <- SACAT:::.sacat_abundance_control(81L)
    failed_search <- list(ok = FALSE, reason = "conditional_present_gradient")
    refined <- SACAT:::.sacat_refine_structural_beta(
        failed_search, x$y, x$N, x$X_eta, gh, control
    )
    expect_true(refined$ok)
    expect_true(refined$structural_search$ok)
    expect_identical(refined$search_refinement_reason, failed_search$reason)
    expect_identical(refined$quadrature_Q, 81L)

    search <- SACAT:::zt_fit_beta(x$y, x$N, x$X_eta, x$gh)
    fit <- with_mocked_bindings(
        SACAT:::.sacat_refine_structural_beta(
            search, x$y, x$N, x$X_eta, gh, control
        ),
        .sacat_polish_mark_root = function(...) {
            list(ok = FALSE, reason = "conditional_present_root_polish_no_descent")
        },
        .package = "SACAT"
    )
    expect_true(fit$ok)
    expect_true(fit$structural_search$ok)
    expect_identical(fit$refinement_reason,
                     "conditional_present_root_polish_no_descent")
    expect_equal(fit$nll(fit$par),
                 SACAT:::zt_beta_nll(fit$par, x$y, x$N, x$X_eta, gh))
})

test_that("shared search respects all-positive and rank-deficient designs", {
    x <- .shared_mark_fixture()
    Y <- cbind(positive = pmax(x$y, 1L), one_group = x$y)
    Y[x$group == 1L, "one_group"] <- 0L
    shared <- with_mocked_bindings(
        SACAT:::.sacat_conditional_present_fits(
            Y, x$N, x$group, x$z, 41L, fit_all_positive = FALSE
        ),
        zt_fit_beta = function(...) stop("No optimization is required."),
        .package = "SACAT"
    )
    expect_null(shared$positive)
    expect_false(shared$one_group$ok)
    result <- SACAT:::zt_count_structural_test(
        Y[, "positive"], x$N, x$group, x$z,
        conditional_present_fit = list(ok = FALSE)
    )
    expect_true(result$tested)
    expect_equal(result$p, 1)
})
