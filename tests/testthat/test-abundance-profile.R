.abundance_profile_fixture <- local({
    cache <- new.env(parent = emptyenv())
    function(with_covariate = FALSE) {
        key <- as.character(with_covariate)
        if (exists(key, cache, inherits = FALSE)) return(cache[[key]])
        set.seed(1809)
        n <- 80L
        group <- rep(0:1, each = n / 2L)
        depth <- rep(c(3000L, 8000L), length.out = n)
        covariate <- rep(seq(-1, 1, length.out = n / 2L), 2L)
        z <- if (with_covariate) matrix(covariate, ncol = 1L) else NULL
        latent <- -4.8 + 0.35 * group + 0.2 * covariate +
            rnorm(n, sd = 0.45)
        y <- rbinom(n, depth, plogis(latent))
        y[c(3L, 9L, 44L, 58L)] <- 0L
        control <- SACAT:::.sacat_abundance_control()
        gh <- SACAT:::.sacat_make_abundance_gh_rule(control$quadrature_Q)
        design <- SACAT:::.sacat_abundance_designs(group, z)
        fit <- SACAT:::.sacat_abundance_fit_taxon(
            y, depth, group, z, gh, gh, control
        )
        cache[[key]] <- list(
            fit = fit, y = y, N = depth, gh = gh,
            X = design$X_eta, X_b = design$X_b
        )
        cache[[key]]
    }
})

.abundance_profile_family <- local({
    fixture <- NULL
    function() {
        if (!is.null(fixture)) return(fixture)
        set.seed(1810)
        n <- 80L
        group <- rep(0:1, each = n / 2L)
        depth <- rep(c(3000L, 8000L), length.out = n)
        effects <- c(0.35, -0.3, 0, 0, 0, 0.1)
        Y <- vapply(seq_along(effects), function(j) {
            latent <- -5 + j / 10 + effects[j] * group +
                rnorm(n, sd = 0.4)
            y <- rbinom(n, depth, plogis(latent))
            y[c(3L, 44L)] <- 0L
            y
        }, integer(n))
        colnames(Y) <- paste0("Taxon_", seq_along(effects))
        control <- SACAT:::.sacat_abundance_control()
        gh <- SACAT:::.sacat_make_abundance_gh_rule(control$quadrature_Q)
        fits <- lapply(seq_len(ncol(Y)), function(j) {
            SACAT:::.sacat_abundance_fit_taxon(
                Y[, j], depth, group, NULL, gh, gh, control
            )
        })
        names(fits) <- colnames(Y)
        corrected <- SACAT:::.sacat_abundance_correct(
            fits, colnames(Y), n, keep_diagnostics = TRUE
        )
        X_b <- SACAT:::.sacat_abundance_designs(group, NULL)$X_b
        tested <- SACAT:::.sacat_abundance_profile_tests(
            corrected, fits, Y, depth, X_b, gh
        )
        fixture <<- list(
            corrected = corrected, tested = tested, fits = fits,
            Y = Y, N = depth, X_b = X_b, gh = gh
        )
        fixture
    }
})

test_that("profile derivatives agree with the fitted conditional likelihood", {
    for (with_covariate in c(FALSE, TRUE)) {
        input <- .abundance_profile_fixture(with_covariate)
        expect_true(input$fit$available)
        beta <- input$fit$theta + 0.015 * seq_along(input$fit$theta)
        step <- 1e-4
        effect_fn <- function(value) {
            SACAT:::.sacat_abundance_mark_effect(value, input$X_b, input$gh)
        }
        nll_fn <- function(value) {
            -sum(SACAT:::zt_beta_loglik_by_sample_inference(
                value, input$y, input$N, input$X, input$gh
            ))
        }
        difference <- function(fn) {
            vapply(seq_along(beta), function(j) {
                plus <- minus <- beta
                plus[j] <- plus[j] + step
                minus[j] <- minus[j] - step
                (fn(plus) - fn(minus)) / (2 * step)
            }, numeric(1))
        }
        effect <- SACAT:::.sacat_abundance_effect_derivatives(
            beta, input$X_b, input$gh
        )
        derivatives <- SACAT:::.sacat_abundance_profile_derivatives(
            beta, input$y, input$N, input$X, input$gh
        )
        expect_equal(effect$value, effect_fn(beta), tolerance = 1e-12)
        expect_equal(effect$gradient, difference(effect_fn), tolerance = 1e-7)
        expect_equal(
            effect$hessian, unname(optimHess(beta, effect_fn)),
            tolerance = 1e-5
        )
        expect_equal(derivatives$nll, nll_fn(beta), tolerance = 1e-11)
        expect_equal(
            unname(derivatives$score), -difference(nll_fn), tolerance = 1e-5
        )
        expect_equal(
            derivatives$information, unname(optimHess(beta, nll_fn)),
            tolerance = 1e-4
        )
        positive <- input$y > 0
        positive_only <- SACAT:::.sacat_abundance_profile_derivatives(
            beta, input$y[positive], input$N[positive],
            input$X[positive, , drop = FALSE], input$gh
        )
        expect_equal(derivatives, positive_only, tolerance = 1e-13)
    }
})

test_that("profile decoding delays the full effect Hessian until the solution", {
    input <- .abundance_profile_fixture(TRUE)
    beta <- input$fit$theta + 0.015 * seq_along(input$fit$theta)
    full <- SACAT:::.sacat_abundance_effect_derivatives(
        beta, input$X_b, input$gh
    )
    first_order <- SACAT:::.sacat_abundance_effect_derivatives(
        beta, input$X_b, input$gh, need_hessian = FALSE
    )
    expect_identical(first_order, full[c("value", "gradient")])

    derivative_function <- SACAT:::.sacat_abundance_effect_derivatives
    full_calls <- first_order_calls <- 0L
    root_calls <- 0L
    X_b <- cbind(Intercept = 1, z = seq(-1, 1, length.out = 20L))
    decode <- function(lambda, target) {
        with_mocked_bindings(
            SACAT:::.sacat_abundance_profile_decode(
                lambda, target, X_b, input$gh
            ),
            .sacat_abundance_effect_derivatives = function(
                    beta, X_b, gh, need_hessian = TRUE) {
                if (need_hessian) {
                    full_calls <<- full_calls + 1L
                } else {
                    first_order_calls <<- first_order_calls + 1L
                }
                derivative_function(beta, X_b, gh, need_hessian)
            },
            uniroot = function(...) {
                root_calls <<- root_calls + 1L
                stats::uniroot(...)
            },
            .package = "SACAT"
        )
    }

    newton <- decode(c(-3, 0.2, log(0.8)), 0.3)
    root_fallback <- decode(c(5, 0, log(0.01)), -1)
    expect_identical(
        newton$effect,
        derivative_function(newton$beta, X_b, input$gh)
    )
    expect_identical(
        root_fallback$effect,
        derivative_function(root_fallback$beta, X_b, input$gh)
    )
    expect_identical(full_calls, 2L)
    expect_identical(root_calls, 1L)
    expect_gt(first_order_calls, 2L)
})

test_that("the fitted contrast is a zero-loss profile constraint", {
    for (with_covariate in c(FALSE, TRUE)) {
        input <- .abundance_profile_fixture(with_covariate)
        fit <- input$fit
        tested <- SACAT:::.sacat_abundance_profile_test(
            fit, input$y, input$N, input$X_b, fit$raw_delta,
            fit$raw_se^2, input$gh
        )
        expect_lte(tested$fit$constraint_error, 1e-8)
        expect_lte(tested$fit$score_residue, 1e-6)
        expect_lte(tested$likelihood_ratio, 1e-7)
        expect_equal(tested$likelihood_ratio, 0)
        expect_equal(tested$statistic, 0)
        expect_equal(tested$z, 0)
        expect_equal(tested$p, 1)
        expect_equal(tested$fit$theta, unname(fit$theta), tolerance = 1e-5)
    }
})

test_that("profile tests use the constrained contrast and signed likelihood root", {
    input <- .abundance_profile_fixture(TRUE)
    fit <- input$fit
    for (target in c(0, fit$raw_delta - 0.2, fit$raw_delta + 0.2)) {
        variance <- 1.4 * fit$raw_se^2
        tested <- SACAT:::.sacat_abundance_profile_test(
            fit, input$y, input$N, input$X_b, target, variance, input$gh
        )
        expected_curvature <- sum(fit$solver_diagnostics$effect_gradient *
            solve(fit$solver_diagnostics$information,
                  fit$solver_diagnostics$effect_gradient))
        expect_equal(tested$fit$effect, target, tolerance = 1e-8)
        expect_lte(tested$fit$score_residue, 1e-6)
        expect_gt(tested$likelihood_ratio, 0)
        expect_equal(tested$curvature_variance, expected_curvature)
        expect_equal(
            tested$statistic,
            tested$likelihood_ratio * expected_curvature / variance
        )
        expect_identical(sign(tested$z), sign(fit$raw_delta - target))
        expect_equal(tested$z^2, tested$statistic)
        expect_equal(tested$p, 2 * pnorm(-abs(tested$z)), tolerance = 1e-14)
        if (target == 0) expect_equal(tested$fit$theta[2L], 0)
    }
})

test_that("successful profile fits retain the analytic fitting path", {
    local_mocked_bindings(
        .sacat_abundance_profile_numeric_fit = function(...) {
            stop("numerical fallback was called")
        },
        .package = "SACAT"
    )
    for (with_covariate in c(FALSE, TRUE)) {
        input <- .abundance_profile_fixture(with_covariate)
        for (target in c(0, input$fit$raw_delta)) {
            restricted <- SACAT:::.sacat_abundance_profile_fit(
                input$fit, input$y, input$N, input$X_b, target, input$gh
            )
            expect_lte(restricted$score_residue, 1e-6)
            expect_lte(restricted$constraint_error, 1e-8)
        }
    }
})

test_that("profile fitting converges for dispersed positive counts", {
    y <- c(1L, 1L, 1L, 1L, 3L, 1L, 1L, 1L, 22L, 10L, 2L, 22L,
           18L, 16L, 0L, 7L, 165L, 470L, 34L, 370L, 418L, 305L, 0L,
           294L, 131L, 224L, 28L, 473L, 249L, 218L, 0L, 208L)
    group <- c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0,
               1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1)
    z <- cbind(rep(c(1, 0), 16), rep(c(0, 0, 1, 1), 8),
               rep(c(rep(1, 4), rep(0, 4)), 4),
               rep(c(rep(0, 8), rep(1, 8)), 2))
    depth <- rep(1000L, length(y))
    control <- SACAT:::.sacat_abundance_control()
    gh <- SACAT:::.sacat_make_abundance_gh_rule(control$quadrature_Q)
    design <- SACAT:::.sacat_abundance_designs(group, z)
    fit <- SACAT:::.sacat_abundance_fit_taxon(
        y, depth, group, z, gh, gh, control
    )
    expect_true(fit$available)
    numeric_fit <- SACAT:::.sacat_abundance_profile_numeric_fit
    cached_numeric <- numeric_fit(
        fit, y, depth, design$X_b, -0.45, gh
    )
    uncached_numeric <- with_mocked_bindings(
        numeric_fit(fit, y, depth, design$X_b, -0.45, gh),
        .sacat_exact_memoise = function(fn, capacity = 1L) fn,
        .package = "SACAT"
    )
    expect_identical(cached_numeric, uncached_numeric)

    fallback_calls <- 0L
    local_mocked_bindings(
        .sacat_abundance_profile_numeric_fit = function(...) {
            fallback_calls <<- fallback_calls + 1L
            numeric_fit(...)
        },
        .package = "SACAT"
    )
    nuisance <- setdiff(seq_along(fit$theta), 2L)
    for (target in c(0, -0.45)) {
        restricted <- SACAT:::.sacat_abundance_profile_fit(
            fit, y, depth, design$X_b, target, gh
        )
        expect_lte(restricted$score_residue, 1e-6)
        expect_lte(restricted$constraint_error, 1e-8)
        expect_equal(restricted$effect, target, tolerance = 1e-8)
        expect_gt(min(eigen(restricted$information, symmetric = TRUE)$values), 0)
        if (target == 0) expect_equal(restricted$theta[2L], 0)

        objective <- function(lambda) {
            beta <- SACAT:::.sacat_abundance_profile_decode(
                lambda, target, design$X_b, gh
            )$beta
            -sum(SACAT:::zt_beta_loglik_by_sample_inference(
                beta, y, depth, design$X_eta, gh
            ))
        }
        expect_equal(restricted$nll, objective(restricted$theta[nuisance]))
        for (step in c(1e-4, 1e-5)) {
            gradient <- SACAT:::zt_central_derivative_matrix_fixed(
                objective, restricted$theta[nuisance],
                rep(step, length(nuisance))
            )
            expect_lte(max(abs(gradient)) / sum(y > 0), 1e-6)
        }
    }
    expect_gt(fallback_calls, 0L)
})

test_that("profile inference preserves corrected effects and standard errors", {
    input <- .abundance_profile_family()
    corrected <- input$corrected
    tested <- input$tested
    expect_true(all(corrected$formed))
    expect_true(all(tested$formed))
    expect_identical(tested$estimate, corrected$estimate)
    expect_identical(tested$se, corrected$se)
    expect_identical(
        tested$diagnostics$taxon$reference_taxa,
        corrected$diagnostics$taxon$reference_taxa
    )
    expect_equal(tested$p, 2 * pnorm(-abs(tested$z)), tolerance = 1e-14)
    expect_equal(tested$diagnostics$taxon$corrected_p_value, tested$p)
    expect_equal(tested$diagnostics$taxon$profile_statistic, tested$z^2)
    expect_named(tested$diagnostics$profile_fits, colnames(input$Y))
    expect_null(tested$background)

    compact <- corrected
    compact$diagnostics <- NULL
    compact <- SACAT:::.sacat_abundance_profile_tests(
        compact, input$fits, input$Y, input$N, input$X_b, input$gh
    )
    expect_null(compact$diagnostics)
    expect_identical(compact$p, tested$p)
    expect_identical(compact$z, tested$z)
    expect_identical(compact$estimate, tested$estimate)
    expect_identical(compact$se, tested$se)
})

test_that("compact profile tasks preserve serial and parallel results", {
    input <- .abundance_profile_family()
    index <- which(input$corrected$formed)
    full_tests <- lapply(index, function(j) {
        SACAT:::.sacat_abundance_profile_test(
            input$fits[[j]], as.numeric(input$Y[, j]), input$N, input$X_b,
            input$corrected$background[j], input$corrected$se[j]^2, input$gh
        )
    })
    tasks <- lapply(index, function(j) {
        fit <- input$fits[[j]]
        list(
            fit = list(
                theta = fit$theta,
                raw_delta = fit$raw_delta,
                solver_diagnostics = fit$solver_diagnostics[c(
                    "X_eta", "information", "effect_gradient"
                )]
            ),
            y = as.numeric(input$Y[, j]),
            background = input$corrected$background[j],
            variance = input$corrected$se[j]^2
        )
    })
    compact_tests <- lapply(
        tasks, SACAT:::.sacat_abundance_profile_task,
        N = input$N, X_b = input$X_b, gh = input$gh
    )
    expect_identical(compact_tests, full_tests)

    serial <- SACAT:::.sacat_abundance_profile_tests(
        input$corrected, input$fits, input$Y, input$N, input$X_b, input$gh
    )
    expect_identical(serial, input$tested)

    cluster <- tryCatch(
        parallel::makePSOCKcluster(2L),
        error = function(e) NULL
    )
    if (is.null(cluster)) skip("PSOCK workers are unavailable")
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    paths <- .libPaths()
    parallel::clusterCall(cluster, function(paths) {
        .libPaths(paths)
        loadNamespace("SACAT")
        TRUE
    }, paths = paths)
    parallel_result <- SACAT:::.sacat_abundance_profile_tests(
        input$corrected, input$fits, input$Y, input$N, input$X_b, input$gh,
        cluster = cluster
    )
    expect_identical(parallel_result, serial)
})

test_that("a failed profile preserves effects and other targets' references", {
    input <- .abundance_profile_family()
    fits <- input$fits
    fits[[1L]]$solver_diagnostics$information[,] <- 0
    tested <- SACAT:::.sacat_abundance_profile_tests(
        input$corrected, fits, input$Y, input$N, input$X_b, input$gh
    )
    expect_false(tested$formed[1L])
    expect_match(tested$reason[1L], "abundance_profile_failed")
    expect_equal(tested$p[1L], 1)
    expect_true(is.na(tested$z[1L]))
    expect_identical(tested$estimate, input$corrected$estimate)
    expect_identical(tested$se, input$corrected$se)
    expect_identical(
        tested$diagnostics$taxon$reference_taxa,
        input$corrected$diagnostics$taxon$reference_taxa
    )
    expect_identical(tested$p[-1L], input$tested$p[-1L])
    expect_identical(tested$z[-1L], input$tested$z[-1L])
    expect_identical(tested$formed[-1L], input$tested$formed[-1L])
    expect_null(tested$diagnostics$profile_fits[[1L]])
})
