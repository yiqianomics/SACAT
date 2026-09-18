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
