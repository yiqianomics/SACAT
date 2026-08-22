test_that("exact C0 helper classifies both sides and equality", {
    y <- c(1, 0)

    negative <- DASRA:::zt_intercept_alpha_boundary(
        rep(log(1 / 3), length(y)), y
    )
    equality <- DASRA:::zt_intercept_alpha_boundary(
        rep(log(1 / 2), length(y)), y
    )

    epsilon <- 1e-10
    odds <- 1 + epsilon
    tiny_positive_probability <- odds / (1 + odds)
    tiny_positive <- DASRA:::zt_intercept_alpha_boundary(
        rep(log(tiny_positive_probability), length(y)), y
    )

    expect_lt(negative$C0, 0)
    expect_true(negative$at_zero)
    expect_equal(equality$C0, 0, tolerance = 16 * .Machine$double.eps)
    expect_true(equality$at_zero)
    expect_gt(tiny_positive$C0, 0)
    expect_lt(tiny_positive$C0, 1e-8)
    expect_false(tiny_positive$at_zero)
})

test_that("exact C0 structural boundary returns its full contract", {
    y <- c(1, 1, 0, 0, 1, 1, 0, 0)
    N <- rep(100L, length(y))
    g <- rep(c(0, 1), each = length(y) / 2)

    run_case <- function(probability) {
        with_mocked_bindings(
            DASRA:::zt_count_structural_test(
                y = y, N = N, g = g, Q = 3L,
                min_positive_samples = 3L, keep_fit = TRUE
            ),
            make_structural_gh_rule = function(Q) list(Q = Q),
            zt_fit_beta = function(...) {
                list(
                    ok = TRUE,
                    reason = "ok",
                    par = 0,
                    numerical_warnings = character()
                )
            },
            zt_beta_detection_components = function(...) {
                list(log_r = rep(log(probability), length(y)))
            },
            .package = "DASRA"
        )
    }

    negative <- run_case(0.1)
    equality <- run_case(0.5)

    for (result in list(negative, equality)) {
        expect_equal(result$p, 1)
        expect_true(result$tested)
        expect_false(result$regular)
        expect_identical(
            result$reason, "structural_absence_boundary_at_zero"
        )
        expect_true(result$diagnostics$exact_boundary)
        expect_equal(result$gamma, numeric(length(y)))
        expect_equal(result$rho, numeric(length(y)))
        expect_identical(
            result$diagnostics$fit$structural_absence$method,
            "exact_C0_boundary"
        )
    }
    expect_lt(negative$diagnostics$intercept_alpha_C0, 0)
    expect_equal(
        equality$diagnostics$intercept_alpha_C0,
        0,
        tolerance = 16 * .Machine$double.eps
    )
})

test_that("positive C0 proceeds to the finite structural fit", {
    y <- c(1, 1, 0, 0, 1, 0, 0, 0)
    N <- rep(100L, length(y))
    g <- rep(c(0, 1), each = length(y) / 2)
    finite_path_reached <- FALSE

    result <- with_mocked_bindings(
        DASRA:::zt_count_structural_test(
            y = y, N = N, g = g, Q = 3L,
            min_positive_samples = 3L
        ),
        make_structural_gh_rule = function(Q) list(Q = Q),
        zt_fit_beta = function(...) {
            list(
                ok = TRUE,
                reason = "ok",
                par = 0,
                numerical_warnings = character()
            )
        },
        zt_beta_detection_components = function(...) {
            list(log_r = rep(log(0.5), length(y)))
        },
        zt_fit_alpha = function(...) {
            finite_path_reached <<- TRUE
            list(ok = FALSE, reason = "finite_alpha_path_reached")
        },
        .package = "DASRA"
    )

    expect_true(finite_path_reached)
    expect_identical(result$reason, "finite_alpha_path_reached")
    expect_false(identical(
        result$reason, "structural_absence_boundary_at_zero"
    ))
})

test_that("Cauchy omnibus uses only regular components for retained taxa", {
    scenario <- data.frame(
        structural_p = c(0.2, 1.0, 1.0, 0.4),
        structural_formed = c(TRUE, TRUE, TRUE, TRUE),
        structural_regular = c(TRUE, FALSE, FALSE, TRUE),
        abundance_p = c(1.0, 0.3, 1.0, 0.6),
        abundance_formed = c(FALSE, TRUE, FALSE, TRUE),
        row.names = c(
            "structural_only", "abundance_only", "none", "both"
        )
    )

    local_mocked_bindings(
        .dasra_structural_arm = function(
                Y, N, g, z, keep_diagnostics,
                conditional_present_starts = 1L) {
            selected <- scenario[colnames(Y), , drop = FALSE]
            list(
                p = selected$structural_p,
                formed = selected$structural_formed,
                regular = selected$structural_regular,
                reason = rep("controlled_test_state", ncol(Y)),
                score_z = rep(NA_real_, ncol(Y)),
                warning = rep("", ncol(Y)),
                diagnostics = NULL
            )
        },
        .dasra_abundance_arm = function(Y, N, g, z, keep_diagnostics) {
            selected <- scenario[colnames(Y), , drop = FALSE]
            list(
                p = selected$abundance_p,
                formed = selected$abundance_formed,
                reason = rep("controlled_test_state", ncol(Y)),
                estimate = rep(NA_real_, ncol(Y)),
                se = rep(NA_real_, ncol(Y)),
                z = rep(NA_real_, ncol(Y)),
                warning = rep("", ncol(Y)),
                diagnostics = NULL
            )
        },
        .package = "DASRA"
    )

    samples <- paste0("sample_", seq_len(6L))
    counts <- rbind(
        structural_only = rep(1, 6L),
        abundance_only = rep(1, 6L),
        none = rep(1, 6L),
        both = rep(1, 6L),
        not_retained = c(1, 1, 0, 0, 0, 0)
    )
    colnames(counts) <- samples
    metadata <- data.frame(
        group = factor(rep(c("reference", "comparison"), each = 3L)),
        reads = rep(100L, 6L),
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

    expect_true(all(
        fit$diagnostics[rownames(scenario), "retained"]
    ))
    expect_identical(
        unname(fit$results["structural_only", "p_omnibus_cauchy"]),
        0.2
    )
    expect_identical(
        unname(fit$results["abundance_only", "p_omnibus_cauchy"]),
        0.3
    )
    expect_identical(
        unname(fit$results["none", "p_omnibus_cauchy"]),
        1.0
    )
    expect_equal(
        unname(fit$results["both", "p_omnibus_cauchy"]),
        DASRA:::cauchy_combination(c(0.4, 0.6)),
        tolerance = 1e-15
    )
    expect_identical(
        unname(fit$results[rownames(scenario), "components_used_cauchy"]),
        c(
            "structural_absence", "relative_abundance", "none", "both"
        )
    )

    expect_false(fit$diagnostics["not_retained", "retained"])
    expect_true(is.na(
        fit$results["not_retained", "p_omnibus_cauchy"]
    ))
    expect_identical(
        unname(fit$results["not_retained", "components_used_cauchy"]),
        "not_retained"
    )
})
