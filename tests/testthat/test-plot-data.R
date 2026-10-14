test_that("zero-source posteriors use the group-inclusive plotting fit", {
    group <- c(0, 0, 1, 1)
    design <- cbind(`(Intercept)` = 1, Group = group)
    fit <- list(
        tested = TRUE, regular = TRUE,
        diagnostics = list(fit = list(
            conditional_present = list(par = c(0, 0, 0)),
            X_eta = design, quadrature = list(nodes = 0, weights = 1)
        ))
    )
    y <- c(0, 1, 0, 1)
    detection <- 0.8
    alpha <- c(stats::qlogis(0.2), log(2))
    rho <- stats::plogis(as.numeric(design %*% alpha))
    fit_calls <- 0L
    component_calls <- 0L
    profile <- with_mocked_bindings(
        SACAT:::.sacat_structural_plot_profile(y, rep(10, 4), fit),
        zt_beta_detection_components = function(...) {
            component_calls <<- component_calls + 1L
            list(log_r = rep(log(detection), 4))
        },
        zt_fit_alpha = function(..., X_rho) {
            fit_calls <<- fit_calls + 1L
            expect_identical(X_rho, design)
            list(ok = TRUE, par = alpha, nll = function(alpha) 1)
        },
        zt_structural_zero_limit_comparison = function(...) {
            list(dominated = FALSE)
        },
        .package = "SACAT"
    )
    expected <- ifelse(
        y == 0, rho / (rho + (1 - rho) * (1 - detection)), 0
    )
    expect_equal(profile$posterior_sa, expected, tolerance = 1e-12)
    expect_identical(profile$status, "ok")
    expect_identical(fit_calls, 1L)
    expect_identical(component_calls, 1L)
})

test_that("plotting summaries retain observed groups and zero-source shares", {
    Y <- cbind(
        Taxon_1 = c(0, 2, 4, 0, 0, 6),
        Taxon_2 = c(1, 2, 3, 4, 0, 0)
    )
    depth <- c(100, 200, 400, 100, 200, 300)
    group <- rep(0:1, each = 3)
    posterior <- c(0.7, 0, 0, 0.4, 0.8, 0)
    summaries <- with_mocked_bindings(
        SACAT:::.sacat_prepare_plot_data(
            Y = Y, N = depth, group = group,
            results = data.frame(taxon = colnames(Y)),
            structural_fits = list(Taxon_1 = 1L, Taxon_2 = 2L),
            abundance_fits = NULL,
            contrast = c(reference = "Control", comparison = "Case")
        ),
        .sacat_structural_plot_profile = function(y, N, fit) {
            if (fit == 2L) {
                return(SACAT:::.sacat_plot_profile_failure("unavailable"))
            }
            list(reference = 0.2, comparison = 0.3,
                 posterior_sa = posterior, status = "ok")
        },
        .sacat_abundance_plot_profile = function(...) {
            list(reference = 1, comparison = 2, status = "ok")
        },
        .package = "SACAT"
    )
    expect_identical(summaries$version, 2L)
    expect_identical(summaries$zero_sources$group, rep(0:1, 2))
    expect_equal(summaries$zero_sources$n, rep(3, 4))
    first <- summaries$zero_sources[1:2, ]
    expect_equal(first$detected, c(2 / 3, 1 / 3))
    expect_equal(first$estimated_sa, c(0.7 / 3, 1.2 / 3))
    expect_equal(first$estimated_nd, c(0.3 / 3, 0.8 / 3))
    share_columns <- c("detected", "estimated_sa", "estimated_nd")
    expect_equal(unname(rowSums(first[, share_columns])), c(1, 1))
    expect_equal(
        unlist(summaries$zero_sources[3, share_columns]),
        c(detected = 1, estimated_sa = 0, estimated_nd = 0)
    )
    expect_equal(summaries$zero_sources$detected[[4]], 1 / 3)
    expect_true(is.na(summaries$zero_sources$estimated_sa[[4]]))
    expect_true(is.na(summaries$zero_sources$estimated_nd[[4]]))
    expect_identical(summaries$zero_sources$status[[4]], "unavailable")
    expected_positive <- do.call(rbind, lapply(seq_len(ncol(Y)), function(j) {
        keep <- Y[, j] > 0
        data.frame(feature = colnames(Y)[[j]], group = group[keep],
                   abundance = 100 * Y[keep, j] / depth[keep])
    }))
    expect_identical(summaries$positive_abundance, expected_positive)
})

test_that("observed-data plotting summaries preserve covariate-adjusted inference", {
    set.seed(20261014)
    n <- 80L
    p <- 6L
    group <- rep(0:1, each = n / 2)
    covariate <- stats::rnorm(n)
    depth <- rep(10000L, n)
    probability <- sapply(seq_len(p), function(j) {
        present <- stats::runif(n) > 0.15
        present * stats::plogis(
            -5.8 + 0.1 * j + 0.15 * covariate +
                0.15 * group * (j == 1L) + stats::rnorm(n, sd = 0.3)
        )
    })
    counts <- vapply(seq_len(n), function(i) {
        draw <- stats::rmultinom(
            1, depth[[i]], c(probability[i, ], 1 - sum(probability[i, ]))
        )
        draw[seq_len(p), 1]
    }, numeric(p))
    dimnames(counts) <- list(
        paste0("Taxon_", seq_len(p)), paste0("Sample_", seq_len(n))
    )
    metadata <- data.frame(
        group = factor(group), covariate = covariate,
        reads = depth, row.names = colnames(counts)
    )
    args <- list(
        counts = counts, metadata = metadata, formula = ~ group + covariate,
        group = "group", library_size = "reads", full_output = TRUE,
        structural_quadrature_points = 31L, abundance_quadrature_points = 31L
    )
    plain <- do.call(sacat, c(args, list(store_plot_data = FALSE)))
    prepared <- do.call(sacat, c(args, list(store_plot_data = TRUE)))
    expect_identical(prepared$results, plain$results)
    expect_identical(prepared$diagnostics, plain$diagnostics)
    expect_identical(prepared$fits, plain$fits)
    shares <- prepared$plot_data$zero_sources
    available <- shares$status == "ok"
    expect_true(any(available))
    share_columns <- c("detected", "estimated_sa", "estimated_nd")
    expect_equal(
        unname(rowSums(shares[available, share_columns])),
        rep(1, sum(available)), tolerance = 1e-12
    )
    expect_true(all(shares$estimated_sa[available] >= 0))
    expect_true(all(shares$estimated_nd[available] >= 0))
    expect_equal(
        prepared$plot_data$positive_abundance$abundance,
        unlist(lapply(seq_len(p), function(j) {
            selected <- counts[j, ] > 0
            100 * counts[j, selected] / depth[selected]
        }), use.names = FALSE)
    )
    observed <- prepared$plot_data$profiles
    raw <- prepared$fits$relative_abundance$taxon
    keep <- observed$abundance_status == "ok"
    expect_equal(
        log(observed$abundance_comparison[keep] /
                observed$abundance_reference[keep]),
        unname(raw$raw_estimate[keep]), tolerance = 1e-10
    )
})
