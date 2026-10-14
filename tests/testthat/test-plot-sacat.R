.make_plot_fixture <- function(n = 6L) {
    stopifnot(n >= 6L)
    feature <- c(
        "Faecalibacterium_prausnitzii",
        "Bacteroides_uniformis",
        "Roseburia_intestinalis",
        "Akkermansia_muciniphila",
        "Blautia",
        "Prevotella_copri",
        paste0("Taxon_", seq_len(max(0L, n - 6L)))
    )
    feature <- feature[seq_len(n)]
    p_omnibus <- c(2e-5, 0.008, 0.03, 0.12, 0.4, 0.8,
                   rep(0.9, max(0L, n - 6L)))
    p_adjusted <- c(0.0008, 0.009, 0.04, 0.2, 0.6, 0.95,
                    rep(1, max(0L, n - 6L)))
    z_structural <- c(-5.1, 2.5, -1.9, 0.8, -0.4, 0,
                      rep(0, max(0L, n - 6L)))
    z_abundance <- c(2.3, -3.8, 1.2, -0.7, 0.2, 0,
                     rep(0, max(0L, n - 6L)))
    p_structural_raw <- 2 * stats::pnorm(-abs(z_structural))
    p_abundance_raw <- 2 * stats::pnorm(-abs(z_abundance))
    p_structural <- stats::p.adjust(p_structural_raw, method = "BH")
    p_abundance <- stats::p.adjust(p_abundance_raw, method = "BH")

    results <- data.frame(
        taxon = feature,
        p_structural_absence = p_structural_raw,
        p_adj_structural_absence = p_structural,
        z_structural_absence = z_structural,
        p_relative_abundance = p_abundance_raw,
        p_adj_relative_abundance = p_abundance,
        z_relative_abundance = z_abundance,
        p_omnibus = p_omnibus,
        p_adj_omnibus = p_adjusted,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    diagnostics <- data.frame(
        taxon = feature,
        retained = rep(TRUE, n),
        formed_omnibus = rep(TRUE, n),
        stringsAsFactors = FALSE
    )
    profiles <- data.frame(
        feature = feature,
        structural_reference = seq(0.08, 0.28, length.out = n),
        structural_comparison = seq(0.18, 0.12, length.out = n),
        structural_status = rep("ok", n),
        abundance_reference = 10 ^ seq(-1.3, -0.4, length.out = n),
        abundance_comparison = 10 ^ seq(-1.1, -0.5, length.out = n),
        abundance_status = rep("ok", n),
        stringsAsFactors = FALSE
    )
    structure(
        list(
            results = results,
            diagnostics = diagnostics,
            settings = list(
                component = "all",
                contrast = c(reference = "healthy", comparison = "disease"),
                p_adjust_method = "BH"
            ),
            call = quote(sacat()),
            plot_data = list(
                version = 2L,
                contrast = c(reference = "healthy", comparison = "disease"),
                profiles = profiles,
                zero_sources = data.frame(
                    feature = rep(feature, each = 2),
                    group = rep(0:1, n),
                    detected = rep(c(0.65, 0.60), n),
                    estimated_sa = as.vector(t(cbind(
                        profiles$structural_reference,
                        profiles$structural_comparison
                    ))),
                    estimated_nd = 1 - rep(c(0.65, 0.60), n) -
                        as.vector(t(cbind(
                            profiles$structural_reference,
                            profiles$structural_comparison
                        ))),
                    n = 20L, status = "ok"
                ),
                positive_abundance = do.call(rbind, lapply(seq_len(n), function(i) {
                    data.frame(
                        feature = feature[[i]], group = rep(0:1, each = 10),
                        abundance = rep(seq(0.3, 2, length.out = 10), 2) *
                            rep(c(profiles$abundance_reference[[i]],
                                  profiles$abundance_comparison[[i]]), each = 10)
                    )
                }))
            )
        ),
        class = "sacat"
    )
}

.make_plot_analysis_fixture <- function(
        n = 72L, p = 6L, seed = 20260825L) {
    stopifnot(n %% 2L == 0L, p >= 2L)
    set.seed(seed)
    group <- rep(c(0L, 1L), each = n / 2L)
    depth <- rep(10000L, n)
    shift <- c(0.25, -0.20, rep(0, p - 2L))
    latent <- sapply(seq_len(p), function(j) {
        -5.7 + (j - 1L) * 0.12 + shift[[j]] * group +
            stats::rnorm(n, sd = 0.25)
    })
    present <- matrix(stats::runif(n * p) > 0.15, nrow = n)
    probability <- present * stats::plogis(latent)
    count_by_sample <- t(vapply(seq_len(n), function(i) {
        draw <- stats::rmultinom(
            1L, depth[[i]],
            c(probability[i, ], 1 - sum(probability[i, ]))
        )
        draw[seq_len(p), 1L]
    }, numeric(p)))
    counts <- t(count_by_sample)
    rownames(counts) <- paste0("Taxon_", seq_len(p))
    colnames(counts) <- paste0("Sample_", seq_len(n))
    metadata <- data.frame(
        group = factor(
            group, levels = c(0L, 1L),
            labels = c("reference", "comparison")
        ),
        reads = depth,
        row.names = colnames(counts)
    )
    list(counts = counts, metadata = metadata)
}

.plot_group_colors <- c("#4C78A8", "#E45756")
.plot_zero_colors <- c("#6A9DC0", "#E5A36C", "#E9E9E9")

test_that("plot selection and significance annotations are deterministic", {
    fit <- .make_plot_fixture()
    spec <- SACAT:::.sacat_plot_build_spec(
        fit, NULL, "significant", 24L, 0.05, .plot_group_colors, .plot_zero_colors
    )

    expect_identical(
        spec$data$feature,
        fit$results$taxon[c(1L, 2L, 3L)]
    )
    expect_identical(spec$data$stars, c("***", "**", "*"))
    expect_identical(spec$contrast, c("healthy", "disease"))
    expect_equal(
        spec$component_guides$structural$z,
        stats::qnorm(
            0.05 * 2 / (2 * nrow(fit$results)), lower.tail = FALSE
        )
    )
    expect_equal(
        spec$component_guides$abundance$z,
        stats::qnorm(
            0.05 / (2 * nrow(fit$results)), lower.tail = FALSE
        )
    )

    explicit <- SACAT:::.sacat_plot_build_spec(
        fit, rev(fit$results$taxon[c(2L, 5L)]), "significant",
        1L, 0.001, .plot_group_colors, .plot_zero_colors
    )
    expect_identical(
        explicit$data$feature,
        rev(fit$results$taxon[c(2L, 5L)])
    )

    top <- SACAT:::.sacat_plot_build_spec(
        fit, NULL, "top", 2L, 0.05, .plot_group_colors, .plot_zero_colors
    )
    expect_identical(top$data$feature, fit$results$taxon[1:2])
})

test_that("component BH gates use the complete fitted families", {
    fit <- .make_plot_fixture()
    full <- SACAT:::.sacat_plot_build_spec(
        fit, NULL, "all", 24L, 0.05, .plot_group_colors, .plot_zero_colors
    )
    subset <- SACAT:::.sacat_plot_build_spec(
        fit, fit$results$taxon[[1L]], "top", 1L, 0.05,
        .plot_group_colors, .plot_zero_colors
    )

    expect_identical(
        subset$component_guides, full$component_guides
    )
    expect_identical(full$component_guides$structural$discoveries, 2L)
    expect_identical(full$component_guides$abundance$discoveries, 1L)
    expect_identical(full$component_guides$structural$family_size, 6L)
    expect_identical(full$component_guides$abundance$family_size, 6L)

    no_discoveries <- fit
    no_discoveries$results$z_structural_absence[] <- 0
    no_discoveries$results$z_relative_abundance[] <- 0
    no_discoveries$results$p_structural_absence[] <- 1
    no_discoveries$results$p_relative_abundance[] <- 1
    no_discoveries$results$p_adj_structural_absence[] <- 1
    no_discoveries$results$p_adj_relative_abundance[] <- 1
    no_guide <- SACAT:::.sacat_plot_build_spec(
        no_discoveries, NULL, "top", 2L, 0.05, .plot_group_colors, .plot_zero_colors
    )
    expect_identical(
        no_guide$component_guides$structural$discoveries, 0L
    )
    expect_identical(
        no_guide$component_guides$abundance$discoveries, 0L
    )
    expect_true(is.na(no_guide$component_guides$structural$z))
    expect_true(is.na(no_guide$component_guides$abundance$z))

    disabled <- SACAT:::.sacat_plot_build_spec(
        fit, NULL, "top", 2L, 0.05, .plot_group_colors, .plot_zero_colors, FALSE
    )
    expect_false(disabled$component_guides$enabled)
    expect_error(
        SACAT:::.sacat_plot_build_spec(
            fit, NULL, "top", 2L, 0.05, .plot_group_colors, .plot_zero_colors, NA
        ),
        "TRUE or FALSE"
    )

    holm <- fit
    holm$settings$p_adjust_method <- "holm"
    holm$results$p_adj_structural_absence <- stats::p.adjust(
        holm$results$p_structural_absence, method = "holm"
    )
    holm$results$p_adj_relative_abundance <- stats::p.adjust(
        holm$results$p_relative_abundance, method = "holm"
    )
    holm_spec <- SACAT:::.sacat_plot_build_spec(
        holm, NULL, "top", 2L, 0.05, .plot_group_colors, .plot_zero_colors
    )
    expect_false(holm_spec$component_guides$enabled)

    malformed <- fit
    malformed$results$p_structural_absence[[1L]] <- 0.2
    expect_warning(
        malformed_spec <- SACAT:::.sacat_plot_build_spec(
            malformed, NULL, "top", 2L, 0.05, .plot_group_colors, .plot_zero_colors
        ),
        "structural absence"
    )
    expect_true(is.na(malformed_spec$component_guides$structural$z))
    expect_equal(
        malformed_spec$component_guides$abundance$z,
        full$component_guides$abundance$z
    )

    malformed_retained <- fit
    malformed_retained$diagnostics$retained[[1L]] <- NA
    expect_error(
        SACAT:::.sacat_plot_build_spec(
            malformed_retained, NULL, "top", 2L, 0.05, .plot_group_colors, .plot_zero_colors
        ),
        "retained-taxon indicator"
    )
})

test_that("component BH gates tolerate floating-point boundary rounding", {
    raw_p <- c(rep(0.04166666666666667, 5L), 1)
    adjusted_p <- stats::p.adjust(raw_p, method = "BH")
    z <- stats::qnorm(raw_p / 2, lower.tail = FALSE)

    guide <- expect_no_error(
        SACAT:::.sacat_plot_component_bh_guide(
            rep(TRUE, 6L), raw_p, adjusted_p, z, 0.05
        )
    )
    expect_identical(guide$discoveries, 5L)
    expect_equal(guide$p, 0.05 * 5 / 6)
    expect_equal(
        guide$z,
        stats::qnorm((0.05 * 5 / 6) / 2, lower.tail = FALSE)
    )
})

test_that("unverifiable component gates are omitted with one warning", {
    fit <- .make_plot_fixture()
    fit$results$p_structural_absence[[1L]] <- 0
    fit$results$p_adj_structural_absence <- stats::p.adjust(
        fit$results$p_structural_absence, method = "BH"
    )
    fit$results$z_structural_absence[[1L]] <- Inf

    expect_warning(
        spec <- SACAT:::.sacat_plot_build_spec(
            fit, NULL, "top", 2L,
            .Machine$double.xmin * .Machine$double.eps,
            .plot_group_colors, .plot_zero_colors
        ),
        "structural absence"
    )
    expect_true(is.na(spec$component_guides$structural$z))
    expect_true(is.na(spec$component_guides$abundance$z))

    expect_no_warning(
        disabled <- SACAT:::.sacat_plot_build_spec(
            fit, NULL, "top", 2L, 0.05, .plot_group_colors, .plot_zero_colors, FALSE
        )
    )
    expect_false(disabled$component_guides$enabled)
})

test_that("plot colors identify groups and estimated zero sources", {
    fit <- .make_plot_fixture()
    expect_error(
        SACAT:::.sacat_plot_build_spec(
            fit, NULL, "top", 2L, 0.05, "blue", .plot_zero_colors
        ), "two|2 valid colors"
    )
    expect_error(
        SACAT:::.sacat_plot_build_spec(
            fit, NULL, "top", 2L, 0.05, .plot_group_colors, c("blue", "red")
        ), "3 valid colors"
    )
    method <- getS3method("plot", "sacat")
    expect_length(eval(formals(method)$group_colors), 2L)
    expect_length(eval(formals(method)$zero_colors), 3L)
})

test_that("observed abundance boxes are summarized on the plotted log scale", {
    values <- c(0.01, 0.02, 0.03, 0.04, 0.5, 1, 10)
    z <- log10(values)
    q <- quantile(z, c(0.25, 0.5, 0.75), names = FALSE)
    inside <- z[z >= q[1] - 1.5 * diff(q[c(1, 3)]) &
                z <= q[3] + 1.5 * diff(q[c(1, 3)])]
    expect_equal(SACAT:::.sacat_plot_box_summary(values),
                 c(min(inside), q, max(inside)))
    expect_equal(SACAT:::.sacat_plot_box_summary(0.1), rep(-1, 5))
    expect_true(all(is.na(SACAT:::.sacat_plot_box_summary(numeric()))))
})

test_that("abundance log-axis ticks label percentages on a regular log grid", {
    scale <- SACAT:::.sacat_plot_abundance_scale(10 ^ c(-3.8, 0.8), max_ticks = 4L)
    expect_identical(scale$limits, c(-4, 1))
    expect_identical(scale$ticks, c(-4, -2, 0))
    expect_true(length(unique(diff(scale$ticks))) == 1L)
    expect_lte(length(scale$ticks), 4L)
    expect_equal(as.numeric(scale$labels), 10 ^ scale$ticks)

    same_decade <- SACAT:::.sacat_plot_abundance_scale(c(0.12, 0.18))
    expect_identical(same_decade$limits, c(-1, 0))
    expect_identical(same_decade$ticks, c(-1, 0))
})

test_that("plot dispatch is side-effect free and writes vector output", {
    fit <- .make_plot_fixture()
    before <- serialize(fit, NULL)
    set.seed(9917)
    rng_before <- .Random.seed
    path <- tempfile(fileext = ".pdf")

    spec <- plot(
        fit, selection = "top", max_features = 4L,
        file = path, width = 7.2, height = 3.8
    )

    expect_type(spec, "list")
    expect_true(file.exists(path))
    expect_gt(file.info(path)$size, 1000)
    expect_identical(serialize(fit, NULL), before)
    expect_identical(.Random.seed, rng_before)
})

test_that("display group labels do not change the fitted contrast", {
    fit <- .make_plot_fixture()
    original_contrast <- fit$settings$contrast
    path <- tempfile(fileext = ".pdf")

    spec <- plot(
        fit,
        selection = "top",
        max_features = 2L,
        group_labels = c("Control", "Case"),
        file = path,
        width = 7.2,
        height = 3.2
    )

    expect_identical(spec$contrast, c("Control", "Case"))
    expect_identical(fit$settings$contrast, original_contrast)
    expect_error(
        plot(fit, group_labels = "Control"),
        "two non-empty labels"
    )

    positional_path <- tempfile(fileext = ".pdf")
    method <- getS3method("plot", "sacat")
    method(
        fit, NULL, "top", 2L, 0.05,
        eval(formals(method)$group_colors),
        eval(formals(method)$zero_colors),
        positional_path, 7.2, 3.2, 300
    )
    expect_true(file.exists(positional_path))
})

test_that("plot reports incomplete objects and unavailable rows clearly", {
    fit <- .make_plot_fixture()
    missing_payload <- fit
    missing_payload$plot_data <- NULL
    expect_error(plot(missing_payload), "store_plot_data = TRUE")

    unsupported <- fit
    unsupported$plot_data$version <- 3L
    expect_error(plot(unsupported), "required plotting data")

    malformed <- fit
    malformed$plot_data$profiles$abundance_status <- NULL
    expect_error(plot(malformed), "summaries are malformed")

    wrong_contrast <- fit
    wrong_contrast$plot_data$contrast[[2L]] <- "another group"
    expect_error(plot(wrong_contrast), "fitted contrast")

    single_arm <- fit
    single_arm$settings$component <- "relative_abundance"
    expect_error(plot(single_arm), "component = \"all\"")

    unavailable <- fit
    unavailable$diagnostics$formed_omnibus[[1L]] <- FALSE
    expect_error(
        plot(unavailable, features = unavailable$results$taxon[[1L]]),
        "without a formed omnibus result"
    )
    expect_error(
        plot(fit, features = "not_a_feature"),
        "Unknown feature"
    )
})

test_that("plot-summary preparation leaves fitted inference unchanged", {
    input <- .make_plot_analysis_fixture()
    fit_args <- list(
        counts = input$counts,
        metadata = input$metadata,
        formula = ~ group,
        group = "group",
        library_size = "reads",
        component = "all",
        full_output = TRUE,
        structural_quadrature_points = 31L,
        abundance_quadrature_points = 31L,
        workers = 1L
    )

    plain <- do.call(sacat, c(fit_args, list(store_plot_data = FALSE)))
    prepared <- do.call(sacat, c(fit_args, list(store_plot_data = TRUE)))

    expect_identical(prepared$results, plain$results)
    expect_identical(prepared$diagnostics, plain$diagnostics)
    expect_identical(prepared$fits, plain$fits)
    expect_true(any(vapply(
        plain$fits$structural_absence,
        function(fit) isTRUE(fit$diagnostics$quadrature_checked),
        logical(1)
    )))
    expect_true(any(vapply(
        plain$fits$relative_abundance$raw_fits,
        function(fit) isTRUE(fit$quadrature_checked),
        logical(1)
    )))
    prepared_settings <- prepared$settings
    prepared_settings$plot_data_stored <- NULL
    expect_identical(prepared_settings, plain$settings)

    detail <- prepared$fits$relative_abundance$taxon
    profile <- prepared$plot_data$profiles
    formed <- which(
        detail$formed & profile$abundance_status == "ok"
    )
    expect_gt(length(formed), 0L)
    profile_contrast <- log(
        profile$abundance_comparison[formed] /
            profile$abundance_reference[formed]
    )
    expect_equal(
        unname(profile_contrast),
        unname(detail$raw_estimate[formed]),
        tolerance = 1e-10
    )
    expect_identical(
        unname(prepared$results$estimate_relative_abundance[formed]),
        unname(detail$corrected_estimate[formed])
    )
    expect_true(any(abs(
        profile_contrast - detail$corrected_estimate[formed]
    ) > 1e-8))

    compact_args <- fit_args
    compact_args$full_output <- FALSE
    compact_plain <- do.call(
        sacat, c(compact_args, list(store_plot_data = FALSE))
    )
    quadrature_calls <- 0L
    compact_prepared <- with_mocked_bindings(
        do.call(sacat, c(
            compact_args, list(store_plot_data = TRUE)
        )),
        .sacat_higher_order_quadrature_points = function(...) {
            quadrature_calls <<- quadrature_calls + 1L
            stop("unexpected high-order quadrature check")
        },
        .package = "SACAT"
    )
    expect_identical(quadrature_calls, 0L)
    expect_identical(compact_prepared$results, compact_plain$results)
    expect_identical(
        compact_prepared$diagnostics, compact_plain$diagnostics
    )
    expect_false("fits" %in% names(compact_prepared))
    expect_true("plot_data" %in% names(compact_prepared))
    compact_settings <- compact_prepared$settings
    compact_settings$plot_data_stored <- NULL
    expect_identical(compact_settings, compact_plain$settings)
})

test_that("structural plot summaries reuse the conditional-present state", {
    X_eta <- cbind(
        `(Intercept)` = 1,
        Group = c(0, 0, 1, 1)
    )
    fit <- list(
        tested = TRUE,
        regular = TRUE,
        diagnostics = list(fit = list(
            conditional_present = list(par = c(0, 0, 0)),
            X_eta = X_eta,
            quadrature = list(nodes = 0, weights = 1)
        ))
    )
    conditional_present <- list(log_r = rep(log(0.8), 4L))
    component_calls <- 0L
    received_component <- NULL

    profile <- with_mocked_bindings(
        SACAT:::.sacat_structural_plot_profile(
            y = c(0, 1, 0, 1),
            N = rep(10, 4L),
            fit = fit
        ),
        zt_beta_detection_components = function(...) {
            component_calls <<- component_calls + 1L
            conditional_present
        },
        zt_fit_alpha = function(..., conditional_present = NULL) {
            received_component <<- conditional_present
            list(
                ok = TRUE,
                par = c(stats::qlogis(0.2), log(2)),
                nll = function(alpha) 1
            )
        },
        zt_structural_zero_limit_comparison = function(...) {
            list(dominated = FALSE)
        },
        .package = "SACAT"
    )

    expect_identical(component_calls, 1L)
    expect_identical(received_component, conditional_present)
    expect_equal(profile$reference, 0.2)
    expect_equal(profile$comparison, 1 / 3)
    expect_identical(profile$status, "ok")
})

test_that("structural plot summaries preserve companion failure reasons", {
    X_eta <- cbind(`(Intercept)` = 1, Group = c(0, 1))
    fit <- list(
        tested = TRUE,
        regular = TRUE,
        diagnostics = list(fit = list(
            conditional_present = list(par = c(0, 0, 0)),
            X_eta = X_eta,
            quadrature = list(nodes = 0, weights = 1)
        ))
    )
    alpha_calls <- 0L
    unavailable <- with_mocked_bindings(
        SACAT:::.sacat_structural_plot_profile(
            y = c(0, 1), N = c(10, 10), fit = fit
        ),
        zt_beta_detection_components = function(...) NULL,
        zt_fit_alpha = function(...) {
            alpha_calls <<- alpha_calls + 1L
            stop("should not be called")
        },
        .package = "SACAT"
    )
    expect_identical(alpha_calls, 0L)
    expect_identical(
        unavailable$status,
        "structural_absence_no_finite_fit"
    )

    optimization_failure <- expect_no_warning(with_mocked_bindings(
        SACAT:::.sacat_structural_plot_profile(
            y = c(0, 1), N = c(10, 10), fit = fit
        ),
        zt_beta_detection_components = function(...) {
            list(log_r = rep(log(0.8), 2L))
        },
        zt_fit_alpha = function(...) {
            warning("deliberate optimizer warning")
            list(ok = FALSE, reason = "deliberate_fit_failure")
        },
        .package = "SACAT"
    ))
    expect_identical(
        optimization_failure$status,
        "deliberate_fit_failure"
    )

    failed <- with_mocked_bindings(
        SACAT:::.sacat_structural_plot_profile(
            y = c(0, 1), N = c(10, 10), fit = fit
        ),
        zt_beta_detection_components = function(...) {
            stop("deliberate detection failure")
        },
        .package = "SACAT"
    )
    expect_identical(failed$status, "companion_fit_unavailable")
})

test_that("profile output works across supported graphics devices", {
    fit <- .make_plot_fixture()
    for (extension in c(".pdf", ".png")) {
        path <- tempfile(fileext = extension)
        expect_no_warning(
            spec <- plot(
                fit, selection = "top", max_features = 3L,
                file = path, width = 7.2, height = 3.2, dpi = 150
            )
        )
        expect_type(spec, "list")
        expect_true(file.exists(path))
        expect_gt(file.info(path)$size, 1000)
    }

    svg_path <- tempfile(fileext = ".svg")
    if (isTRUE(capabilities("cairo"))) {
        expect_no_warning(
            plot(
                fit, selection = "top", max_features = 3L,
                file = svg_path, width = 7.2, height = 3.2
            )
        )
        expect_true(file.exists(svg_path))
        expect_gt(file.info(svg_path)$size, 1000)
    } else {
        expect_error(
            plot(fit, file = svg_path),
            "Cairo support"
        )
    }
})

test_that("two PSOCK workers preserve inference and plotting summaries", {
    package_path <- find.package("SACAT")
    skip_if_not(
        file.exists(file.path(package_path, "Meta", "package.rds")),
        "PSOCK regression test requires an installed package"
    )
    probe <- tryCatch(
        parallel::makePSOCKcluster(2L),
        error = function(e) NULL
    )
    skip_if(is.null(probe), "Two local PSOCK workers are unavailable")
    parallel::stopCluster(probe)

    original_library_paths <- .libPaths()
    on.exit(.libPaths(original_library_paths), add = TRUE)
    package_library <- normalizePath(
        dirname(package_path), winslash = "/", mustWork = TRUE
    )
    retained_library_paths <- original_library_paths[
        normalizePath(
            original_library_paths, winslash = "/", mustWork = FALSE
        ) != package_library
    ]
    skip_if(
        !length(retained_library_paths),
        "PSOCK regression test requires another package library"
    )
    .libPaths(retained_library_paths)

    input <- .make_plot_analysis_fixture(
        n = 48L, p = 5L, seed = 20260828L
    )
    fit_args <- list(
        counts = input$counts,
        metadata = input$metadata,
        formula = ~ group,
        group = "group",
        library_size = "reads",
        component = "all",
        full_output = FALSE,
        store_plot_data = TRUE,
        structural_quadrature_points = 31L,
        abundance_quadrature_points = 31L
    )

    set.seed(4802)
    rng_before <- .Random.seed
    serial <- suppressWarnings(do.call(
        sacat, c(fit_args, list(workers = 1L))
    ))
    expect_identical(.Random.seed, rng_before)
    parallel_fit <- suppressWarnings(do.call(
        sacat, c(fit_args, list(workers = 2L))
    ))
    expect_identical(.Random.seed, rng_before)

    expect_identical(parallel_fit$results, serial$results)
    expect_identical(parallel_fit$diagnostics, serial$diagnostics)
    expect_identical(parallel_fit$plot_data, serial$plot_data)
    serial$settings$workers_requested <- NULL
    serial$settings$workers_used <- NULL
    parallel_fit$settings$workers_requested <- NULL
    parallel_fit$settings$workers_used <- NULL
    expect_identical(parallel_fit$settings, serial$settings)
})
