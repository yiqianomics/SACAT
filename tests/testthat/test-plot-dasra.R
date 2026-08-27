.make_plot_contract_fixture <- function(n = 6L) {
    stopifnot(n >= 6L)
    feature <- c(
        "Faecalibacterium_prausnitzii",
        "Bacteroides_uniformis",
        "Roseburia_intestinalis",
        "Akkermansia_muciniphila",
        "Blautia",
        "Unicode_æ_feature",
        paste0("Taxon_", seq_len(max(0L, n - 6L)))
    )
    feature <- feature[seq_len(n)]
    p_omnibus <- c(2e-5, 0.008, 0.03, 0.12, 0.4, 0.8,
                   rep(0.9, max(0L, n - 6L)))
    p_adjusted <- c(0.0008, 0.009, 0.04, 0.2, 0.6, 0.95,
                    rep(1, max(0L, n - 6L)))
    p_structural <- c(4.7e-7, 0.004, 0.03, 0.2, 0.7, 1,
                      rep(1, max(0L, n - 6L)))
    p_abundance <- c(0.02, 0.0004, 0.08, 0.5, 0.9, 1,
                     rep(1, max(0L, n - 6L)))
    z_structural <- c(-5.1, 2.5, -1.9, 0.8, -0.4, 0,
                      rep(0, max(0L, n - 6L)))
    z_abundance <- c(2.3, -3.8, 1.2, -0.7, 0.2, 0,
                     rep(0, max(0L, n - 6L)))

    results <- data.frame(
        taxon = feature,
        p_structural_absence = p_structural,
        p_adj_structural_absence = p_structural,
        z_structural_absence = z_structural,
        p_relative_abundance = p_abundance,
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
            call = quote(dasra()),
            plot_data = list(
                version = 1L,
                contrast = c(reference = "healthy", comparison = "disease"),
                profiles = profiles
            )
        ),
        class = "dasra"
    )
}

.plot_contract_palette <- c(
    "#6090c1", "#acd2e5", "#fef9b7",
    "#fee395", "#f2724d", "#d7312d"
)

test_that("plot selection and significance annotations are deterministic", {
    fit <- .make_plot_contract_fixture()
    spec <- DASRA:::.dasra_plot_build_spec(
        fit, NULL, "significant", 24L, 0.05, "adaptive",
        .plot_contract_palette, c("#6090c1", "#f28e4b")
    )

    expect_identical(
        spec$data$feature,
        fit$results$taxon[c(1L, 2L, 3L)]
    )
    expect_identical(spec$data$stars, c("***", "**", "*"))
    expect_identical(spec$contrast, c("healthy", "disease"))

    explicit <- DASRA:::.dasra_plot_build_spec(
        fit, rev(fit$results$taxon[c(2L, 5L)]), "significant",
        1L, 0.001, "adaptive", .plot_contract_palette,
        c("#6090c1", "#f28e4b")
    )
    expect_identical(
        explicit$data$feature,
        rev(fit$results$taxon[c(2L, 5L)])
    )

    top <- DASRA:::.dasra_plot_build_spec(
        fit, NULL, "top", 2L, 0.05, "adaptive",
        .plot_contract_palette, c("#6090c1", "#f28e4b")
    )
    expect_identical(top$data$feature, fit$results$taxon[1:2])
})

test_that("adaptive colors use one shared displayed-component scale", {
    fit <- .make_plot_contract_fixture()
    spec <- DASRA:::.dasra_plot_build_spec(
        fit, fit$results$taxon[[1L]], "significant", 24L, 0.05,
        "adaptive", .plot_contract_palette,
        c("#6090c1", "#f28e4b")
    )
    expect_identical(
        unname(spec$color_limits), c(1e-7, 1)
    )

    manual <- DASRA:::.dasra_plot_build_spec(
        fit, fit$results$taxon[[1L]], "significant", 24L, 0.05,
        c(1e-8, 0.5), .plot_contract_palette,
        c("#6090c1", "#f28e4b")
    )
    expect_identical(unname(manual$color_limits), c(1e-8, 0.5))
    expect_error(
        DASRA:::.dasra_plot_build_spec(
            fit, NULL, "top", 2L, 0.05, c(1, 1e-4),
            .plot_contract_palette, c("#6090c1", "#f28e4b")
        ),
        "0 < lower < upper <= 1"
    )

    all_one <- fit
    all_one$results$p_adj_structural_absence[] <- 1
    all_one$results$p_adj_relative_abundance[] <- 1
    one_data <- data.frame(
        p_structural = all_one$results$p_adj_structural_absence,
        z_structural = all_one$results$z_structural_absence,
        p_abundance = all_one$results$p_adj_relative_abundance,
        z_abundance = all_one$results$z_relative_abundance
    )
    one_limits <- DASRA:::.dasra_plot_color_limits(
        one_data, "adaptive"
    )
    expect_identical(unname(one_limits), c(0.1, 1))
    one_colors <- DASRA:::.dasra_plot_palette_function(
        .plot_contract_palette, one_limits
    )(rep(1, nrow(all_one$results)))
    expect_true(length(unique(one_colors)) == 1L)

    subnormal <- one_data[1L, , drop = FALSE]
    subnormal$p_structural <- 5e-324
    subnormal$p_abundance <- 5e-324
    tiny_limits <- DASRA:::.dasra_plot_color_limits(
        subnormal, "adaptive"
    )
    expect_identical(tiny_limits[["lower"]], .Machine$double.xmin)
})

test_that("default evidence colors are sequential and distinct from groups", {
    method <- getS3method("plot", "dasra")
    evidence <- eval(formals(method)$evidence_palette)
    groups <- eval(formals(method)$group_colors)
    lightness <- grDevices::convertColor(
        t(grDevices::col2rgb(evidence)) / 255,
        from = "sRGB", to = "Lab"
    )[, "L"]

    expect_true(all(diff(lightness) < 0))
    expect_length(intersect(tolower(evidence), tolower(groups)), 0L)
})

test_that("abundance log-axis ticks are regular and carry percent units", {
    scale <- DASRA:::.dasra_plot_abundance_scale(10 ^ c(-3.8, 0.8))
    expect_identical(scale$limits, c(-4, 1))
    expect_identical(scale$ticks, c(-4, -2, 0))
    expect_true(length(unique(diff(scale$ticks))) == 1L)
    expect_lte(length(scale$ticks), 4L)
    expect_true(all(grepl("%", scale$labels, fixed = TRUE)))

    same_decade <- DASRA:::.dasra_plot_abundance_scale(c(0.12, 0.18))
    expect_identical(same_decade$limits, c(-1, 0))
    expect_identical(same_decade$ticks, c(-1, 0))
})

test_that("plot dispatch is side-effect free and writes vector output", {
    fit <- .make_plot_contract_fixture()
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
    fit <- .make_plot_contract_fixture()
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

    legacy_path <- tempfile(fileext = ".pdf")
    method <- getS3method("plot", "dasra")
    method(
        fit, NULL, "top", 2L, 0.05, "adaptive",
        eval(formals(method)$evidence_palette),
        eval(formals(method)$group_colors),
        legacy_path, 7.2, 3.2, 300
    )
    expect_true(file.exists(legacy_path))
})

test_that("plot reports incomplete objects and unavailable rows clearly", {
    fit <- .make_plot_contract_fixture()
    missing_payload <- fit
    missing_payload$plot_data <- NULL
    expect_error(plot(missing_payload), "store_plot_data = TRUE")

    unsupported <- fit
    unsupported$plot_data$version <- 2L
    expect_error(plot(unsupported), "version is not supported")

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
    set.seed(20260825)
    n <- 72L
    p <- 6L
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
            1L,
            depth[[i]],
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
    fit_args <- list(
        counts = counts,
        metadata = metadata,
        formula = ~ group,
        group = "group",
        library_size = "reads",
        component = "all",
        full_output = TRUE,
        structural_quadrature_points = 31L,
        abundance_quadrature_points = 31L,
        workers = 1L
    )

    plain <- do.call(dasra, c(fit_args, list(store_plot_data = FALSE)))
    prepared <- do.call(dasra, c(fit_args, list(store_plot_data = TRUE)))

    expect_identical(prepared$results, plain$results)
    expect_identical(prepared$diagnostics, plain$diagnostics)
    expect_identical(prepared$fits, plain$fits)
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
})
