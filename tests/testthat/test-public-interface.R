.capture_public_group_encoding <- function(group, reference = NULL) {
    sample_names <- paste0("Sample_", seq_along(group))
    counts <- rbind(
        Taxon_1 = rep(2L, length(group)),
        Taxon_2 = rep(3L, length(group))
    )
    colnames(counts) <- sample_names
    metadata <- data.frame(
        group = group,
        reads = rep(100L, length(group)),
        row.names = sample_names
    )
    encoded <- NULL

    fit <- with_mocked_bindings(
        dasra(
            counts = counts,
            metadata = metadata,
            formula = ~ group,
            group = "group",
            library_size = "reads",
            reference = reference,
            component = "structural_absence"
        ),
        .dasra_structural_arm = function(
                Y, N, g, z, keep_diagnostics,
                conditional_present_starts = 1L, ...) {
            encoded <<- g
            list(
                p = rep(0.5, ncol(Y)),
                formed = rep(TRUE, ncol(Y)),
                regular = rep(TRUE, ncol(Y)),
                reason = rep("ok", ncol(Y)),
                score_z = rep(0, ncol(Y)),
                warning = rep("", ncol(Y)),
                diagnostics = NULL
            )
        },
        .package = "DASRA"
    )

    list(g = encoded, contrast = fit$settings$contrast)
}

.run_public_interface_fixture <- function(component = "all") {
    scenario <- data.frame(
        structural_p = c(0.01, 0.20, 1.00, 0.04),
        structural_formed = rep(TRUE, 4L),
        structural_regular = c(TRUE, FALSE, TRUE, FALSE),
        abundance_p = c(0.05, 0.10, 0.30, 1.00),
        abundance_formed = c(TRUE, TRUE, FALSE, TRUE),
        row.names = paste0("Taxon_", 1:4)
    )
    structural_warning <- c(
        "expanded_numerical_bounds", "", "", "weak_identification"
    )
    abundance_warning <- c(
        "", "ill_conditioned_information", "", ""
    )
    names(structural_warning) <- names(abundance_warning) <-
        rownames(scenario)

    local_mocked_bindings(
        .dasra_structural_arm = function(
                Y, N, g, z, keep_diagnostics,
                conditional_present_starts = 1L, ...) {
            selected <- scenario[colnames(Y), , drop = FALSE]
            list(
                p = selected$structural_p,
                formed = selected$structural_formed,
                regular = selected$structural_regular,
                reason = rep("controlled_test_state", ncol(Y)),
                score_z = rep(NA_real_, ncol(Y)),
                warning = unname(structural_warning[colnames(Y)]),
                diagnostics = NULL
            )
        },
        .dasra_abundance_arm = function(
                Y, N, g, z, keep_diagnostics, ...) {
            selected <- scenario[colnames(Y), , drop = FALSE]
            list(
                p = selected$abundance_p,
                formed = selected$abundance_formed,
                reason = rep("controlled_test_state", ncol(Y)),
                estimate = rep(NA_real_, ncol(Y)),
                se = rep(NA_real_, ncol(Y)),
                z = rep(NA_real_, ncol(Y)),
                warning = unname(abundance_warning[colnames(Y)]),
                diagnostics = NULL
            )
        },
        .package = "DASRA"
    )

    samples <- paste0("Sample_", seq_len(6L))
    counts <- rbind(
        Taxon_1 = rep(1L, 6L),
        Taxon_2 = rep(1L, 6L),
        Taxon_3 = rep(1L, 6L),
        Taxon_4 = rep(1L, 6L),
        not_retained = c(1L, 1L, 0L, 0L, 0L, 0L)
    )
    colnames(counts) <- samples
    metadata <- data.frame(
        group = factor(rep(c("reference", "comparison"), each = 3L)),
        reads = rep(100L, 6L),
        row.names = samples
    )

    dasra(
        counts = counts,
        metadata = metadata,
        formula = ~ group,
        group = "group",
        library_size = "reads",
        reference = "reference",
        p_adjust_method = "holm",
        component = component
    )
}

.make_warning_abundance_fit <- function() {
    fit <- DASRA:::.dasra_abundance_empty_fit("ok", 6L)
    fit$available <- TRUE
    fit$status <- "ok"
    fit$raw_delta <- 0.1
    fit$raw_se <- 0.2
    fit$raw_p <- 2 * pnorm(-0.5)
    fit$phi <- c(-0.12, -0.04, -0.01, 0.03, 0.05, 0.09)
    fit$numerical_warning <- "expanded_numerical_bounds"
    fit
}

test_that("group encodings preserve the comparison-minus-reference direction", {
    labels <- rep(c("case", "control"), 3L)
    expected <- as.integer(labels == "case")

    factor_result <- .capture_public_group_encoding(
        factor(labels, levels = c("control", "case"))
    )
    character_result <- .capture_public_group_encoding(
        labels, reference = "control"
    )
    numeric_result <- .capture_public_group_encoding(
        ifelse(labels == "case", 2, 4), reference = 4
    )

    expect_identical(factor_result$g, expected)
    expect_identical(character_result$g, expected)
    expect_identical(numeric_result$g, expected)
    expect_identical(
        unname(factor_result$contrast), c("control", "case")
    )
    expect_identical(
        unname(character_result$contrast), c("control", "case")
    )
    expect_identical(unname(numeric_result$contrast), c("4", "2"))
})

test_that("adjusted p-value columns have a method-neutral schema", {
    fit <- .run_public_interface_fixture()
    retained <- fit$diagnostics$retained
    adjusted_columns <- c(
        "p_adj_structural_absence",
        "p_adj_relative_abundance",
        "p_adj_omnibus",
        "p_adj_omnibus_cauchy"
    )

    expect_true(all(adjusted_columns %in% names(fit$results)))
    expect_false(any(startsWith(names(fit$results), "q_")))
    for (column in adjusted_columns) {
        p_column <- sub("^p_adj_", "p_", column)
        expect_equal(
            fit$results[retained, column],
            p.adjust(fit$results[retained, p_column], method = "holm")
        )
        expect_true(all(is.na(fit$results[!retained, column])))
    }
    expect_identical(fit$settings$p_adjust_method, "holm")
})

test_that("lightweight warnings and nonregular status are returned by default", {
    fit <- .run_public_interface_fixture()

    expect_identical(
        fit$diagnostics$warning_structural_absence,
        c(
            "expanded_numerical_bounds", "", "", "weak_identification", ""
        )
    )
    expect_identical(
        fit$diagnostics$warning_relative_abundance,
        c("", "ill_conditioned_information", "", "", "")
    )
    expect_identical(
        fit$diagnostics$nonregular_structural_absence,
        c(FALSE, TRUE, FALSE, TRUE, FALSE)
    )
    expect_false("fits" %in% names(fit))
})

test_that("lightweight diagnostics follow the requested components", {
    structural <- .run_public_interface_fixture("structural_absence")
    abundance <- .run_public_interface_fixture("relative_abundance")

    expect_true(all(c(
        "warning_structural_absence",
        "nonregular_structural_absence"
    ) %in% names(structural$diagnostics)))
    expect_false(
        "warning_relative_abundance" %in% names(structural$diagnostics)
    )
    expect_true(
        "warning_relative_abundance" %in% names(abundance$diagnostics)
    )
    expect_false(any(c(
        "warning_structural_absence",
        "nonregular_structural_absence"
    ) %in% names(abundance$diagnostics)))
})

test_that("component arms retain warning codes without full output", {
    structural <- NULL
    expect_warning(
        structural <- with_mocked_bindings(
            DASRA:::.dasra_structural_arm(
                Y = matrix(
                    c(1, 0, 2, 1, 0, 2),
                    ncol = 1L,
                    dimnames = list(NULL, "Taxon_1")
                ),
                N = rep(100L, 6L),
                g = rep(c(0, 1), each = 3L),
                z = NULL,
                keep_diagnostics = FALSE
            ),
            zt_count_structural_test = function(...) {
                list(
                    p = 0.4,
                    tested = TRUE,
                    regular = TRUE,
                    reason = "ok",
                    diagnostics = list(
                        score_z = 0.2,
                        numerical_warnings = c(
                            "expanded_numerical_bounds",
                            "weak_identification"
                        )
                    )
                )
            },
            .package = "DASRA"
        ),
        "warning_structural_absence"
    )
    expect_identical(
        unname(structural$warning),
        "expanded_numerical_bounds;weak_identification"
    )
    expect_null(structural$diagnostics)

    abundance_fit <- .make_warning_abundance_fit()
    abundance <- NULL
    expect_warning(
        abundance <- with_mocked_bindings(
            DASRA:::.dasra_abundance_arm(
                Y = matrix(
                    c(1, 0, 2, 1, 0, 2),
                    ncol = 1L,
                    dimnames = list(NULL, "Taxon_1")
                ),
                N = rep(100L, 6L),
                g = rep(c(0, 1), each = 3L),
                z = NULL,
                keep_diagnostics = FALSE
            ),
            .dasra_abundance_fit_taxon = function(...) abundance_fit,
            .package = "DASRA"
        ),
        "warning_relative_abundance"
    )
    expect_identical(
        unname(abundance$warning), "expanded_numerical_bounds"
    )
    expect_null(abundance$diagnostics)
})

test_that("print separates regular and conservative structural results", {
    fit <- .run_public_interface_fixture()
    output <- capture.output(print(fit))

    expect_true(any(grepl(
        "Structural results returned: 4/4", output, fixed = TRUE
    )))
    expect_true(any(grepl(
        "Regular structural score tests: 2/4", output, fixed = TRUE
    )))
    expect_true(any(grepl(
        "Conservative nonregular results: 2/4", output, fixed = TRUE
    )))
    expect_true(any(grepl("holm-adjusted p", output, fixed = TRUE)))
})

test_that("public numerical and formation controls have stable defaults", {
    defaults <- formals(DASRA::dasra)

    expect_false("conditional_present_starts" %in% names(defaults))
    expect_identical(
        names(defaults)[11L],
        "structural_conditional_present_starts"
    )
    expect_identical(
        eval(defaults$structural_conditional_present_starts),
        c("adaptive", "full")
    )
    expect_identical(eval(defaults$min_positive_samples), 3L)
    expect_identical(eval(defaults$min_reference_taxa), 4L)
    expect_identical(eval(defaults$structural_quadrature_points), 1001L)
    expect_identical(eval(defaults$abundance_quadrature_points), 41L)
    expect_identical(eval(defaults$workers), 1L)
    expect_identical(eval(defaults$verbose), FALSE)
    expect_identical(eval(defaults$store_plot_data), FALSE)
    expect_identical(getNamespaceExports("DASRA"), "dasra")
})

test_that("plot summaries are an explicit dual-component opt-in", {
    samples <- paste0("Sample_", seq_len(8L))
    counts <- matrix(
        2L,
        nrow = 5L,
        ncol = length(samples),
        dimnames = list(paste0("Taxon_", 1:5), samples)
    )
    metadata <- data.frame(
        group = factor(rep(c("reference", "comparison"), each = 4L)),
        reads = rep(100L, length(samples)),
        row.names = samples
    )

    expect_error(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            component = "structural_absence", store_plot_data = TRUE
        ),
        "requires `component = \"all\"`"
    )
    expect_error(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            store_plot_data = 1L
        ),
        "store_plot_data"
    )
})

test_that("public controls are checked and recorded", {
    samples <- paste0("Sample_", seq_len(8L))
    counts <- matrix(
        2L,
        nrow = 5L,
        ncol = length(samples),
        dimnames = list(paste0("Taxon_", 1:5), samples)
    )
    metadata <- data.frame(
        group = factor(rep(c("reference", "comparison"), each = 4L)),
        reads = rep(100L, length(samples)),
        row.names = samples
    )
    captured <- NULL

    fit <- with_mocked_bindings(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            component = "structural_absence",
            structural_conditional_present_starts = "full",
            min_positive_samples = 2L,
            min_reference_taxa = 3L,
            structural_quadrature_points = 31L
        ),
        .dasra_structural_arm = function(
                Y, N, g, z, keep_diagnostics,
                conditional_present_starts, min_positive_samples,
                quadrature_points, cluster, verbose,
                check_quadrature) {
            captured <<- list(
                starts = conditional_present_starts,
                minimum = min_positive_samples,
                Q = quadrature_points,
                cluster = cluster,
                verbose = verbose,
                check_quadrature = check_quadrature
            )
            list(
                p = rep(0.5, ncol(Y)),
                formed = rep(TRUE, ncol(Y)),
                regular = rep(TRUE, ncol(Y)),
                reason = rep("ok", ncol(Y)),
                score_z = rep(0, ncol(Y)),
                warning = rep("", ncol(Y)),
                diagnostics = NULL
            )
        },
        .package = "DASRA"
    )

    expect_identical(captured$starts, 5L)
    expect_identical(captured$minimum, 2L)
    expect_identical(captured$Q, 31L)
    expect_null(captured$cluster)
    expect_false(captured$verbose)
    expect_false(captured$check_quadrature)
    expect_identical(
        fit$settings$structural_conditional_present_starts, "full"
    )
    expect_false("conditional_present_starts" %in% names(fit$settings))
    expect_identical(fit$settings$min_positive_samples_retained, 2L)
    expect_identical(fit$settings$min_reference_taxa, 3L)
    expect_identical(fit$settings$structural_quadrature_Q, 31L)
    expect_identical(fit$settings$workers_used, 1L)

    expect_error(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            min_positive_samples = 0L
        ),
        "min_positive_samples"
    )
    expect_error(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            min_reference_taxa = 2L
        ),
        "min_reference_taxa"
    )
    expect_error(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            structural_quadrature_points = 2L
        ),
        "structural_quadrature_points"
    )
    expect_error(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            abundance_quadrature_points = 2L
        ),
        "abundance_quadrature_points"
    )
    expect_error(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            conditional_present_starts = "full"
        ),
        "unused argument.*conditional_present_starts"
    )
})

test_that("the abundance quadrature control is passed and recorded", {
    samples <- paste0("Sample_", seq_len(8L))
    counts <- matrix(
        2L,
        nrow = 5L,
        ncol = length(samples),
        dimnames = list(paste0("Taxon_", 1:5), samples)
    )
    metadata <- data.frame(
        group = factor(rep(c("reference", "comparison"), each = 4L)),
        reads = rep(100L, length(samples)),
        row.names = samples
    )
    captured_Q <- NULL

    fit <- with_mocked_bindings(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            component = "relative_abundance",
            abundance_quadrature_points = 31L
        ),
        .dasra_abundance_arm = function(
                Y, N, g, z, keep_diagnostics,
                min_positive_samples, min_reference_taxa,
                quadrature_points, cluster, verbose,
                check_quadrature) {
            captured_Q <<- quadrature_points
            expect_false(check_quadrature)
            list(
                p = rep(0.5, ncol(Y)),
                formed = rep(TRUE, ncol(Y)),
                reason = rep("ok", ncol(Y)),
                estimate = rep(0, ncol(Y)),
                se = rep(1, ncol(Y)),
                z = rep(0, ncol(Y)),
                warning = rep("", ncol(Y)),
                diagnostics = NULL
            )
        },
        .package = "DASRA"
    )

    expect_identical(captured_Q, 31L)
    expect_identical(fit$settings$abundance_quadrature_Q, 31L)
})

test_that("only metadata used by the analysis must be complete", {
    samples <- paste0("Sample_", seq_len(8L))
    counts <- matrix(
        2L,
        nrow = 5L,
        ncol = length(samples),
        dimnames = list(paste0("Taxon_", 1:5), samples)
    )
    metadata <- data.frame(
        group = factor(rep(c("reference", "comparison"), each = 4L)),
        age = seq(20, 55, by = 5),
        unused = c(NA, rep(1, 7L)),
        reads = rep(100L, length(samples)),
        row.names = samples
    )

    expect_no_error(with_mocked_bindings(
        dasra(
            counts, metadata, ~ group, "group", "reads",
            component = "structural_absence"
        ),
        .dasra_structural_arm = function(Y, ...) {
            list(
                p = rep(0.5, ncol(Y)),
                formed = rep(TRUE, ncol(Y)),
                regular = rep(TRUE, ncol(Y)),
                reason = rep("ok", ncol(Y)),
                score_z = rep(0, ncol(Y)),
                warning = rep("", ncol(Y)),
                diagnostics = NULL
            )
        },
        .package = "DASRA"
    ))

    metadata$age[2L] <- NA
    expect_error(
        dasra(
            counts, metadata, ~ group + age, "group", "reads",
            component = "structural_absence"
        ),
        "Could not construct the model"
    )
})
