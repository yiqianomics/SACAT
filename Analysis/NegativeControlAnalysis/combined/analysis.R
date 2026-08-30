#!/usr/bin/env Rscript

# Summarize negative-control results across all public datasets.

options(stringsAsFactors = FALSE, warn = 1)

scenario_order <- c("balanced", "fourfold")
scenario_labels <- c(
    balanced = "Balanced depth",
    fourfold = "Fourfold depth difference"
)
method_order <- c("DASRA", "ZINQ", "MaAsLin3")
required_dasra_version <- "0.6.0"
dataset_order <- c(
    "NogueraJulianHIV",
    "BaxterE_2016",
    "HMPV35Throat",
    "Qiita13631ASD",
    "ArtPrize_2015_Forehead",
    "MehtaRS_2018",
    "VatanenT_2016",
    "LeChatelierE_2013",
    "ORIGINS_2022_Healthy_Plaque",
    "NHANESOral_2011_2012",
    "NielsenHB_2014",
    "NHANESOral_2009_2010",
    "JieZ_2017",
    "Qiita11993Colombia",
    "LiJ_2014",
    "Atlas1006",
    "SchirmerM_2016",
    "ZeeviD_2015",
    "VilaAV_2018",
    "QinJ_2012"
)
dataset_count <- length(dataset_order)
alpha <- 0.05
replicates <- 100L
taxa_per_replicate <- 30L

# Locate the combined output directory from the executed script.
script_directory <- function() {
    arguments <- commandArgs(trailingOnly = FALSE)
    file_argument <- grep("^--file=", arguments, value = TRUE)
    if (!length(file_argument)) return(normalizePath(getwd()))
    normalizePath(dirname(sub("^--file=", "", file_argument[[1L]])))
}

# Write a complete table without exposing a partially written file.
atomic_write_csv <- function(table, path) {
    temporary <- tempfile(pattern = basename(path), tmpdir = dirname(path))
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    utils::write.csv(table, temporary, row.names = FALSE, na = "")
    if (!file.rename(temporary, path)) {
        if (!file.copy(temporary, path, overwrite = TRUE)) {
            stop("Could not write output: ", path, call. = FALSE)
        }
        unlink(temporary)
    }
    invisible(path)
}

# Read the ordered dataset metadata from the stored analysis inputs.
read_dataset_catalog <- function(root) {
    rows <- lapply(seq_along(dataset_order), function(index) {
        dataset <- dataset_order[[index]]
        path <- file.path(root, dataset, "analysis_input.rds")
        if (!file.exists(path)) {
            stop("Missing analysis input for ", dataset, ".",
                 call. = FALSE)
        }
        input <- readRDS(path)
        if (!is.list(input) ||
            !identical(as.character(input$dataset_id), dataset) ||
            length(input$label) != 1L || !nzchar(input$label) ||
            length(input$dataset_index) != 1L ||
            !is.finite(input$dataset_index)) {
            stop("The analysis input metadata is malformed for ", dataset,
                 ".", call. = FALSE)
        }
        data.frame(
            study_order = index,
            dataset_index = as.integer(input$dataset_index),
            dataset = dataset,
            label = as.character(input$label),
            stringsAsFactors = FALSE
        )
    })
    catalog <- do.call(rbind, rows)
    if (nrow(catalog) != dataset_count || anyDuplicated(catalog$dataset) ||
        anyDuplicated(catalog$label) || anyDuplicated(catalog$dataset_index)) {
        stop("The analysis inputs must have unique dataset metadata.",
             call. = FALSE)
    }
    rownames(catalog) <- NULL
    catalog
}

# Read one complete taxon-level result table.
read_dataset_result <- function(dataset, label, root) {
    path <- file.path(root, dataset, "results", "taxon_pvalues.csv")
    if (!file.exists(path)) {
        stop("Missing result table for ", dataset, ".", call. = FALSE)
    }
    result <- utils::read.csv(path, check.names = FALSE)
    expected_names <- c(
        "dataset", "dataset_label", "replicate", "scenario", "method",
        "taxon", "p_value", "native_p_value", "available", "status",
        "reason"
    )
    if (!identical(names(result), expected_names) ||
        nrow(result) != replicates * length(scenario_order) *
            length(method_order) * taxa_per_replicate ||
        !identical(unique(result$dataset), dataset) ||
        !identical(unique(result$dataset_label), label)) {
        stop("The result table is malformed for ", dataset, ".",
             call. = FALSE)
    }
    result
}

# Check the complete result grid and the treatment of unavailable tests.
validate_results <- function(result, catalog) {
    expected_rows <- nrow(catalog) * replicates * length(scenario_order) *
        length(method_order) * taxa_per_replicate
    if (nrow(result) != expected_rows ||
        !setequal(unique(result$dataset), catalog$dataset) ||
        !setequal(unique(result$scenario), scenario_order) ||
        !setequal(unique(result$method), method_order)) {
        stop("The result grid is incomplete.", call. = FALSE)
    }
    replicate_value <- suppressWarnings(as.numeric(result$replicate))
    if (any(!is.finite(replicate_value)) ||
        any(replicate_value != as.integer(replicate_value)) ||
        !identical(sort(unique(as.integer(replicate_value))),
                   seq_len(replicates))) {
        stop("Replicate identifiers must be 1 through 100.", call. = FALSE)
    }
    if (anyNA(result$taxon) || any(!nzchar(result$taxon)) ||
        any(!is.finite(result$p_value)) ||
        any(result$p_value < 0 | result$p_value > 1) ||
        !is.logical(result$available) || anyNA(result$available)) {
        stop("Taxon identifiers, p-values, or availability are invalid.",
             call. = FALSE)
    }
    expected_status <- ifelse(result$available, "available", "unavailable")
    if (anyNA(result$status) || any(result$status != expected_status) ||
        anyNA(result$reason) || any(!nzchar(result$reason)) ||
        any(result$p_value[!result$available] != 1) ||
        any(!is.finite(result$native_p_value[result$available])) ||
        any(abs(result$p_value[result$available] -
                result$native_p_value[result$available]) > 1e-12)) {
        stop("Result status and p-value fields are inconsistent.",
             call. = FALSE)
    }
    full_key <- paste(
        result$dataset, result$replicate, result$scenario, result$method,
        result$taxon, sep = "\r"
    )
    cell_key <- paste(
        result$dataset, result$replicate, result$scenario, result$method,
        sep = "\r"
    )
    cell_counts <- table(cell_key)
    dataset_taxa <- vapply(
        split(result$taxon, result$dataset),
        function(value) length(unique(value)), integer(1)
    )
    if (anyDuplicated(full_key) ||
        length(cell_counts) != nrow(catalog) * replicates *
            length(scenario_order) * length(method_order) ||
        any(cell_counts != taxa_per_replicate) ||
        any(dataset_taxa != taxa_per_replicate)) {
        stop("The result grid contains duplicate or incomplete cells.",
             call. = FALSE)
    }
    invisible(TRUE)
}

# Check the stored input, fixed panel, and positive-count summaries.
validate_dataset_outputs <- function(dataset, result, catalog_entry, root) {
    input_path <- file.path(root, dataset, "analysis_input.rds")
    diagnostics_path <- file.path(root, dataset, "results", "diagnostics.csv")
    settings_path <- file.path(
        root, dataset, "results", "analysis_settings.csv"
    )
    if (!all(file.exists(c(
        input_path, diagnostics_path, settings_path
    )))) {
        stop("Required analysis files are missing for ", dataset, ".",
             call. = FALSE)
    }

    input <- readRDS(input_path)
    tolerance <- 1e-12
    if (!is.list(input) || !is.matrix(input$probability) ||
        !identical(as.character(input$dataset_id), dataset) ||
        !identical(as.character(input$label), catalog_entry$label) ||
        !identical(as.integer(input$dataset_index),
                   as.integer(catalog_entry$dataset_index)) ||
        nrow(input$probability) != 200L ||
        length(input$evaluation_taxa) != taxa_per_replicate ||
        anyDuplicated(input$evaluation_taxa) ||
        !setequal(input$evaluation_taxa, unique(result$taxon))) {
        stop("The taxon panel is malformed for ", dataset, ".",
             call. = FALSE)
    }

    probability <- input$probability
    if (any(!is.finite(probability)) || any(probability < 0) ||
        any(abs(rowSums(probability) - 1) > 1e-10) ||
        is.null(rownames(probability)) || is.null(colnames(probability)) ||
        anyDuplicated(rownames(probability)) ||
        anyDuplicated(colnames(probability))) {
        stop("The probability matrix is invalid for ", dataset, ".",
             call. = FALSE)
    }
    candidate_taxa <- setdiff(colnames(probability), "Other_unmodeled")
    source_prevalence <- colMeans(
        probability[, candidate_taxa, drop = FALSE] > 0
    )
    mean_abundance <- colMeans(probability[, candidate_taxa, drop = FALSE])
    eligible <- candidate_taxa[
        source_prevalence >= 0.10 & mean_abundance >= 1e-5
    ]
    eligible <- eligible[order(
        -mean_abundance[eligible], eligible, method = "radix"
    )]
    if (length(eligible) < taxa_per_replicate) {
        stop("Fewer than 30 source-eligible taxa remain for ", dataset, ".",
             call. = FALSE)
    }
    selected <- eligible[seq_len(taxa_per_replicate)]
    if (!identical(selected, as.character(input$evaluation_taxa))) {
        stop("The fixed panel is not source-defined for ", dataset, ".",
             call. = FALSE)
    }
    selection <- input$taxon_selection
    if (!is.list(selection) ||
        !identical(as.numeric(selection$source_prevalence_minimum), 0.10) ||
        !identical(as.numeric(selection$mean_abundance_minimum), 1e-5) ||
        !identical(as.integer(selection$target_taxa), taxa_per_replicate)) {
        stop("The panel definition is incomplete for ", dataset, ".",
             call. = FALSE)
    }

    diagnostics <- utils::read.csv(diagnostics_path, check.names = FALSE)
    required_diagnostics <- c(
        "dataset", "replicate", "scenario", "taxon",
        "positive_samples_H", "positive_samples_Case", "dasra_retained",
        "dasra_formed_structural_absence",
        "dasra_regular_structural_absence",
        "dasra_nonregular_structural_absence",
        "dasra_formed_relative_abundance", "dasra_formed_omnibus",
        "dasra_warning_structural_absence",
        "dasra_warning_relative_abundance"
    )
    if (!all(required_diagnostics %in% names(diagnostics)) ||
        nrow(diagnostics) != replicates * length(scenario_order) *
            taxa_per_replicate ||
        !identical(unique(diagnostics$dataset), dataset) ||
        !setequal(unique(diagnostics$scenario), scenario_order)) {
        stop("Diagnostics are incomplete for ", dataset, ".",
             call. = FALSE)
    }
    dasra_result <- result[result$method == "DASRA", , drop = FALSE]
    diagnostic_key <- paste(
        diagnostics$replicate, diagnostics$scenario, diagnostics$taxon,
        sep = "\r"
    )
    result_key <- paste(
        dasra_result$replicate, dasra_result$scenario, dasra_result$taxon,
        sep = "\r"
    )
    diagnostic_index <- match(result_key, diagnostic_key)
    logical_fields <- c(
        "dasra_retained", "dasra_formed_structural_absence",
        "dasra_regular_structural_absence",
        "dasra_nonregular_structural_absence",
        "dasra_formed_relative_abundance", "dasra_formed_omnibus"
    )
    if (anyDuplicated(diagnostic_key) || anyNA(diagnostic_index) ||
        any(!vapply(
            diagnostics[, logical_fields, drop = FALSE],
            is.logical, logical(1)
        )) || anyNA(as.matrix(diagnostics[, logical_fields, drop = FALSE]))) {
        stop("DASRA diagnostics are invalid for ", dataset, ".",
             call. = FALSE)
    }
    expected_omnibus <- diagnostics$dasra_formed_structural_absence |
        diagnostics$dasra_formed_relative_abundance
    expected_nonregular <- diagnostics$dasra_formed_structural_absence &
        !diagnostics$dasra_regular_structural_absence
    expected_available <-
        diagnostics$dasra_retained[diagnostic_index] &
        diagnostics$dasra_formed_omnibus[diagnostic_index]
    if (any(diagnostics$dasra_formed_omnibus != expected_omnibus) ||
        any(diagnostics$dasra_nonregular_structural_absence !=
                expected_nonregular) ||
        any(dasra_result$available != expected_available)) {
        stop("DASRA availability is inconsistent for ", dataset, ".",
             call. = FALSE)
    }

    fourfold <- diagnostics[
        diagnostics$scenario == "fourfold", , drop = FALSE
    ]
    fourfold$positive_samples_H <- as.integer(fourfold$positive_samples_H)
    fourfold$positive_samples_Case <- as.integer(
        fourfold$positive_samples_Case
    )
    fourfold$minimum_both <- pmin(
        fourfold$positive_samples_H, fourfold$positive_samples_Case
    )
    rows <- lapply(input$evaluation_taxa, function(taxon) {
        subset <- fourfold[fourfold$taxon == taxon, , drop = FALSE]
        if (nrow(subset) != replicates ||
            !identical(sort(as.integer(subset$replicate)),
                       seq_len(replicates))) {
            stop("Fourfold support records are incomplete for ", dataset,
                 ".", call. = FALSE)
        }
        data.frame(
            dataset = dataset,
            taxon = taxon,
            minimum_positive_H = min(subset$positive_samples_H),
            minimum_positive_Case = min(subset$positive_samples_Case),
            minimum_positive_both = min(subset$minimum_both),
            mean_positive_H = mean(subset$positive_samples_H),
            mean_positive_Case = mean(subset$positive_samples_Case),
            mean_positive_both = mean(subset$minimum_both),
            replicates_at_least_20 = sum(subset$minimum_both >= 20L),
            replicates_at_least_18 = sum(subset$minimum_both >= 18L),
            replicates_below_18 = sum(subset$minimum_both < 18L),
            stringsAsFactors = FALSE
        )
    })
    fourfold_summary <- do.call(rbind, rows)
    rownames(fourfold_summary) <- NULL

    settings <- utils::read.csv(settings_path, check.names = FALSE)
    required_settings <- c(
        "dataset", "dataset_label", "n_samples", "group_size", "n_taxa",
        "replicates", "scenarios",
        "fourfold_minimum_positive_H",
        "fourfold_minimum_positive_Case",
        "fourfold_minimum_positive_both",
        "fourfold_fraction_records_at_least_20",
        "fourfold_fraction_records_at_least_18",
        "fourfold_max_taxon_replicates_below_18",
        "DASRA_omnibus_definition", "DASRA_version"
    )
    minimum_positive_both <- fourfold$minimum_both
    below_18_by_taxon <- tapply(
        minimum_positive_both < 18L, fourfold$taxon, sum
    )
    settings_numbers <- suppressWarnings(as.numeric(c(
        settings$n_samples, settings$group_size, settings$n_taxa,
        settings$replicates,
        settings$fourfold_minimum_positive_H,
        settings$fourfold_minimum_positive_Case,
        settings$fourfold_minimum_positive_both,
        settings$fourfold_fraction_records_at_least_20,
        settings$fourfold_fraction_records_at_least_18,
        settings$fourfold_max_taxon_replicates_below_18
    )))
    expected_settings_numbers <- c(
        200, 100, 30, 100,
        min(fourfold$positive_samples_H), min(fourfold$positive_samples_Case),
        min(minimum_positive_both), mean(minimum_positive_both >= 20L),
        mean(minimum_positive_both >= 18L), max(below_18_by_taxon)
    )
    if (!all(required_settings %in% names(settings)) ||
        nrow(settings) != 1L ||
        !identical(as.character(settings$dataset), dataset) ||
        !identical(as.character(settings$dataset_label),
                   unique(result$dataset_label)) ||
        any(!is.finite(settings_numbers)) ||
        any(abs(settings_numbers - expected_settings_numbers) > tolerance) ||
        !identical(as.character(settings$scenarios), "balanced | fourfold") ||
        !identical(as.character(settings$DASRA_omnibus_definition),
                   paste(
                       "Bonferroni omnibus p-value when at least one",
                       "component is formed"
                   )) ||
        !identical(as.character(settings$DASRA_version),
                   required_dasra_version)) {
        stop("The analysis settings are inconsistent for ", dataset, ".",
             call. = FALSE)
    }
    fourfold_summary
}

# Apply catalog order consistently to a dataset-method-scenario summary.
ordered_summary <- function(summary, catalog) {
    summary$dataset_order <- match(summary$dataset, catalog$dataset)
    summary$scenario_order <- match(summary$scenario, scenario_order)
    summary$method_order <- match(summary$method, method_order)
    summary <- summary[order(
        summary$dataset_order, summary$method_order, summary$scenario_order
    ), , drop = FALSE]
    summary$dataset_order <- NULL
    summary$scenario_order <- NULL
    summary$method_order <- NULL
    rownames(summary) <- NULL
    summary
}

# Collapse the fixed 30-taxon family to one result per randomization cell.
make_replicate_summary <- function(result, catalog) {
    cell <- interaction(
        result$dataset, result$replicate, result$scenario, result$method,
        drop = TRUE, lex.order = TRUE
    )
    rows <- lapply(split(seq_len(nrow(result)), cell), function(index) {
        subset <- result[index, , drop = FALSE]
        subset$bh_adjusted_p <- stats::p.adjust(subset$p_value, method = "BH")
        subset$rejected <- subset$bh_adjusted_p <= alpha
        data.frame(
            dataset = subset$dataset[[1L]],
            dataset_label = subset$dataset_label[[1L]],
            replicate = subset$replicate[[1L]],
            scenario = subset$scenario[[1L]],
            method = subset$method[[1L]],
            n_taxa = nrow(subset),
            n_available = sum(subset$available),
            availability_rate = mean(subset$available),
            n_bh_rejections = sum(subset$rejected),
            familywise_type_i_error = as.integer(any(subset$rejected)),
            stringsAsFactors = FALSE
        )
    })
    summary <- do.call(rbind, rows)
    summary$dataset_order <- match(summary$dataset, catalog$dataset)
    summary$scenario_order <- match(summary$scenario, scenario_order)
    summary$method_order <- match(summary$method, method_order)
    summary <- summary[order(
        summary$dataset_order, summary$replicate,
        summary$method_order, summary$scenario_order
    ), , drop = FALSE]
    summary$dataset_order <- NULL
    summary$scenario_order <- NULL
    summary$method_order <- NULL
    rownames(summary) <- NULL
    summary
}

# Summarize family-wise Type I error and Monte Carlo uncertainty by dataset.
make_dataset_summary <- function(replicate_summary, catalog) {
    cell <- interaction(
        replicate_summary$dataset, replicate_summary$scenario,
        replicate_summary$method, drop = TRUE, lex.order = TRUE
    )
    rows <- lapply(split(seq_len(nrow(replicate_summary)), cell),
                   function(index) {
        subset <- replicate_summary[index, , drop = FALSE]
        mean_error <- mean(subset$familywise_type_i_error)
        monte_carlo_se <- stats::sd(subset$familywise_type_i_error) /
            sqrt(nrow(subset))
        margin <- stats::qnorm(0.975) * monte_carlo_se
        data.frame(
            dataset = subset$dataset[[1L]],
            dataset_label = subset$dataset_label[[1L]],
            scenario = subset$scenario[[1L]],
            method = subset$method[[1L]],
            n_replicates = nrow(subset),
            n_taxa_per_replicate = unique(subset$n_taxa),
            empirical_familywise_type_i_error = mean_error,
            monte_carlo_se = monte_carlo_se,
            mc_ci95_lower = max(0, mean_error - margin),
            mc_ci95_upper = min(1, mean_error + margin),
            stringsAsFactors = FALSE
        )
    })
    ordered_summary(do.call(rbind, rows), catalog)
}

# Summarize method availability across all fixed-family tests.
make_availability_summary <- function(result, catalog) {
    cell <- interaction(
        result$dataset, result$scenario, result$method,
        drop = TRUE, lex.order = TRUE
    )
    rows <- lapply(split(seq_len(nrow(result)), cell), function(index) {
        subset <- result[index, , drop = FALSE]
        data.frame(
            dataset = subset$dataset[[1L]],
            dataset_label = subset$dataset_label[[1L]],
            scenario = subset$scenario[[1L]],
            method = subset$method[[1L]],
            n_tests = nrow(subset),
            n_available = sum(subset$available),
            n_unavailable = sum(!subset$available),
            availability_rate = mean(subset$available),
            stringsAsFactors = FALSE
        )
    })
    ordered_summary(do.call(rbind, rows), catalog)
}

# Record whether every dataset-method-scenario cell is complete.
make_completeness_summary <- function(result, catalog) {
    cell <- interaction(
        result$dataset, result$scenario, result$method,
        drop = TRUE, lex.order = TRUE
    )
    rows <- lapply(split(seq_len(nrow(result)), cell), function(index) {
        subset <- result[index, , drop = FALSE]
        records_by_replicate <- table(subset$replicate)
        complete <- nrow(subset) == replicates * taxa_per_replicate &&
            length(records_by_replicate) == replicates &&
            all(records_by_replicate == taxa_per_replicate) &&
            length(unique(subset$taxon)) == taxa_per_replicate
        data.frame(
            dataset = subset$dataset[[1L]],
            dataset_label = subset$dataset_label[[1L]],
            scenario = subset$scenario[[1L]],
            method = subset$method[[1L]],
            n_records = nrow(subset),
            n_replicates = length(records_by_replicate),
            n_taxa = length(unique(subset$taxon)),
            expected_records = replicates * taxa_per_replicate,
            complete = complete,
            stringsAsFactors = FALSE
        )
    })
    ordered_summary(do.call(rbind, rows), catalog)
}

# Draw the dataset-level family-wise Type I error comparison.
draw_figure <- function(summary, catalog, output_directory) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The ggplot2 package is required to draw the figure.",
             call. = FALSE)
    }
    plot_data <- summary
    plot_data$dataset_label <- factor(
        plot_data$dataset_label, levels = rev(catalog$label)
    )
    plot_data$scenario <- factor(
        scenario_labels[plot_data$scenario], levels = unname(scenario_labels)
    )
    plot_data$method <- factor(plot_data$method, levels = method_order)
    dodge <- ggplot2::position_dodge(width = 0.60)
    plot_upper <- max(alpha * 1.25, plot_data$mc_ci95_upper,
                      na.rm = TRUE) * 1.04

    figure <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = empirical_familywise_type_i_error, y = dataset_label,
            color = scenario, shape = scenario
        )
    ) +
        ggplot2::geom_vline(
            xintercept = alpha, linetype = "dashed", linewidth = 0.45,
            color = "#777777"
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(xmin = mc_ci95_lower, xmax = mc_ci95_upper),
            orientation = "y", width = 0.18, linewidth = 0.55,
            position = dodge
        ) +
        ggplot2::geom_point(size = 2.45, stroke = 0.9, position = dodge) +
        ggplot2::facet_grid(
            . ~ method, axes = "all_y", axis.labels = "margins"
        ) +
        ggplot2::scale_color_manual(
            values = c(
                "Balanced depth" = "#667786",
                "Fourfold depth difference" = "#7A4E8E"
            ),
            drop = FALSE
        ) +
        ggplot2::scale_shape_manual(
            values = c("Balanced depth" = 1,
                       "Fourfold depth difference" = 16),
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            breaks = pretty(c(0, plot_upper), n = 4),
            labels = function(value) {
                text <- format(
                    round(value, 3), trim = TRUE, scientific = FALSE
                )
                positive <- !is.na(value) & value > 0
                text[positive] <- sub("^0\\.", ".", text[positive])
                text
            },
            expand = ggplot2::expansion(mult = c(0.04, 0.08))
        ) +
        ggplot2::coord_cartesian(xlim = c(0, plot_upper), clip = "on") +
        ggplot2::labs(
            x = "Family-wise Type I error after BH adjustment",
            y = NULL, color = NULL, shape = NULL
        ) +
        ggplot2::theme_classic(base_size = 10.5, base_family = "sans") +
        ggplot2::theme(
            legend.position = "top",
            legend.justification = "center",
            legend.direction = "horizontal",
            legend.box.spacing = grid::unit(0.15, "lines"),
            legend.spacing.x = grid::unit(0.35, "lines"),
            legend.key.width = grid::unit(1.25, "lines"),
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(face = "bold", size = 11),
            panel.spacing.x = grid::unit(1.15, "lines"),
            axis.title.x = ggplot2::element_text(
                size = 11.5, margin = ggplot2::margin(t = 8)
            ),
            axis.text.x = ggplot2::element_text(
                size = 9.5, color = "#222222"
            ),
            axis.text.y = ggplot2::element_text(
                size = 9.7, color = "#222222"
            ),
            axis.ticks = ggplot2::element_line(
                linewidth = 0.4, color = "#333333"
            ),
            axis.line = ggplot2::element_line(
                linewidth = 0.5, color = "#333333"
            ),
            plot.margin = ggplot2::margin(6, 10, 6, 6)
        )

    ggplot2::ggsave(
        file.path(output_directory, "negative_control_type1_error.pdf"),
        figure, device = grDevices::cairo_pdf,
        width = 10.4, height = 9.8, units = "in", bg = "white"
    )
    invisible(figure)
}

# Validate all inputs, write summaries, and render the comparison figure.
main <- function() {
    output_directory <- script_directory()
    root <- normalizePath(file.path(output_directory, ".."))
    catalog <- read_dataset_catalog(root)
    result <- do.call(rbind, Map(
        read_dataset_result,
        dataset = catalog$dataset,
        label = catalog$label,
        MoreArgs = list(root = root)
    ))
    rownames(result) <- NULL
    validate_results(result, catalog)
    fourfold_support <- do.call(rbind, Map(
        function(dataset) {
            validate_dataset_outputs(
                dataset,
                result[result$dataset == dataset, , drop = FALSE],
                catalog[catalog$dataset == dataset, , drop = FALSE],
                root
            )
        },
        catalog$dataset
    ))
    rownames(fourfold_support) <- NULL

    replicate_summary <- make_replicate_summary(result, catalog)
    dataset_summary <- make_dataset_summary(replicate_summary, catalog)
    availability_summary <- make_availability_summary(result, catalog)
    completeness_summary <- make_completeness_summary(result, catalog)
    expected_dataset_summaries <-
        dataset_count * length(scenario_order) * length(method_order)
    if (nrow(fourfold_support) != dataset_count * taxa_per_replicate ||
        nrow(replicate_summary) != expected_dataset_summaries * replicates ||
        nrow(dataset_summary) != expected_dataset_summaries ||
        nrow(availability_summary) != expected_dataset_summaries ||
        nrow(completeness_summary) != expected_dataset_summaries ||
        any(!completeness_summary$complete)) {
        stop("The combined summaries are incomplete.", call. = FALSE)
    }

    atomic_write_csv(
        fourfold_support,
        file.path(output_directory, "fourfold_positive_count_support.csv")
    )
    atomic_write_csv(
        replicate_summary,
        file.path(output_directory, "replicate_type1_error.csv")
    )
    atomic_write_csv(
        dataset_summary,
        file.path(output_directory, "dataset_type1_error.csv")
    )
    atomic_write_csv(
        availability_summary,
        file.path(output_directory, "availability_summary.csv")
    )
    atomic_write_csv(
        completeness_summary,
        file.path(output_directory, "completeness.csv")
    )
    draw_figure(dataset_summary, catalog, output_directory)
    message("All ", dataset_count, " datasets were validated and summarized.")
}

tryCatch(
    main(),
    error = function(error) {
        message("Error: ", conditionMessage(error))
        quit(save = "no", status = 1L)
    }
)
