# Summarize completed real-data analyses without rerunning any method
#
# A dataset is included after its result, input-summary, and method-status
# tables have been checked for completeness. The `significant` column is the
# within-result-set BH decision written by analysis.R.

options(stringsAsFactors = FALSE)

locate_real_data_directory <- function() {
    file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                          value = TRUE)
    if (length(file_argument) == 1L) {
        return(dirname(normalizePath(sub("^--file=", "", file_argument))))
    }
    file.path(getwd(), "Analysis", "RealDataAnalysis")
}

real_data_directory <- locate_real_data_directory()
output_directory <- file.path(real_data_directory, "summary")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

expected_methods <- c(
    "DASRA structural absence",
    "DASRA present-conditional abundance",
    "DASRA combined",
    "MaAsLin3 prevalence",
    "MaAsLin3 abundance",
    "MaAsLin3 combined",
    "ZINQ prevalence",
    "ZINQ abundance",
    "ZINQ combined",
    "ANCOM-BC2",
    "LinDA",
    "corncob",
    "edgeR",
    "DESeq2",
    "metagenomeSeq"
)
expected_families <- c(
    "DASRA", "MaAsLin3", "ZINQ", "ANCOM-BC2", "LinDA", "corncob",
    "edgeR", "DESeq2", "metagenomeSeq"
)
expected_dataset_ids <- c(
    "crc_baxter",
    "cdi_schubert",
    "gems_pediatric_diarrhea",
    "korean_hypertension",
    "microbiomehd_zupancic_obesity",
    "qiita_1939_pediatric_crohn",
    "ravel_vaginal_ethnicity"
)
dataset_directories <- file.path(real_data_directory, expected_dataset_ids)
if (any(!dir.exists(dataset_directories))) {
    stop(sprintf(
        "Expected dataset directories are missing: %s.",
        paste(expected_dataset_ids[!dir.exists(dataset_directories)],
              collapse = ", ")
    ), call. = FALSE)
}

read_completed_dataset <- function(dataset_directory, expected_dataset_id) {
    table_directory <- file.path(dataset_directory, "table")
    results_file <- list.files(
        table_directory,
        pattern = "_method_results_all_taxa[.]csv$",
        full.names = TRUE
    )
    input_summary_file <- list.files(
        table_directory,
        pattern = "_analysis_input_summary[.]csv$",
        full.names = TRUE
    )
    status_file <- list.files(
        table_directory,
        pattern = "_method_status[.]csv$",
        full.names = TRUE
    )
    if (length(results_file) != 1L || length(input_summary_file) != 1L ||
        length(status_file) != 1L) {
        stop(sprintf(
            paste(
                "Expected exactly one result, input-summary, and method-status",
                "table for %s."
            ),
            expected_dataset_id
        ), call. = FALSE)
    }

    results <- utils::read.csv(
        results_file, check.names = FALSE, stringsAsFactors = FALSE
    )
    input_summary <- utils::read.csv(
        input_summary_file, check.names = FALSE, stringsAsFactors = FALSE
    )
    status <- utils::read.csv(
        status_file, check.names = FALSE, stringsAsFactors = FALSE
    )

    required_results_columns <- c(
        "dataset", "taxon", "method", "family", "component",
        "available", "reason", "significant", "components_used"
    )
    required_status_columns <- c("dataset", "family", "status")
    if (!all(required_results_columns %in% names(results)) ||
        !all(c("item", "value") %in% names(input_summary)) ||
        !all(required_status_columns %in% names(status))) {
        stop(sprintf(
            "Completed tables for %s do not follow the real-data output schema.",
            basename(dataset_directory)
        ), call. = FALSE)
    }

    summary_values <- stats::setNames(
        input_summary$value, input_summary$item
    )
    dataset_id <- unname(summary_values[["dataset"]])
    tested_taxa <- suppressWarnings(as.integer(
        unname(summary_values[["tested taxa"]])
    ))
    valid_grid <- length(tested_taxa) == 1L && is.finite(tested_taxa) &&
        tested_taxa > 0L &&
        nrow(results) == tested_taxa * length(expected_methods) &&
        length(unique(results$taxon)) == tested_taxa &&
        !anyDuplicated(results[, c("taxon", "method")]) &&
        setequal(unique(results$method), expected_methods) &&
        all(table(results$method) == tested_taxa)
    valid_status <- nrow(status) == length(expected_families) &&
        !anyDuplicated(status$family) &&
        setequal(status$family, expected_families) &&
        all(status$status %in% c(
            "completed", "completed with unavailable taxa"
        ))
    if (!identical(dataset_id, expected_dataset_id) ||
        !identical(unique(results$dataset), dataset_id) ||
        !identical(unique(status$dataset), dataset_id) ||
        !valid_grid || !valid_status) {
        stop(sprintf(
            "Completed outputs are incomplete or inconsistent for %s.",
            basename(dataset_directory)
        ), call. = FALSE)
    }

    results$available <- as.logical(results$available)
    results$significant <- as.logical(results$significant)
    if (anyNA(results$available) || anyNA(results$significant) ||
        any(results$significant & !results$available) ||
        any(stats::aggregate(
            available ~ method, results, sum
        )$available == 0L)) {
        stop(sprintf(
            "Availability and discovery fields are incomplete for %s.",
            basename(dataset_directory)
        ), call. = FALSE)
    }

    fixed_labels <- c(
        crc_baxter = "Baxter colorectal cancer",
        cdi_schubert = "Schubert C. difficile infection"
    )
    if (dataset_id %in% names(fixed_labels)) {
        dataset_label <- unname(fixed_labels[[dataset_id]])
    } else {
        readme_title <- readLines(
            file.path(dataset_directory, "README.md"), n = 1L, warn = FALSE
        )
        dataset_label <- trimws(sub("^#[[:space:]]*", "", readme_title))
        if (!length(dataset_label) || !nzchar(dataset_label)) {
            dataset_label <- dataset_id
        }
    }

    list(
        dataset = dataset_id,
        dataset_label = dataset_label,
        reference_samples = as.integer(summary_values[["reference samples"]]),
        comparison_samples = as.integer(summary_values[["comparison samples"]]),
        tested_taxa = as.integer(summary_values[["tested taxa"]]),
        results = results
    )
}

completed_datasets <- Map(
    read_completed_dataset, dataset_directories, expected_dataset_ids
)
all_results <- do.call(rbind, lapply(
    completed_datasets, function(dataset_input) dataset_input$results
))
row.names(all_results) <- NULL

humanize_unavailable_reason <- function(reason) {
    descriptions <- c(
        conditional_present_gradient =
            "Conditional-present fit retained a non-negligible gradient",
        conditional_present_information_nonpositive =
            "Conditional-present fit had non-positive information",
        conditional_present_nonconvergence =
            "Conditional-present fit did not converge",
        conditional_present_persistent_boundary =
            "Conditional-present fit remained on a parameter boundary",
        conditional_present_root_polish_no_descent =
            "Conditional-present fit root refinement did not improve the objective",
        conditional_present_root_polish_not_closed =
            "Conditional-present fit root refinement did not close",
        inference_derivative_unstable =
            "Structural-absence inference derivative was unstable",
        no_observed_zeros =
            "No observed zeros for the structural-absence model",
        positive_part_design_rank_deficient =
            "Conditional-present fit design was rank deficient",
        rank_deficient_positive_mark_design =
            "Conditional-present fit design was rank deficient",
        structural_absence_nonoptimal_nuisance_fit =
            "Structural-absence nuisance fit was not optimal",
        structural_absence_persistent_boundary =
            "Structural-absence fit remained on a parameter boundary",
        `one or both components were unavailable` =
            "One or both DASRA components were unavailable"
    )
    translated <- unname(descriptions[reason])
    missing_translation <- is.na(translated)
    translated[missing_translation] <- gsub(
        "_", " ", reason[missing_translation], fixed = TRUE
    )
    translated
}

dasra_methods <- c(
    "DASRA structural absence",
    "DASRA present-conditional abundance",
    "DASRA combined"
)
dasra_result_labels <- c(
    "DASRA structural absence" = "Structural absence component",
    "DASRA present-conditional abundance" =
        "Present-conditional abundance component",
    "DASRA combined" = "Combined omnibus test"
)

availability_rows <- list()
signal_category_rows <- list()
comparison_overlap_rows <- list()
signal_overlap_rows <- list()

signal_categories <- c(
    "Structural absence only",
    "Present-conditional abundance only",
    "Both components",
    "Omnibus only"
)
comparison_overlap_categories <- c(
    "Significant in at least one comparison method",
    "Not significant in any comparison method"
)
# MaAsLin3 and ZINQ each contribute their package-provided combined result.
# This gives every comparison family one discovery set.
comparison_methods <- c(
    "MaAsLin3 combined", "ZINQ combined", "ANCOM-BC2", "LinDA",
    "corncob", "edgeR", "DESeq2", "metagenomeSeq"
)

# Availability is compared across the same prespecified taxon--dataset grid.
# DASRA contributes its omnibus result; MaAsLin3 and ZINQ contribute their
# package-provided combined results; each remaining family contributes its
# primary differential-abundance result.
availability_result_labels <- c(
    "DASRA combined" = "DASRA",
    "MaAsLin3 combined" = "MaAsLin3",
    "ZINQ combined" = "ZINQ",
    "ANCOM-BC2" = "ANCOM-BC2",
    "LinDA" = "LinDA",
    "corncob" = "corncob",
    "edgeR" = "edgeR",
    "DESeq2" = "DESeq2",
    "metagenomeSeq" = "metagenomeSeq"
)
dasra_availability_result_sets <- dasra_methods

for (dataset_input in completed_datasets) {
    dataset_id <- dataset_input$dataset
    results <- dataset_input$results

    for (method_name in dasra_methods) {
        method_rows <- results[results$method == method_name, , drop = FALSE]
        unavailable_rows <- method_rows[!method_rows$available, , drop = FALSE]
        unavailable_reasons <- table(humanize_unavailable_reason(
            unavailable_rows$reason
        ))
        if (!length(unavailable_reasons)) {
            unavailable_reasons <- stats::setNames(0L, "None")
        }

        for (reason_name in names(unavailable_reasons)) {
            unavailable_count <- nrow(unavailable_rows)
            reason_count <- as.integer(unavailable_reasons[[reason_name]])
            availability_rows[[length(availability_rows) + 1L]] <- data.frame(
                dataset = dataset_id,
                dataset_label = dataset_input$dataset_label,
                reference_samples = dataset_input$reference_samples,
                comparison_samples = dataset_input$comparison_samples,
                tested_taxa = dataset_input$tested_taxa,
                dasra_result = unname(dasra_result_labels[[method_name]]),
                available_taxa = sum(method_rows$available),
                unavailable_taxa = unavailable_count,
                availability_percent =
                    100 * mean(method_rows$available),
                unavailable_reason_description = reason_name,
                taxa_with_reason = reason_count,
                percent_of_unavailable_taxa = if (unavailable_count > 0L) {
                    100 * reason_count / unavailable_count
                } else {
                    NA_real_
                },
                stringsAsFactors = FALSE
            )
        }
    }

    structural_rows <- results[
        results$method == "DASRA structural absence", , drop = FALSE
    ]
    abundance_rows <- results[
        results$method == "DASRA present-conditional abundance",
        , drop = FALSE
    ]
    combined_discoveries <- results[
        results$method == "DASRA combined" &
            results$available & results$significant,
        , drop = FALSE
    ]

    structural_significant <- stats::setNames(
        structural_rows$available & structural_rows$significant,
        structural_rows$taxon
    )
    abundance_significant <- stats::setNames(
        abundance_rows$available & abundance_rows$significant,
        abundance_rows$taxon
    )
    total_combined_discoveries <- nrow(combined_discoveries)
    two_component_discoveries <- combined_discoveries[
        combined_discoveries$components_used == "both",
        , drop = FALSE
    ]
    discovered_taxa <- two_component_discoveries$taxon
    two_component_discovery_count <- length(discovered_taxa)
    one_component_discovery_count <-
        total_combined_discoveries - two_component_discovery_count
    structural_discovery <- unname(structural_significant[discovered_taxa])
    abundance_discovery <- unname(abundance_significant[discovered_taxa])
    discovery_category <- ifelse(
        structural_discovery & abundance_discovery, "Both components",
        ifelse(
            structural_discovery, "Structural absence only",
            ifelse(
                abundance_discovery,
                "Present-conditional abundance only",
                "Omnibus only"
            )
        )
    )
    category_counts <- table(factor(
        discovery_category, levels = signal_categories
    ))

    for (category_name in signal_categories) {
        category_count <- as.integer(category_counts[[category_name]])
        signal_category_rows[[length(signal_category_rows) + 1L]] <-
            data.frame(
                dataset = dataset_id,
                dataset_label = dataset_input$dataset_label,
                dasra_combined_bh_discoveries =
                    total_combined_discoveries,
                discoveries_with_both_components_available =
                    two_component_discovery_count,
                discoveries_with_one_component_available =
                    one_component_discovery_count,
                signal_category = category_name,
                taxa = category_count,
                percent_of_two_component_discoveries =
                    if (two_component_discovery_count > 0L) {
                        100 * category_count /
                            two_component_discovery_count
                    } else {
                        NA_real_
                    },
                stringsAsFactors = FALSE
            )
    }

    comparison_significant_taxa <- unique(results$taxon[
        results$method %in% comparison_methods &
            results$available & results$significant
    ])
    all_discovered_taxa <- combined_discoveries$taxon
    combined_significant_in_comparison_method <-
        all_discovered_taxa %in% comparison_significant_taxa
    significant_in_comparison_method <-
        discovered_taxa %in% comparison_significant_taxa
    overlap_count <- c(
        "Significant in at least one comparison method" =
            sum(combined_significant_in_comparison_method),
        "Not significant in any comparison method" =
            sum(!combined_significant_in_comparison_method)
    )

    for (overlap_name in comparison_overlap_categories) {
        count <- unname(overlap_count[[overlap_name]])
        comparison_overlap_rows[[length(comparison_overlap_rows) + 1L]] <-
            data.frame(
                dataset = dataset_id,
                dataset_label = dataset_input$dataset_label,
                dasra_combined_bh_discoveries =
                    total_combined_discoveries,
                comparison_overlap_status = overlap_name,
                taxa = count,
                percent_of_dasra_combined_discoveries =
                    if (total_combined_discoveries > 0L) {
                        100 * count / total_combined_discoveries
                    } else {
                        NA_real_
                    },
                stringsAsFactors = FALSE
            )
    }

    for (category_name in signal_categories) {
        in_category <- discovery_category == category_name
        category_total <- sum(in_category)
        signal_overlap_rows[[length(signal_overlap_rows) + 1L]] <-
            data.frame(
                dataset = dataset_id,
                dataset_label = dataset_input$dataset_label,
                signal_category = category_name,
                dasra_discoveries_in_category = category_total,
                significant_in_at_least_one_comparison_method = sum(
                    in_category & significant_in_comparison_method
                ),
                not_significant_in_any_comparison_method = sum(
                    in_category & !significant_in_comparison_method
                ),
                percent_in_category_not_significant_in_any_comparison_method =
                    if (category_total > 0L) {
                        100 * sum(
                            in_category &
                                !significant_in_comparison_method
                        ) / category_total
                    } else {
                        NA_real_
                    },
                stringsAsFactors = FALSE
            )
    }
}

availability_summary <- do.call(rbind, availability_rows)
signal_category_summary <- do.call(rbind, signal_category_rows)
comparison_overlap_summary <- do.call(rbind, comparison_overlap_rows)
signal_overlap_summary <- do.call(rbind, signal_overlap_rows)

expected_availability_results <- sum(vapply(
    completed_datasets, function(dataset_input) dataset_input$tested_taxa,
    integer(1)
))
expected_results_by_dataset <- stats::setNames(
    vapply(
        completed_datasets, function(dataset_input) dataset_input$tested_taxa,
        integer(1)
    ),
    expected_dataset_ids
)

method_availability_rows <- lapply(
    seq_along(availability_result_labels),
    function(index) {
        result_set <- names(availability_result_labels)[[index]]
        result_rows <- all_results[
            all_results$method == result_set, , drop = FALSE
        ]
        observed_by_dataset <- table(factor(
            result_rows$dataset, levels = expected_dataset_ids
        ))
        complete_grid <- nrow(result_rows) == expected_availability_results &&
            !anyDuplicated(result_rows[, c("dataset", "taxon")]) &&
            identical(
                as.integer(observed_by_dataset),
                unname(expected_results_by_dataset)
            )
        if (!complete_grid || anyNA(result_rows$available)) {
            stop(sprintf(
                "Availability rows are incomplete for %s.", result_set
            ), call. = FALSE)
        }

        data.frame(
            display_order = index,
            method = unname(availability_result_labels[[result_set]]),
            result_set = result_set,
            taxon_dataset_results = nrow(result_rows),
            available_results = sum(result_rows$available),
            unavailable_results = sum(!result_rows$available),
            availability_percent = 100 * mean(result_rows$available),
            stringsAsFactors = FALSE
        )
    }
)
method_availability_summary <- do.call(rbind, method_availability_rows)

dasra_unavailability_rows <- lapply(
    dasra_availability_result_sets,
    function(result_set) {
        result_rows <- all_results[
            all_results$method == result_set, , drop = FALSE
        ]
        unavailable_rows <- result_rows[
            !result_rows$available, , drop = FALSE
        ]
        reasons <- trimws(unavailable_rows$reason)
        if (!nrow(unavailable_rows)) {
            return(data.frame(
                dasra_result = character(0),
                result_set = character(0),
                taxon_dataset_results = integer(0),
                unavailable_results = integer(0),
                unavailable_reason = character(0),
                results_with_reason = integer(0),
                percent_of_all_results = numeric(0),
                percent_of_unavailable_results = numeric(0),
                stringsAsFactors = FALSE
            ))
        }
        if (anyNA(reasons) || any(!nzchar(reasons)) ||
            any(reasons == "Available")) {
            stop(sprintf(
                "Unavailable DASRA reasons are incomplete for %s.", result_set
            ), call. = FALSE)
        }
        reason_counts <- sort(table(reasons), decreasing = TRUE)

        data.frame(
            dasra_result = unname(
                dasra_result_labels[[result_set]]
            ),
            result_set = result_set,
            taxon_dataset_results = nrow(result_rows),
            unavailable_results = nrow(unavailable_rows),
            unavailable_reason = names(reason_counts),
            results_with_reason = as.integer(reason_counts),
            percent_of_all_results =
                100 * as.integer(reason_counts) / nrow(result_rows),
            percent_of_unavailable_results =
                100 * as.integer(reason_counts) / nrow(unavailable_rows),
            stringsAsFactors = FALSE
        )
    }
)
dasra_unavailability_summary <- do.call(rbind, dasra_unavailability_rows)

reason_totals_by_result <- stats::setNames(
    integer(length(dasra_availability_result_sets)),
    dasra_availability_result_sets
)
if (nrow(dasra_unavailability_summary)) {
    observed_reason_totals <- tapply(
        dasra_unavailability_summary$results_with_reason,
        dasra_unavailability_summary$result_set,
        sum
    )
    reason_totals_by_result[names(observed_reason_totals)] <- as.integer(
        observed_reason_totals
    )
}
expected_unavailable_by_result <- stats::setNames(
    vapply(
        dasra_availability_result_sets,
        function(result_set) {
            sum(
                all_results$method == result_set &
                    !all_results$available
            )
        },
        integer(1)
    ),
    dasra_availability_result_sets
)
if (!identical(
    unname(reason_totals_by_result),
    unname(expected_unavailable_by_result)
)) {
    stop("DASRA unavailability reasons do not match availability totals.",
         call. = FALSE)
}

equal_dataset_rows <- lapply(
    list(
        "DASRA signal category" = list(
            data = signal_category_summary,
            label_column = "signal_category",
            percentage_column = "percent_of_two_component_discoveries"
        ),
        "Overlap with comparison methods" = list(
            data = comparison_overlap_summary,
            label_column = "comparison_overlap_status",
            percentage_column = "percent_of_dasra_combined_discoveries"
        )
    ),
    function(section) {
        data <- section$data
        labels <- unique(data[[section$label_column]])
        do.call(rbind, lapply(labels, function(label) {
            percentages <- data[[section$percentage_column]][
                data[[section$label_column]] == label &
                    !is.na(data[[section$percentage_column]])
            ]
            data.frame(
                classification = label,
                datasets_contributing = length(percentages),
                equal_dataset_mean_percent = mean(percentages),
                median_dataset_percent = stats::median(percentages),
                minimum_dataset_percent = min(percentages),
                maximum_dataset_percent = max(percentages),
                calculation = paste(
                    "Arithmetic summary of within-dataset percentages;",
                    "each dataset receives equal weight"
                ),
                stringsAsFactors = FALSE
            )
        }))
    }
)
equal_dataset_summary <- do.call(rbind, Map(function(name, rows) {
    rows$summary_type <- name
    rows[, c(
        "summary_type", "classification", "datasets_contributing",
        "equal_dataset_mean_percent", "median_dataset_percent",
        "minimum_dataset_percent", "maximum_dataset_percent", "calculation"
    )]
}, names(equal_dataset_rows), equal_dataset_rows))

overall_signal_rows <- lapply(signal_categories, function(category_name) {
    category_rows <- signal_category_summary[
        signal_category_summary$signal_category == category_name,
        , drop = FALSE
    ]
    dataset_discovery_counts <- unique(signal_category_summary[c(
        "dataset", "dasra_combined_bh_discoveries",
        "discoveries_with_both_components_available",
        "discoveries_with_one_component_available"
    )])
    total_combined_discoveries <- sum(
        dataset_discovery_counts$dasra_combined_bh_discoveries
    )
    total_two_component_discoveries <- sum(
        dataset_discovery_counts$discoveries_with_both_components_available
    )
    total_one_component_discoveries <- sum(
        dataset_discovery_counts$discoveries_with_one_component_available
    )
    category_count <- sum(category_rows$taxa)
    category_overlap_rows <- signal_overlap_summary[
        signal_overlap_summary$signal_category == category_name,
        , drop = FALSE
    ]
    significant_in_comparison_count <- sum(
        category_overlap_rows$significant_in_at_least_one_comparison_method
    )
    not_significant_in_comparison_count <- sum(
        category_overlap_rows$not_significant_in_any_comparison_method
    )
    equal_dataset_percent <- equal_dataset_summary$equal_dataset_mean_percent[
        equal_dataset_summary$summary_type == "DASRA signal category" &
            equal_dataset_summary$classification == category_name
    ]
    equal_dataset_n <- equal_dataset_summary$datasets_contributing[
        equal_dataset_summary$summary_type == "DASRA signal category" &
            equal_dataset_summary$classification == category_name
    ]
    data.frame(
        signal_category = category_name,
        datasets_included = length(unique(signal_category_summary$dataset)),
        all_dasra_combined_bh_discoveries = total_combined_discoveries,
        discoveries_with_both_components_available =
            total_two_component_discoveries,
        discoveries_with_one_component_available =
            total_one_component_discoveries,
        taxa_across_datasets = category_count,
        percent_of_two_component_discoveries =
            if (total_two_component_discoveries > 0L) {
                100 * category_count / total_two_component_discoveries
            } else {
                NA_real_
            },
        significant_in_at_least_one_comparison_method =
            significant_in_comparison_count,
        not_significant_in_any_comparison_method =
            not_significant_in_comparison_count,
        percent_in_category_not_significant_in_any_comparison_method =
            if (category_count > 0L) {
                100 * not_significant_in_comparison_count / category_count
            } else {
                NA_real_
            },
        equal_dataset_descriptive_percent = equal_dataset_percent,
        datasets_contributing_to_equal_dataset_summary = equal_dataset_n,
        calculation_note = paste(
            "Mechanism categories use combined discoveries for which both",
            "component results were available; pooled summaries are",
            "descriptive, and dataset-averaged percentages average the",
            "within-dataset percentages among contributing datasets"
        ),
        stringsAsFactors = FALSE
    )
})
overall_signal_summary <- do.call(rbind, overall_signal_rows)

availability_summary$availability_percent <- round(
    availability_summary$availability_percent, 2
)
availability_summary$percent_of_unavailable_taxa <- round(
    availability_summary$percent_of_unavailable_taxa, 2
)
signal_category_summary$percent_of_two_component_discoveries <- round(
    signal_category_summary$percent_of_two_component_discoveries, 2
)
comparison_overlap_summary$percent_of_dasra_combined_discoveries <- round(
    comparison_overlap_summary$percent_of_dasra_combined_discoveries, 2
)
signal_overlap_summary$
    percent_in_category_not_significant_in_any_comparison_method <- round(
        signal_overlap_summary$
            percent_in_category_not_significant_in_any_comparison_method,
        2
    )
equal_dataset_percentage_columns <- c(
    "equal_dataset_mean_percent", "median_dataset_percent",
    "minimum_dataset_percent", "maximum_dataset_percent"
)
equal_dataset_summary[equal_dataset_percentage_columns] <- lapply(
    equal_dataset_summary[equal_dataset_percentage_columns], round, digits = 2
)
overall_signal_summary$percent_of_two_component_discoveries <- round(
    overall_signal_summary$percent_of_two_component_discoveries, 2
)
overall_signal_summary$equal_dataset_descriptive_percent <- round(
    overall_signal_summary$equal_dataset_descriptive_percent, 2
)
overall_signal_summary$
    percent_in_category_not_significant_in_any_comparison_method <- round(
        overall_signal_summary$
            percent_in_category_not_significant_in_any_comparison_method,
        2
    )
method_availability_output <- method_availability_summary[, c(
    "method", "result_set", "taxon_dataset_results", "available_results",
    "unavailable_results", "availability_percent"
), drop = FALSE]
method_availability_output$availability_percent <- round(
    method_availability_output$availability_percent, 2
)
dasra_unavailability_output <- dasra_unavailability_summary
dasra_unavailability_output$percent_of_all_results <- round(
    dasra_unavailability_output$percent_of_all_results, 2
)
dasra_unavailability_output$percent_of_unavailable_results <- round(
    dasra_unavailability_output$percent_of_unavailable_results, 2
)

utils::write.csv(
    availability_summary,
    file.path(output_directory, "dasra_availability_by_dataset.csv"),
    row.names = FALSE, quote = TRUE, na = ""
)
utils::write.csv(
    signal_category_summary,
    file.path(output_directory, "dasra_signal_categories_by_dataset.csv"),
    row.names = FALSE, quote = TRUE, na = ""
)
utils::write.csv(
    comparison_overlap_summary,
    file.path(
        output_directory,
        "dasra_combined_discovery_comparison_overlap_by_dataset.csv"
    ),
    row.names = FALSE, quote = TRUE, na = ""
)
utils::write.csv(
    signal_overlap_summary,
    file.path(
        output_directory,
        "dasra_signal_category_comparison_overlap_by_dataset.csv"
    ),
    row.names = FALSE, quote = TRUE, na = ""
)
utils::write.csv(
    equal_dataset_summary,
    file.path(output_directory, "equal_dataset_descriptive_percentages.csv"),
    row.names = FALSE, quote = TRUE, na = ""
)
utils::write.csv(
    overall_signal_summary,
    file.path(
        output_directory, "dasra_signal_overall_descriptive_summary.csv"
    ),
    row.names = FALSE, quote = TRUE, na = ""
)
utils::write.csv(
    method_availability_output,
    file.path(output_directory, "method_availability.csv"),
    row.names = FALSE, quote = TRUE, na = ""
)
utils::write.csv(
    dasra_unavailability_output,
    file.path(output_directory, "dasra_unavailability_reasons.csv"),
    row.names = FALSE, quote = TRUE, na = ""
)

required_plot_packages <- c("ggplot2", "patchwork")
missing_plot_packages <- required_plot_packages[!vapply(
    required_plot_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_plot_packages)) {
    stop(sprintf(
        "Required plotting packages are unavailable: %s.",
        paste(missing_plot_packages, collapse = ", ")
    ), call. = FALSE)
}

plot_data <- signal_category_summary
plot_data$signal_category <- factor(
    plot_data$signal_category, levels = signal_categories
)
dataset_totals <- unique(plot_data[c(
    "dataset", "dataset_label",
    "discoveries_with_both_components_available"
)])
dataset_totals <- dataset_totals[order(
    dataset_totals$discoveries_with_both_components_available,
    decreasing = TRUE
), , drop = FALSE]
dataset_totals$plot_dataset_label <- sprintf(
    "%s (n=%d)",
    dataset_totals$dataset_label,
    dataset_totals$discoveries_with_both_components_available
)
plot_data$plot_dataset_label <- dataset_totals$plot_dataset_label[
    match(plot_data$dataset, dataset_totals$dataset)
]
plot_data$plot_dataset_label <- factor(
    plot_data$plot_dataset_label,
    levels = rev(dataset_totals$plot_dataset_label)
)
plot_data$plot_percent <- ifelse(
    plot_data$discoveries_with_both_components_available == 0L,
    0,
    100 * plot_data$taxa /
        plot_data$discoveries_with_both_components_available
)

signal_colors <- c(
    "Structural absence only" = "#0072B2",
    "Present-conditional abundance only" = "#D55E00",
    "Both components" = "#009E73",
    "Omnibus only" = "#CC79A7"
)

signal_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
        x = plot_dataset_label, y = plot_percent, fill = signal_category
    )
) +
    ggplot2::geom_col(
        width = 0.72, position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::coord_flip(ylim = c(0, 100)) +
    ggplot2::scale_fill_manual(values = signal_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
        breaks = seq(0, 100, by = 20),
        labels = function(value) paste0(value, "%"),
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
        x = NULL,
        y = "Percentage of classified discoveries",
        fill = NULL
    ) +
    ggplot2::guides(
        fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::theme_classic(base_size = 10.5) +
    ggplot2::theme(
        legend.position = "top",
        legend.text = ggplot2::element_text(size = 9.4),
        legend.box.just = "center",
        panel.grid.major.x = ggplot2::element_line(
            color = "grey88", linewidth = 0.35
        ),
        axis.text.y = ggplot2::element_text(color = "grey15"),
        plot.margin = ggplot2::margin(7, 16, 7, 7)
    )

signal_figure_file <- file.path(
    output_directory, "dasra_signal_categories_by_dataset.pdf"
)
grDevices::pdf(
    signal_figure_file,
    width = 174 / 25.4,
    height = max(3.8, 1.55 + 0.36 * nrow(dataset_totals)),
    family = "Helvetica",
    useDingbats = FALSE,
    version = "1.5"
)
print(signal_plot)
grDevices::dev.off()

method_plot_data <- method_availability_summary
method_plot_data$method <- factor(
    method_plot_data$method,
    levels = rev(unname(availability_result_labels))
)
method_plot_data$plot_group <- "Comparison method"
method_plot_data$plot_group[
    method_plot_data$result_set == "DASRA combined"
] <- "DASRA"

availability_colors <- c(
    "DASRA" = "#0072B2",
    "Comparison method" = "#606060"
)
method_availability_plot <- ggplot2::ggplot(
    method_plot_data,
    ggplot2::aes(
        x = method, y = availability_percent, fill = plot_group
    )
) +
    ggplot2::geom_col(width = 0.64) +
    ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.1f%%", availability_percent)),
        hjust = -0.18,
        color = "grey10",
        size = 3.15
    ) +
    ggplot2::coord_flip(ylim = c(0, 106), clip = "off") +
    ggplot2::scale_fill_manual(values = availability_colors) +
    ggplot2::scale_y_continuous(
        breaks = seq(0, 100, by = 25),
        labels = function(value) paste0(value, "%"),
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
        x = NULL,
        y = "Availability",
        title = "Availability across methods",
        subtitle = sprintf(
            "%d taxon-by-dataset results per method across seven datasets",
            expected_availability_results
        )
    ) +
    ggplot2::theme_classic(base_size = 10.5) +
    ggplot2::theme(
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_line(
            color = "grey88", linewidth = 0.35
        ),
        axis.text.y = ggplot2::element_text(color = "grey15", size = 9.2),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        plot.subtitle = ggplot2::element_text(color = "grey30", size = 9.2),
        plot.margin = ggplot2::margin(7, 28, 4, 7)
    )

reason_totals <- stats::aggregate(
    results_with_reason ~ unavailable_reason,
    dasra_unavailability_summary,
    sum
)
reason_totals <- reason_totals[order(
    -reason_totals$results_with_reason,
    reason_totals$unavailable_reason,
    method = "radix"
), , drop = FALSE]
reason_order <- reason_totals$unavailable_reason
wrapped_reason_labels <- stats::setNames(vapply(
    reason_order,
    function(reason) paste(strwrap(reason, width = 49), collapse = "\n"),
    character(1)
), reason_order)

dasra_result_plot_labels <- stats::setNames(vapply(
    dasra_availability_result_sets,
    function(result_set) {
        result_rows <- all_results[
            all_results$method == result_set, , drop = FALSE
        ]
        result_name <- if (result_set == "DASRA structural absence") {
            "Structural-absence\ncomponent"
        } else if (
            result_set == "DASRA present-conditional abundance"
        ) {
            "Present-conditional\nabundance component"
        } else {
            "Omnibus\nresult"
        }
        sprintf(
            "%s\n%d/%d unavailable",
            result_name,
            sum(!result_rows$available),
            nrow(result_rows)
        )
    },
    character(1)
), dasra_availability_result_sets)

reason_plot_data <- dasra_unavailability_summary
reason_plot_data$reason_label <- factor(
    unname(wrapped_reason_labels[reason_plot_data$unavailable_reason]),
    levels = rev(unname(wrapped_reason_labels))
)
reason_plot_data$result_label <- factor(
    unname(dasra_result_plot_labels[reason_plot_data$result_set]),
    levels = unname(dasra_result_plot_labels)
)
reason_plot_data$cell_label <- sprintf(
    "%d (%.1f%%)",
    reason_plot_data$results_with_reason,
    reason_plot_data$percent_of_unavailable_results
)

reason_grid <- expand.grid(
    reason_label = levels(reason_plot_data$reason_label),
    result_label = levels(reason_plot_data$result_label),
    stringsAsFactors = FALSE
)
reason_grid$reason_label <- factor(
    reason_grid$reason_label,
    levels = levels(reason_plot_data$reason_label)
)
reason_grid$result_label <- factor(
    reason_grid$result_label,
    levels = levels(reason_plot_data$result_label)
)

dasra_result_cell_colors <- stats::setNames(
    c("#DCEEF8", "#F8E2D5", "#DCEFE8"),
    unname(dasra_result_plot_labels)
)
dasra_reason_plot <- ggplot2::ggplot(
    reason_plot_data,
    ggplot2::aes(x = result_label, y = reason_label)
) +
    ggplot2::geom_tile(
        data = reason_grid,
        ggplot2::aes(x = result_label, y = reason_label),
        inherit.aes = FALSE,
        fill = "grey97",
        color = "white",
        linewidth = 0.7,
        width = 0.94,
        height = 0.88
    ) +
    ggplot2::geom_tile(
        ggplot2::aes(fill = result_label),
        color = "white",
        linewidth = 0.7,
        width = 0.94,
        height = 0.88
    ) +
    ggplot2::geom_text(
        ggplot2::aes(label = cell_label),
        size = 3.05,
        color = "grey10"
    ) +
    ggplot2::scale_fill_manual(values = dasra_result_cell_colors) +
    ggplot2::scale_x_discrete(position = "top", drop = FALSE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::labs(
        x = NULL,
        y = NULL,
        title = "DASRA unavailability reasons",
        subtitle = paste(
            "Cells show count (percentage of unavailable results",
            "within each result set)"
        )
    ) +
    ggplot2::theme_minimal(base_size = 10.5) +
    ggplot2::theme(
        legend.position = "none",
        panel.grid = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
            color = "grey15", face = "bold", size = 9.1,
            lineheight = 0.95, margin = ggplot2::margin(b = 5)
        ),
        axis.text.y = ggplot2::element_text(
            color = "grey15", size = 8.6, lineheight = 0.95
        ),
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        plot.subtitle = ggplot2::element_text(color = "grey30", size = 9.2),
        plot.margin = ggplot2::margin(5, 10, 7, 7)
    )

method_availability_figure <- patchwork::wrap_plots(
    patchwork::free(method_availability_plot, side = "l"),
    dasra_reason_plot,
    ncol = 1,
    heights = c(0.82, 1.35)
) +
    patchwork::plot_annotation(
        tag_levels = "a",
        theme = ggplot2::theme(
            plot.tag = ggplot2::element_text(face = "bold", size = 11)
        )
    )

method_availability_figure_file <- file.path(
    output_directory, "method_availability.pdf"
)
grDevices::pdf(
    method_availability_figure_file,
    width = 190 / 25.4,
    height = 185 / 25.4,
    family = "Helvetica",
    useDingbats = FALSE,
    version = "1.5",
    timestamp = FALSE
)
print(method_availability_figure)
grDevices::dev.off()

message(sprintf(
    "Summarized %d completed real-data datasets in %s",
    length(completed_datasets), output_directory
))
