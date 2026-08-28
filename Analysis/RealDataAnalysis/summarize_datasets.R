# Summarize completed real-data analyses without rerunning any method
#
# A dataset is included only when one method-results table and one input-summary
# table are present in its table directory. The method-results table is assumed
# to contain one row per taxon and method. Its `significant` column is the
# formal within-result-set BH decision written by analysis.R.

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

dataset_directories <- list.dirs(
    real_data_directory, full.names = TRUE, recursive = FALSE
)

read_completed_dataset <- function(dataset_directory) {
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
    if (length(results_file) != 1L || length(input_summary_file) != 1L) {
        return(NULL)
    }

    results <- utils::read.csv(
        results_file, check.names = FALSE, stringsAsFactors = FALSE
    )
    input_summary <- utils::read.csv(
        input_summary_file, check.names = FALSE, stringsAsFactors = FALSE
    )

    required_results_columns <- c(
        "dataset", "taxon", "method", "family", "component",
        "available", "reason", "significant"
    )
    if (!all(required_results_columns %in% names(results)) ||
        !all(c("item", "value") %in% names(input_summary))) {
        stop(sprintf(
            "Completed tables for %s do not follow the real-data output schema.",
            basename(dataset_directory)
        ), call. = FALSE)
    }

    summary_values <- stats::setNames(
        input_summary$value, input_summary$item
    )
    dataset_id <- unname(summary_values[["dataset"]])
    if (!identical(unique(results$dataset), dataset_id)) {
        stop(sprintf(
            "Dataset identifiers disagree between the two completed tables for %s.",
            basename(dataset_directory)
        ), call. = FALSE)
    }

    results$available <- as.logical(results$available)
    results$significant <- as.logical(results$significant)

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

completed_datasets <- Filter(
    Negate(is.null), lapply(dataset_directories, read_completed_dataset)
)
if (!length(completed_datasets)) {
    stop("No completed real-data result pairs were found.", call. = FALSE)
}

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
selection_rows <- list()
signal_selection_rows <- list()

signal_categories <- c(
    "Structural absence only",
    "Present-conditional abundance only",
    "Both components",
    "Omnibus only"
)
selection_categories <- c(
    "Selected by at least one non-DASRA method",
    "Not selected by any non-DASRA method"
)
# MaAsLin3 and ZINQ each contribute their formal combined result so that every
# non-DASRA method family contributes one selection set.
non_dasra_methods <- c(
    "MaAsLin3 combined", "ZINQ combined", "ANCOM-BC2", "LinDA",
    "corncob", "edgeR", "DESeq2", "metagenomeSeq"
)

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

    structural_selected <- stats::setNames(
        structural_rows$available & structural_rows$significant,
        structural_rows$taxon
    )
    abundance_selected <- stats::setNames(
        abundance_rows$available & abundance_rows$significant,
        abundance_rows$taxon
    )
    discovered_taxa <- combined_discoveries$taxon
    structural_discovery <- unname(structural_selected[discovered_taxa])
    abundance_discovery <- unname(abundance_selected[discovered_taxa])
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
    total_discoveries <- length(discovered_taxa)

    for (category_name in signal_categories) {
        category_count <- as.integer(category_counts[[category_name]])
        signal_category_rows[[length(signal_category_rows) + 1L]] <-
            data.frame(
                dataset = dataset_id,
                dataset_label = dataset_input$dataset_label,
                dasra_combined_bh_discoveries = total_discoveries,
                signal_category = category_name,
                taxa = category_count,
                percent_of_dasra_combined_discoveries =
                    if (total_discoveries > 0L) {
                        100 * category_count / total_discoveries
                    } else {
                        NA_real_
                    },
                stringsAsFactors = FALSE
            )
    }

    non_dasra_selected_taxa <- unique(results$taxon[
        results$method %in% non_dasra_methods &
            results$available & results$significant
    ])
    selected_by_other_method <- discovered_taxa %in% non_dasra_selected_taxa
    selection_count <- c(
        "Selected by at least one non-DASRA method" =
            sum(selected_by_other_method),
        "Not selected by any non-DASRA method" =
            sum(!selected_by_other_method)
    )

    for (selection_name in selection_categories) {
        count <- unname(selection_count[[selection_name]])
        selection_rows[[length(selection_rows) + 1L]] <- data.frame(
            dataset = dataset_id,
            dataset_label = dataset_input$dataset_label,
            dasra_combined_bh_discoveries = total_discoveries,
            non_dasra_selection_status = selection_name,
            taxa = count,
            percent_of_dasra_combined_discoveries =
                if (total_discoveries > 0L) {
                    100 * count / total_discoveries
                } else {
                    NA_real_
                },
            stringsAsFactors = FALSE
        )
    }

    for (category_name in signal_categories) {
        in_category <- discovery_category == category_name
        category_total <- sum(in_category)
        signal_selection_rows[[length(signal_selection_rows) + 1L]] <-
            data.frame(
                dataset = dataset_id,
                dataset_label = dataset_input$dataset_label,
                signal_category = category_name,
                dasra_discoveries_in_category = category_total,
                selected_by_at_least_one_non_dasra_method = sum(
                    in_category & selected_by_other_method
                ),
                not_selected_by_any_non_dasra_method = sum(
                    in_category & !selected_by_other_method
                ),
                percent_in_category_not_selected_by_any_non_dasra =
                    if (category_total > 0L) {
                        100 * sum(in_category & !selected_by_other_method) /
                            category_total
                    } else {
                        NA_real_
                    },
                stringsAsFactors = FALSE
            )
    }
}

availability_summary <- do.call(rbind, availability_rows)
signal_category_summary <- do.call(rbind, signal_category_rows)
selection_summary <- do.call(rbind, selection_rows)
signal_selection_summary <- do.call(rbind, signal_selection_rows)

equal_dataset_rows <- lapply(
    list(
        "DASRA signal category" = list(
            data = signal_category_summary,
            label_column = "signal_category"
        ),
        "Selection by non-DASRA methods" = list(
            data = selection_summary,
            label_column = "non_dasra_selection_status"
        )
    ),
    function(section) {
        data <- section$data
        labels <- unique(data[[section$label_column]])
        do.call(rbind, lapply(labels, function(label) {
            percentages <- data$percent_of_dasra_combined_discoveries[
                data[[section$label_column]] == label &
                    !is.na(data$percent_of_dasra_combined_discoveries)
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
    total_discoveries <- sum(unique(
        signal_category_summary[c(
            "dataset", "dasra_combined_bh_discoveries"
        )]
    )$dasra_combined_bh_discoveries)
    category_count <- sum(category_rows$taxa)
    category_selection_rows <- signal_selection_summary[
        signal_selection_summary$signal_category == category_name,
        , drop = FALSE
    ]
    selected_by_other_count <- sum(
        category_selection_rows$selected_by_at_least_one_non_dasra_method
    )
    not_selected_by_other_count <- sum(
        category_selection_rows$not_selected_by_any_non_dasra_method
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
        all_dasra_combined_bh_discoveries = total_discoveries,
        taxa_across_datasets = category_count,
        percent_of_all_dasra_combined_discoveries =
            if (total_discoveries > 0L) {
                100 * category_count / total_discoveries
            } else {
                NA_real_
            },
        selected_by_at_least_one_non_dasra_method = selected_by_other_count,
        not_selected_by_any_non_dasra_method = not_selected_by_other_count,
        percent_in_category_not_selected_by_any_non_dasra =
            if (category_count > 0L) {
                100 * not_selected_by_other_count / category_count
            } else {
                NA_real_
            },
        equal_dataset_descriptive_percent = equal_dataset_percent,
        datasets_contributing_to_equal_dataset_summary = equal_dataset_n,
        calculation_note = paste(
            "Pooled counts and percentages are descriptive; the equal-dataset",
            "percentage is the arithmetic mean of within-dataset percentages",
            "among datasets with at least one combined discovery; no",
            "meta-analysis was performed"
        ),
        stringsAsFactors = FALSE
    )
})
overall_signal_summary <- do.call(rbind, overall_signal_rows)

dataset_availability <- unique(availability_summary[c(
    "dataset", "dasra_result", "available_taxa", "unavailable_taxa",
    "availability_percent"
)])
minimum_dataset_availability <- 70

summarize_workability <- function(analysis_set, dataset_ids) {
    component_labels <- c(
        structural = "Structural absence component",
        abundance = "Present-conditional abundance component",
        combined = "Combined omnibus test"
    )
    selected <- dataset_availability[
        dataset_availability$dataset %in% dataset_ids, , drop = FALSE
    ]
    if (!length(dataset_ids)) {
        return(data.frame(
            analysis_set = analysis_set,
            suitability_criterion = paste(
                "Broadly non-working if any DASRA component availability",
                sprintf("is below %d%%", minimum_dataset_availability)
            ),
            datasets_analyzed = 0L,
            broadly_non_working_datasets = 0L,
            broadly_non_working_percent = 0,
            structural_median_dataset_availability_percent = NA_real_,
            structural_minimum_dataset_availability_percent = NA_real_,
            abundance_median_dataset_availability_percent = NA_real_,
            abundance_minimum_dataset_availability_percent = NA_real_,
            combined_median_dataset_availability_percent = NA_real_,
            combined_minimum_dataset_availability_percent = NA_real_,
            pooled_available_component_taxon_results = 0L,
            pooled_unavailable_component_taxon_results = 0L,
            pooled_component_taxon_results = 0L,
            pooled_available_percent = 0,
            pooled_unavailable_percent = 0,
            stringsAsFactors = FALSE
        ))
    }

    minimum_by_dataset <- tapply(
        selected$availability_percent, selected$dataset, min
    )
    component_availability <- lapply(component_labels, function(label) {
        selected$availability_percent[selected$dasra_result == label]
    })
    pooled_available <- sum(selected$available_taxa)
    pooled_unavailable <- sum(selected$unavailable_taxa)
    pooled_total <- pooled_available + pooled_unavailable

    data.frame(
        analysis_set = analysis_set,
        suitability_criterion = paste(
            "Broadly non-working if any DASRA component availability",
            sprintf("is below %d%%", minimum_dataset_availability)
        ),
        datasets_analyzed = length(dataset_ids),
        broadly_non_working_datasets = sum(
            minimum_by_dataset < minimum_dataset_availability
        ),
        broadly_non_working_percent =
            100 * mean(minimum_by_dataset < minimum_dataset_availability),
        structural_median_dataset_availability_percent =
            stats::median(component_availability$structural),
        structural_minimum_dataset_availability_percent =
            min(component_availability$structural),
        abundance_median_dataset_availability_percent =
            stats::median(component_availability$abundance),
        abundance_minimum_dataset_availability_percent =
            min(component_availability$abundance),
        combined_median_dataset_availability_percent =
            stats::median(component_availability$combined),
        combined_minimum_dataset_availability_percent =
            min(component_availability$combined),
        pooled_available_component_taxon_results = pooled_available,
        pooled_unavailable_component_taxon_results = pooled_unavailable,
        pooled_component_taxon_results = pooled_total,
        pooled_available_percent = 100 * pooled_available / pooled_total,
        pooled_unavailable_percent = 100 * pooled_unavailable / pooled_total,
        stringsAsFactors = FALSE
    )
}

all_completed_dataset_ids <- unique(dataset_availability$dataset)
additional_dataset_ids <- setdiff(
    all_completed_dataset_ids, c("crc_baxter", "cdi_schubert")
)
workability_summary <- rbind(
    summarize_workability(
        "Additional datasets excluding CRC and CDI",
        additional_dataset_ids
    ),
    summarize_workability(
        "All completed datasets", all_completed_dataset_ids
    )
)

candidate_screen_file <- file.path(
    real_data_directory, "screening", "dasra_candidate_screen.csv"
)
candidate_screening_summary <- NULL
if (file.exists(candidate_screen_file)) {
    candidate_screen <- utils::read.csv(
        candidate_screen_file, stringsAsFactors = FALSE, check.names = FALSE
    )
    logical_columns <- c(
        "broadly_non_working", "common_background_incompatible",
        "screen_passed"
    )
    candidate_screen[logical_columns] <- lapply(
        candidate_screen[logical_columns], as.logical
    )
    screened_n <- nrow(candidate_screen)
    broadly_non_working_n <- sum(candidate_screen$broadly_non_working)
    background_incompatible_n <- sum(
        candidate_screen$common_background_incompatible
    )
    rejected_n <- sum(!candidate_screen$screen_passed)
    passed_n <- sum(candidate_screen$screen_passed)
    candidate_screening_summary <- data.frame(
        candidates_reaching_dasra_screen = screened_n,
        broadly_non_working_candidates = broadly_non_working_n,
        broadly_non_working_percent = 100 * broadly_non_working_n / screened_n,
        common_background_incompatible_candidates = background_incompatible_n,
        common_background_incompatible_percent =
            100 * background_incompatible_n / screened_n,
        candidates_rejected_by_either_rule = rejected_n,
        candidates_rejected_percent = 100 * rejected_n / screened_n,
        candidates_passing_screen = passed_n,
        candidates_passing_screen_percent = 100 * passed_n / screened_n,
        denominator_note = paste(
            "The denominator includes every candidate that reached DASRA",
            "screening, including candidates not advanced to the final",
            "formal analysis"
        ),
        stringsAsFactors = FALSE
    )
    percentage_columns <- grep(
        "_percent$", names(candidate_screening_summary), value = TRUE
    )
    candidate_screening_summary[percentage_columns] <- lapply(
        candidate_screening_summary[percentage_columns], round, digits = 2
    )
}

availability_summary$availability_percent <- round(
    availability_summary$availability_percent, 2
)
availability_summary$percent_of_unavailable_taxa <- round(
    availability_summary$percent_of_unavailable_taxa, 2
)
signal_category_summary$percent_of_dasra_combined_discoveries <- round(
    signal_category_summary$percent_of_dasra_combined_discoveries, 2
)
selection_summary$percent_of_dasra_combined_discoveries <- round(
    selection_summary$percent_of_dasra_combined_discoveries, 2
)
signal_selection_summary$
    percent_in_category_not_selected_by_any_non_dasra <- round(
        signal_selection_summary$
            percent_in_category_not_selected_by_any_non_dasra,
        2
    )
equal_dataset_percentage_columns <- c(
    "equal_dataset_mean_percent", "median_dataset_percent",
    "minimum_dataset_percent", "maximum_dataset_percent"
)
equal_dataset_summary[equal_dataset_percentage_columns] <- lapply(
    equal_dataset_summary[equal_dataset_percentage_columns], round, digits = 2
)
overall_signal_summary$percent_of_all_dasra_combined_discoveries <- round(
    overall_signal_summary$percent_of_all_dasra_combined_discoveries, 2
)
overall_signal_summary$equal_dataset_descriptive_percent <- round(
    overall_signal_summary$equal_dataset_descriptive_percent, 2
)
overall_signal_summary$
    percent_in_category_not_selected_by_any_non_dasra <- round(
        overall_signal_summary$
            percent_in_category_not_selected_by_any_non_dasra,
        2
    )
workability_percentage_columns <- c(
    "broadly_non_working_percent",
    "structural_median_dataset_availability_percent",
    "structural_minimum_dataset_availability_percent",
    "abundance_median_dataset_availability_percent",
    "abundance_minimum_dataset_availability_percent",
    "combined_median_dataset_availability_percent",
    "combined_minimum_dataset_availability_percent",
    "pooled_available_percent", "pooled_unavailable_percent"
)
workability_summary[workability_percentage_columns] <- lapply(
    workability_summary[workability_percentage_columns], round, digits = 2
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
    selection_summary,
    file.path(
        output_directory,
        "dasra_combined_discovery_selection_by_dataset.csv"
    ),
    row.names = FALSE, quote = TRUE, na = ""
)
utils::write.csv(
    signal_selection_summary,
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
    workability_summary,
    file.path(output_directory, "dasra_workability_summary.csv"),
    row.names = FALSE, quote = TRUE, na = ""
)
if (!is.null(candidate_screening_summary)) {
    utils::write.csv(
        candidate_screening_summary,
        file.path(output_directory, "dasra_candidate_screening_summary.csv"),
        row.names = FALSE, quote = TRUE, na = ""
    )
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to draw the summary figure.",
         call. = FALSE)
}

plot_data <- signal_category_summary
plot_data$signal_category <- factor(
    plot_data$signal_category, levels = signal_categories
)
dataset_totals <- unique(plot_data[c(
    "dataset", "dataset_label", "dasra_combined_bh_discoveries"
)])
dataset_totals <- dataset_totals[order(
    dataset_totals$dasra_combined_bh_discoveries,
    decreasing = TRUE
), , drop = FALSE]
dataset_totals$plot_dataset_label <- sprintf(
    "%s (n=%d)",
    dataset_totals$dataset_label,
    dataset_totals$dasra_combined_bh_discoveries
)
plot_data$plot_dataset_label <- dataset_totals$plot_dataset_label[
    match(plot_data$dataset, dataset_totals$dataset)
]
plot_data$plot_dataset_label <- factor(
    plot_data$plot_dataset_label,
    levels = rev(dataset_totals$plot_dataset_label)
)
plot_data$plot_percent <- ifelse(
    plot_data$dasra_combined_bh_discoveries == 0L,
    0,
    100 * plot_data$taxa / plot_data$dasra_combined_bh_discoveries
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
        title = "Composition of formal DASRA combined discoveries",
        subtitle = paste(
            "Within-dataset percentages based on component-level BH decisions",
            "n denotes the number of combined discoveries",
            sep = "\n"
        ),
        x = NULL,
        y = "DASRA combined BH discoveries",
        fill = "Signal category"
    ) +
    ggplot2::guides(
        fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::theme_classic(base_size = 10.5) +
    ggplot2::theme(
        legend.position = "top",
        legend.title = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        panel.grid.major.x = ggplot2::element_line(
            color = "grey88", linewidth = 0.35
        ),
        axis.text.y = ggplot2::element_text(color = "grey15"),
        plot.margin = ggplot2::margin(7, 16, 7, 7)
    )

ggplot2::ggsave(
    filename = file.path(
        output_directory, "dasra_signal_categories_by_dataset.pdf"
    ),
    plot = signal_plot,
    width = 9,
    height = max(4.8, 2.2 + 0.32 * nrow(dataset_totals)),
    units = "in"
)

message(sprintf(
    "Summarized %d completed real-data datasets in %s",
    length(completed_datasets), output_directory
))
