# Candidate-screening rules were fixed before candidate-level DASRA results
# were inspected. Their definitions and rationale are recorded in
# Analysis/RealDataAnalysis/DATASET_SELECTION.md.

options(stringsAsFactors = FALSE, warn = 1)

script_argument <- grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1L]
script_path <- sub("^--file=", "", script_argument)
analysis_directory <- dirname(normalizePath(script_path))
selection_document <- file.path(
    analysis_directory, "DATASET_SELECTION.md"
)

local_library <- file.path(analysis_directory, "R_lib")
if (dir.exists(local_library)) {
    .libPaths(c(normalizePath(local_library), .libPaths()))
}

if (!file.exists(selection_document)) {
    stop("DATASET_SELECTION.md was not found.", call. = FALSE)
}
if (!requireNamespace("DASRA", quietly = TRUE)) {
    stop("DASRA is not available.", call. = FALSE)
}
dasra_version <- as.character(utils::packageVersion("DASRA"))
if (!identical(dasra_version, "0.6.0")) {
    stop("This candidate screen was specified for DASRA 0.6.0.",
         call. = FALSE)
}

candidate_datasets <- data.frame(
    dataset_id = c(
        "gems_pediatric_diarrhea",
        "korean_hypertension",
        "microbiomehd_zupancic_obesity",
        "qiita_1939_pediatric_crohn",
        "ravel_vaginal_ethnicity"
    ),
    input_file = c(
        "gems_pediatric_diarrhea/processed/gems_pediatric_diarrhea_dasra_input.rds",
        "korean_hypertension/processed/korean_hypertension_dasra_input.rds",
        "microbiomehd_zupancic_obesity/processed/microbiomehd_zupancic_obesity_dasra_input.rds",
        "qiita_1939_pediatric_crohn/processed/qiita_1939_pediatric_crohn_dasra_input.rds",
        "ravel_vaginal_ethnicity/processed/ravel_vaginal_ethnicity_dasra_input.rds"
    ),
    stringsAsFactors = FALSE
)
candidate_datasets$dataset_label <- vapply(
    candidate_datasets$dataset_id,
    function(dataset_id) {
        title <- readLines(
            file.path(analysis_directory, dataset_id, "README.md"),
            n = 1L,
            warn = FALSE
        )
        trimws(sub("^#\\s*", "", title))
    },
    character(1)
)

humanize_reason <- function(reason) {
    reason <- as.character(reason)
    missing_reason <- is.na(reason) | !nzchar(trimws(reason))
    reason[missing_reason] <- "result unavailable without a recorded reason"
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
            "Conditional-present root refinement did not improve the objective",
        conditional_present_root_polish_not_closed =
            "Conditional-present root refinement did not close",
        inference_derivative_unstable =
            "Structural-absence inference derivative was unstable",
        positive_part_design_rank_deficient =
            "Conditional-present fit design was rank deficient",
        rank_deficient_positive_mark_design =
            "Conditional-present fit design was rank deficient",
        structural_absence_nonoptimal_nuisance_fit =
            "Structural-absence nuisance fit was not optimal",
        structural_absence_persistent_boundary =
            "Structural-absence fit remained on a parameter boundary"
    )
    translated <- unname(descriptions[reason])
    untranslated <- is.na(translated)
    translated[untranslated] <- gsub(
        "_", " ", reason[untranslated], fixed = TRUE
    )
    translated
}

count_unavailability_reasons <- function(dataset_id, dataset_label,
                                         component, available, reason) {
    unavailable_reason <- humanize_reason(reason[!available])
    if (!length(unavailable_reason)) {
        return(data.frame(
            dataset_id = character(),
            dataset_label = character(),
            component = character(),
            reason = character(),
            taxa_count = integer(),
            stringsAsFactors = FALSE
        ))
    }
    reason_counts <- sort(table(unavailable_reason), decreasing = TRUE)
    data.frame(
        dataset_id = dataset_id,
        dataset_label = dataset_label,
        component = component,
        reason = names(reason_counts),
        taxa_count = as.integer(reason_counts),
        stringsAsFactors = FALSE
    )
}

screen_one_dataset <- function(candidate) {
    input <- readRDS(file.path(
        analysis_directory, candidate$input_file
    ))
    prespecified_taxa <- rownames(input$counts)
    counts <- input$counts[prespecified_taxa, , drop = FALSE]
    metadata <- input$metadata[colnames(counts), , drop = FALSE]
    preprocessing <- input$preprocessing
    planned_covariates <- preprocessing$planned_adjustment_variables
    reference_group <- preprocessing$reference_group
    comparison_group <- preprocessing$comparison_group
    metadata$group <- factor(
        as.character(metadata$group),
        levels = c(reference_group, comparison_group)
    )
    model_formula <- stats::reformulate(c(planned_covariates, "group"))

    started <- proc.time()[["elapsed"]]
    fit <- DASRA::dasra(
        counts = counts,
        metadata = metadata,
        formula = model_formula,
        group = "group",
        library_size = "library_size",
        taxa_are_rows = TRUE,
        reference = reference_group,
        p_adjust_method = "BH",
        component = "all",
        full_output = FALSE,
        structural_conditional_present_starts = 1L,
        workers = 1L,
        verbose = FALSE,
        store_plot_data = FALSE
    )
    elapsed_seconds <- proc.time()[["elapsed"]] - started

    result_index <- match(prespecified_taxa, fit$results$taxon)
    diagnostic_index <- match(prespecified_taxa, fit$diagnostics$taxon)
    if (anyNA(result_index) || anyNA(diagnostic_index)) {
        stop("DASRA did not return every prespecified taxon.",
             call. = FALSE)
    }
    results <- fit$results[result_index, , drop = FALSE]
    diagnostics <- fit$diagnostics[diagnostic_index, , drop = FALSE]

    structural_available <- diagnostics$retained &
        diagnostics$formed_structural_absence &
        is.finite(results$p_structural_absence)
    abundance_available <- diagnostics$retained &
        diagnostics$formed_relative_abundance &
        is.finite(results$p_relative_abundance)
    combined_available <- structural_available & abundance_available &
        is.finite(results$p_omnibus)

    n_taxa <- length(prespecified_taxa)
    structural_rate <- mean(structural_available)
    abundance_rate <- mean(abundance_available)
    combined_rate <- mean(combined_available)

    abundance_p_for_bh <- ifelse(
        abundance_available, results$p_relative_abundance, 1
    )
    abundance_q <- stats::p.adjust(abundance_p_for_bh, method = "BH")
    abundance_discovery <- abundance_available & abundance_q <= 0.05
    abundance_discovery_n <- sum(abundance_discovery)
    abundance_discovery_rate <- abundance_discovery_n / n_taxa

    broadly_non_working <- any(c(
        structural_rate, abundance_rate, combined_rate
    ) < 0.70)
    common_background_incompatible <- abundance_discovery_rate >= 0.50

    structural_reason <- ifelse(
        diagnostics$retained,
        diagnostics$reason_structural_absence,
        "taxon not retained"
    )
    abundance_reason <- ifelse(
        diagnostics$retained,
        diagnostics$reason_relative_abundance,
        "taxon not retained"
    )
    combined_reason <- ifelse(
        !structural_available & !abundance_available,
        "both component results unavailable",
        ifelse(
            !structural_available,
            "structural-absence result unavailable",
            ifelse(
                !abundance_available,
                "present-conditional abundance result unavailable",
                "combined p-value unavailable"
            )
        )
    )

    reason_table <- rbind(
        count_unavailability_reasons(
            candidate$dataset_id, candidate$dataset_label,
            "structural absence",
            structural_available, structural_reason
        ),
        count_unavailability_reasons(
            candidate$dataset_id, candidate$dataset_label,
            "present-conditional abundance",
            abundance_available, abundance_reason
        ),
        count_unavailability_reasons(
            candidate$dataset_id, candidate$dataset_label, "combined",
            combined_available, combined_reason
        )
    )

    summary_row <- data.frame(
        dataset_id = candidate$dataset_id,
        dataset_label = candidate$dataset_label,
        reference_group = reference_group,
        comparison_group = comparison_group,
        planned_covariates = if (length(planned_covariates)) {
            paste(planned_covariates, collapse = "; ")
        } else {
            "None"
        },
        model_formula = paste(deparse(model_formula), collapse = " "),
        n_samples = ncol(counts),
        n_taxa = n_taxa,
        structural_absence_available_n = sum(structural_available),
        structural_absence_available_rate = structural_rate,
        present_conditional_abundance_available_n =
            sum(abundance_available),
        present_conditional_abundance_available_rate = abundance_rate,
        combined_available_n = sum(combined_available),
        combined_available_rate = combined_rate,
        ra_bh_discovery_n = abundance_discovery_n,
        ra_bh_discovery_rate = abundance_discovery_rate,
        broadly_non_working = broadly_non_working,
        common_background_incompatible =
            common_background_incompatible,
        screen_passed = !broadly_non_working &
            !common_background_incompatible,
        dasra_version = dasra_version,
        elapsed_seconds = round(elapsed_seconds, 2),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    list(summary = summary_row, reasons = reason_table)
}

outer_cores <- min(3L, nrow(candidate_datasets))
screened <- parallel::mclapply(
    seq_len(nrow(candidate_datasets)),
    function(index) {
        screen_one_dataset(candidate_datasets[index, , drop = FALSE])
    },
    mc.cores = outer_cores,
    mc.preschedule = FALSE
)

screen_summary <- do.call(rbind, lapply(screened, `[[`, "summary"))
unavailability_reasons <- do.call(
    rbind, lapply(screened, `[[`, "reasons")
)

screening_directory <- file.path(analysis_directory, "screening")
dir.create(screening_directory, showWarnings = FALSE)
summary_file <- file.path(
    screening_directory, "dasra_candidate_screen.csv"
)
reason_file <- file.path(
    screening_directory, "dasra_unavailability_reasons.csv"
)

# Candidates removed after screening no longer have local input data, but they
# remain part of the prespecified screening denominator. Retain those historical
# rows while replacing results for every currently configured candidate.
if (file.exists(summary_file)) {
    previous_summary <- utils::read.csv(
        summary_file, stringsAsFactors = FALSE, check.names = FALSE
    )
    historical_summary <- previous_summary[
        !previous_summary$dataset_id %in% candidate_datasets$dataset_id,
        , drop = FALSE
    ]
    screen_summary <- rbind(historical_summary, screen_summary)
}
if (file.exists(reason_file)) {
    previous_reasons <- utils::read.csv(
        reason_file, stringsAsFactors = FALSE, check.names = FALSE
    )
    historical_reasons <- previous_reasons[
        !previous_reasons$dataset_id %in% candidate_datasets$dataset_id,
        , drop = FALSE
    ]
    unavailability_reasons <- rbind(
        historical_reasons, unavailability_reasons
    )
}
utils::write.csv(
    screen_summary,
    summary_file,
    row.names = FALSE,
    na = ""
)
utils::write.csv(
    unavailability_reasons,
    reason_file,
    row.names = FALSE,
    na = ""
)
