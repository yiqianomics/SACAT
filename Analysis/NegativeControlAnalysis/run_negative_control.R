options(stringsAsFactors = FALSE, warn = 1)

# The depth settings differ only in the group-specific median library size.
scenario_table <- data.frame(
    scenario = c("balanced", "fourfold"),
    median_H = c(4000, 1500),
    median_Case = c(4000, 6000),
    scenario_index = seq_len(2L),
    stringsAsFactors = FALSE
)

method_names <- c("DASRA", "ZINQ", "MaAsLin3")
reference_group <- "H"
comparison_group <- "Case"
depth_sdlog <- 0.45
minimum_depth <- 300L
maximum_depth <- 30000L
minimum_source_prevalence <- 0.10
minimum_mean_abundance <- 1e-5
required_dasra_version <- "0.6.0"
checkpoint_schema_version <- 5L
checkpoint_contract <- "balanced_fourfold_public_or_v3_dasra_0.6.0"

`%||%` <- function(x, y) {
    if (is.null(x) || !length(x)) y else x
}

validate_analysis_input <- function(object, input_path) {
    if (!is.list(object)) {
        stop("analysis_input.rds must contain a named list.", call. = FALSE)
    }
    probability <- object$probability %||% object$probabilities
    taxa <- object$evaluation_taxa %||% object$selected_taxa
    dataset_id <- as.character(object$dataset_id)
    dataset_label <- as.character(object$label %||% dataset_id)
    dataset_index <- suppressWarnings(as.integer(object$dataset_index))

    if (length(dataset_id) != 1L || !grepl("^[A-Za-z0-9._-]+$", dataset_id)) {
        stop("The analysis input must have one safe dataset_id.", call. = FALSE)
    }
    if (length(dataset_label) != 1L || !nzchar(dataset_label)) {
        stop("The analysis input must have one nonempty label.", call. = FALSE)
    }
    if (length(dataset_index) != 1L || !is.finite(dataset_index) ||
        dataset_index < 1L) {
        stop("The analysis input must have one positive dataset_index.",
             call. = FALSE)
    }
    if (is.data.frame(probability)) probability <- as.matrix(probability)
    if (!is.matrix(probability) || !is.numeric(probability) ||
        is.null(rownames(probability)) || is.null(colnames(probability))) {
        stop("probability must be a named numeric sample-by-taxon matrix.",
             call. = FALSE)
    }
    taxa <- as.character(taxa)
    if (length(taxa) != 30L || anyNA(taxa) || any(!nzchar(taxa)) ||
        anyDuplicated(taxa)) {
        stop("evaluation_taxa must contain exactly 30 unique taxa.",
             call. = FALSE)
    }
    if (!all(taxa %in% colnames(probability)) &&
        all(taxa %in% rownames(probability))) {
        probability <- t(probability)
    }
    if (!all(taxa %in% colnames(probability))) {
        stop("Every evaluation taxon must be a probability column.",
             call. = FALSE)
    }
    if ("Other_unmodeled" %in% taxa) {
        stop("Other_unmodeled cannot be an evaluation taxon.", call. = FALSE)
    }
    if (nrow(probability) != 200L) {
        stop("Each analysis requires exactly 200 samples.", call. = FALSE)
    }
    if (anyNA(probability) || any(!is.finite(probability)) ||
        any(probability < 0)) {
        stop("All probability entries must be finite and nonnegative.",
             call. = FALSE)
    }
    if (anyDuplicated(rownames(probability)) ||
        anyDuplicated(colnames(probability))) {
        stop("Probability row and column names must be unique.",
             call. = FALSE)
    }
    row_total <- rowSums(probability)
    if (any(!is.finite(row_total)) || any(row_total <= 0)) {
        stop("Every probability row must have positive mass.", call. = FALSE)
    }
    probability <- sweep(probability, 1L, row_total, "/")

    candidate_taxa <- setdiff(colnames(probability), "Other_unmodeled")
    source_prevalence <- colMeans(
        probability[, candidate_taxa, drop = FALSE] > 0
    )
    mean_abundance <- colMeans(
        probability[, candidate_taxa, drop = FALSE]
    )
    eligible <- candidate_taxa[
        source_prevalence >= minimum_source_prevalence &
        mean_abundance >= minimum_mean_abundance
    ]
    eligible <- eligible[order(
        -mean_abundance[eligible], eligible, method = "radix"
    )]
    if (length(eligible) < 30L || !identical(taxa, eligible[seq_len(30L)])) {
        stop("evaluation_taxa must be the source-defined top-30 panel.",
             call. = FALSE)
    }
    selection <- object$taxon_selection
    if (!is.list(selection) ||
        !identical(as.numeric(selection$source_prevalence_minimum),
                   minimum_source_prevalence) ||
        !identical(as.numeric(selection$mean_abundance_minimum),
                   minimum_mean_abundance) ||
        !identical(as.integer(selection$target_taxa), 30L)) {
        stop("The source panel-selection settings are incomplete.",
             call. = FALSE)
    }

    list(
        probability = probability,
        evaluation_taxa = taxa,
        dataset_id = dataset_id,
        label = dataset_label,
        dataset_index = dataset_index,
        taxon_selection = selection,
        input_path = normalizePath(input_path)
    )
}

make_seed <- function(analysis_input, replicate, stream, scenario_index = 0L) {
    modulus <- 2147483646
    value <- 900000 + analysis_input$dataset_index * 1000003 +
        as.integer(replicate) * 1009 + as.integer(stream) * 97 +
        as.integer(scenario_index) * 100003
    as.integer((value %% modulus) + 1)
}

draw_balanced_group <- function(analysis_input, replicate) {
    set.seed(make_seed(analysis_input, replicate, stream = 1L))
    labels <- rep(c(reference_group, comparison_group), each = 100L)
    factor(sample(labels, replace = FALSE),
           levels = c(reference_group, comparison_group))
}

draw_scenario_counts <- function(analysis_input, group, scenario, replicate) {
    setting <- scenario_table[
        scenario_table$scenario == scenario, , drop = FALSE
    ]
    if (nrow(setting) != 1L) {
        stop("Unknown depth scenario: ", scenario, call. = FALSE)
    }
    set.seed(make_seed(
        analysis_input, replicate, stream = 2L,
        scenario_index = setting$scenario_index
    ))
    target_median <- ifelse(
        group == reference_group, setting$median_H, setting$median_Case
    )
    library_size <- round(stats::rlnorm(
        length(group), meanlog = log(target_median), sdlog = depth_sdlog
    ))
    library_size <- as.integer(pmin(
        pmax(library_size, minimum_depth), maximum_depth
    ))

    tested_probability <- analysis_input$probability[,
        analysis_input$evaluation_taxa, drop = FALSE]
    other_probability <- pmax(0, 1 - rowSums(tested_probability))
    sampling_probability <- cbind(
        tested_probability,
        Other_unmodeled = other_probability
    )
    sampling_probability <- sweep(
        sampling_probability, 1L, rowSums(sampling_probability), "/"
    )
    counts <- vapply(seq_len(nrow(sampling_probability)), function(index) {
        as.integer(stats::rmultinom(
            1L, size = library_size[[index]],
            prob = sampling_probability[index, ]
        ))
    }, integer(ncol(sampling_probability)))
    rownames(counts) <- colnames(sampling_probability)
    colnames(counts) <- rownames(analysis_input$probability)
    storage.mode(counts) <- "integer"
    if (!identical(as.integer(colSums(counts)), library_size)) {
        stop("Generated counts do not match the drawn library sizes.",
             call. = FALSE)
    }

    metadata <- data.frame(
        group = factor(group, levels = c(reference_group, comparison_group)),
        library_size = library_size,
        stringsAsFactors = FALSE,
        row.names = colnames(counts)
    )
    metadata$log_library_size_z <- as.numeric(
        scale(log(metadata$library_size))
    )
    if (any(!is.finite(metadata$log_library_size_z))) {
        stop("Standardized log library size is not finite.", call. = FALSE)
    }
    list(counts = counts, metadata = metadata)
}

failure_result <- function(taxa, reason) {
    data.frame(
        taxon = taxa,
        native_p_value = NA_real_,
        available = FALSE,
        reason = reason,
        stringsAsFactors = FALSE
    )
}

empty_dasra_diagnostics <- function(taxa, reason) {
    data.frame(
        taxon = taxa,
        dasra_retained = FALSE,
        dasra_formed_structural_absence = FALSE,
        dasra_regular_structural_absence = FALSE,
        dasra_nonregular_structural_absence = FALSE,
        dasra_formed_relative_abundance = FALSE,
        dasra_formed_omnibus = FALSE,
        dasra_both_components_formed = FALSE,
        dasra_reason_structural_absence = reason,
        dasra_reason_relative_abundance = reason,
        dasra_warning_structural_absence = "",
        dasra_warning_relative_abundance = "",
        stringsAsFactors = FALSE
    )
}

sanitize_result <- function(result, taxa, method) {
    required <- c("taxon", "native_p_value", "available", "reason")
    if (!is.data.frame(result) || !all(required %in% names(result))) {
        result <- failure_result(taxa, "method returned an invalid result table")
    }
    original_taxa <- result$taxon
    result <- result[match(taxa, original_taxa), , drop = FALSE]
    result$taxon <- taxa
    missing_row <- is.na(match(taxa, original_taxa))
    result$native_p_value <- suppressWarnings(as.numeric(
        result$native_p_value
    ))
    result$available <- as.logical(result$available)
    result$available[is.na(result$available)] <- FALSE
    invalid <- missing_row | !is.finite(result$native_p_value) |
        result$native_p_value < 0 | result$native_p_value > 1
    result$available[invalid] <- FALSE
    result$reason <- as.character(result$reason)
    result$reason[missing_row] <- "taxon not returned"
    result$reason[invalid & !missing_row] <- "invalid p-value"
    result$reason[result$available &
                      (is.na(result$reason) | !nzchar(result$reason))] <-
        "available"
    result$p_value <- ifelse(result$available, result$native_p_value, 1)
    result$method <- method
    result[, c(
        "taxon", "method", "p_value", "native_p_value", "available",
        "reason"
    )]
}

run_with_warnings <- function(expression) {
    warnings <- character()
    value <- tryCatch(
        withCallingHandlers(
            expression,
            warning = function(warning) {
                warnings <<- c(warnings, conditionMessage(warning))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(error) error
    )
    list(value = value, warnings = unique(warnings))
}

run_dasra <- function(counts, metadata, taxa) {
    outcome <- run_with_warnings(DASRA::dasra(
        counts = counts,
        metadata = metadata,
        formula = ~ group,
        group = "group",
        library_size = "library_size",
        taxa_are_rows = TRUE,
        reference = reference_group,
        p_adjust_method = "BH",
        component = "all",
        full_output = FALSE,
        structural_conditional_present_starts = 1L
    ))
    if (inherits(outcome$value, "error")) {
        failure_reason <- paste(
            "method error:", conditionMessage(outcome$value)
        )
        return(list(
            result = failure_result(taxa, failure_reason),
            diagnostics = empty_dasra_diagnostics(taxa, failure_reason),
            warnings = outcome$warnings
        ))
    }
    fit <- outcome$value
    if (!is.data.frame(fit$results) || !is.data.frame(fit$diagnostics)) {
        return(list(
            result = failure_result(taxa, "DASRA returned incomplete output"),
            diagnostics = empty_dasra_diagnostics(
                taxa, "DASRA returned incomplete output"
            ),
            warnings = outcome$warnings
        ))
    }

    result_index <- match(taxa, fit$results$taxon)
    diagnostic_index <- match(taxa, fit$diagnostics$taxon)
    native_p <- rep(NA_real_, length(taxa))
    available <- rep(FALSE, length(taxa))
    reason <- rep("taxon not returned", length(taxa))
    diagnostics <- empty_dasra_diagnostics(taxa, "taxon not returned")
    has_result <- !is.na(result_index)
    has_diagnostic <- !is.na(diagnostic_index)
    if ("p_omnibus" %in% names(fit$results)) {
        native_p[has_result] <- suppressWarnings(as.numeric(
            fit$results$p_omnibus[result_index[has_result]]
        ))
    }
    required_diagnostics <- c("retained", "formed_omnibus")
    if (all(required_diagnostics %in% names(fit$diagnostics))) {
        diagnostic_ok <- rep(FALSE, length(taxa))
        matched_diagnostics <- fit$diagnostics[
            diagnostic_index[has_diagnostic], , drop = FALSE
        ]
        retained_values <- as.logical(matched_diagnostics$retained)
        omnibus_values <- as.logical(matched_diagnostics$formed_omnibus)
        retained_values[is.na(retained_values)] <- FALSE
        omnibus_values[is.na(omnibus_values)] <- FALSE
        diagnostic_ok[has_diagnostic] <-
            retained_values & omnibus_values
        available <- has_result & has_diagnostic & diagnostic_ok &
            is.finite(native_p)
        reason[has_result & has_diagnostic & diagnostic_ok] <- "available"
        retained <- rep(FALSE, length(taxa))
        formed_omnibus <- rep(FALSE, length(taxa))
        retained[has_diagnostic] <- retained_values
        formed_omnibus[has_diagnostic] <- omnibus_values
        reason[has_diagnostic & !retained] <- "taxon was not retained"
        reason[has_diagnostic & retained & !formed_omnibus] <-
            "no DASRA component was formed"
        reason[has_result & has_diagnostic & diagnostic_ok &
                   !is.finite(native_p)] <- "invalid omnibus p-value"
    } else {
        reason[has_result] <- "DASRA diagnostics were incomplete"
    }

    if (any(has_diagnostic)) {
        source <- fit$diagnostics[
            diagnostic_index[has_diagnostic], , drop = FALSE
        ]
        target <- which(has_diagnostic)
        logical_fields <- c(
            dasra_retained = "retained",
            dasra_formed_structural_absence = "formed_structural_absence",
            dasra_regular_structural_absence =
                "regular_structural_absence",
            dasra_nonregular_structural_absence =
                "nonregular_structural_absence",
            dasra_formed_relative_abundance = "formed_relative_abundance",
            dasra_formed_omnibus = "formed_omnibus"
        )
        for (field in names(logical_fields)) {
            source_field <- logical_fields[[field]]
            if (source_field %in% names(source)) {
                values <- as.logical(source[[source_field]])
                values[is.na(values)] <- FALSE
                diagnostics[[field]][target] <- values
            }
        }
        diagnostics$dasra_both_components_formed <-
            diagnostics$dasra_retained &
            diagnostics$dasra_formed_structural_absence &
            diagnostics$dasra_formed_relative_abundance
        reason_fields <- c(
            dasra_reason_structural_absence = "reason_structural_absence",
            dasra_reason_relative_abundance = "reason_relative_abundance",
            dasra_warning_structural_absence =
                "warning_structural_absence",
            dasra_warning_relative_abundance =
                "warning_relative_abundance"
        )
        for (field in names(reason_fields)) {
            source_field <- reason_fields[[field]]
            if (source_field %in% names(source)) {
                diagnostics[[field]][target] <- as.character(
                    source[[source_field]]
                )
            }
        }
    }
    list(
        result = data.frame(
            taxon = taxa,
            native_p_value = native_p,
            available = available,
            reason = reason,
            stringsAsFactors = FALSE
        ),
        diagnostics = diagnostics,
        warnings = outcome$warnings
    )
}

run_zinq <- function(counts, metadata, taxa) {
    relative_abundance <- sweep(
        t(counts[taxa, , drop = FALSE]), 1L, metadata$library_size, "/"
    )
    group_binary <- as.integer(metadata$group == comparison_group)
    covariate_matrix <- stats::model.matrix(
        ~ log_library_size_z, metadata
    )[, -1L, drop = FALSE]
    colnames(covariate_matrix) <- "covariate1"
    model_formula <- y ~ group + covariate1
    taus <- c(0.25, 0.50, 0.75)
    native_p <- rep(NA_real_, length(taxa))
    available <- rep(FALSE, length(taxa))
    reason <- rep("model did not return a valid p-value", length(taxa))
    warnings <- character()

    for (index in seq_along(taxa)) {
        taxon <- taxa[[index]]
        model_data <- data.frame(
            y = as.numeric(relative_abundance[, taxon]),
            group = group_binary,
            covariate1 = covariate_matrix[, 1L],
            stringsAsFactors = FALSE
        )
        fit_outcome <- run_with_warnings(ZINQ::ZINQ_tests(
            formula.logistic = model_formula,
            formula.quantile = model_formula,
            C = "group",
            y_CorD = "C",
            data = model_data,
            taus = taus,
            seed = 2026L + index
        ))
        warnings <- c(warnings, fit_outcome$warnings)
        if (inherits(fit_outcome$value, "error")) {
            reason[[index]] <- paste(
                "model error:", conditionMessage(fit_outcome$value)
            )
            next
        }
        combination <- run_with_warnings(ZINQ::ZINQ_combination(
            fit_outcome$value, method = "Cauchy", taus = taus
        ))
        warnings <- c(warnings, combination$warnings)
        if (inherits(combination$value, "error")) {
            reason[[index]] <- paste(
                "combination error:", conditionMessage(combination$value)
            )
            next
        }
        value <- suppressWarnings(as.numeric(combination$value))
        if (length(value) == 1L && is.finite(value) &&
            value >= 0 && value <= 1) {
            native_p[[index]] <- value
            available[[index]] <- TRUE
            reason[[index]] <- "available"
        } else {
            reason[[index]] <- "invalid native Cauchy p-value"
        }
    }
    list(
        result = data.frame(
            taxon = taxa,
            native_p_value = native_p,
            available = available,
            reason = reason,
            stringsAsFactors = FALSE
        ),
        warnings = unique(warnings)
    )
}

extract_maaslin_joint <- function(fit, taxa) {
    extract_table <- function(table) {
        required <- c(
            "feature", "metadata", "value", "pval_individual", "pval_joint"
        )
        if (!is.data.frame(table) || !all(required %in% names(table))) {
            return(list(
                individual = rep(NA_real_, length(taxa)),
                joint = rep(NA_real_, length(taxa)),
                unique = rep(FALSE, length(taxa))
            ))
        }
        table <- table[
            table$metadata == "group" & table$value == comparison_group,
            , drop = FALSE
        ]
        individual <- joint <- rep(NA_real_, length(taxa))
        unique <- rep(FALSE, length(taxa))
        for (index in seq_along(taxa)) {
            rows <- which(table$feature == taxa[[index]])
            if (length(rows) == 1L) {
                individual[[index]] <- suppressWarnings(as.numeric(
                    table$pval_individual[rows]
                ))
                joint[[index]] <- suppressWarnings(as.numeric(
                    table$pval_joint[rows]
                ))
                unique[[index]] <- TRUE
            }
        }
        list(individual = individual, joint = joint, unique = unique)
    }

    prevalence <- extract_table(fit$fit_data_prevalence$results)
    abundance <- extract_table(fit$fit_data_abundance$results)
    agreement <- is.finite(prevalence$joint) & is.finite(abundance$joint) &
        abs(prevalence$joint - abundance$joint) <= 1e-12
    available <- prevalence$unique & abundance$unique &
        is.finite(prevalence$individual) & is.finite(abundance$individual) &
        agreement
    data.frame(
        taxon = taxa,
        native_p_value = prevalence$joint,
        available = available,
        reason = ifelse(
            available,
            "available",
            "both component models and a consistent joint p-value were required"
        ),
        stringsAsFactors = FALSE
    )
}

run_maaslin3 <- function(counts, metadata, taxa, work_directory) {
    dir.create(work_directory, recursive = TRUE, showWarnings = FALSE)
    output_directory <- file.path(work_directory, "maaslin3")
    if (dir.exists(output_directory)) {
        unlink(output_directory, recursive = TRUE, force = TRUE)
    }
    on.exit(try(maaslin3::maaslin_log_reset(), silent = TRUE), add = TRUE)
    on.exit(unlink(work_directory, recursive = TRUE, force = TRUE), add = TRUE)
    relative_abundance <- sweep(
        t(counts[taxa, , drop = FALSE]),
        1L,
        metadata$library_size,
        "/"
    )
    outcome <- run_with_warnings(maaslin3::maaslin3(
        input_data = as.data.frame(relative_abundance, check.names = FALSE),
        input_metadata = metadata,
        output = output_directory,
        formula = ~ log_library_size_z + group,
        min_abundance = 0,
        min_prevalence = 0,
        max_prevalence = 1.01,
        zero_threshold = 0,
        min_variance = 0,
        normalization = "NONE",
        transform = "LOG",
        correction = "BH",
        standardize = TRUE,
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        warn_prevalence = TRUE,
        augment = TRUE,
        cores = 1,
        max_significance = 1,
        plot_summary_plot = FALSE,
        plot_associations = FALSE,
        save_models = FALSE,
        save_plots_rds = FALSE,
        verbosity = "WARN",
        reference = paste0("group,", reference_group)
    ))
    try(maaslin3::maaslin_log_reset(), silent = TRUE)
    if (inherits(outcome$value, "error")) {
        result <- failure_result(
            taxa, paste("method error:", conditionMessage(outcome$value))
        )
    } else {
        result <- tryCatch(
            extract_maaslin_joint(outcome$value, taxa),
            error = function(error) failure_result(
                taxa, paste("result extraction error:", conditionMessage(error))
            )
        )
    }
    list(result = result, warnings = outcome$warnings)
}

run_timed_method <- function(method, runner, taxa) {
    started <- proc.time()[["elapsed"]]
    outcome <- tryCatch(
        runner(),
        error = function(error) list(
            result = failure_result(
                taxa, paste("method error:", conditionMessage(error))
            ),
            warnings = character()
        )
    )
    list(
        result = sanitize_result(outcome$result, taxa, method),
        diagnostics = outcome$diagnostics %||% NULL,
        elapsed_seconds = proc.time()[["elapsed"]] - started,
        warnings = paste(unique(outcome$warnings), collapse = " | ")
    )
}

make_taxon_diagnostics <- function(counts, metadata, taxa) {
    tested <- counts[taxa, , drop = FALSE]
    is_H <- metadata$group == reference_group
    is_Case <- metadata$group == comparison_group
    data.frame(
        taxon = taxa,
        observed_prevalence = rowMeans(tested > 0),
        prevalence_H = rowMeans(tested[, is_H, drop = FALSE] > 0),
        prevalence_Case = rowMeans(tested[, is_Case, drop = FALSE] > 0),
        positive_samples_H = rowSums(tested[, is_H, drop = FALSE] > 0),
        positive_samples_Case = rowSums(
            tested[, is_Case, drop = FALSE] > 0
        ),
        total_count = rowSums(tested),
        stringsAsFactors = FALSE
    )
}

checkpoint_path <- function(result_directory, scenario, replicate) {
    file.path(
        result_directory, ".checkpoints", scenario,
        sprintf("rep_%04d.rds", replicate)
    )
}

atomic_save_rds <- function(object, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(pattern = basename(path), tmpdir = dirname(path))
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(object, temporary, version = 3)
    if (!file.rename(temporary, path)) {
        if (!file.copy(temporary, path, overwrite = TRUE)) {
            stop("Could not write checkpoint: ", path, call. = FALSE)
        }
        unlink(temporary)
    }
    invisible(path)
}

run_unit <- function(analysis_input, scenario, replicate, result_directory,
                     work_root) {
    unit_started <- proc.time()[["elapsed"]]
    group <- draw_balanced_group(analysis_input, replicate)
    generated <- draw_scenario_counts(analysis_input, group, scenario, replicate)
    taxa <- analysis_input$evaluation_taxa
    counts <- generated$counts[taxa, , drop = FALSE]
    metadata <- generated$metadata
    work_directory <- file.path(
        work_root, scenario, sprintf("rep_%04d", replicate)
    )

    dasra <- run_timed_method(
        "DASRA", function() run_dasra(counts, metadata, taxa), taxa
    )
    zinq <- run_timed_method(
        "ZINQ", function() run_zinq(counts, metadata, taxa), taxa
    )
    maaslin <- run_timed_method(
        "MaAsLin3",
        function() run_maaslin3(counts, metadata, taxa, work_directory),
        taxa
    )
    method_results <- list(dasra, zinq, maaslin)
    taxon_pvalues <- do.call(rbind, lapply(method_results, `[[`, "result"))
    taxon_pvalues$dataset <- analysis_input$dataset_id
    taxon_pvalues$dataset_label <- analysis_input$label
    taxon_pvalues$replicate <- as.integer(replicate)
    taxon_pvalues$scenario <- scenario
    taxon_pvalues$status <- ifelse(
        taxon_pvalues$available, "available", "unavailable"
    )
    taxon_pvalues <- taxon_pvalues[, c(
        "dataset", "dataset_label", "replicate", "scenario", "method",
        "taxon", "p_value", "native_p_value", "available", "status",
        "reason"
    )]

    replicate_metrics <- do.call(rbind, lapply(method_results, function(item) {
        result <- item$result
        data.frame(
            dataset = analysis_input$dataset_id,
            dataset_label = analysis_input$label,
            replicate = as.integer(replicate),
            scenario = scenario,
            method = unique(result$method),
            n_taxa = length(taxa),
            n_available = sum(result$available),
            availability_rate = mean(result$available),
            n_rejections = sum(result$p_value <= 0.05),
            type1_error = mean(result$p_value <= 0.05),
            stringsAsFactors = FALSE
        )
    }))
    timing <- do.call(rbind, Map(function(item, method) {
        data.frame(
            dataset = analysis_input$dataset_id,
            dataset_label = analysis_input$label,
            replicate = as.integer(replicate),
            scenario = scenario,
            method = method,
            elapsed_seconds = item$elapsed_seconds,
            warnings = item$warnings,
            stringsAsFactors = FALSE
        )
    }, method_results, method_names))

    diagnostics <- make_taxon_diagnostics(counts, metadata, taxa)
    dasra_diagnostics <- dasra$diagnostics %||%
        empty_dasra_diagnostics(taxa, "DASRA diagnostics were unavailable")
    dasra_diagnostics <- dasra_diagnostics[
        match(taxa, dasra_diagnostics$taxon), , drop = FALSE
    ]
    if (anyNA(dasra_diagnostics$taxon)) {
        dasra_diagnostics <- empty_dasra_diagnostics(
            taxa, "DASRA diagnostics were incomplete"
        )
    }
    diagnostics <- cbind(
        diagnostics,
        dasra_diagnostics[, setdiff(names(dasra_diagnostics), "taxon"),
                          drop = FALSE]
    )
    diagnostics$dataset <- analysis_input$dataset_id
    diagnostics$dataset_label <- analysis_input$label
    diagnostics$replicate <- as.integer(replicate)
    diagnostics$scenario <- scenario
    diagnostics <- diagnostics[, c(
        "dataset", "dataset_label", "replicate", "scenario", "taxon",
        "observed_prevalence", "prevalence_H", "prevalence_Case",
        "positive_samples_H", "positive_samples_Case", "total_count",
        "dasra_retained", "dasra_formed_structural_absence",
        "dasra_regular_structural_absence",
        "dasra_nonregular_structural_absence",
        "dasra_formed_relative_abundance", "dasra_formed_omnibus",
        "dasra_both_components_formed",
        "dasra_reason_structural_absence",
        "dasra_reason_relative_abundance",
        "dasra_warning_structural_absence",
        "dasra_warning_relative_abundance"
    )]
    depth_diagnostics <- data.frame(
        dataset = analysis_input$dataset_id,
        dataset_label = analysis_input$label,
        replicate = as.integer(replicate),
        scenario = scenario,
        n_H = sum(metadata$group == reference_group),
        n_Case = sum(metadata$group == comparison_group),
        median_depth_H = stats::median(
            metadata$library_size[metadata$group == reference_group]
        ),
        median_depth_Case = stats::median(
            metadata$library_size[metadata$group == comparison_group]
        ),
        depth_ratio_Case_to_H = stats::median(
            metadata$library_size[metadata$group == comparison_group]
        ) / stats::median(
            metadata$library_size[metadata$group == reference_group]
        ),
        minimum_depth = min(metadata$library_size),
        maximum_depth = max(metadata$library_size),
        group_seed = make_seed(analysis_input, replicate, stream = 1L),
        count_seed = make_seed(
            analysis_input, replicate, stream = 2L,
            scenario_index = scenario_table$scenario_index[
                match(scenario, scenario_table$scenario)
            ]
        ),
        elapsed_seconds = proc.time()[["elapsed"]] - unit_started,
        stringsAsFactors = FALSE
    )

    atomic_save_rds(list(
        schema_version = checkpoint_schema_version,
        analysis_contract = checkpoint_contract,
        DASRA_version = as.character(utils::packageVersion("DASRA")),
        dataset_id = analysis_input$dataset_id,
        evaluation_taxa = taxa,
        replicate = as.integer(replicate),
        scenario = scenario,
        taxon_pvalues = taxon_pvalues,
        replicate_metrics = replicate_metrics,
        timing = timing,
        diagnostics = diagnostics,
        depth_diagnostics = depth_diagnostics
    ), checkpoint_path(result_directory, scenario, replicate))
    TRUE
}

checkpoint_is_current <- function(path, analysis_input, scenario, replicate) {
    if (!file.exists(path)) return(FALSE)
    object <- tryCatch(readRDS(path), error = function(error) NULL)
    !is.null(object) &&
        identical(object$schema_version, checkpoint_schema_version) &&
        identical(object$analysis_contract, checkpoint_contract) &&
        identical(object$DASRA_version, required_dasra_version) &&
        identical(object$dataset_id, analysis_input$dataset_id) &&
        identical(object$evaluation_taxa, analysis_input$evaluation_taxa) &&
        identical(object$scenario, scenario) &&
        identical(object$replicate, as.integer(replicate)) &&
        is.data.frame(object$taxon_pvalues) &&
        nrow(object$taxon_pvalues) == 90L &&
        is.data.frame(object$replicate_metrics) &&
        nrow(object$replicate_metrics) == 3L &&
        is.data.frame(object$diagnostics) &&
        nrow(object$diagnostics) == 30L &&
        is.data.frame(object$depth_diagnostics) &&
        nrow(object$depth_diagnostics) == 1L
}

atomic_write_csv <- function(table, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
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

aggregate_checkpoints <- function(result_directory, analysis_input,
                                  expected_replicates = 100L) {
    units <- expand.grid(
        replicate = seq_len(expected_replicates),
        scenario = scenario_table$scenario,
        stringsAsFactors = FALSE
    )
    paths <- mapply(
        checkpoint_path,
        scenario = units$scenario,
        replicate = units$replicate,
        MoreArgs = list(result_directory = result_directory),
        USE.NAMES = FALSE
    )
    current <- mapply(
        checkpoint_is_current,
        path = paths,
        scenario = units$scenario,
        replicate = units$replicate,
        MoreArgs = list(analysis_input = analysis_input),
        USE.NAMES = FALSE
    )
    if (!all(current)) {
        stop(sum(!current), " scenario units are incomplete.", call. = FALSE)
    }
    checkpoints <- lapply(paths, readRDS)
    formal_diagnostics <- do.call(rbind, lapply(
        checkpoints, `[[`, "diagnostics"
    ))
    fourfold_diagnostics <- formal_diagnostics[
        formal_diagnostics$scenario == "fourfold", , drop = FALSE
    ]
    minimum_both <- pmin(
        fourfold_diagnostics$positive_samples_H,
        fourfold_diagnostics$positive_samples_Case
    )
    below_18_by_taxon <- tapply(
        minimum_both < 18L, fourfold_diagnostics$taxon, sum
    )
    components <- c(
        taxon_pvalues = "taxon_pvalues.csv",
        replicate_metrics = "replicate_metrics.csv",
        diagnostics = "diagnostics.csv",
        depth_diagnostics = "depth_diagnostics.csv"
    )
    for (component in names(components)) {
        table <- do.call(rbind, lapply(checkpoints, `[[`, component))
        if (component == "depth_diagnostics") {
            table$elapsed_seconds <- NULL
        }
        ordering <- intersect(
            c("replicate", "scenario", "method", "taxon"), names(table)
        )
        order_values <- lapply(ordering, function(variable) {
            values <- table[[variable]]
            if (variable == "scenario") {
                match(values, scenario_table$scenario)
            } else if (variable == "method") {
                match(values, method_names)
            } else {
                values
            }
        })
        table <- table[do.call(order, order_values), , drop = FALSE]
        rownames(table) <- NULL
        atomic_write_csv(table, file.path(
            result_directory, components[[component]]
        ))
    }

    manifest <- data.frame(
        dataset = analysis_input$dataset_id,
        dataset_label = analysis_input$label,
        n_samples = nrow(analysis_input$probability),
        group_size = nrow(analysis_input$probability) / 2L,
        n_taxa = length(analysis_input$evaluation_taxa),
        replicates = expected_replicates,
        scenarios = paste(scenario_table$scenario, collapse = " | "),
        depth_medians_H = paste(scenario_table$median_H, collapse = " | "),
        depth_medians_Case = paste(
            scenario_table$median_Case, collapse = " | "
        ),
        depth_sdlog = depth_sdlog,
        panel_source_prevalence_minimum = minimum_source_prevalence,
        panel_mean_abundance_minimum = minimum_mean_abundance,
        formal_fourfold_minimum_positive_H = min(
            fourfold_diagnostics$positive_samples_H
        ),
        formal_fourfold_minimum_positive_Case = min(
            fourfold_diagnostics$positive_samples_Case
        ),
        formal_fourfold_minimum_positive_both = min(minimum_both),
        formal_fourfold_fraction_records_at_least_20 =
            mean(minimum_both >= 20L),
        formal_fourfold_fraction_records_at_least_18 =
            mean(minimum_both >= 18L),
        formal_fourfold_max_taxon_replicates_below_18 =
            max(below_18_by_taxon),
        sampling_bin = "Other_unmodeled excluded from method inputs",
        DASRA_omnibus_contract = "p_omnibus with formed_omnibus",
        MaAsLin3_input = "taxon relative abundance with normalization NONE",
        DASRA_version = as.character(utils::packageVersion("DASRA")),
        ZINQ_version = as.character(utils::packageVersion("ZINQ")),
        maaslin3_version = as.character(utils::packageVersion("maaslin3")),
        R_version = paste(R.version$major, R.version$minor, sep = "."),
        stringsAsFactors = FALSE
    )
    atomic_write_csv(manifest, file.path(result_directory, "run_manifest.csv"))
    invisible(manifest)
}

run_negative_control <- function(dataset_directory, workers = 7L,
                                 replicates = 100L, force = FALSE) {
    workers <- as.integer(workers)
    replicates <- as.integer(replicates)
    if (!identical(workers, 7L)) {
        stop("This analysis is configured for seven PSOCK workers.",
             call. = FALSE)
    }
    if (!identical(replicates, 100L)) {
        stop("This analysis requires exactly 100 randomizations.",
             call. = FALSE)
    }
    required_packages <- c("DASRA", "ZINQ", "maaslin3")
    missing_packages <- required_packages[!vapply(
        required_packages, requireNamespace, logical(1), quietly = TRUE
    )]
    if (length(missing_packages)) {
        stop(
            "Required R packages are unavailable: ",
            paste(missing_packages, collapse = ", "),
            call. = FALSE
        )
    }
    package_versions <- vapply(
        required_packages,
        function(package) as.character(getNamespaceVersion(package)),
        character(1)
    )
    if (!identical(
        package_versions[["DASRA"]],
        required_dasra_version
    )) {
        stop(
            "This analysis requires DASRA ", required_dasra_version, ".",
            call. = FALSE
        )
    }

    dataset_directory <- normalizePath(dataset_directory)
    runner_path <- normalizePath(file.path(dataset_directory, "..",
                                           "run_negative_control.R"))
    input_path <- file.path(dataset_directory, "analysis_input.rds")
    if (!file.exists(input_path)) {
        stop("Missing analysis_input.rds in ", dataset_directory,
             call. = FALSE)
    }
    analysis_input <- validate_analysis_input(readRDS(input_path), input_path)
    if (!identical(analysis_input$dataset_id, basename(dataset_directory))) {
        stop("The folder name and analysis input dataset_id must match.",
             call. = FALSE)
    }
    result_directory <- file.path(dataset_directory, "results")
    work_root <- file.path(dataset_directory, ".work")
    dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
    dir.create(work_root, recursive = TRUE, showWarnings = FALSE)

    units <- expand.grid(
        replicate = seq_len(replicates),
        scenario = scenario_table$scenario,
        stringsAsFactors = FALSE
    )
    units$scenario <- factor(units$scenario, levels = scenario_table$scenario)
    units <- units[order(units$replicate, units$scenario), , drop = FALSE]
    units$scenario <- as.character(units$scenario)
    paths <- mapply(
        checkpoint_path,
        scenario = units$scenario,
        replicate = units$replicate,
        MoreArgs = list(result_directory = result_directory),
        USE.NAMES = FALSE
    )
    if (force) {
        current <- rep(FALSE, nrow(units))
    } else {
        current <- mapply(
            checkpoint_is_current,
            path = paths,
            scenario = units$scenario,
            replicate = units$replicate,
            MoreArgs = list(analysis_input = analysis_input),
            USE.NAMES = FALSE
        )
    }
    remaining <- units[!current, , drop = FALSE]

    message(sprintf(
        "%s: 200 samples (100 per group), 30 taxa, %d units remaining.",
        analysis_input$label, nrow(remaining)
    ))
    thread_environment <- c(
        OMP_NUM_THREADS = "1",
        OMP_THREAD_LIMIT = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1"
    )
    do.call(Sys.setenv, as.list(thread_environment))
    worker_library_paths <- normalizePath(
        .libPaths(), winslash = "/", mustWork = TRUE
    )

    if (nrow(remaining)) {
        cluster <- parallel::makePSOCKcluster(
            min(workers, nrow(remaining)),
            outfile = "",
            rscript = file.path(R.home("bin"), "Rscript")
        )
        on.exit(try(parallel::stopCluster(cluster), silent = TRUE), add = TRUE)
        initialized <- parallel::clusterCall(
            cluster,
            function(runner_path, input_path, result_directory, work_root,
                     thread_environment, library_paths, package_versions) {
                .libPaths(unique(c(library_paths, .libPaths())))
                do.call(Sys.setenv, as.list(thread_environment))
                worker_versions <- vapply(
                    names(package_versions),
                    function(package) {
                        if (!requireNamespace(package, quietly = TRUE)) {
                            stop("A required package is unavailable on a worker.",
                                 call. = FALSE)
                        }
                        as.character(getNamespaceVersion(package))
                    },
                    character(1)
                )
                if (!identical(worker_versions, package_versions)) {
                    stop("Worker package versions differ from the main process.",
                         call. = FALSE)
                }
                sys.source(runner_path, envir = .GlobalEnv)
                validator <- get(
                    "validate_analysis_input", envir = .GlobalEnv
                )
                assign(
                    ".negative_control_input",
                    validator(readRDS(input_path), input_path),
                    envir = .GlobalEnv
                )
                assign(
                    ".negative_control_results", result_directory,
                    envir = .GlobalEnv
                )
                assign(
                    ".negative_control_work", work_root,
                    envir = .GlobalEnv
                )
                worker_versions
            },
            runner_path = runner_path,
            input_path = input_path,
            result_directory = result_directory,
            work_root = work_root,
            thread_environment = thread_environment,
            library_paths = worker_library_paths,
            package_versions = package_versions
        )
        if (!all(vapply(
            initialized, identical, logical(1), package_versions
        ))) {
            stop("One or more workers could not be initialized.",
                 call. = FALSE)
        }
        unit_list <- lapply(seq_len(nrow(remaining)), function(index) {
            list(
                replicate = as.integer(remaining$replicate[[index]]),
                scenario = as.character(remaining$scenario[[index]])
            )
        })
        outcomes <- parallel::parLapplyLB(cluster, unit_list, function(unit) {
            tryCatch({
                unit_runner <- get("run_unit", envir = .GlobalEnv)
                unit_runner(
                    analysis_input = get(
                        ".negative_control_input", envir = .GlobalEnv
                    ),
                    scenario = unit$scenario,
                    replicate = unit$replicate,
                    result_directory = get(
                        ".negative_control_results", envir = .GlobalEnv
                    ),
                    work_root = get(
                        ".negative_control_work", envir = .GlobalEnv
                    )
                )
                list(ok = TRUE, error = "")
            }, error = function(error) {
                list(ok = FALSE, error = conditionMessage(error))
            })
        })
        parallel::stopCluster(cluster)
        cluster <- NULL
        failed <- !vapply(outcomes, `[[`, logical(1), "ok")
        if (any(failed)) {
            errors <- unique(vapply(
                outcomes[failed], `[[`, character(1), "error"
            ))
            stop(
                sum(failed), " units failed before checkpointing: ",
                paste(errors, collapse = " | "),
                call. = FALSE
            )
        }
    }

    manifest <- aggregate_checkpoints(
        result_directory,
        analysis_input,
        expected_replicates = replicates
    )
    taxon_results <- utils::read.csv(
        file.path(result_directory, "taxon_pvalues.csv"),
        check.names = FALSE
    )
    replicate_results <- utils::read.csv(
        file.path(result_directory, "replicate_metrics.csv"),
        check.names = FALSE
    )
    diagnostics <- utils::read.csv(
        file.path(result_directory, "diagnostics.csv"),
        check.names = FALSE
    )
    depth_results <- utils::read.csv(
        file.path(result_directory, "depth_diagnostics.csv"),
        check.names = FALSE
    )
    expected_units <- replicates * nrow(scenario_table)
    if (nrow(taxon_results) !=
            expected_units * length(method_names) *
                length(analysis_input$evaluation_taxa) ||
        nrow(replicate_results) != expected_units * length(method_names) ||
        nrow(diagnostics) !=
            expected_units * length(analysis_input$evaluation_taxa) ||
        nrow(depth_results) != expected_units ||
        !setequal(unique(taxon_results$scenario), scenario_table$scenario) ||
        !setequal(unique(replicate_results$scenario),
                  scenario_table$scenario) ||
        !setequal(unique(replicate_results$method), method_names)) {
        stop("The completed result tables have unexpected dimensions.",
             call. = FALSE)
    }
    required_diagnostics <- c(
        "dasra_regular_structural_absence",
        "dasra_nonregular_structural_absence",
        "dasra_warning_structural_absence",
        "dasra_warning_relative_abundance"
    )
    if (!all(required_diagnostics %in% names(diagnostics))) {
        stop("The completed DASRA diagnostics are incomplete.",
             call. = FALSE)
    }
    unlink(file.path(result_directory, ".checkpoints"),
           recursive = TRUE, force = TRUE)
    unlink(work_root, recursive = TRUE, force = TRUE)
    message(analysis_input$label, ": complete and validated.")
    invisible(manifest)
}
