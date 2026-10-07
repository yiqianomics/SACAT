#!/usr/bin/env Rscript

# metaSPARSim Gamma abundances and read sampling with known component effects.
# Usage: Rscript run_simulation.R design | preflight | smoke | replicate ID | summarize

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

.metaspar_script_dir <- local({
    source_files <- vapply(sys.frames(), function(frame) {
        if (is.null(frame$ofile)) "" else as.character(frame$ofile)
    }, character(1L))
    source_files <- source_files[nzchar(source_files)]
    if (length(source_files)) {
        dirname(normalizePath(tail(source_files, 1L), mustWork = TRUE))
    } else {
        file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
        if (length(file_argument)) {
            dirname(normalizePath(sub("^--file=", "", file_argument[1L]),
                                  mustWork = TRUE))
        } else getwd()
    }
})

CONFIG <- list(
    root = Sys.getenv("METASPAR_ROOT",
                      file.path(dirname(.metaspar_script_dir), "output")),
    n_taxa = 100L,
    replications = 100L,
    null_replications = 500L,
    setting_blocks = 2L,
    base_seed = 202609160L,
    simulation_seed = 202611460L,
    alpha = 0.05,
    zinq_taus = c(0.25, 0.50, 0.75),
    maaslin3_warn_prevalence = TRUE,
    corncob_robust = TRUE,
    metagenomeseq_maxit = 10L
)

REQUIRED_PACKAGES <- c(
    "metaSPARSim", "SACAT", "ZINQ", "maaslin3", "edgeR", "DESeq2",
    "ANCOMBC", "MicrobiomeStat", "corncob", "metagenomeSeq", "Biobase",
    "limma", "data.table", "TreeSummarizedExperiment",
    "SummarizedExperiment", "S4Vectors", "statmod"
)

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) y else x
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
messagef <- function(fmt, ...) message(sprintf(fmt, ...))

clamp <- function(x, lower, upper) pmin(pmax(x, lower), upper)

safe_dir_create <- function(path) {
    if (!dir.exists(path)) {
        ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (!ok && !dir.exists(path)) stopf("Could not create directory: %s", path)
    }
    invisible(path)
}
collapse_text <- function(x) {
    x <- unique(trimws(as.character(x)))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) NA_character_ else paste(x, collapse = " | ")
}

safe_capture <- function(expr) {
    warnings <- character()
    started <- proc.time()[["elapsed"]]
    value <- tryCatch(
        withCallingHandlers(
            expr,
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) structure(
            list(message = conditionMessage(e)), class = "sacat_sim_error"
        )
    )
    elapsed <- proc.time()[["elapsed"]] - started
    list(value = value, warnings = unique(warnings), elapsed = elapsed)
}

is_sim_error <- function(x) inherits(x, "sacat_sim_error")
cauchy_combine <- function(p_values, weights = NULL) {
    p_values <- as.numeric(p_values)
    keep <- is.finite(p_values) & p_values >= 0 & p_values <= 1
    if (!any(keep)) return(NA_real_)
    p_values <- clamp(p_values[keep], 1e-15, 1 - 1e-15)
    if (is.null(weights)) {
        weights <- rep(1 / length(p_values), length(p_values))
    } else {
        weights <- as.numeric(weights)[keep]
        if (length(weights) != length(p_values) || any(!is.finite(weights)) ||
            any(weights < 0) || sum(weights) <= 0) {
            return(NA_real_)
        }
        weights <- weights / sum(weights)
    }
    statistic <- sum(weights * tan((0.5 - p_values) * pi))
    clamp(1 - stats::pcauchy(statistic), 0, 1)
}

first_existing_column <- function(data, candidates) {
    candidates <- candidates[candidates %in% names(data)]
    if (!length(candidates)) return(NULL)
    candidates[[1L]]
}

as_numeric_column <- function(data, candidates) {
    column <- first_existing_column(data, candidates)
    if (is.null(column)) return(rep(NA_real_, nrow(data)))
    suppressWarnings(as.numeric(data[[column]]))
}

method_result_template <- function(taxa, method, component, status,
                                   elapsed = NA_real_) {
    data.frame(
        taxon = taxa,
        method = method,
        component = component,
        p_value = NA_real_,
        p_value_unfiltered = NA_real_,
        package_q_value = NA_real_,
        estimate = NA_real_,
        standard_error = NA_real_,
        statistic = NA_real_,
        available = FALSE,
        status = status,
        warning = NA_character_,
        method_detail = NA_character_,
        abundance_induced_flag = NA,
        passed_sensitivity = NA,
        runtime_seconds = elapsed,
        runtime_seconds_total_setting = NA_real_,
        p_tau_025 = NA_real_,
        p_tau_050 = NA_real_,
        p_tau_075 = NA_real_,
        stringsAsFactors = FALSE
    )
}

valid_p <- function(x) is.finite(x) & x >= 0 & x <= 1

design_formula <- function(metadata, include_depth = FALSE, response = NULL,
                           group = "group", group_last = FALSE) {
    terms <- c(group, if ("z" %in% names(metadata)) "z")
    if (group_last) terms <- rev(terms)
    if (include_depth) terms <- c(terms, "log_depth")
    stats::reformulate(terms, response = response)
}

run_sacat <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture(
        SACAT::sacat(
            counts = simulation$counts[taxa, , drop = FALSE],
            metadata = simulation$metadata,
            formula = design_formula(simulation$metadata),
            group = "group",
            library_size = "reads",
            taxa_are_rows = TRUE,
            reference = "control",
            p_adjust_method = "BH",
            component = "all",
            full_output = FALSE
        )
    )
    if (is_sim_error(captured$value)) {
        status <- paste0("error: ", captured$value$message)
        return(rbind(
            method_result_template(
                taxa, "SACAT", "structural_absence", status, captured$elapsed
            ),
            method_result_template(
                taxa, "SACAT", "present_conditional_abundance", status,
                captured$elapsed
            )
        ))
    }

    fit <- captured$value
    result <- fit$results
    diagnostics <- fit$diagnostics
    result <- result[match(taxa, result$taxon), , drop = FALSE]
    diagnostics <- diagnostics[match(taxa, diagnostics$taxon), , drop = FALSE]

    structural <- method_result_template(
        taxa, "SACAT", "structural_absence", "unavailable", captured$elapsed
    )
    structural$p_value <- suppressWarnings(as.numeric(
        result$p_structural_absence
    ))
    structural$package_q_value <- suppressWarnings(as.numeric(
        result$p_adj_structural_absence
    ))
    structural$statistic <- suppressWarnings(as.numeric(
        result$z_structural_absence
    ))
    structural$available <- valid_p(structural$p_value) &
        as.logical(diagnostics$formed_structural_absence)
    structural$status <- ifelse(
        structural$available &
            as.logical(diagnostics$regular_structural_absence),
        "ok",
        as.character(diagnostics$reason_structural_absence)
    )

    abundance <- method_result_template(
        taxa, "SACAT", "present_conditional_abundance", "unavailable",
        captured$elapsed
    )
    abundance$p_value <- suppressWarnings(as.numeric(
        result$p_relative_abundance
    ))
    abundance$package_q_value <- suppressWarnings(as.numeric(
        result$p_adj_relative_abundance
    ))
    abundance$estimate <- suppressWarnings(as.numeric(
        result$estimate_relative_abundance
    ))
    abundance$standard_error <- suppressWarnings(as.numeric(
        result$se_relative_abundance
    ))
    abundance$statistic <- suppressWarnings(as.numeric(
        result$z_relative_abundance
    ))
    abundance$available <- valid_p(abundance$p_value) &
        as.logical(diagnostics$formed_relative_abundance)
    abundance$status <- ifelse(
        abundance$available,
        "ok",
        as.character(diagnostics$reason_relative_abundance)
    )

    structural$warning <- as.character(
        diagnostics$warning_structural_absence
    )
    abundance$warning <- as.character(
        diagnostics$warning_relative_abundance
    )
    rbind(structural, abundance)
}

run_zinq <- function(simulation, replication_id) {
    taxa <- simulation$evaluation_taxa
    count <- simulation$counts[taxa, , drop = FALSE]
    relative_abundance <- sweep(count, 2L, simulation$depth, "/")
    prevalence_rows <- vector("list", length(taxa))
    abundance_rows <- vector("list", length(taxa))
    tau_weights <- pmin(CONFIG$zinq_taus, 1 - CONFIG$zinq_taus)
    tau_weights <- tau_weights / sum(tau_weights)

    started_all <- proc.time()[["elapsed"]]
    for (j in seq_along(taxa)) {
        data <- data.frame(
            y = as.numeric(relative_abundance[j, ]),
            group_num = simulation$metadata$group_num
        )
        if ("z" %in% names(simulation$metadata)) data$z <- simulation$metadata$z
        data$log_depth <- simulation$metadata$log_depth
        set.seed(
            CONFIG$base_seed + replication_id * 100003L +
                simulation$setting$setting_index * 1009L + j
        )
        captured <- safe_capture(
            ZINQ::ZINQ_tests(
                formula.logistic = design_formula(
                    data, include_depth = TRUE, response = "y", group = "group_num"
                ),
                formula.quantile = design_formula(
                    data, include_depth = TRUE, response = "y", group = "group_num"
                ),
                C = "group_num",
                y_CorD = "C",
                data = data,
                taus = CONFIG$zinq_taus,
                seed = CONFIG$base_seed + replication_id * 1009L + j
            )
        )

        prevalence <- method_result_template(
            taxa[j], "ZINQ", "observed_prevalence", "unavailable",
            captured$elapsed
        )
        abundance <- method_result_template(
            taxa[j], "ZINQ", "detected_quantile_abundance", "unavailable",
            captured$elapsed
        )
        warning_text <- collapse_text(captured$warnings)
        prevalence$warning <- warning_text
        abundance$warning <- warning_text

        if (is_sim_error(captured$value)) {
            status <- paste0("error: ", captured$value$message)
            prevalence$status <- status
            abundance$status <- status
        } else {
            fit <- captured$value
            prevalence$p_value <- suppressWarnings(as.numeric(
                fit$pvalue.logistic
            )[1L])
            prevalence$available <- valid_p(prevalence$p_value)
            prevalence$status <- ifelse(
                prevalence$available, "ok", "nonfinite_prevalence_p_value"
            )

            quantile_p <- suppressWarnings(as.numeric(fit$pvalue.quantile))
            if (length(quantile_p) == length(CONFIG$zinq_taus)) {
                abundance$p_tau_025 <- quantile_p[1L]
                abundance$p_tau_050 <- quantile_p[2L]
                abundance$p_tau_075 <- quantile_p[3L]
                if (all(valid_p(quantile_p))) {
                    abundance$p_value <- cauchy_combine(
                        quantile_p, tau_weights
                    )
                }
            }
            abundance$available <- valid_p(abundance$p_value)
            abundance$status <- ifelse(
                abundance$available, "ok", "quantile_component_unavailable"
            )
        }
        prevalence_rows[[j]] <- prevalence
        abundance_rows[[j]] <- abundance
    }
    elapsed_all <- proc.time()[["elapsed"]] - started_all
    prevalence <- do.call(rbind, prevalence_rows)
    abundance <- do.call(rbind, abundance_rows)
    prevalence$runtime_seconds_total_setting <- elapsed_all
    abundance$runtime_seconds_total_setting <- elapsed_all
    rbind(prevalence, abundance)
}

maaslin_result_frame <- function(fit, component, output_dir) {
    candidate <- if (component == "abundance") {
        fit$fit_data_abundance %||% fit$abundance
    } else {
        fit$fit_data_prevalence %||% fit$prevalence
    }
    if (is.data.frame(candidate)) return(candidate)
    if (is.list(candidate) && is.data.frame(candidate$results)) {
        return(candidate$results)
    }
    if (is.list(fit) && is.data.frame(fit$results)) return(fit$results)

    files <- list.files(
        output_dir,
        pattern = "(all_results|results).*\\.(tsv|txt|csv)$",
        recursive = TRUE, full.names = TRUE, ignore.case = TRUE
    )
    for (file in files) {
        data <- tryCatch(
            if (grepl("\\.csv$", file, ignore.case = TRUE)) {
                utils::read.csv(file, check.names = FALSE)
            } else {
                utils::read.delim(file, check.names = FALSE)
            },
            error = function(e) NULL
        )
        if (is.data.frame(data) && nrow(data)) return(data)
    }
    NULL
}

select_maaslin_row <- function(data, taxon, component) {
    if (is.null(data) || !is.data.frame(data) || !nrow(data)) return(NULL)
    feature_column <- first_existing_column(
        data, c("feature", "taxon", "name_feature")
    )
    candidates <- data
    if (!is.null(feature_column)) {
        candidates <- candidates[
            as.character(candidates[[feature_column]]) == taxon,
            , drop = FALSE
        ]
    }
    if (!nrow(candidates)) return(NULL)

    if ("model" %in% names(candidates)) {
        model <- tolower(trimws(as.character(candidates$model)))
        accepted <- if (component == "abundance") {
            c("linear", "abundance")
        } else {
            c("logistic", "prevalence")
        }
        candidates <- candidates[model %in% accepted, , drop = FALSE]
    }
    if (!nrow(candidates)) return(NULL)

    if ("metadata" %in% names(candidates)) {
        keep <- as.character(candidates$metadata) == "group"
        if (any(keep, na.rm = TRUE)) candidates <- candidates[keep, , drop = FALSE]
    }
    if (nrow(candidates) > 1L && "value" %in% names(candidates)) {
        keep <- tolower(as.character(candidates$value)) == "case"
        if (any(keep, na.rm = TRUE)) candidates <- candidates[keep, , drop = FALSE]
    }
    if (nrow(candidates) > 1L && "name" %in% names(candidates)) {
        normalized <- gsub("[^a-z0-9]", "", tolower(as.character(candidates$name)))
        keep <- grepl("group.*case|case.*group", normalized)
        if (any(keep, na.rm = TRUE)) candidates <- candidates[keep, , drop = FALSE]
    }
    if (!nrow(candidates)) NULL else candidates[1L, , drop = FALSE]
}

extract_maaslin_component <- function(data, taxa, component, elapsed, warnings) {
    component_label <- if (component == "abundance") {
        "detected_log_abundance"
    } else {
        "observed_prevalence"
    }
    output <- method_result_template(
        taxa, "MaAsLin 3", component_label, "row_unavailable", elapsed
    )
    output$warning <- collapse_text(warnings)
    for (j in seq_along(taxa)) {
        row <- select_maaslin_row(data, taxa[j], component)
        if (is.null(row)) next
        output$p_value[j] <- as_numeric_column(
            row, c("pval_individual", "pval", "p_value", "pvalue")
        )[1L]
        output$package_q_value[j] <- as_numeric_column(
            row, c("qval_individual", "qval", "q_value", "qvalue")
        )[1L]
        output$estimate[j] <- as_numeric_column(
            row, c("coef", "coefficient", "estimate")
        )[1L]
        output$standard_error[j] <- as_numeric_column(
            row, c("stderr", "std_error", "se")
        )[1L]
        output$statistic[j] <- as_numeric_column(
            row, c("stat", "test_statistic", "z")
        )[1L]
        output$available[j] <- valid_p(output$p_value[j])
        error_column <- first_existing_column(row, c("error", "warning"))
        error_text <- if (is.null(error_column)) NA_character_ else {
            collapse_text(row[[error_column]])
        }
        output$method_detail[j] <- error_text
        if (component == "prevalence") {
            output$abundance_induced_flag[j] <- !is.na(error_text) &&
                grepl("abundance|induced", error_text, ignore.case = TRUE)
        }
        output$status[j] <- if (output$available[j]) {
            "ok"
        } else if (!is.na(error_text)) {
            error_text
        } else {
            "nonfinite_p_value"
        }
    }
    output
}

run_maaslin3 <- function(simulation, replication_id, temporary_root) {
    taxa <- simulation$evaluation_taxa
    output_dir <- file.path(
        temporary_root,
        paste0("maaslin3_", simulation$setting$setting_id)
    )
    unlink(output_dir, recursive = TRUE, force = TRUE)
    safe_dir_create(output_dir)
    on.exit({
        try(maaslin3::maaslin_log_reset(), silent = TRUE)
        unlink(output_dir, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    input_data <- as.data.frame(t(simulation$counts))
    metadata <- simulation$metadata[, intersect(
        c("group", "z", "log_depth"), names(simulation$metadata)
    ), drop = FALSE]
    arguments <- list(
        input_data = input_data,
        input_metadata = metadata,
        output = output_dir,
        formula = paste(deparse(design_formula(
            metadata, include_depth = TRUE
        )), collapse = " "),
        min_abundance = 0,
        min_prevalence = 0,
        max_prevalence = 1.01,
        zero_threshold = 0,
        min_variance = 0,
        normalization = "TSS",
        transform = "LOG",
        correction = "BH",
        standardize = TRUE,
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        warn_prevalence = isTRUE(CONFIG$maaslin3_warn_prevalence),
        augment = TRUE,
        evaluate_only = NULL,
        plot_summary_plot = FALSE,
        plot_associations = FALSE,
        max_pngs = 1,
        cores = 1,
        save_models = FALSE,
        save_plots_rds = FALSE,
        verbosity = "ERROR",
        reference = "group,control"
    )
    available <- names(formals(maaslin3::maaslin3))
    arguments <- arguments[names(arguments) %in% available]

    set.seed(
        CONFIG$base_seed + replication_id * 100003L +
            simulation$setting$setting_index * 4001L
    )
    try(maaslin3::maaslin_log_reset(), silent = TRUE)
    captured <- safe_capture({
        invisible(capture.output(
            fit_value <- suppressMessages(do.call(maaslin3::maaslin3, arguments))
        ))
        fit_value
    })
    if (is_sim_error(captured$value)) {
        status <- paste0("error: ", captured$value$message)
        return(rbind(
            method_result_template(
                taxa, "MaAsLin 3", "observed_prevalence", status,
                captured$elapsed
            ),
            method_result_template(
                taxa, "MaAsLin 3", "detected_log_abundance", status,
                captured$elapsed
            )
        ))
    }

    abundance_data <- maaslin_result_frame(
        captured$value, "abundance", output_dir
    )
    prevalence_data <- maaslin_result_frame(
        captured$value, "prevalence", output_dir
    )
    rbind(
        extract_maaslin_component(
            prevalence_data, taxa, "prevalence", captured$elapsed,
            captured$warnings
        ),
        extract_maaslin_component(
            abundance_data, taxa, "abundance", captured$elapsed,
            captured$warnings
        )
    )
}

run_edger <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        design <- stats::model.matrix(
            design_formula(simulation$metadata), data = simulation$metadata
        )
        coefficient <- grep("^group", colnames(design))
        if (length(coefficient) != 1L) {
            stop("Could not identify the edgeR group coefficient.")
        }
        dge <- edgeR::DGEList(counts = simulation$counts)
        keep <- edgeR::filterByExpr(dge, design = design)
        if (!any(keep)) stop("edgeR filterByExpr retained no taxa.")
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        dge <- edgeR::estimateDisp(dge, design = design, robust = TRUE)
        fit <- edgeR::glmQLFit(dge, design = design, robust = TRUE)
        test <- edgeR::glmQLFTest(fit, coef = coefficient)
        table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
        list(table = table, keep = keep)
    })
    output <- method_result_template(
        taxa, "edgeR", "general_abundance", "filtered_or_unavailable",
        captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }
    table <- captured$value$table
    matched <- match(taxa, rownames(table))
    present <- !is.na(matched)
    output$p_value[present] <- suppressWarnings(as.numeric(
        table$PValue[matched[present]]
    ))
    output$package_q_value[present] <- suppressWarnings(as.numeric(
        table$FDR[matched[present]]
    ))
    output$estimate[present] <- suppressWarnings(as.numeric(
        table$logFC[matched[present]]
    ))
    output$statistic[present] <- suppressWarnings(as.numeric(
        table$F[matched[present]]
    ))
    output$available <- valid_p(output$p_value)
    output$status <- ifelse(output$available, "ok", "filterByExpr_excluded")
    output
}

run_deseq2 <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        col_data <- simulation$metadata[, intersect(
            c("group", "z"), names(simulation$metadata)
        ), drop = FALSE]
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = round(simulation$counts),
            colData = col_data,
            design = design_formula(col_data, group_last = TRUE)
        )
        dds <- DESeq2::DESeq(
            dds,
            test = "Wald",
            fitType = "parametric",
            sfType = "poscounts",
            betaPrior = FALSE,
            quiet = TRUE
        )
        coefficient <- grep(
            "group.*case.*control|group_case_vs_control",
            DESeq2::resultsNames(dds), value = TRUE, ignore.case = TRUE
        )
        if (length(coefficient) != 1L) {
            stop("Could not identify the DESeq2 group coefficient.")
        }
        result <- DESeq2::results(
            dds,
            name = coefficient,
            alpha = CONFIG$alpha,
            independentFiltering = TRUE,
            cooksCutoff = TRUE
        )
        as.data.frame(result)
    })
    output <- method_result_template(
        taxa, "DESeq2", "general_abundance", "unavailable", captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }
    table <- captured$value
    matched <- match(taxa, rownames(table))
    present <- !is.na(matched)
    output$p_value[present] <- suppressWarnings(as.numeric(
        table$pvalue[matched[present]]
    ))
    output$package_q_value[present] <- suppressWarnings(as.numeric(
        table$padj[matched[present]]
    ))
    output$estimate[present] <- suppressWarnings(as.numeric(
        table$log2FoldChange[matched[present]]
    ))
    output$standard_error[present] <- suppressWarnings(as.numeric(
        table$lfcSE[matched[present]]
    ))
    output$statistic[present] <- suppressWarnings(as.numeric(
        table$stat[matched[present]]
    ))
    output$available <- valid_p(output$p_value)
    output$status <- ifelse(
        output$available, "ok", "independent_filter_or_nonfinite_p"
    )
    output
}

find_ancom_column <- function(names_vector, prefix, term_pattern) {
    normalized <- tolower(names_vector)
    prefix_pattern <- paste0("^", tolower(prefix), "[_\\.]")
    candidates <- which(
        grepl(prefix_pattern, normalized) & grepl(term_pattern, normalized)
    )
    if (!length(candidates)) return(NA_integer_)
    candidates[[1L]]
}

run_ancombc2 <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        arguments <- list(
            data = simulation$counts,
            taxa_are_rows = TRUE,
            assay_name = "counts",
            tax_level = NULL,
            aggregate_data = NULL,
            meta_data = simulation$metadata,
            fix_formula = paste(attr(stats::terms(design_formula(
                simulation$metadata
            )), "term.labels"), collapse = " + "),
            rand_formula = NULL,
            p_adj_method = "BH",
            pseudo = 0,
            pseudo_sens = TRUE,
            prv_cut = 0,
            lib_cut = 0,
            s0_perc = 0.05,
            group = "group",
            struc_zero = TRUE,
            neg_lb = TRUE,
            alpha = CONFIG$alpha,
            n_cl = 1,
            verbose = FALSE,
            global = FALSE,
            pairwise = FALSE,
            dunnet = FALSE,
            trend = FALSE
        )
        available_arguments <- names(formals(ANCOMBC::ancombc2))
        arguments <- arguments[names(arguments) %in% available_arguments]
        fit <- suppressMessages(do.call(ANCOMBC::ancombc2, arguments))
        result <- fit$res
        if (!is.data.frame(result)) result <- as.data.frame(result)
        result
    })
    output <- method_result_template(
        taxa, "ANCOM-BC2", "general_abundance", "unavailable",
        captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }

    table <- captured$value
    taxon_column <- first_existing_column(
        table, c("taxon", "feature", "otu")
    )
    table_taxa <- if (is.null(taxon_column)) {
        rownames(table)
    } else {
        as.character(table[[taxon_column]])
    }
    matched <- match(taxa, table_taxa)
    present <- !is.na(matched)
    term_pattern <- "group.*case|case.*group"
    p_col <- find_ancom_column(names(table), "p", term_pattern)
    q_col <- find_ancom_column(names(table), "q", term_pattern)
    lfc_col <- find_ancom_column(names(table), "lfc", term_pattern)
    se_col <- find_ancom_column(names(table), "se", term_pattern)
    w_col <- find_ancom_column(names(table), "w", term_pattern)
    sensitivity_col <- find_ancom_column(
        names(table), "passed_ss", term_pattern
    )

    if (is.na(p_col)) {
        output$status <- "group_p_value_column_unavailable"
        return(output)
    }

    raw_p <- rep(NA_real_, length(taxa))
    raw_p[present] <- suppressWarnings(as.numeric(
        table[matched[present], p_col]
    ))
    output$p_value_unfiltered <- raw_p
    output$p_value <- raw_p

    if (!is.na(sensitivity_col)) {
        passed <- rep(NA, length(taxa))
        passed[present] <- as.logical(
            table[matched[present], sensitivity_col]
        )
        output$passed_sensitivity <- passed
        failed <- !is.na(passed) & !passed & valid_p(raw_p)
        output$p_value[failed] <- 1
        output$method_detail[failed] <-
            "pseudocount_sensitivity_failed; p-value set to one"
    } else {
        output$method_detail[present] <-
            "pseudocount_sensitivity_column_not_returned"
    }

    if (!is.na(q_col)) {
        output$package_q_value[present] <- suppressWarnings(as.numeric(
            table[matched[present], q_col]
        ))
    }
    if (!is.na(lfc_col)) {
        output$estimate[present] <- suppressWarnings(as.numeric(
            table[matched[present], lfc_col]
        ))
    }
    if (!is.na(se_col)) {
        output$standard_error[present] <- suppressWarnings(as.numeric(
            table[matched[present], se_col]
        ))
    }
    if (!is.na(w_col)) {
        output$statistic[present] <- suppressWarnings(as.numeric(
            table[matched[present], w_col]
        ))
    }

    output$available <- valid_p(output$p_value)
    output$status <- ifelse(
        output$available,
        ifelse(
            !is.na(output$passed_sensitivity) &
                !output$passed_sensitivity,
            "ok_after_sensitivity_exclusion",
            "ok"
        ),
        "nonfinite_p_value"
    )
    output
}
run_linda <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    linda_formals <- names(formals(MicrobiomeStat::linda))

    arguments <- if ("feature.dat" %in% linda_formals) {
        list(
            feature.dat = simulation$counts,
            meta.dat = simulation$metadata[, intersect(
                c("group", "z"), names(simulation$metadata)
            ), drop = FALSE],
            formula = paste(deparse(design_formula(simulation$metadata)),
                            collapse = " "),
            feature.dat.type = "count",
            prev.filter = 0,
            mean.abund.filter = 0,
            max.abund.filter = 0,
            adaptive = TRUE,
            p.adj.method = "BH",
            alpha = CONFIG$alpha,
            n.cores = 1,
            verbose = FALSE
        )
    } else {
        list(
            otu.tab = simulation$counts,
            meta = simulation$metadata[, intersect(
                c("group", "z"), names(simulation$metadata)
            ), drop = FALSE],
            formula = paste(deparse(design_formula(simulation$metadata)),
                            collapse = " "),
            type = "count",
            adaptive = TRUE,
            p.adj.method = "BH",
            alpha = CONFIG$alpha,
            prev.cut = 0,
            lib.cut = 0,
            n.cores = 1
        )
    }
    arguments <- arguments[names(arguments) %in% linda_formals]

    captured <- safe_capture(
        suppressMessages(do.call(MicrobiomeStat::linda, arguments))
    )
    output <- method_result_template(
        taxa, "LinDA", "general_abundance", "unavailable",
        captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }

    fit <- captured$value
    if (is.null(fit$output) || !is.list(fit$output)) {
        output$status <- "output_list_unavailable"
        return(output)
    }
    output_names <- names(fit$output)
    coefficient_name <- output_names[
        grepl("group.*case|case.*group", output_names, ignore.case = TRUE)
    ]
    if (length(coefficient_name) > 1L) {
        normalized <- gsub(
            "[^a-z0-9]", "", tolower(coefficient_name)
        )
        exact <- which(normalized %in% c("groupcase", "groupcasevscontrol"))
        if (length(exact)) coefficient_name <- coefficient_name[exact[1L]]
    }
    if (length(coefficient_name) != 1L) {
        output$status <- "group_output_unavailable"
        return(output)
    }

    table <- as.data.frame(fit$output[[coefficient_name]])
    matched <- match(taxa, rownames(table))
    present <- !is.na(matched)
    selected <- table[matched[present], , drop = FALSE]
    output$p_value[present] <- as_numeric_column(
        selected, c("pvalue", "p_value", "p.value")
    )
    output$package_q_value[present] <- as_numeric_column(
        selected, c("padj", "qvalue", "q_value", "adj.p")
    )
    output$estimate[present] <- as_numeric_column(
        selected, c("log2FoldChange", "coef", "estimate")
    )
    output$standard_error[present] <- as_numeric_column(
        selected, c("lfcSE", "se", "stderr", "standard_error")
    )
    output$statistic[present] <- as_numeric_column(
        selected, c("stat", "t", "z")
    )
    output$available <- valid_p(output$p_value)
    output$status <- ifelse(
        output$available, "ok", "nonfinite_p_value"
    )
    output
}
run_corncob <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    rows <- vector("list", length(taxa))
    started <- proc.time()[["elapsed"]]
    for (j in seq_along(taxa)) {
        data <- data.frame(
            W = as.numeric(simulation$counts[taxa[j], ]),
            M = as.numeric(simulation$depth),
            group = simulation$metadata$group
        )
        if ("z" %in% names(simulation$metadata)) data$z <- simulation$metadata$z
        captured <- safe_capture(
            corncob::bbdml(
                formula = design_formula(data, response = "cbind(W, M - W)"),
                phi.formula = design_formula(data),
                data = data,
                method = "trust",
                robust = isTRUE(CONFIG$corncob_robust)
            )
        )
        row <- method_result_template(
            taxa[j], "corncob", "general_abundance", "unavailable",
            captured$elapsed
        )
        row$warning <- collapse_text(captured$warnings)
        if (is_sim_error(captured$value)) {
            row$status <- paste0("error: ", captured$value$message)
            rows[[j]] <- row
            next
        }
        fit <- captured$value
        summary_fit <- tryCatch(summary(fit), error = function(e) NULL)
        coefficient_table <- if (is.null(summary_fit)) NULL else summary_fit$coefficients
        if (is.null(coefficient_table) || !is.matrix(coefficient_table)) {
            row$status <- "coefficient_table_unavailable"
            rows[[j]] <- row
            next
        }
        abundance_table <- coefficient_table[
            seq_len(min(fit$np.mu, nrow(coefficient_table))), , drop = FALSE
        ]
        coefficient_row <- grep(
            "group.*case|case.*group", rownames(abundance_table),
            ignore.case = TRUE
        )
        if (length(coefficient_row) != 1L) {
            row$status <- "group_coefficient_unavailable"
            rows[[j]] <- row
            next
        }
        columns <- colnames(abundance_table)
        estimate_col <- grep("estimate", columns, ignore.case = TRUE)[1L]
        se_col <- grep("std.*error|standard.*error", columns, ignore.case = TRUE)[1L]
        stat_col <- grep("t value|z value|stat", columns, ignore.case = TRUE)[1L]
        p_col <- grep("pr\\(|p.value|p value", columns, ignore.case = TRUE)[1L]
        if (is.na(p_col)) {
            row$status <- "group_p_value_unavailable"
            rows[[j]] <- row
            next
        }
        row$p_value <- suppressWarnings(as.numeric(
            abundance_table[coefficient_row, p_col]
        ))
        if (!is.na(estimate_col)) {
            row$estimate <- suppressWarnings(as.numeric(
                abundance_table[coefficient_row, estimate_col]
            ))
        }
        if (!is.na(se_col)) {
            row$standard_error <- suppressWarnings(as.numeric(
                abundance_table[coefficient_row, se_col]
            ))
        }
        if (!is.na(stat_col)) {
            row$statistic <- suppressWarnings(as.numeric(
                abundance_table[coefficient_row, stat_col]
            ))
        }
        row$available <- valid_p(row$p_value)
        row$status <- ifelse(row$available, "ok", "nonfinite_p_value")
        rows[[j]] <- row
    }
    output <- do.call(rbind, rows)
    output$runtime_seconds_total_setting <- proc.time()[["elapsed"]] - started
    output
}

run_metagenomeseq <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        pheno <- Biobase::AnnotatedDataFrame(
            simulation$metadata[, intersect(
                c("group", "z"), names(simulation$metadata)
            ), drop = FALSE]
        )
        object <- metagenomeSeq::newMRexperiment(
            counts = simulation$counts,
            phenoData = pheno
        )
        percentile <- metagenomeSeq::cumNormStatFast(
            object, pFlag = FALSE
        )
        object <- metagenomeSeq::cumNorm(object, p = percentile)
        design <- stats::model.matrix(
            design_formula(simulation$metadata), data = simulation$metadata
        )
        coefficient <- grep("^group", colnames(design))
        if (length(coefficient) != 1L) {
            stop("Could not identify the metagenomeSeq group coefficient.")
        }
        coefficient_name <- colnames(design)[coefficient]
        control <- metagenomeSeq::zigControl(
            maxit = CONFIG$metagenomeseq_maxit,
            verbose = FALSE
        )
        fit <- metagenomeSeq::fitZig(
            obj = object, mod = design, control = control
        )

        table <- tryCatch(
            metagenomeSeq::MRcoefs(
                fit,
                by = coefficient_name,
                coef = coefficient_name,
                number = Inf,
                group = 4,
                adjustMethod = "BH",
                eff = 0,
                counts = 0
            ),
            error = function(e) NULL
        )
        if (is.null(table)) {
            fit_slot <- methods::slot(fit, "fit")
            table <- limma::topTable(
                fit_slot,
                coef = coefficient,
                number = Inf,
                sort.by = "none",
                adjust.method = "BH"
            )
        }
        list(
            table = as.data.frame(table),
            coefficient_name = coefficient_name
        )
    })
    output <- method_result_template(
        taxa, "metagenomeSeq", "general_abundance", "unavailable",
        captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }

    table <- captured$value$table
    coefficient_name <- captured$value$coefficient_name
    taxon_column <- first_existing_column(
        table, c("taxon", "feature", "otu")
    )
    table_taxa <- if (is.null(taxon_column)) {
        rownames(table)
    } else {
        as.character(table[[taxon_column]])
    }
    matched <- match(taxa, table_taxa)
    present <- !is.na(matched)
    selected <- table[matched[present], , drop = FALSE]

    output$p_value[present] <- as_numeric_column(
        selected, c("pvalues", "P.Value", "pvalue", "p_value")
    )
    output$package_q_value[present] <- as_numeric_column(
        selected, c("adjPvalues", "adj.P.Val", "padj", "q_value")
    )
    estimate_candidates <- c(
        "logFC", "coef", "coefficient", "estimate"
    )
    if (!is.null(first_existing_column(selected, coefficient_name))) {
        estimate_candidates <- c(coefficient_name, estimate_candidates)
    }
    output$estimate[present] <- as_numeric_column(
        selected, estimate_candidates
    )
    output$standard_error[present] <- as_numeric_column(
        selected, c("se", "SE", "stderr", "standard_error")
    )
    output$statistic[present] <- as_numeric_column(
        selected, c("t", "stat", "z")
    )
    output$available <- valid_p(output$p_value)
    output$status <- ifelse(
        output$available, "ok", "nonfinite_p_value"
    )
    output
}
run_generic_abundance_methods <- function(simulation) {
    rbind(
        run_edger(simulation),
        run_deseq2(simulation),
        run_ancombc2(simulation),
        run_linda(simulation),
        run_corncob(simulation),
        run_metagenomeseq(simulation)
    )
}

apply_common_multiplicity <- function(results) {
    data <- data.table::as.data.table(results)
    data[, p_value_for_multiplicity := ifelse(valid_p(p_value), p_value, 1)]
    data[, q_value_bh := stats::p.adjust(
        p_value_for_multiplicity, method = "BH"
    ), by = .(method, component)]
    data[, reject_bh := q_value_bh < CONFIG$alpha]
    as.data.frame(data)
}

set_evaluation_truth <- function(result) {
    structural <- result$component == "structural_absence"
    abundance <- result$component %in% c(
        "present_conditional_abundance", "detected_quantile_abundance",
        "detected_log_abundance"
    )
    result$evaluation_target <- ifelse(
        structural, "structural_absence",
        ifelse(abundance, "conditional_abundance", "perturbed_taxa")
    )
    result$truth_for_component <- ifelse(
        structural, result$truth_structural,
        ifelse(abundance, result$truth_abundance, result$truth_any)
    )
    result
}

analyze_setting <- function(simulation, replication_id, temporary_root) {
    setting <- simulation$setting
    method_seed <- CONFIG$base_seed + replication_id * 100003L +
        setting$setting_index * 10007L
    set.seed(method_seed + 1L)
    sacat_result <- run_sacat(simulation)
    set.seed(method_seed + 2L)
    zinq_result <- run_zinq(simulation, replication_id)
    set.seed(method_seed + 3L)
    maaslin_result <- run_maaslin3(simulation, replication_id, temporary_root)
    set.seed(method_seed + 101L)
    general_result <- run_generic_abundance_methods(simulation)
    result <- apply_common_multiplicity(rbind(
        sacat_result, zinq_result, maaslin_result, general_result
    ))
    result$setting_id <- setting$setting_id
    result$replication <- as.integer(replication_id)
    truth <- simulation$truth[match(result$taxon, simulation$truth$taxon), ]
    for (column in setdiff(names(truth), "taxon")) {
        result[[column]] <- truth[[column]]
    }
    rownames(result) <- NULL
    set_evaluation_truth(result)
}

make_metasparsim_grid <- function() {
    sa_effect <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.40)
    ra_effect <- c(0.20, 0.40, 0.50, 0.60, 0.80, 1.20)

    null <- data.frame(
        include_covariate = c(FALSE, TRUE),
        scenario = "global_null", effect_index = 0L,
        structural_difference = 0, abundance_log_effect = 0
    )

    alternative <- expand.grid(
        include_covariate = c(FALSE, TRUE),
        scenario = c("structural_only", "abundance_only", "joint"),
        effect_index = seq_along(sa_effect),
        stringsAsFactors = FALSE
    )
    alternative$structural_difference <- ifelse(
        alternative$scenario == "abundance_only", 0,
        sa_effect[alternative$effect_index]
    )
    alternative$abundance_log_effect <- ifelse(
        alternative$scenario == "structural_only", 0,
        ra_effect[alternative$effect_index]
    )
    grid <- rbind(null, alternative)
    rownames(grid) <- NULL
    grid$setting_index <- seq_len(nrow(grid))
    grid$setting_id <- sprintf("setting_%04d", grid$setting_index)
    grid$n_control <- 152L
    grid$n_case <- 91L
    grid$total_sample_size <- grid$n_control + grid$n_case
    grid$n_taxa <- CONFIG$n_taxa
    grid$n_signal_target <- ifelse(grid$scenario == "global_null", 0L, 20L)
    grid$signal_fraction <- grid$n_signal_target / grid$n_taxa
    grid$depth_control <- 5855.5
    grid$depth_case <- 3678L
    grid$gamma_dispersion <- 0.75
    grid$focal_mass <- 0.10
    grid$n_replications <- ifelse(
        grid$scenario == "global_null",
        CONFIG$null_replications, CONFIG$replications
    )
    stopifnot(
        nrow(grid) == 38L, sum(grid$scenario == "global_null") == 2L,
        !anyDuplicated(grid[c("include_covariate", "scenario", "effect_index")])
    )
    grid
}

make_metasparsim_template <- function() {
    probability <- exp(seq(log(1e-4), log(2.8e-3), length.out = 50L))
    candidates <- round(seq(1L, 50L, length.out = 20L))
    candidates <- candidates[c(seq(1L, 19L, 2L), seq(2L, 20L, 2L))]
    set.seed(CONFIG$base_seed + 1701L)
    order <- c(candidates, sample(setdiff(seq_len(50L), candidates)))
    gamma_mean <- probability[order]
    gamma_mean[seq_len(20L)] <- exp(seq(log(5e-4), log(2e-3), length.out = 20L))
    profiles <- data.frame(
        gamma_mean = gamma_mean,
        baseline_rho = 0.35,
        abundance_covariate = stats::runif(50L, -0.16, 0.16),
        structural_covariate = stats::runif(50L, -0.22, 0.22)
    )
    template <- profiles[rep(seq_len(nrow(profiles)), 2L), ]
    template$taxon <- sprintf("Taxon_%03d", seq_len(nrow(template)))
    rownames(template) <- NULL
    stopifnot(nrow(template) == CONFIG$n_taxa)
    template
}

simulate_metasparsim <- function(replication_id, setting,
                                template = make_metasparsim_template()) {
    J <- CONFIG$n_taxa
    n_control <- as.integer(setting$n_control)
    n_case <- as.integer(setting$n_case)
    n <- n_control + n_case
    taxa <- template$taxon
    samples <- sprintf("Sample_%03d", seq_len(n))
    group_num <- rep(0:1, times = c(n_control, n_case))
    seed <- CONFIG$simulation_seed + replication_id * 1000000L

    set.seed(seed + 1L)
    z <- if (setting$include_covariate) {
        as.numeric(scale(stats::rnorm(n)))
    } else rep(0, n)

    set.seed(seed + 2L)
    control_multiplier <- exp(stats::rnorm(n_control, sd = 0.45))
    case_multiplier <- exp(stats::rnorm(n_case, sd = 0.45))
    control_depth <- as.integer(round(setting$depth_control *
        control_multiplier / median(control_multiplier)))
    case_depth <- as.integer(round(setting$depth_case *
        case_multiplier / median(case_multiplier)))
    depth <- as.integer(c(control_depth, case_depth))
    metadata <- data.frame(
        group = factor(ifelse(group_num == 0, "control", "case"),
                       levels = c("control", "case")),
        group_num = group_num, reads = depth,
        log_depth = as.numeric(scale(log(depth))), row.names = samples
    )
    if (setting$include_covariate) metadata$z <- z

    signal <- seq_len(J) <= setting$n_signal_target
    sa_sign <- rep(c(-1, 1), length.out = J)
    ra_sign <- if (setting$scenario == "joint") {
        rep(c(-1, 1, 1, -1), length.out = J)
    } else sa_sign
    half_shift <- sa_sign * signal * setting$structural_difference / 2
    absence_control <- template$baseline_rho - half_shift
    absence_case <- template$baseline_rho + half_shift
    delta <- qlogis(absence_case) - qlogis(absence_control)
    zeta <- ra_sign * signal * setting$abundance_log_effect

    set.seed(seed + 3L)
    gamma <- getFromNamespace("simulate_biovar", "metaSPARSim")(
        avgAbund = stats::setNames(template$gamma_mean, taxa),
        phi = rep(setting$gamma_dispersion, J), nrep = n
    )
    positive_baseline <- sweep(gamma, 2L, colSums(gamma), "/") * setting$focal_mass
    positive_probability <- positive_baseline * exp(
        outer(template$abundance_covariate, z) + outer(zeta, group_num - 0.5)
    )
    rho_control <- plogis(
        outer(template$structural_covariate, z) + qlogis(absence_control)
    )
    rho_case <- plogis(
        outer(template$structural_covariate, z) + qlogis(absence_case)
    )
    rho <- rho_control
    rho[, group_num == 1L] <- rho_case[, group_num == 1L]

    set.seed(seed + 4L)
    profile_size <- J %/% 2L
    uniform <- rbind(
        matrix(runif(profile_size * n), profile_size, n),
        matrix(runif(profile_size * n), profile_size, n)
    )
    present <- uniform > rho
    # The background preserves the relative-abundance truth of other focal taxa.
    probability <- rbind(
        positive_probability * present,
        Other_unmodeled = 1 - colSums(positive_probability * present)
    )
    dimnames(probability) <- list(c(taxa, "Other_unmodeled"), samples)
    dimnames(positive_baseline) <- dimnames(positive_probability) <-
        dimnames(present) <- dimnames(rho) <- list(taxa, samples)
    stopifnot(
        all(probability >= 0), all(probability[J + 1L, ] > 0),
        max(abs(colSums(probability) - 1)) < 1e-12,
        abs(median(depth[group_num == 0L]) - setting$depth_control) <= 0.5,
        abs(median(depth[group_num == 1L]) - setting$depth_case) <= 0.5
    )

    # Gamma variation has already been generated by the package.
    parameters <- lapply(seq_len(n), function(i) {
        list(intensity = probability[, i], variability = rep(0, J + 1L),
             lib_size = depth[i])
    })
    names(parameters) <- samples
    set.seed(seed + 5L)
    invisible(utils::capture.output(native <- withCallingHandlers(
        metaSPARSim::metaSPARSim(parameters),
        warning = function(w) {
            if (identical(conditionMessage(w),
                          "no non-missing arguments to min; returning Inf")) {
                invokeRestart("muffleWarning")
            }
        }
    )))
    counts <- native$counts
    fragments <- native$gamma
    dimnames(counts) <- dimnames(fragments) <- dimnames(probability)
    storage.mode(counts) <- "integer"
    stopifnot(
        all(counts >= 0L), identical(as.integer(colSums(counts)), depth),
        all(counts[seq_len(J), , drop = FALSE][!present] == 0L)
    )
    fragment_probability <- sweep(fragments, 2L, colSums(fragments), "/")

    truth <- data.frame(
        taxon = taxa, truth_structural = delta != 0,
        truth_abundance = zeta != 0, truth_any = delta != 0 | zeta != 0,
        absence_control_at_z_zero = absence_control,
        absence_case_at_z_zero = absence_case,
        structural_effect = absence_case - absence_control,
        structural_effect_standardized = rowMeans(rho_case - rho_control),
        structural_log_odds_effect = delta,
        abundance_effect = zeta,
        joint_relationship = ifelse(
            !signal | setting$scenario != "joint", "none",
            ifelse(sa_sign == ra_sign, "opposing", "reinforcing")
        )
    )
    diagnostics <- data.frame(
        depth_median_control = median(depth[group_num == 0L]),
        depth_median_case = median(depth[group_num == 1L]),
        covariate_mean_difference = mean(z[group_num == 1L]) - mean(z[group_num == 0L]),
        zero_fraction = mean(counts[seq_len(J), , drop = FALSE] == 0L),
        biological_absence_fraction = mean(!present),
        sampling_zero_given_presence = mean(
            counts[seq_len(J), , drop = FALSE][present] == 0L
        ),
        rounding_zero_given_presence = mean(
            fragments[seq_len(J), , drop = FALSE][present] == 0
        ),
        max_fragment_probability_error = max(abs(fragment_probability - probability)),
        fragment_pool = unique(colSums(fragments)),
        minimum_background_probability = min(probability[J + 1L, ])
    )
    list(
        counts = counts, evaluation_taxa = taxa, depth = depth,
        metadata = metadata, setting = setting, truth = truth,
        diagnostics = diagnostics,
        latent = list(
            positive_baseline = positive_baseline,
            positive_probability = positive_probability, present = present,
            absence_probability = rho, probability = probability,
            fragments = fragments
        )
    )
}

save_result <- function(object, path) {
    safe_dir_create(dirname(path))
    temporary <- tempfile("result_", tmpdir = dirname(path))
    saveRDS(object, temporary, compress = "gzip")
    if (!file.rename(temporary, path)) stop("Could not save ", path)
}

software_versions <- function() {
    data.frame(
        package = c("R", REQUIRED_PACKAGES),
        version = c(as.character(getRversion()), vapply(
            REQUIRED_PACKAGES, function(p) as.character(utils::packageVersion(p)),
            character(1L)
        ))
    )
}

check_packages <- function() {
    missing <- REQUIRED_PACKAGES[!vapply(
        REQUIRED_PACKAGES, requireNamespace, logical(1L), quietly = TRUE
    )]
    if (length(missing)) stop(
        "Packages not available in this R environment: ",
        paste(missing, collapse = ", "), call. = FALSE
    )
    if (utils::packageVersion("SACAT") < "0.7.0") stop(
        "SACAT 0.7.0 or later is required.", call. = FALSE
    )
    random_uniform <- getFromNamespace("random_unif", "metaSPARSim")
    set.seed(901L)
    first <- random_uniform(100L, 7L)
    second <- random_uniform(100L, 7L)
    set.seed(901L)
    repeated <- random_uniform(100L, 7L)
    if (!identical(first, repeated) || identical(first, second)) stop(
        "metaSPARSim must use R-seeded random sampling.", call. = FALSE
    )
    invisible(TRUE)
}

add_omnibus <- function(result) {
    sa <- result[result$method == "SACAT" &
                     result$component == "structural_absence", , drop = FALSE]
    ra <- result[result$method == "SACAT" &
                     result$component == "present_conditional_abundance", ,
                 drop = FALSE]
    ra <- ra[match(sa$taxon, ra$taxon), , drop = FALSE]
    omnibus <- sa
    omnibus$component <- "primary_omnibus"
    p_sa <- ifelse(valid_p(sa$p_value), sa$p_value, 1)
    p_ra <- ifelse(valid_p(ra$p_value), ra$p_value, 1)
    omnibus$p_value <- pmin(1, 2 * pmin(p_sa, p_ra))
    omnibus$p_value_unfiltered <- omnibus$p_value
    omnibus$package_q_value <- NA_real_
    omnibus$estimate <- omnibus$standard_error <- omnibus$statistic <- NA_real_
    omnibus$available <- sa$available | ra$available
    omnibus$status <- ifelse(omnibus$available, "ok", "both_components_unavailable")
    omnibus$warning <- NA_character_
    omnibus$truth_for_component <- omnibus$truth_structural | omnibus$truth_abundance
    omnibus <- apply_common_multiplicity(omnibus)
    set_evaluation_truth(rbind(result, omnibus))
}

validate_analysis <- function(result) {
    stopifnot(nrow(result) == 13L * CONFIG$n_taxa,
              !anyNA(result$truth_structural),
              !anyNA(result$truth_abundance),
              all(valid_p(result$q_value_bh)))
    families <- split(result, interaction(result$method, result$component,
                                          drop = TRUE))
    stopifnot(length(families) == 13L,
              all(vapply(families, nrow, integer(1L)) == CONFIG$n_taxa))
    for (x in families) {
        if (all(grepl("^error:", x$status))) stop(
            x$method[1L], " / ", x$component[1L],
            " failed for every taxon: ", x$status[1L], call. = FALSE
        )
        expected <- p.adjust(ifelse(valid_p(x$p_value), x$p_value, 1), "BH")
        stopifnot(isTRUE(all.equal(x$q_value_bh, expected)))
    }
    invisible(TRUE)
}

mean_selected <- function(x, selected) {
    if (any(selected)) mean(x[selected]) else NA_real_
}

replication_metrics <- function(result) {
    groups <- split(result, interaction(result$method, result$component, drop = TRUE))
    do.call(rbind, lapply(groups, function(x) {
        truth <- x$truth_for_component
        rejected <- x$reject_bh
        n_discoveries <- sum(rejected)
        n_false <- sum(rejected & !truth)
        out <- data.frame(
            setting_id = x$setting_id[1L], replication = x$replication[1L],
            method = x$method[1L], component = x$component[1L],
            evaluation_target = x$evaluation_target[1L],
            n_taxa = nrow(x), n_signals = sum(truth), n_null = sum(!truth),
            discoveries = n_discoveries,
            true_discoveries = sum(rejected & truth),
            false_discoveries = n_false,
            power = mean_selected(rejected, truth),
            fdp = n_false / max(n_discoveries, 1L),
            type1_error = mean_selected(rejected, !truth),
            availability = mean(x$available),
            signal_availability = mean_selected(x$available, truth),
            null_availability = mean_selected(x$available, !truth),
            recovery_structural_targets = mean_selected(rejected, x$truth_structural),
            recovery_abundance_targets = mean_selected(rejected, x$truth_abundance),
            recovery_any_targets = mean_selected(rejected, x$truth_any),
            non_target_component_rejection = mean_selected(
                rejected, x$truth_any & !truth
            ),
            error_fraction = mean(grepl("^error:", x$status))
        )
        if (x$component[1L] == "present_conditional_abundance") {
            usable <- truth & x$available & is.finite(x$estimate) &
                is.finite(x$standard_error) & x$standard_error > 0
            difference <- x$estimate - x$abundance_effect
            out$effect_bias <- mean_selected(difference, usable)
            out$effect_rmse <- if (any(usable)) {
                sqrt(mean(difference[usable]^2))
            } else NA_real_
            out$interval_coverage <- mean_selected(
                abs(difference) <= qnorm(0.975) * x$standard_error, usable
            )
        } else {
            out$effect_bias <- out$effect_rmse <- out$interval_coverage <- NA_real_
        }
        out
    }))
}

summarize_metrics <- function(metrics) {
    fields <- c(
        "power", "fdp", "type1_error", "availability", "signal_availability",
        "null_availability", "recovery_structural_targets",
        "recovery_abundance_targets", "recovery_any_targets",
        "non_target_component_rejection", "error_fraction",
        "effect_bias", "effect_rmse", "interval_coverage",
        "discoveries", "true_discoveries", "false_discoveries"
    )
    groups <- split(metrics, interaction(metrics$method, metrics$component, drop = TRUE))
    do.call(rbind, lapply(groups, function(x) {
        out <- x[1L, c(
            "setting_id", "method", "component", "evaluation_target",
            "n_taxa", "n_signals", "n_null"
        ), drop = FALSE]
        out$n_replications <- nrow(x)
        for (field in fields) {
            values <- x[[field]][is.finite(x[[field]])]
            name <- if (field == "fdp") "empirical_fdr" else field
            out[[name]] <- if (length(values)) mean(values) else NA_real_
            out[[paste0(name, "_mcse")]] <- if (length(values) > 1L) {
                sd(values) / sqrt(length(values))
            } else NA_real_
        }
        out
    }))
}

result_path <- function(setting, replication) {
    file.path(CONFIG$root, "results", setting$setting_id,
              sprintf("rep_%03d.rds", replication))
}

write_design <- function(grid, template) {
    safe_dir_create(CONFIG$root)
    utils::write.csv(grid, file.path(CONFIG$root, "design.csv"), row.names = FALSE)
    utils::write.csv(template, file.path(CONFIG$root, "taxon_template.csv"),
                     row.names = FALSE)
}

run_cell <- function(setting, replication, template, temporary_root, progress = NULL) {
    path <- result_path(setting, replication)
    progress_label <- if (is.null(progress)) "" else paste0(" | ", progress)
    if (file.exists(path)) {
        previous <- readRDS(path)
        stopifnot(identical(previous$setting, setting),
                  identical(previous$replication, as.integer(replication)))
        message(sprintf("%s | replication %d%s | already completed",
                        setting$setting_id, replication, progress_label))
        return(invisible(path))
    }
    started <- proc.time()[["elapsed"]]
    simulation <- simulate_metasparsim(replication, setting, template)
    result <- analyze_setting(simulation, replication, temporary_root)
    result <- add_omnibus(result)
    validate_analysis(result)
    metrics <- replication_metrics(result)
    save_result(list(
        setting = setting, replication = as.integer(replication),
        simulation = simulation, results = result, metrics = metrics,
        software = software_versions(),
        elapsed_seconds = proc.time()[["elapsed"]] - started
    ), path)
    message(sprintf("%s | replication %d%s | %.1f minutes",
                    setting$setting_id, replication, progress_label,
                    (proc.time()[["elapsed"]] - started) / 60))
    invisible(path)
}

run_replication <- function(task, grid, template) {
    if (length(task) != 1L || is.na(task) || task < 1L ||
        task > CONFIG$replications * CONFIG$setting_blocks) {
        stop("Array task must be between 1 and 200.")
    }
    check_packages()
    replication_index <- (task - 1L) %% CONFIG$replications + 1L
    block <- (task - 1L) %/% CONFIG$replications + 1L
    temporary_root <- file.path(CONFIG$root, "work", sprintf("task_%03d", task))
    safe_dir_create(temporary_root)
    if (task == 1L) write_design(grid, template)
    indices <- seq.int(block, nrow(grid), by = CONFIG$setting_blocks)
    total <- sum(vapply(indices, function(i) {
        length(seq.int(replication_index, grid$n_replications[i],
                       by = CONFIG$replications))
    }, integer(1L)))
    completed <- 0L
    for (i in indices) {
        setting <- grid[i, , drop = FALSE]
        n_reps <- if (setting$scenario == "global_null") {
            CONFIG$null_replications
        } else CONFIG$replications
        ids <- seq.int(replication_index, n_reps, by = CONFIG$replications)
        for (replication in ids) {
            completed <- completed + 1L
            run_cell(setting, replication, template, temporary_root,
                     progress = paste0(completed, "/", total))
        }
    }
    message("Array task ", task, " completed.")
}

stream_csv <- function(data, connection, header) {
    utils::write.table(data, connection, sep = ",", row.names = FALSE,
                       col.names = header, quote = TRUE, na = "")
}

archive_replications <- function(paths, output) {
    archive <- file.path(normalizePath(output), "replication_objects.tar")
    previous_directory <- getwd()
    on.exit(setwd(previous_directory), add = TRUE)
    setwd(file.path(CONFIG$root, "results"))
    utils::tar(archive, files = paths, compression = "none", tar = "internal")
}

summarize_results <- function(grid, template) {
    for (i in seq_len(nrow(grid))) {
        setting <- grid[i, , drop = FALSE]
        n_reps <- if (setting$scenario == "global_null") {
            CONFIG$null_replications
        } else CONFIG$replications
        missing <- which(!file.exists(vapply(seq_len(n_reps), function(r) {
            result_path(setting, r)
        }, character(1L))))
        if (length(missing)) stop(
            "Results still missing for ", setting$setting_id, ": replication ",
            paste(missing, collapse = ", "), call. = FALSE
        )
    }
    output <- file.path(CONFIG$root, "summary")
    safe_dir_create(output)
    streams <- lapply(c("taxon_results", "replication_metrics",
                       "generator_diagnostics", "availability_reasons"), function(name) {
        gzfile(file.path(output, paste0(name, ".csv.gz")), "wt")
    })
    names(streams) <- c("taxon", "replication", "generator", "availability")
    on.exit(lapply(streams, close), add = TRUE)
    summary_rows <- vector("list", nrow(grid))
    replication_index <- vector("list", sum(grid$n_replications))
    completed <- 0L
    for (i in seq_len(nrow(grid))) {
        setting <- grid[i, , drop = FALSE]
        n_reps <- if (setting$scenario == "global_null") {
            CONFIG$null_replications
        } else CONFIG$replications
        cell_metrics <- vector("list", n_reps)
        for (replication in seq_len(n_reps)) {
            saved <- readRDS(result_path(setting, replication))
            stopifnot(identical(saved$setting, setting),
                      identical(saved$replication, as.integer(replication)))
            completed <- completed + 1L
            replication_index[[completed]] <- data.frame(
                setting_id = setting$setting_id,
                replication = replication,
                result_file = file.path(setting$setting_id,
                                        sprintf("rep_%03d.rds", replication)),
                simulation_seed = CONFIG$simulation_seed + replication * 1000000L,
                method_seed = CONFIG$base_seed + replication * 100003L +
                    setting$setting_index * 10007L,
                elapsed_seconds = saved$elapsed_seconds
            )
            header <- i == 1L && replication == 1L
            stream_csv(saved$results, streams$taxon, header)
            stream_csv(saved$metrics, streams$replication, header)
            diagnostic <- saved$simulation$diagnostics
            diagnostic$setting_id <- setting$setting_id
            diagnostic$replication <- replication
            stream_csv(diagnostic, streams$generator, header)
            status <- data.table::as.data.table(saved$results)
            counts <- status[, .(n_taxa = .N),
                             by = .(setting_id, replication, method,
                                    component, available, status)]
            stream_csv(counts, streams$availability, header)
            cell_metrics[[replication]] <- saved$metrics
        }
        summary_rows[[i]] <- summarize_metrics(do.call(rbind, cell_metrics))
        if (i %% 20L == 0L) message("Summarized ", i, " / ", nrow(grid), " settings")
    }
    summary <- do.call(rbind, summary_rows)
    design_columns <- setdiff(names(grid), names(summary))
    summary <- merge(summary, grid[, c("setting_id", design_columns), drop = FALSE],
                     by = "setting_id", sort = FALSE)
    utils::write.csv(summary, file.path(output, "setting_summary.csv"),
                     row.names = FALSE, na = "")
    utils::write.csv(grid, file.path(output, "design.csv"), row.names = FALSE)
    utils::write.csv(template, file.path(output, "taxon_template.csv"),
                     row.names = FALSE)
    first <- readRDS(result_path(grid[1L, , drop = FALSE], 1L))
    utils::write.csv(first$software, file.path(output, "software_versions.csv"),
                     row.names = FALSE)
    index <- do.call(rbind, replication_index)
    utils::write.csv(index, file.path(output, "replication_index.csv"),
                     row.names = FALSE)
    saveRDS(CONFIG, file.path(output, "configuration.rds"))
    message("Archiving ", nrow(index), " complete replication objects")
    archive_replications(index$result_file, output)
    message("Summary saved to ", output)
}

run_smoke <- function(grid, template) {
    check_packages()
    chosen <- c(
        which(grid$scenario == "global_null"),
        which(grid$scenario == "joint" & grid$effect_index == 3L)
    )
    stopifnot(length(chosen) == 4L)
    output <- file.path(CONFIG$root, "smoke")
    safe_dir_create(output)
    for (i in chosen) {
        setting <- grid[i, , drop = FALSE]
        simulation <- simulate_metasparsim(1L, setting, template)
        repeated <- simulate_metasparsim(1L, setting, template)
        stopifnot(identical(simulation, repeated))
        result <- add_omnibus(analyze_setting(simulation, 1L, output))
        validate_analysis(result)
        save_result(list(
            setting = setting, simulation = simulation,
            results = result, metrics = replication_metrics(result)
        ), file.path(output, paste0(setting$setting_id, ".rds")))
    }
    message("Smoke test completed: generator, all methods, and BH summaries.")
}

main <- function() {
    arguments <- commandArgs(trailingOnly = TRUE)
    mode <- if (length(arguments)) arguments[1L] else "design"
    grid <- make_metasparsim_grid()
    template <- make_metasparsim_template()
    if (mode == "design") {
        write_design(grid, template)
        print(with(grid, table(scenario, include_covariate)))
        message("38 settings: 2 global nulls with 500 replications each; ",
                "36 alternatives with 100 replications each.")
    } else if (mode == "preflight") {
        check_packages()
        print(software_versions(), row.names = FALSE)
        message("Package checks passed.")
    } else if (mode == "smoke") {
        run_smoke(grid, template)
    } else if (mode == "replicate") {
        if (length(arguments) != 2L) stop("Supply the array task number.")
        run_replication(as.integer(arguments[2L]), grid, template)
    } else if (mode == "summarize") {
        summarize_results(grid, template)
    } else stop("Mode must be design, preflight, smoke, replicate, or summarize.")
}

if (sys.nframe() == 0L) main()
