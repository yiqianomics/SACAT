#!/usr/bin/env Rscript

# Joint SA and RA effects under the internal logit-normal model.
# Rscript run_joint.R replicate <id> | summarize | design | check

script_file <- local({
    paths <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
    if (length(paths)) tail(paths, 1L)[[1L]] else
        sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
})
study_dir <- dirname(dirname(normalizePath(script_file)))
options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1")

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) y else x
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

messagef <- function(fmt, ...) message(sprintf(fmt, ...))

clamp <- function(x, lower, upper) pmin(pmax(x, lower), upper)

expit <- function(x) stats::plogis(x)

logit <- function(p) stats::qlogis(clamp(p, 1e-12, 1 - 1e-12))

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

atomic_save_rds <- function(object, path, compress = "xz") {
    safe_dir_create(dirname(path))
    temporary <- paste0(
        path, ".tmp_", Sys.getpid(), "_", sample.int(.Machine$integer.max, 1L)
    )
    saveRDS(object, temporary, compress = compress)
    if (!file.rename(temporary, path)) {
        unlink(temporary, force = TRUE)
        stopf("Could not atomically write %s", path)
    }
    invisible(path)
}

package_version_or_na <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package))
}

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

normal_gh_rule <- local({
    cache <- new.env(parent = emptyenv())
    function(order = CONFIG$gh_order_truth) {
        key <- as.character(order)
        if (exists(key, envir = cache, inherits = FALSE)) {
            return(get(key, envir = cache, inherits = FALSE))
        }
        rule <- statmod::gauss.quad(order, kind = "hermite")
        value <- list(
            z = sqrt(2) * as.numeric(rule$nodes),
            w = as.numeric(rule$weights) / sqrt(pi)
        )
        assign(key, value, envir = cache)
        value
    }
})

expected_nondetection <- function(eta, sigma, depth,
                                  gh = normal_gh_rule()) {
    eta <- as.numeric(eta)
    depth <- as.numeric(depth)
    if (length(eta) == 1L) eta <- rep(eta, length(depth))
    if (length(eta) != length(depth)) {
        stop("eta and depth must have compatible lengths.", call. = FALSE)
    }
    vapply(seq_along(depth), function(i) {
        probability <- expit(eta[i] + sigma * gh$z)
        sum(gh$w * exp(depth[i] * log1p(-probability)))
    }, numeric(1))
}

expected_detection <- function(eta, sigma, depth,
                               gh = normal_gh_rule()) {
    clamp(1 - expected_nondetection(eta, sigma, depth, gh), 0, 1)
}

expected_log_probability <- function(eta, sigma,
                                     gh = normal_gh_rule()) {
    eta <- as.numeric(eta)
    vapply(eta, function(mu) {
        sum(gh$w * log(expit(mu + sigma * gh$z)))
    }, numeric(1))
}

solve_scalar_root <- function(objective, lower, upper,
                              tolerance = 1e-10) {
    f_lower <- objective(lower)
    f_upper <- objective(upper)
    if (!is.finite(f_lower) || !is.finite(f_upper) ||
        f_lower * f_upper > 0) {
        return(NA_real_)
    }
    stats::uniroot(
        objective, interval = c(lower, upper), tol = tolerance
    )$root
}

solve_structural_shift <- function(
        target_difference, baseline_rho, rho_z, standardized_z) {
    rho_zero <- expit(logit(baseline_rho) + rho_z * standardized_z)
    objective <- function(delta) {
        rho_one <- expit(
            logit(baseline_rho) + rho_z * standardized_z + delta
        )
        mean(rho_zero - rho_one) - target_difference
    }
    solve_scalar_root(objective, lower = -12, upper = 0)
}

solve_present_log_effect_shift <- function(
        target_effect, eta_base, sigma) {
    objective <- function(shift) {
        mean(
            expected_log_probability(eta_base + shift, sigma) -
                expected_log_probability(eta_base, sigma)
        ) - target_effect
    }
    solve_scalar_root(objective, lower = -8, upper = 8)
}

draw_multinomial_counts <- function(probability, depth, other_probability,
                                    seed) {
    n <- nrow(probability)
    J <- ncol(probability)
    set.seed(seed)
    count_by_sample <- vapply(seq_len(n), function(i) {
        as.numeric(stats::rmultinom(
            1L,
            size = depth[i],
            prob = c(probability[i, ], other_probability[i])
        )[, 1L])
    }, numeric(J + 1L))
    count_by_sample
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

analysis_terms <- function(simulation, group = "group", depth = FALSE) {
    c(
        group,
        if (isTRUE(simulation$setting$with_covariate)) "z",
        if (depth) "log_depth"
    )
}

analysis_formula <- function(simulation, response = NULL,
                             group = "group", depth = FALSE) {
    stats::reformulate(
        analysis_terms(simulation, group, depth),
        response = response
    )
}

analysis_formula_text <- function(simulation, depth = FALSE,
                                  include_tilde = TRUE) {
    paste0(
        if (include_tilde) "~ " else "",
        paste(analysis_terms(simulation, depth = depth), collapse = " + ")
    )
}

run_sacat <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture(
        SACAT::sacat(
            counts = simulation$counts[taxa, , drop = FALSE],
            metadata = simulation$metadata,
            formula = analysis_formula(simulation),
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
            ),
            method_result_template(
                taxa, "SACAT", "omnibus", status, captured$elapsed
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
    omnibus <- method_result_template(
        taxa, "SACAT", "omnibus", "unavailable", captured$elapsed
    )
    omnibus$p_value <- as.numeric(result$p_omnibus)
    omnibus$package_q_value <- as.numeric(result$p_adj_omnibus)
    omnibus$available <- valid_p(omnibus$p_value) &
        (structural$available | abundance$available)
    omnibus$status <- ifelse(
        omnibus$available, "ok",
        ifelse(
            !structural$available & !abundance$available,
            "both_components_unavailable", "nonfinite_p_value"
        )
    )
    rbind(structural, abundance, omnibus)
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

run_edger <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        design <- stats::model.matrix(
            analysis_formula(simulation), data = simulation$metadata
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
        col_data <- simulation$metadata[, analysis_terms(simulation), drop = FALSE]
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = round(simulation$counts),
            colData = col_data,
            design = if (isTRUE(simulation$setting$with_covariate)) {
                ~ z + group
            } else {
                ~ group
            }
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
            fix_formula = analysis_formula_text(simulation, include_tilde = FALSE),
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
            meta.dat = simulation$metadata[, analysis_terms(simulation), drop = FALSE],
            formula = analysis_formula_text(simulation),
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
            meta = simulation$metadata[, analysis_terms(simulation), drop = FALSE],
            formula = analysis_formula_text(simulation),
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
            group = simulation$metadata$group,
            z = design_covariate(simulation$metadata)
        )
        captured <- safe_capture(
            corncob::bbdml(
                formula = analysis_formula(simulation, "cbind(W, M - W)"),
                phi.formula = analysis_formula(simulation),
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
            simulation$metadata[, analysis_terms(simulation), drop = FALSE]
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
            analysis_formula(simulation), data = simulation$metadata
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

save_native_joint_output <- function(method, simulation, replication_id, value) {
    path <- file.path(
        CONFIG$root, "native_package_outputs",
        sprintf("replicate_%04d", replication_id),
        simulation$setting$setting_id, paste0(method, ".rds")
    )
    safe_dir_create(dirname(path))
    temporary <- paste0(path, ".tmp_", Sys.getpid())
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(list(
        method = method, replication = as.integer(replication_id),
        setting = simulation$setting, evaluation_taxa = simulation$evaluation_taxa,
        package_version = as.character(utils::packageVersion(
            if (method == "ZINQ") "ZINQ" else "maaslin3"
        )), output = value
    ), temporary, compress = "gzip")
    if (!file.rename(temporary, path)) stopf("Could not save native output: %s", path)
    invisible(path)
}

run_zinq <- function(simulation, replication_id) {
    taxa <- simulation$evaluation_taxa
    count <- simulation$counts[taxa, , drop = FALSE]
    relative_abundance <- sweep(count, 2L, simulation$depth, "/")
    prevalence_rows <- vector("list", length(taxa))
    abundance_rows <- vector("list", length(taxa))
    joint_rows <- vector("list", length(taxa))
    package_outputs <- stats::setNames(vector("list", length(taxa)), taxa)
    tau_weights <- pmin(CONFIG$zinq_taus, 1 - CONFIG$zinq_taus)
    tau_weights <- tau_weights / sum(tau_weights)

    started_all <- proc.time()[["elapsed"]]
    for (j in seq_along(taxa)) {
        data <- data.frame(
            y = as.numeric(relative_abundance[j, ]),
            group_num = simulation$metadata$group_num,
            z = design_covariate(simulation$metadata),
            log_depth = simulation$metadata$log_depth
        )
        set.seed(
            CONFIG$base_seed + replication_id * 100003L +
                simulation$setting$setting_index * 1009L + j
        )
        captured <- safe_capture(
            ZINQ::ZINQ_tests(
                formula.logistic = analysis_formula(
                    simulation, "y", "group_num", depth = TRUE
                ),
                formula.quantile = analysis_formula(
                    simulation, "y", "group_num", depth = TRUE
                ),
                C = "group_num",
                y_CorD = "C",
                data = data,
                taus = CONFIG$zinq_taus,
                seed = CONFIG$base_seed + replication_id * 1009L + j
            )
        )

        prevalence <- method_result_template(
            taxa[j], "ZINQ", "observed_prevalence", "unavailable", captured$elapsed
        )
        abundance <- method_result_template(
            taxa[j], "ZINQ", "detected_quantile_abundance", "unavailable", captured$elapsed
        )
        joint <- method_result_template(
            taxa[j], "ZINQ", "native_joint", "unavailable", captured$elapsed
        )
        warning_text <- collapse_text(captured$warnings)
        prevalence$warning <- abundance$warning <- joint$warning <- warning_text
        native_joint <- NULL

        if (is_sim_error(captured$value)) {
            status <- paste0("error: ", captured$value$message)
            prevalence$status <- abundance$status <- joint$status <- status
        } else {
            fit <- captured$value
            prevalence$p_value <- suppressWarnings(as.numeric(fit$pvalue.logistic)[1L])
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
                    abundance$p_value <- cauchy_combine(quantile_p, tau_weights)
                }
            }
            abundance$available <- valid_p(abundance$p_value)
            abundance$status <- ifelse(
                abundance$available, "ok", "quantile_component_unavailable"
            )

            # Combine prevalence and quantile evidence using the package's test.
            native_joint <- safe_capture(ZINQ::ZINQ_combination(
                fit, method = "Cauchy", taus = CONFIG$zinq_taus
            ))
            joint$runtime_seconds <- captured$elapsed + native_joint$elapsed
            joint$warning <- collapse_text(c(captured$warnings, native_joint$warnings))
            if (is_sim_error(native_joint$value)) {
                joint$status <- paste0("error: ", native_joint$value$message)
            } else {
                joint$p_value <- as.numeric(native_joint$value)[1L]
                joint$available <- valid_p(joint$p_value)
                joint$status <- ifelse(joint$available, "ok", "nonfinite_native_joint_p_value")
            }
        }
        package_outputs[[j]] <- list(fit = captured$value,
                                     warnings = captured$warnings,
                                     native_joint = native_joint)
        prevalence_rows[[j]] <- prevalence
        abundance_rows[[j]] <- abundance
        joint_rows[[j]] <- joint
    }
    output <- rbind(do.call(rbind, prevalence_rows), do.call(rbind, abundance_rows),
                    do.call(rbind, joint_rows))
    output$runtime_seconds_total_setting <- proc.time()[["elapsed"]] - started_all
    save_native_joint_output("ZINQ", simulation, replication_id, package_outputs)
    output
}

extract_maaslin_native_joint <- function(abundance_data, prevalence_data, taxa,
                                         elapsed, warnings) {
    output <- method_result_template(
        taxa, "MaAsLin 3", "native_joint", "row_unavailable", elapsed
    )
    output$warning <- collapse_text(warnings)
    for (j in seq_along(taxa)) {
        abundance <- select_maaslin_row(abundance_data, taxa[j], "abundance")
        prevalence <- select_maaslin_row(prevalence_data, taxa[j], "prevalence")
        candidates <- Filter(function(row) {
            !is.null(row) && "pval_joint" %in% names(row)
        }, list(abundance, prevalence))
        if (!length(candidates)) {
            output$status[j] <- "package_native_joint_row_unavailable"
            next
        }
        # Both component tables carry the package's joint p-value.
        if (length(candidates) == 2L &&
            !identical(as.numeric(candidates[[1L]]$pval_joint),
                       as.numeric(candidates[[2L]]$pval_joint))) {
            output$status[j] <- "package_native_joint_component_mismatch"
            next
        }
        row <- candidates[[1L]]
        output$p_value[j] <- as.numeric(row$pval_joint)[1L]
        if ("qval_joint" %in% names(row)) {
            output$package_q_value[j] <- as.numeric(row$qval_joint)[1L]
        }
        errors <- unlist(lapply(list(abundance, prevalence), function(part) {
            if (!is.null(part) && "error" %in% names(part)) as.character(part$error)
        }), use.names = FALSE)
        output$method_detail[j] <- collapse_text(errors)
        output$available[j] <- valid_p(output$p_value[j])
        output$status[j] <- if (output$available[j]) "ok" else "nonfinite_native_joint_p_value"
    }
    output
}

run_maaslin3 <- function(simulation, replication_id, temporary_root) {
    taxa <- simulation$evaluation_taxa
    output_dir <- file.path(temporary_root, paste0("maaslin3_", simulation$setting$setting_id))
    unlink(output_dir, recursive = TRUE, force = TRUE)
    safe_dir_create(output_dir)
    on.exit({
        try(maaslin3::maaslin_log_reset(), silent = TRUE)
        unlink(output_dir, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    input_data <- as.data.frame(t(simulation$counts))
    metadata <- simulation$metadata[, analysis_terms(simulation, depth = TRUE), drop = FALSE]
    arguments <- list(
        input_data = input_data, input_metadata = metadata, output = output_dir,
        formula = analysis_formula_text(simulation, depth = TRUE),
        min_abundance = 0, min_prevalence = 0, max_prevalence = 1.01,
        zero_threshold = 0, min_variance = 0, normalization = "TSS", transform = "LOG",
        correction = "BH", standardize = TRUE, median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        warn_prevalence = isTRUE(CONFIG$maaslin3_warn_prevalence), augment = TRUE,
        evaluate_only = NULL, plot_summary_plot = FALSE, plot_associations = FALSE,
        max_pngs = 1, cores = 1, save_models = FALSE, save_plots_rds = FALSE,
        verbosity = "ERROR", reference = "group,control"
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
        save_native_joint_output("MaAsLin3", simulation, replication_id,
                                 list(error = captured$value, warnings = captured$warnings))
        status <- paste0("error: ", captured$value$message)
        return(rbind(
            method_result_template(taxa, "MaAsLin 3", "observed_prevalence", status, captured$elapsed),
            method_result_template(taxa, "MaAsLin 3", "detected_log_abundance", status, captured$elapsed),
            method_result_template(taxa, "MaAsLin 3", "native_joint", status, captured$elapsed)
        ))
    }
    abundance_data <- maaslin_result_frame(captured$value, "abundance", output_dir)
    prevalence_data <- maaslin_result_frame(captured$value, "prevalence", output_dir)
    # Keep full result tables before temporary files are removed.
    save_native_joint_output("MaAsLin3", simulation, replication_id, list(
        abundance = abundance_data, prevalence = prevalence_data,
        warnings = captured$warnings
    ))
    rbind(
        extract_maaslin_component(prevalence_data, taxa, "prevalence",
                                 captured$elapsed, captured$warnings),
        extract_maaslin_component(abundance_data, taxa, "abundance",
                                 captured$elapsed, captured$warnings),
        extract_maaslin_native_joint(abundance_data, prevalence_data, taxa,
                                    captured$elapsed, captured$warnings)
    )
}

CONFIG <- list(
    root = Sys.getenv("SACAT_JOINT_ROOT", file.path(study_dir, "joint_output")),
    n_taxa = 50L,
    n_signal = 10L,
    n_replications = 100L,
    sample_sizes_per_group = c(60L, 80L, 120L),
    covariate_levels = c("without_covariate", "with_covariate"),
    structural_effects = c(0, 0.10, 0.15, 0.20, 0.25, 0.35, 0.50),
    abundance_effects = c(0, 0.20, 0.40, 0.50, 0.60, 0.80, 1.40),
    alpha = 0.05,
    base_seed = 202608190L,
    depth_median = 8000,
    depth_sdlog = 0.45,
    depth_min = 1500L,
    depth_max = 40000L,
    gh_order_truth = 81L,
    zinq_taus = c(0.25, 0.50, 0.75),
    maaslin3_warn_prevalence = TRUE,
    corncob_robust = TRUE,
    metagenomeseq_maxit = 10L
)

REQUIRED_PACKAGES <- c(
    "SACAT", "ZINQ", "maaslin3", "edgeR", "DESeq2", "ANCOMBC",
    "MicrobiomeStat", "corncob", "metagenomeSeq", "Biobase", "limma",
    "statmod", "data.table"
)

check_packages <- function() {
    missing <- REQUIRED_PACKAGES[!vapply(
        REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE
    )]
    if (length(missing)) stopf("Packages unavailable: %s", paste(missing, collapse = ", "))
    if (utils::packageVersion("SACAT") < "0.7.1") {
        stop("SACAT 0.7.1 or later is required.", call. = FALSE)
    }
    invisible(TRUE)
}

make_setting_grid <- function() {
    grid <- expand.grid(
        effect_index = seq_along(CONFIG$structural_effects) - 1L,
        n_per_group = CONFIG$sample_sizes_per_group,
        covariate = CONFIG$covariate_levels,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    grid$setting_index <- seq_len(nrow(grid))
    grid$setting_id <- sprintf("setting_%04d", grid$setting_index)
    grid$scenario <- "joint_effect"
    grid$structural_difference <- CONFIG$structural_effects[grid$effect_index + 1L]
    grid$abundance_log_effect <- CONFIG$abundance_effects[grid$effect_index + 1L]
    grid$with_covariate <- grid$covariate == "with_covariate"
    grid$total_sample_size <- 2L * grid$n_per_group
    grid$n_taxa <- CONFIG$n_taxa
    grid$n_signal_target <- ifelse(grid$effect_index == 0L, 0L, CONFIG$n_signal)
    grid$n_replications <- CONFIG$n_replications
    grid
}

make_joint_template <- function() {
    set.seed(CONFIG$base_seed + 1701L)
    J <- CONFIG$n_taxa
    probability <- exp(stats::runif(J, log(1e-4), log(2.8e-3)))
    template <- data.frame(
        taxon = sprintf("Taxon_%03d", seq_len(J)),
        baseline_probability = probability,
        baseline_rho = stats::runif(J, 0.15, 0.45),
        sigma = stats::runif(J, 0.50, 0.90),
        eta_z = stats::runif(J, -0.16, 0.16),
        rho_z = stats::runif(J, -0.22, 0.22)
    )
    pool <- seq_len(20L)
    template$baseline_probability[pool] <- exp(seq(
        log(4e-4), log(1.5e-3), length.out = length(pool)
    ))
    template$baseline_rho[pool] <- seq(0.55, 0.68, length.out = length(pool))
    template$sigma[pool] <- seq(0.52, 0.82, length.out = length(pool))
    template$baseline_eta <- logit(template$baseline_probability)
    template$designated_taxon <- seq_len(J) <= CONFIG$n_signal
    template$abundance_direction <- 0L
    template$abundance_direction[template$designated_taxon] <- rep(c(1L, -1L), 5L)
    template
}

make_replication_design <- function(replication_id, n_per_group, covariate) {
    n <- 2L * n_per_group
    set.seed(CONFIG$base_seed + replication_id * 100003L + n_per_group * 1009L + 11L)
    group_num <- rep(0:1, each = n_per_group)
    z_draw <- stats::rnorm(n)
    z <- if (covariate == "with_covariate") as.numeric(scale(z_draw)) else rep(0, n)
    depth_control <- as.integer(clamp(round(exp(stats::rnorm(
        n_per_group, log(CONFIG$depth_median), CONFIG$depth_sdlog
    ))), CONFIG$depth_min, CONFIG$depth_max))
    depth_case <- sample(depth_control, n_per_group, replace = FALSE)
    depth <- c(depth_control, depth_case)
    metadata <- data.frame(
        group = factor(ifelse(group_num == 0L, "control", "case"),
                       levels = c("control", "case")),
        group_num = group_num,
        z = z,
        reads = depth,
        log_depth = as.numeric(scale(log(depth))),
        row.names = sprintf("Sample_%03d", seq_len(n))
    )
    if (covariate == "without_covariate") metadata$z <- NULL
    metadata
}

design_covariate <- function(metadata) {
    if ("z" %in% names(metadata)) metadata$z else rep(0, nrow(metadata))
}

simulate_joint_effect <- function(replication_id, setting, template = make_joint_template()) {
    metadata <- make_replication_design(replication_id, setting$n_per_group, setting$covariate)
    z <- design_covariate(metadata)
    depth <- metadata$reads
    n <- nrow(metadata)
    J <- nrow(template)
    delta <- zeta <- numeric(J)
    signal <- template$designated_taxon & setting$effect_index > 0L
    for (j in which(signal)) {
        delta[j] <- solve_structural_shift(
            setting$structural_difference, template$baseline_rho[j], template$rho_z[j], z
        )
        zeta[j] <- solve_present_log_effect_shift(
            template$abundance_direction[j] * setting$abundance_log_effect,
            template$baseline_eta[j] + template$eta_z[j] * z, template$sigma[j]
        )
    }
    stopifnot(all(is.finite(delta)), all(is.finite(zeta)))
    eta0 <- outer(rep(1, n), template$baseline_eta) + outer(z, template$eta_z)
    rho0_linear <- outer(rep(1, n), logit(template$baseline_rho)) + outer(z, template$rho_z)
    eta <- eta0 + outer(metadata$group_num, zeta)
    rho <- expit(rho0_linear + outer(metadata$group_num, delta))
    set.seed(CONFIG$base_seed + replication_id * 100003L + 9000027L)
    absent <- matrix(stats::runif(n * J), n, J) < rho
    latent_normal <- matrix(stats::rnorm(n * J), n, J)
    latent_x <- eta + sweep(latent_normal, 2L, template$sigma, "*")
    probability <- expit(latent_x)
    probability[absent] <- 0
    background <- 1 - rowSums(probability)
    stopifnot(all(is.finite(probability)), all(background > 0))
    counts <- draw_multinomial_counts(
        probability, depth, background,
        CONFIG$base_seed + replication_id * 100003L + setting$setting_index * 1009L + 97L
    )
    dimnames(counts) <- list(c(template$taxon, "Other_unmodeled"), rownames(metadata))
    dimnames(absent) <- dimnames(latent_x) <- dimnames(probability) <-
        list(rownames(metadata), template$taxon)
    structural_effect <- abundance_effect <- prevalence_effect <- numeric(J)
    for (j in seq_len(J)) {
        rho0 <- expit(rho0_linear[, j])
        rho1 <- expit(rho0_linear[, j] + delta[j])
        structural_effect[j] <- mean(rho1 - rho0)
        abundance_effect[j] <- mean(
            expected_log_probability(eta0[, j] + zeta[j], template$sigma[j]) -
                expected_log_probability(eta0[, j], template$sigma[j])
        )
        prevalence_effect[j] <- mean(
            (1 - rho1) * expected_detection(eta0[, j] + zeta[j], template$sigma[j], depth) -
                (1 - rho0) * expected_detection(eta0[, j], template$sigma[j], depth)
        )
    }
    stopifnot(
        max(abs(structural_effect[signal] + setting$structural_difference), 0) < 1e-8,
        max(abs(abundance_effect[signal] -
                    template$abundance_direction[signal] * setting$abundance_log_effect), 0) < 1e-8,
        all(colSums(counts) == depth),
        all(t(counts[template$taxon, , drop = FALSE])[absent] == 0)
    )
    truth <- data.frame(
        taxon = template$taxon,
        designated_taxon = template$designated_taxon,
        truth_structural = signal,
        truth_abundance = signal,
        truth_any = signal,
        truth_observed_prevalence = abs(prevalence_effect) > 1e-8,
        structural_log_odds_effect = delta,
        conditional_logit_abundance_effect = zeta,
        standardized_structural_effect = structural_effect,
        standardized_abundance_log_effect = abundance_effect,
        standardized_prevalence_effect = prevalence_effect,
        joint_relationship = ifelse(!signal, "null", ifelse(
            structural_effect * abundance_effect < 0, "reinforcing", "opposing"
        )),
        baseline_probability = template$baseline_probability,
        baseline_structural_absence = template$baseline_rho,
        sigma = template$sigma
    )
    list(
        setting = setting, counts = counts, metadata = metadata, depth = depth,
        evaluation_taxa = template$taxon, truth = truth,
        structural_absence = absent, latent_x = latent_x,
        latent_probability = probability
    )
}

apply_common_multiplicity <- function(results) {
    data <- data.table::as.data.table(results)
    data[, available := !is.na(available) & available & valid_p(p_value)]
    data[, p_value_for_multiplicity := ifelse(available, p_value, 1)]
    data[, q_value_bh := stats::p.adjust(p_value_for_multiplicity, "BH"),
         by = .(method, component)]
    data[, reject_raw := p_value_for_multiplicity < CONFIG$alpha]
    data[, reject_bh := q_value_bh < CONFIG$alpha]
    as.data.frame(data)
}

analyze_setting <- function(simulation, replication_id, temporary_root) {
    setting <- simulation$setting
    seed <- CONFIG$base_seed + replication_id * 100003L + setting$setting_index * 10007L
    set.seed(seed + 1L)
    sacat <- run_sacat(simulation)
    set.seed(seed + 2L)
    zinq <- run_zinq(simulation, replication_id)
    set.seed(seed + 3L)
    maaslin <- run_maaslin3(simulation, replication_id, temporary_root)
    set.seed(seed + 101L)
    other <- run_generic_abundance_methods(simulation)
    results <- apply_common_multiplicity(rbind(sacat, zinq, maaslin, other))
    truth <- simulation$truth[match(results$taxon, simulation$truth$taxon), ]
    for (column in setdiff(names(truth), "taxon")) results[[column]] <- truth[[column]]
    results$truth_for_component <- ifelse(
        results$component %in% c("omnibus", "native_joint", "general_abundance"),
        results$truth_any,
        ifelse(results$component == "structural_absence", results$truth_structural,
               ifelse(results$component == "observed_prevalence",
                      results$truth_observed_prevalence, results$truth_abundance))
    )
    results$primary_comparison <- results$component %in%
        c("omnibus", "native_joint", "general_abundance")
    results$replication <- replication_id
    for (column in names(setting)) results[[column]] <- setting[[column]][1L]
    data <- data.table::as.data.table(results)
    failed_methods <- data[, .(failed = all(grepl("^error:", status))), by = method]
    if (any(failed_methods$failed)) {
        stopf("Method fitting failed: %s", paste(failed_methods$method[failed_methods$failed], collapse = ", "))
    }
    native_failures <- data[component == "native_joint", .(failed = all(
        status %in% c("package_native_joint_row_unavailable", "package_native_joint_component_mismatch")
    )), by = method]
    if (any(native_failures$failed)) {
        stopf("Native joint output unavailable: %s", paste(native_failures$method[native_failures$failed], collapse = ", "))
    }
    stopifnot(
        nrow(data) == 15L * CONFIG$n_taxa,
        !anyDuplicated(data[, .(method, component, taxon)]),
        all(data[, .N, by = .(method, component)]$N == CONFIG$n_taxa),
        data[primary_comparison == TRUE, data.table::uniqueN(method)] == 9L,
        !anyNA(data$truth_for_component)
    )
    results
}

replication_metrics <- function(results) {
    data <- data.table::as.data.table(results)
    data[, .(
        n_taxa = .N, n_signals = sum(truth_for_component),
        n_null = sum(!truth_for_component), n_available = sum(available),
        availability = mean(available),
        signal_availability = if (any(truth_for_component)) mean(available[truth_for_component]) else NA_real_,
        raw_true_discoveries = sum(reject_raw & truth_for_component),
        raw_false_discoveries = sum(reject_raw & !truth_for_component),
        raw_power = if (any(truth_for_component)) mean(reject_raw[truth_for_component]) else NA_real_,
        raw_type1_error = mean(reject_raw[!truth_for_component]),
        discoveries = sum(reject_bh),
        true_discoveries = sum(reject_bh & truth_for_component),
        false_discoveries = sum(reject_bh & !truth_for_component),
        power = if (any(truth_for_component)) mean(reject_bh[truth_for_component]) else NA_real_,
        fdp = sum(reject_bh & !truth_for_component) / max(sum(reject_bh), 1L),
        bh_null_rejection_rate = mean(reject_bh[!truth_for_component])
    ), by = .(setting_id, replication, method, component, primary_comparison)]
}

run_replication <- function(replication_id) {
    check_packages()
    stopifnot(replication_id >= 1L, replication_id <= CONFIG$n_replications)
    grid <- make_setting_grid()
    selection <- Sys.getenv("SACAT_JOINT_SETTINGS", "")
    if (nzchar(selection)) {
        requested <- strsplit(selection, ",", fixed = TRUE)[[1L]]
        grid <- grid[grid$setting_id %in% requested | as.character(grid$setting_index) %in% requested, ]
        if (!nrow(grid)) stop("SACAT_JOINT_SETTINGS selected no settings.", call. = FALSE)
    }
    template <- make_joint_template()
    scratch_parent <- file.path(CONFIG$root, "tmp")
    safe_dir_create(scratch_parent)
    temporary_root <- tempfile(sprintf("rep_%03d_", replication_id), tmpdir = scratch_parent)
    safe_dir_create(temporary_root)
    on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)
    for (i in seq_len(nrow(grid))) {
        setting <- grid[i, , drop = FALSE]
        path <- file.path(CONFIG$root, "results", setting$setting_id,
                          sprintf("rep_%03d.rds", replication_id))
        if (file.exists(path)) {
            saved <- readRDS(path)
            config_fields <- setdiff(names(CONFIG), "root")
            if (isTRUE(saved$success) && identical(saved$setting, setting) &&
                identical(saved$configuration[config_fields], CONFIG[config_fields])) {
                messagef("%s | replication %d | already complete", setting$setting_id, replication_id)
                next
            }
            stopf("An incomplete or different result exists: %s", path)
        }
        started <- proc.time()[["elapsed"]]
        simulation <- simulate_joint_effect(replication_id, setting, template)
        results <- analyze_setting(simulation, replication_id, temporary_root)
        versions <- data.frame(
            package = c("R", REQUIRED_PACKAGES),
            version = c(as.character(getRversion()), vapply(
                REQUIRED_PACKAGES, package_version_or_na, character(1)
            ))
        )
        atomic_save_rds(list(
            success = TRUE, replication = replication_id, setting = setting,
            configuration = CONFIG, simulation = simulation,
            method_results = results, metrics = replication_metrics(results),
            software_versions = versions,
            elapsed_seconds = proc.time()[["elapsed"]] - started
        ), path, compress = "gzip")
        messagef("%s | replication %d | %.1f minutes", setting$setting_id,
                 replication_id, (proc.time()[["elapsed"]] - started) / 60)
    }
    invisible(NULL)
}

mean_if_present <- function(x) {
    if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
}

mcse <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
}

summarize_results <- function() {
    if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required.")
    data.table::setDTthreads(1L)
    grid <- make_setting_grid()
    output <- file.path(CONFIG$root, "summary")
    safe_dir_create(output)
    collected <- truth_rows <- metrics <- versions <- list()
    status <- expand.grid(setting_id = grid$setting_id,
                          replication = seq_len(CONFIG$n_replications),
                          stringsAsFactors = FALSE)
    status$completed <- FALSE
    for (i in seq_len(nrow(status))) {
        path <- file.path(CONFIG$root, "results", status$setting_id[i],
                          sprintf("rep_%03d.rds", status$replication[i]))
        if (!file.exists(path)) next
        saved <- readRDS(path)
        stopifnot(isTRUE(saved$success), saved$replication == status$replication[i],
                  saved$setting$setting_id == status$setting_id[i])
        status$completed[i] <- TRUE
        collected[[length(collected) + 1L]] <- saved$method_results
        metrics[[length(metrics) + 1L]] <- saved$metrics
        truth <- saved$simulation$truth
        truth$setting_id <- status$setting_id[i]
        truth$replication <- status$replication[i]
        truth_rows[[length(truth_rows) + 1L]] <- truth
        version <- saved$software_versions
        version$setting_id <- status$setting_id[i]
        version$replication <- status$replication[i]
        versions[[length(versions) + 1L]] <- version
    }
    data.table::fwrite(status, file.path(output, "run_status.csv"))
    if (!length(collected)) stop("No completed results found.", call. = FALSE)
    data <- data.table::rbindlist(collected, use.names = TRUE)
    metric <- data.table::rbindlist(metrics, use.names = TRUE)
    stopifnot(!anyDuplicated(data[, .(setting_id, replication, method, component, taxon)]))
    summary <- metric[, .(
        completed_replications = .N,
        n_taxa = unique(n_taxa), n_signals = unique(n_signals),
        power = mean_if_present(power), power_mcse = mcse(power),
        empirical_fdr = mean(fdp), fdr_mcse = mcse(fdp),
        raw_power = mean_if_present(raw_power), raw_power_mcse = mcse(raw_power),
        raw_type1_error = mean(raw_type1_error), raw_type1_mcse = mcse(raw_type1_error),
        bh_null_rejection_rate = mean(bh_null_rejection_rate),
        availability = mean(availability), availability_mcse = mcse(availability),
        signal_availability = mean_if_present(signal_availability),
        mean_discoveries = mean(discoveries),
        mean_true_discoveries = mean(true_discoveries),
        mean_false_discoveries = mean(false_discoveries)
    ), by = .(setting_id, method, component, primary_comparison)]
    summary <- merge(summary, grid[, setdiff(names(grid), "n_taxa")],
                     by = "setting_id", sort = FALSE)
    by_relationship <- data[joint_relationship != "null" & truth_for_component == TRUE, .(
        n_targets = .N, power = mean(reject_bh), raw_power = mean(reject_raw),
        availability = mean(available)
    ), by = .(setting_id, replication, method, component, joint_relationship)]
    relationship_summary <- by_relationship[, .(
        completed_replications = .N, n_targets = unique(n_targets),
        power = mean(power), power_mcse = mcse(power),
        raw_power = mean(raw_power), availability = mean(availability)
    ), by = .(setting_id, method, component, joint_relationship)]
    availability <- data[, .N, by = .(setting_id, method, component, available, status)]
    availability[, proportion := N / sum(N), by = .(setting_id, method, component)]
    data.table::fwrite(grid, file.path(output, "design.csv"))
    data.table::fwrite(make_joint_template(), file.path(output, "taxon_template.csv"))
    data.table::fwrite(data, file.path(output, "taxon_results.csv.gz"))
    data.table::fwrite(metric, file.path(output, "replication_metrics.csv.gz"))
    data.table::fwrite(summary, file.path(output, "setting_summary.csv"))
    data.table::fwrite(summary[primary_comparison == TRUE], file.path(output, "joint_comparison.csv"))
    data.table::fwrite(summary[effect_index == 0L], file.path(output, "global_null_type1.csv"))
    data.table::fwrite(relationship_summary, file.path(output, "power_by_relationship.csv"))
    data.table::fwrite(data.table::rbindlist(truth_rows), file.path(output, "simulation_truth.csv.gz"))
    data.table::fwrite(availability, file.path(output, "availability_by_reason.csv"))
    data.table::fwrite(data.table::rbindlist(versions), file.path(output, "software_versions.csv.gz"))
    messagef("Summarized %d/%d setting-replication results in %s",
             sum(status$completed), nrow(status), output)
    invisible(summary)
}

main <- function() {
    arguments <- commandArgs(trailingOnly = TRUE)
    if (!length(arguments)) stop(
        "Usage: Rscript run_internal_joint.R replicate <id> | summarize | design | check",
        call. = FALSE
    )
    mode <- arguments[1L]
    if (mode == "replicate") {
        if (length(arguments) != 2L) stop("Supply a replication ID.", call. = FALSE)
        id <- suppressWarnings(as.integer(arguments[2L]))
        if (is.na(id)) stop("The replication ID must be an integer.", call. = FALSE)
        run_replication(id)
    } else if (mode == "summarize") {
        summarize_results()
    } else if (mode == "design") {
        print(make_setting_grid(), row.names = FALSE)
    } else if (mode == "check") {
        check_packages()
        print(data.frame(package = REQUIRED_PACKAGES,
                         version = vapply(REQUIRED_PACKAGES, package_version_or_na, character(1))))
    } else stopf("Unknown mode: %s", mode)
}

if (sys.nframe() == 0L) main()
