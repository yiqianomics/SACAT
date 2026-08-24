#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

EXPECTED_REPLICATIONS <- 1:100
EXPECTED_SETTINGS <- c("NULL_BALANCED", "NULL_BALANCED_LOW", "NULL_FOURFOLD")
ALPHA <- 0.05
TAIL_THRESHOLDS <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.10)

EXPECTED_FAMILIES <- data.frame(
    method = c("DASRA", "DASRA", "ZINQ", "ZINQ", "MaAsLin 3", "MaAsLin 3"),
    component = c(
        "structural_absence", "present_conditional_abundance",
        "observed_prevalence", "detected_quantile_abundance",
        "observed_prevalence", "detected_log_abundance"
    ),
    domain = c("prevalence", "abundance", "prevalence", "abundance",
               "prevalence", "abundance"),
    stringsAsFactors = FALSE
)

assert <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
    invisible(TRUE)
}

valid_p <- function(p) is.finite(p) & p >= 0 & p <= 1

cluster_t <- function(values, bounds = c(0, 1)) {
    values <- as.numeric(values)
    assert(length(values) == length(EXPECTED_REPLICATIONS) && all(is.finite(values)),
           "A cluster-t summary does not contain one finite value for each replication 1:100.")
    estimate <- mean(values)
    mcse <- stats::sd(values) / sqrt(length(values))
    half <- stats::qt(0.975, df = length(values) - 1L) * mcse
    c(
        n = length(values), estimate = estimate, mcse = mcse,
        ci95_lower = max(bounds[1L], estimate - half),
        ci95_upper = min(bounds[2L], estimate + half)
    )
}

wilson <- function(events) {
    events <- as.integer(events)
    assert(length(events) == length(EXPECTED_REPLICATIONS) &&
               all(events %in% c(0L, 1L)),
           "A Wilson summary does not contain 100 binary replication-level events.")
    n <- length(events)
    count <- sum(events)
    estimate <- count / n
    z <- stats::qnorm(0.975)
    denominator <- 1 + z^2 / n
    center <- (estimate + z^2 / (2 * n)) / denominator
    half <- z * sqrt(estimate * (1 - estimate) / n + z^2 / (4 * n^2)) /
        denominator
    c(
        n = n, events = count, estimate = estimate,
        mcse = stats::sd(events) / sqrt(n),
        ci95_lower = max(0, center - half),
        ci95_upper = min(1, center + half)
    )
}

prefix_values <- function(values, prefix) {
    stats::setNames(values, paste0(prefix, "_", names(values)))
}

atomic_write_csv <- function(object, path) {
    temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    utils::write.csv(object, temporary, row.names = FALSE, na = "")
    assert(file.rename(temporary, path), sprintf("Could not write %s atomically.", path))
    invisible(path)
}

calibration_role <- function(setting_id) {
    switch(
        setting_id,
        NULL_BALANCED = "primary exact global null; balanced standard depth",
        NULL_BALANCED_LOW = "exact global null; balanced low-depth support stress",
        NULL_FOURFOLD = "no-direct-effect conditional depth-imbalance stress",
        stop("Unknown null setting.", call. = FALSE)
    )
}

interpretation_scope <- function(setting_id) {
    if (setting_id == "NULL_FOURFOLD") {
        "conditional_depth_imbalance_stress"
    } else {
        "exact_global_null"
    }
}

family_key <- function(method, component) paste(method, component, sep = "\r")

display_method <- function(method) {
    ifelse(method == "MaAsLin 3", "MaAsLin3", method)
}

summarize_main <- function(module_root, output_dir) {
    module_root <- normalizePath(module_root, winslash = "/", mustWork = TRUE)
    design_dir <- file.path(module_root, "design")
    raw_dir <- file.path(module_root, "results", "raw")
    contract_path <- file.path(design_dir, "contract.rds")
    assert(file.exists(contract_path), "design/contract.rds is missing.")
    assert(dir.exists(raw_dir), "results/raw is missing.")
    if (!dir.exists(output_dir)) {
        assert(dir.create(output_dir, recursive = TRUE, showWarnings = FALSE),
               "Could not create the requested output directory.")
    }
    output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

    contract <- readRDS(contract_path)
    settings <- contract$settings
    panel <- as.character(contract$panel$taxon)
    assert(identical(as.character(settings$setting_id), EXPECTED_SETTINGS),
           "The n=50 design does not contain the three prespecified null settings in order.")
    assert(all(settings$n_per_group == 50L) && all(settings$total_sample_size == 100L),
           "The design is not exactly 50 samples per group.")
    assert(all(settings$signal_count == 0L) &&
               !any(settings$spike_structural) && !any(settings$spike_abundance),
           "A requested setting is not a true no-direct-spike setting.")
    assert(!any(settings$generic_abundance),
           "The n=50 null run contains methods outside the prespecified three-method comparison.")
    assert(length(panel) == 50L && !anyDuplicated(panel),
           "The frozen evaluation panel is not 50 unique taxa.")
    dasra_version <- as.character(
        contract$package_provenance$version[
            contract$package_provenance$package == "DASRA"
        ]
    )
    assert(identical(dasra_version, "0.4.1"),
           "The frozen n=50 contract is not for DASRA 0.4.1.")

    expected_names <- sprintf("replicate_%04d.rds", EXPECTED_REPLICATIONS)
    observed_names <- sort(list.files(raw_dir, pattern = "^replicate_[0-9]{4}\\.rds$"))
    assert(identical(observed_names, expected_names),
           "The raw cohort must be exactly replicate_0001.rds through replicate_0100.rds.")

    result_parts <- vector("list", length(EXPECTED_REPLICATIONS))
    support_parts <- vector("list", length(EXPECTED_REPLICATIONS) * length(EXPECTED_SETTINGS))
    support_index <- 0L
    decision_replay_count <- 0L
    generator_check_count <- 0L

    expected_pair_keys <- sort(family_key(EXPECTED_FAMILIES$method,
                                          EXPECTED_FAMILIES$component))
    required_result_columns <- c(
        "setting_id", "scenario", "design_id", "n_per_group", "taxon",
        "method", "component", "p_value", "available", "status",
        "q_value_bh", "q_value_by", "reject_raw", "reject_bh", "reject_by",
        "replication", "truth_structural", "truth_observed_prevalence",
        "truth_abundance", "truth_for_component"
    )

    for (replication in EXPECTED_REPLICATIONS) {
        record <- readRDS(file.path(raw_dir, sprintf("replicate_%04d.rds", replication)))
        assert(is.list(record) && isTRUE(record$success),
               sprintf("Replication %d is not a successful aggregate.", replication))
        assert(identical(record$analysis_signature, contract$analysis_signature),
               sprintf("Replication %d has the wrong design signature.", replication))
        assert(identical(record$replication, as.integer(replication)),
               sprintf("Replication %d has the wrong identifier.", replication))
        assert(identical(as.character(record$setting_ids), EXPECTED_SETTINGS) &&
                   identical(names(record$datasets), EXPECTED_SETTINGS),
               sprintf("Replication %d has the wrong setting family or order.", replication))

        results <- record$method_results
        assert(is.data.frame(results) && all(required_result_columns %in% names(results)),
               sprintf("Replication %d has an incomplete method-result schema.", replication))
        assert(nrow(results) == 3L * 6L * 50L,
               sprintf("Replication %d does not contain 900 method-result rows.", replication))
        assert(all(results$replication == replication) &&
                   identical(sort(unique(as.character(results$setting_id))),
                             sort(EXPECTED_SETTINGS)),
               sprintf("Replication %d has inconsistent row identifiers.", replication))
        assert(!anyNA(results[, c(
            "setting_id", "taxon", "method", "component", "available", "status",
            "q_value_bh", "q_value_by", "reject_raw", "reject_bh", "reject_by",
            "truth_structural", "truth_observed_prevalence", "truth_abundance",
            "truth_for_component"
        )]), sprintf("Replication %d has missing required method fields.", replication))
        assert(!any(as.logical(results$truth_structural)) &&
                   !any(as.logical(results$truth_observed_prevalence)) &&
                   !any(as.logical(results$truth_abundance)) &&
                   !any(as.logical(results$truth_for_component)),
               sprintf("Replication %d is not a global no-direct-spike null.", replication))

        record_key <- paste(results$setting_id, results$method, results$component,
                            results$taxon, sep = "\r")
        assert(!anyDuplicated(record_key),
               sprintf("Replication %d contains duplicate result keys.", replication))

        for (setting_id in EXPECTED_SETTINGS) {
            setting <- settings[settings$setting_id == setting_id, , drop = FALSE]
            setting_results <- results[results$setting_id == setting_id, , drop = FALSE]
            observed_pairs <- sort(unique(family_key(setting_results$method,
                                                     setting_results$component)))
            assert(identical(observed_pairs, expected_pair_keys),
                   sprintf("Replication %d setting %s has an unexpected method family.",
                           replication, setting_id))
            assert(all(setting_results$n_per_group == 50L) &&
                       all(setting_results$scenario == "global_null") &&
                       all(setting_results$design_id == setting$design_id),
                   sprintf("Replication %d setting %s has inconsistent design context.",
                           replication, setting_id))

            for (family_index in seq_len(nrow(EXPECTED_FAMILIES))) {
                method <- EXPECTED_FAMILIES$method[family_index]
                component <- EXPECTED_FAMILIES$component[family_index]
                keep <- setting_results$method == method &
                    setting_results$component == component
                family <- setting_results[keep, , drop = FALSE]
                assert(nrow(family) == 50L && !anyDuplicated(family$taxon) &&
                           setequal(as.character(family$taxon), panel),
                       sprintf("Replication %d setting %s family %s/%s changed taxa.",
                               replication, setting_id, method, component))
                p <- as.numeric(family$p_value)
                available <- as.logical(family$available)
                p_is_valid <- valid_p(p)
                assert(all(!available | p_is_valid),
                       sprintf("Replication %d setting %s family %s/%s marks an invalid p-value available.",
                               replication, setting_id, method, component))
                p_for_multiplicity <- ifelse(p_is_valid, p, 1)
                q_bh <- stats::p.adjust(p_for_multiplicity, method = "BH")
                q_by <- stats::p.adjust(p_for_multiplicity, method = "BY")
                assert(max(abs(as.numeric(family$q_value_bh) - q_bh)) < 1e-14 &&
                           max(abs(as.numeric(family$q_value_by) - q_by)) < 1e-14,
                       sprintf("Replication %d setting %s family %s/%s failed BH/BY replay.",
                               replication, setting_id, method, component))
                assert(identical(as.logical(family$reject_raw),
                                 available & p_is_valid & p < ALPHA) &&
                           identical(as.logical(family$reject_bh), q_bh < ALPHA) &&
                           identical(as.logical(family$reject_by), q_by < ALPHA) &&
                           !any(!available & (as.logical(family$reject_raw) |
                                              as.logical(family$reject_bh) |
                                              as.logical(family$reject_by))),
                       sprintf("Replication %d setting %s family %s/%s failed decision replay.",
                               replication, setting_id, method, component))
                decision_replay_count <- decision_replay_count + 1L
            }

            dataset <- record$datasets[[setting_id]]
            assert(is.list(dataset) && identical(as.character(dataset$evaluation_taxa), panel),
                   sprintf("Replication %d setting %s has the wrong frozen panel.",
                           replication, setting_id))
            groups <- table(dataset$metadata$group)
            assert(identical(as.integer(groups[c("control", "case")]), c(50L, 50L)),
                   sprintf("Replication %d setting %s is not 50/50.", replication, setting_id))
            assert(ncol(dataset$counts) == 100L &&
                       identical(as.integer(colSums(dataset$counts)),
                                 as.integer(dataset$depth)),
                   sprintf("Replication %d setting %s has inconsistent library totals.",
                           replication, setting_id))
            assert(is.data.frame(dataset$active_signals) && nrow(dataset$active_signals) == 0L,
                   sprintf("Replication %d setting %s contains active signals.",
                           replication, setting_id))
            generator <- dataset$generator_contract
            assert(identical(as.integer(generator$replication), as.integer(replication)) &&
                       identical(as.character(generator$setting_id), setting_id) &&
                       generator$max_probability_sum_error_before_adjustment < 1e-12 &&
                       identical(as.numeric(generator$max_inactive_tested_probability_error), 0),
                   sprintf("Replication %d setting %s failed generator invariants.",
                           replication, setting_id))
            generator_check_count <- generator_check_count + 1L

            support <- dataset$positive_counts
            assert(is.data.frame(support) &&
                       identical(as.character(support$taxon), panel) &&
                       all(support$positive_control >= 0 & support$positive_control <= 50) &&
                       all(support$positive_case >= 0 & support$positive_case <= 50),
                   sprintf("Replication %d setting %s has invalid positive support.",
                           replication, setting_id))
            support_index <- support_index + 1L
            support_parts[[support_index]] <- data.frame(
                replication = replication, setting_id = setting_id,
                design_id = setting$design_id, taxon = as.character(support$taxon),
                positive_control = as.integer(support$positive_control),
                positive_case = as.integer(support$positive_case),
                stringsAsFactors = FALSE
            )
        }
        result_parts[[replication]] <- results[, required_result_columns, drop = FALSE]
    }

    results <- do.call(rbind, result_parts)
    support <- do.call(rbind, support_parts)
    rownames(results) <- rownames(support) <- NULL
    assert(nrow(results) == 100L * 3L * 6L * 50L,
           "The combined method-result table is not 100 x 3 x 6 x 50.")
    assert(nrow(support) == 100L * 3L * 50L,
           "The combined support table is not 100 x 3 x 50.")

    family_metric_rows <- vector("list", 100L * 3L * 6L)
    tail_rows <- vector("list", 100L * 3L * 6L * length(TAIL_THRESHOLDS))
    metric_index <- 0L
    tail_index <- 0L
    for (replication in EXPECTED_REPLICATIONS) {
        for (setting_id in EXPECTED_SETTINGS) {
            for (family_index in seq_len(nrow(EXPECTED_FAMILIES))) {
                method <- EXPECTED_FAMILIES$method[family_index]
                component <- EXPECTED_FAMILIES$component[family_index]
                keep <- results$replication == replication &
                    results$setting_id == setting_id &
                    results$method == method & results$component == component
                family <- results[keep, , drop = FALSE]
                assert(nrow(family) == 50L, "A replication-level family is incomplete.")
                available <- as.logical(family$available)
                assert(any(available),
                       "Available-conditional calibration is undefined for a zero-availability family.")
                p <- as.numeric(family$p_value)
                p_operational <- ifelse(available & valid_p(p), p, 1)
                metric_index <- metric_index + 1L
                family_metric_rows[[metric_index]] <- data.frame(
                    replication = replication, setting_id = setting_id,
                    method = method, component = component,
                    domain = EXPECTED_FAMILIES$domain[family_index],
                    availability = mean(available),
                    raw_rejection = mean(as.logical(family$reject_raw)),
                    available_conditional_raw = mean(as.logical(family$reject_raw)[available]),
                    bh_fwer_event = as.integer(any(as.logical(family$reject_bh))),
                    by_fwer_event = as.integer(any(as.logical(family$reject_by))),
                    stringsAsFactors = FALSE
                )
                for (threshold in TAIL_THRESHOLDS) {
                    tail_index <- tail_index + 1L
                    tail_rows[[tail_index]] <- data.frame(
                        replication = replication, setting_id = setting_id,
                        method = method, component = component,
                        domain = EXPECTED_FAMILIES$domain[family_index],
                        threshold = threshold,
                        operational_tail = mean(p_operational <= threshold),
                        available_conditional_tail = mean(p[available] <= threshold),
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }
    family_metrics <- do.call(rbind, family_metric_rows)
    tails <- do.call(rbind, tail_rows)

    null_rows <- vector("list", 3L * 6L)
    null_index <- 0L
    for (setting_id in EXPECTED_SETTINGS) {
        setting <- settings[settings$setting_id == setting_id, , drop = FALSE]
        for (family_index in seq_len(nrow(EXPECTED_FAMILIES))) {
            method <- EXPECTED_FAMILIES$method[family_index]
            component <- EXPECTED_FAMILIES$component[family_index]
            values <- family_metrics[
                family_metrics$setting_id == setting_id &
                    family_metrics$method == method &
                    family_metrics$component == component, , drop = FALSE
            ]
            values <- values[order(values$replication), , drop = FALSE]
            assert(identical(as.integer(values$replication), EXPECTED_REPLICATIONS),
                   "A null-calibration family is not complete for replications 1:100.")
            null_index <- null_index + 1L
            null_rows[[null_index]] <- data.frame(
                setting_id = setting_id, scenario = interpretation_scope(setting_id),
                depth_design = as.character(setting$depth_design),
                design_id = as.character(setting$design_id),
                calibration_role = calibration_role(setting_id),
                method = display_method(method), component = component,
                domain = EXPECTED_FAMILIES$domain[family_index], n_taxa = 50L,
                as.list(prefix_values(cluster_t(values$availability), "availability")),
                as.list(prefix_values(cluster_t(values$raw_rejection), "raw_rejection")),
                as.list(prefix_values(
                    cluster_t(values$available_conditional_raw),
                    "available_conditional_raw"
                )),
                as.list(prefix_values(wilson(values$bh_fwer_event), "bh_fwer")),
                as.list(prefix_values(wilson(values$by_fwer_event), "by_fwer")),
                stringsAsFactors = FALSE, check.names = FALSE
            )
        }
    }
    null_calibration <- do.call(rbind, null_rows)

    tail_summary_rows <- vector("list", 3L * 6L * length(TAIL_THRESHOLDS))
    summary_index <- 0L
    for (setting_id in EXPECTED_SETTINGS) {
        setting <- settings[settings$setting_id == setting_id, , drop = FALSE]
        for (family_index in seq_len(nrow(EXPECTED_FAMILIES))) {
            method <- EXPECTED_FAMILIES$method[family_index]
            component <- EXPECTED_FAMILIES$component[family_index]
            for (threshold in TAIL_THRESHOLDS) {
                values <- tails[
                    tails$setting_id == setting_id & tails$method == method &
                        tails$component == component & tails$threshold == threshold,
                    , drop = FALSE
                ]
                values <- values[order(values$replication), , drop = FALSE]
                assert(identical(as.integer(values$replication), EXPECTED_REPLICATIONS),
                       "A p-value-tail family is not complete for replications 1:100.")
                summary_index <- summary_index + 1L
                tail_summary_rows[[summary_index]] <- data.frame(
                    setting_id = setting_id,
                    depth_design = as.character(setting$depth_design),
                    design_id = as.character(setting$design_id),
                    calibration_role = calibration_role(setting_id),
                    method = display_method(method), component = component,
                    domain = EXPECTED_FAMILIES$domain[family_index],
                    threshold = threshold, n_taxa = 50L,
                    as.list(prefix_values(
                        cluster_t(values$operational_tail), "operational_tail"
                    )),
                    as.list(prefix_values(
                        cluster_t(values$available_conditional_tail),
                        "available_conditional_tail"
                    )),
                    stringsAsFactors = FALSE, check.names = FALSE
                )
            }
        }
    }
    pvalue_tail <- do.call(rbind, tail_summary_rows)

    support_rows <- vector("list", 3L * 50L)
    support_summary_index <- 0L
    for (setting_id in EXPECTED_SETTINGS) {
        setting <- settings[settings$setting_id == setting_id, , drop = FALSE]
        for (taxon in panel) {
            values <- support[support$setting_id == setting_id & support$taxon == taxon,
                              , drop = FALSE]
            values <- values[order(values$replication), , drop = FALSE]
            assert(identical(as.integer(values$replication), EXPECTED_REPLICATIONS),
                   "A positive-support taxon is not complete for replications 1:100.")
            control <- values$positive_control
            case <- values$positive_case
            control_mean <- cluster_t(control, c(0, 50))
            case_mean <- cluster_t(case, c(0, 50))
            support_summary_index <- support_summary_index + 1L
            support_rows[[support_summary_index]] <- data.frame(
                setting_id = setting_id,
                depth_design = as.character(setting$depth_design),
                design_id = as.character(setting$design_id), taxon = taxon,
                n_replications = 100L, n_per_group = 50L,
                mean_positive_control = unname(control_mean["estimate"]),
                mean_positive_control_mcse = unname(control_mean["mcse"]),
                mean_positive_control_ci95_lower = unname(control_mean["ci95_lower"]),
                mean_positive_control_ci95_upper = unname(control_mean["ci95_upper"]),
                min_positive_control = min(control),
                q05_positive_control = as.numeric(stats::quantile(control, 0.05, names = FALSE)),
                q25_positive_control = as.numeric(stats::quantile(control, 0.25, names = FALSE)),
                median_positive_control = stats::median(control),
                q75_positive_control = as.numeric(stats::quantile(control, 0.75, names = FALSE)),
                q95_positive_control = as.numeric(stats::quantile(control, 0.95, names = FALSE)),
                max_positive_control = max(control),
                mean_positive_case = unname(case_mean["estimate"]),
                mean_positive_case_mcse = unname(case_mean["mcse"]),
                mean_positive_case_ci95_lower = unname(case_mean["ci95_lower"]),
                mean_positive_case_ci95_upper = unname(case_mean["ci95_upper"]),
                min_positive_case = min(case),
                q05_positive_case = as.numeric(stats::quantile(case, 0.05, names = FALSE)),
                q25_positive_case = as.numeric(stats::quantile(case, 0.25, names = FALSE)),
                median_positive_case = stats::median(case),
                q75_positive_case = as.numeric(stats::quantile(case, 0.75, names = FALSE)),
                q95_positive_case = as.numeric(stats::quantile(case, 0.95, names = FALSE)),
                max_positive_case = max(case),
                probability_control_below_10 = mean(control < 10),
                probability_case_below_10 = mean(case < 10),
                probability_either_group_below_10 = mean(control < 10 | case < 10),
                probability_control_below_20 = mean(control < 20),
                probability_case_below_20 = mean(case < 20),
                probability_either_group_below_20 = mean(control < 20 | case < 20),
                stringsAsFactors = FALSE
            )
        }
    }
    positive_support <- do.call(rbind, support_rows)

    comparison_pairs <- do.call(rbind, lapply(c("prevalence", "abundance"), function(domain) {
        families <- EXPECTED_FAMILIES[EXPECTED_FAMILIES$domain == domain, , drop = FALSE]
        pair_indices <- utils::combn(seq_len(nrow(families)), 2L)
        do.call(rbind, lapply(seq_len(ncol(pair_indices)), function(index) {
            data.frame(
                domain = domain,
                first_method = families$method[pair_indices[1L, index]],
                first_component = families$component[pair_indices[1L, index]],
                second_method = families$method[pair_indices[2L, index]],
                second_component = families$component[pair_indices[2L, index]],
                stringsAsFactors = FALSE
            )
        }))
    }))
    paired_metrics <- c(
        "availability", "raw_rejection", "available_conditional_raw",
        "bh_fwer_event", "by_fwer_event"
    )
    paired_rows <- vector(
        "list", length(EXPECTED_SETTINGS) * nrow(comparison_pairs) * length(paired_metrics)
    )
    paired_index <- 0L
    for (setting_id in EXPECTED_SETTINGS) {
        setting <- settings[settings$setting_id == setting_id, , drop = FALSE]
        for (pair_index in seq_len(nrow(comparison_pairs))) {
            pair <- comparison_pairs[pair_index, , drop = FALSE]
            first <- family_metrics[
                family_metrics$setting_id == setting_id &
                    family_metrics$method == pair$first_method &
                    family_metrics$component == pair$first_component, , drop = FALSE
            ]
            second <- family_metrics[
                family_metrics$setting_id == setting_id &
                    family_metrics$method == pair$second_method &
                    family_metrics$component == pair$second_component, , drop = FALSE
            ]
            first <- first[order(first$replication), , drop = FALSE]
            second <- second[order(second$replication), , drop = FALSE]
            assert(identical(as.integer(first$replication), EXPECTED_REPLICATIONS) &&
                       identical(as.integer(second$replication), EXPECTED_REPLICATIONS),
                   "A paired comparison does not contain replications 1:100.")
            for (metric in paired_metrics) {
                difference <- first[[metric]] - second[[metric]]
                summary <- cluster_t(difference, c(-1, 1))
                paired_index <- paired_index + 1L
                paired_rows[[paired_index]] <- data.frame(
                    setting_id = setting_id,
                    depth_design = as.character(setting$depth_design),
                    design_id = as.character(setting$design_id),
                    domain = pair$domain,
                    metric = metric,
                    first_method = display_method(pair$first_method),
                    first_component = pair$first_component,
                    second_method = display_method(pair$second_method),
                    second_component = pair$second_component,
                    n_paired_replications = unname(summary["n"]),
                    first_estimate = mean(first[[metric]]),
                    second_estimate = mean(second[[metric]]),
                    difference_first_minus_second = unname(summary["estimate"]),
                    difference_mcse = unname(summary["mcse"]),
                    difference_ci95_lower = unname(summary["ci95_lower"]),
                    difference_ci95_upper = unname(summary["ci95_upper"]),
                    interpretation = paste(
                        "paired null-calibration contrast only;",
                        "scientific component estimands remain method-specific"
                    ),
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    paired_differences <- do.call(rbind, paired_rows)

    rownames(null_calibration) <- rownames(pvalue_tail) <- NULL
    rownames(positive_support) <- rownames(paired_differences) <- NULL
    assert(nrow(null_calibration) == 18L, "null_calibration.csv would not have 18 rows.")
    assert(nrow(pvalue_tail) == 108L, "pvalue_tail.csv would not have 108 rows.")
    assert(nrow(positive_support) == 150L,
           "positive_support.csv would not have 150 rows.")
    assert(nrow(paired_differences) == 90L,
           "paired_differences.csv would not have 90 rows.")

    output_paths <- c(
        null_calibration = file.path(output_dir, "null_calibration.csv"),
        pvalue_tail = file.path(output_dir, "pvalue_tail.csv"),
        positive_support = file.path(output_dir, "positive_support.csv"),
        paired_differences = file.path(output_dir, "paired_differences.csv")
    )
    atomic_write_csv(null_calibration, output_paths[["null_calibration"]])
    atomic_write_csv(pvalue_tail, output_paths[["pvalue_tail"]])
    atomic_write_csv(positive_support, output_paths[["positive_support"]])
    atomic_write_csv(paired_differences, output_paths[["paired_differences"]])
    expected_output_rows <- c(18L, 108L, 150L, 90L)
    observed_output_rows <- vapply(output_paths, function(path) {
        nrow(utils::read.csv(path, check.names = FALSE))
    }, integer(1))
    assert(identical(unname(observed_output_rows), expected_output_rows),
           "A written summary failed CSV read-back validation.")

    audit_checks <- data.frame(
        check = c(
            "raw cohort", "design signature", "package contract", "null settings",
            "sample size", "fixed panel", "method-component family", "global-null truth",
            "decision replay", "generator invariants", "positive support",
            "Monte Carlo unit", "multiplicity families", "missing-result policy",
            "paired comparison", "output readback"
        ),
        status = "PASS",
        detail = c(
            "exactly replicate_0001.rds through replicate_0100.rds",
            contract$analysis_signature,
            "frozen design records DASRA 0.4.1",
            "balanced standard, balanced low depth, and fourfold no-direct-spike settings",
            "exactly 50 control and 50 case samples in every frozen dataset",
            "the same 50 prespecified taxa in every method-component family; no support filtering",
            "DASRA, ZINQ, and MaAsLin3; two separately reported components each",
            "all structural, observed-prevalence, and abundance truth flags are false",
            sprintf("raw, common BH, and common BY decisions replayed for %d families",
                    decision_replay_count),
            sprintf("probability and inactive-taxon invariants checked for %d datasets",
                    generator_check_count),
            "100 x 3 x 50 taxon-level support records; no taxa selected by realized support",
            "replication; cluster-t intervals use 100 independent replication summaries",
            "BH and BY are recalculated within each fixed 50-taxon method-component family",
            "unavailable or invalid p-values are non-rejections and p=1 in operational tail summaries",
            "method differences are paired within replication and restricted to the same broad domain",
            paste(sprintf("%s=%d rows", names(observed_output_rows), observed_output_rows),
                  collapse = "; ")
        ),
        stringsAsFactors = FALSE
    )
    audit_path <- file.path(output_dir, "audit_checks.csv")
    atomic_write_csv(audit_checks, audit_path)
    assert(nrow(utils::read.csv(audit_path, check.names = FALSE)) == nrow(audit_checks),
           "audit_checks.csv failed read-back validation.")

    cat(sprintf(
        "n=50 null summary complete: 100 replications, 3 settings, 6 families; output %s\n",
        output_dir
    ))
    invisible(list(
        null_calibration = null_calibration,
        pvalue_tail = pvalue_tail,
        positive_support = positive_support,
        paired_differences = paired_differences,
        audit_checks = audit_checks
    ))
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
    stop(
        "Usage: Rscript summarize_results.R <n50_module_root> <output_directory>",
        call. = FALSE
    )
}
summarize_main(arguments[1L], arguments[2L])
