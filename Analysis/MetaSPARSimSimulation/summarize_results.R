#!/usr/bin/env Rscript

# Recompute simulation summaries from complete taxon-level p-value families.
# Requires data.table. Inputs and tables are relative to this script.
# Run with Rscript; no simulation or model fitting is performed.

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install the data.table R package.", call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript.", call. = FALSE)
root <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
input <- file.path(root, "data")
output <- file.path(root, "tables")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
design <- fread(file.path(input, "design.csv"))
tests <- fread(file.path(input, "tests.csv"))
family_keys <- c("setting_index", "replication", "test_id")
setting_keys <- c("setting_index", "test_id")
alpha <- .05

calibrate_p <- function(p, replication, is_null, family_size) {
    reference <- p[is_null]
    reference_replication <- replication[is_null]
    reference_counts <- table(reference_replication)
    stopifnot(length(reference_counts) == 500L, all(reference_counts == family_size))
    sorted_reference <- sort(reference)
    held_out <- lapply(split(reference, reference_replication), sort)
    rows_by_replication <- split(seq_along(p), replication)
    calibrated <- numeric(length(p))
    for (replication_id in names(rows_by_replication)) {
        rows <- rows_by_replication[[replication_id]]
        count <- findInterval(p[rows], sorted_reference) -
            findInterval(p[rows], held_out[[replication_id]])
        calibrated[rows] <- (1 + count) / (length(reference) - family_size + 1)
    }
    calibrated
}

input_files <- c("raw_sacat.csv.gz", "raw_components.csv.gz",
                 "raw_abundance.csv.gz", "raw_native_joint.csv.gz")
replication_tables <- vector("list", length(input_files))
calibration_tables <- vector("list", length(input_files))
for (i in seq_along(input_files)) {
    message("Summarizing ", input_files[i])
    data <- fread(file.path(input, input_files[i]), colClasses = c(p_value = "character"))
    data[, p_value := as.numeric(p_value)]
    required <- c(family_keys, "taxon_index", "p_value", "available",
                  "truth_for_component", "truth_any")
    stopifnot(all(required %in% names(data)), !anyNA(data[, ..family_keys]),
              !anyDuplicated(data[, c(family_keys, "taxon_index"), with = FALSE]),
              all(data$test_id %in% tests$test_id),
              all(data$available %in% 0:1), all(data$truth_for_component %in% 0:1),
              all(data$truth_any %in% 0:1), all(data$taxon_index %in% seq_len(100L)))
    data <- merge(data, design[, .(setting_index, scenario, include_covariate,
                                  n_replications, n_taxa)], by = "setting_index", sort = FALSE)
    stopifnot(!anyNA(data$scenario), all(data$n_taxa == 100L))
    coverage <- data[, .(rows = .N, replications = uniqueN(replication),
                         expected = unique(n_replications)), by = setting_keys]
    expected_settings <- if (input_files[i] == "raw_native_joint.csv.gz") {
        design[scenario %in% c("global_null", "joint"), setting_index]
    } else design$setting_index
    stopifnot(setequal(data$setting_index, expected_settings),
              nrow(coverage) == length(expected_settings) * uniqueN(data$test_id),
              all(coverage$rows == coverage$expected * 100L),
              all(coverage$replications == coverage$expected),
              all(data[, .N, by = family_keys]$N == 100L),
              all(data$replication >= 1 & data$replication <= data$n_replications))
    data[, `:=`(available = as.logical(available),
                 truth_for_component = as.logical(truth_for_component),
                 truth_any = as.logical(truth_any), is_null = scenario == "global_null")]
    stopifnot(!any(data[is_null == TRUE]$truth_for_component))
    data[, decision_p := fifelse(available & is.finite(p_value) & p_value >= 0 &
                                 p_value <= 1, p_value, 1)]
    data[, calibrated_p := calibrate_p(decision_p, replication, is_null, 100L),
         by = .(include_covariate, test_id)]
    data[, original_bh_q := p.adjust(decision_p, "BH"), by = family_keys]
    data[, calibrated_bh_q := p.adjust(calibrated_p, "BH"), by = family_keys]
    data[, `:=`(raw_reject = decision_p < alpha,
                 calibrated_raw_reject = calibrated_p < alpha,
                 original_bh_reject = original_bh_q < alpha,
                 calibrated_bh_reject = calibrated_bh_q < alpha)]
    stopifnot(all(is.finite(data$calibrated_p)),
              all(data$calibrated_p >= 1 / 49901 & data$calibrated_p <= 1),
              !any(data[available == FALSE]$original_bh_reject),
              !any(data[available == FALSE]$calibrated_bh_reject))
    replication_tables[[i]] <- data[, .(
        n_taxa = .N, n_signals = sum(truth_for_component), n_designated = sum(truth_any),
        availability = mean(available),
        signal_availability = if (any(truth_for_component)) mean(available[truth_for_component]) else NA_real_,
        null_availability = mean(available[!truth_for_component]),
        raw_type1 = mean(raw_reject[!truth_for_component]),
        calibrated_raw_type1 = mean(calibrated_raw_reject[!truth_for_component]),
        original_bh_power = if (any(truth_for_component)) mean(original_bh_reject[truth_for_component]) else NA_real_,
        calibrated_bh_power = if (any(truth_for_component)) mean(calibrated_bh_reject[truth_for_component]) else NA_real_,
        original_bh_fdr = sum(original_bh_reject & !truth_for_component) / max(sum(original_bh_reject), 1L),
        calibrated_bh_fdr = sum(calibrated_bh_reject & !truth_for_component) / max(sum(calibrated_bh_reject), 1L),
        original_bh_type1 = mean(original_bh_reject[!truth_for_component]),
        calibrated_bh_type1 = mean(calibrated_bh_reject[!truth_for_component]),
        original_bh_discoveries = sum(original_bh_reject),
        calibrated_bh_discoveries = sum(calibrated_bh_reject),
        original_bh_designated_rejection = if (any(truth_any)) mean(original_bh_reject[truth_any]) else NA_real_,
        calibrated_bh_designated_rejection = if (any(truth_any)) mean(calibrated_bh_reject[truth_any]) else NA_real_
    ), by = family_keys]
    calibration_tables[[i]] <- rbindlist(lapply(c(.001, .005, .01, .05), function(level) {
        data[is_null == TRUE, .(level = level, null_replications = uniqueN(replication),
            training_tests = 49900L, original_raw_type1 = mean(decision_p < level),
            calibrated_raw_type1 = mean(calibrated_p < level)),
            by = .(include_covariate, test_id)]
    }))
    rm(data)
    invisible(gc())
}

replications <- rbindlist(replication_tables)
metric_columns <- setdiff(names(replications), c(family_keys, "n_taxa", "n_signals", "n_designated"))
summary <- replications[, c(list(n_replications = .N, n_taxa = unique(n_taxa),
                                n_signals = unique(n_signals), n_designated = unique(n_designated)),
                           lapply(.SD, mean)), by = setting_keys, .SDcols = metric_columns]
summary <- merge(summary, tests, by = "test_id", sort = FALSE)
summary <- merge(summary, design[, .(setting_index, setting_id, scenario, include_covariate,
                                     effect_index, structural_difference, abundance_log_effect)],
                 by = "setting_index", sort = FALSE)
setorderv(summary, c("setting_index", "test_id"))
identity_columns <- c("setting_id", "scenario", "include_covariate", "effect_index",
                      "structural_difference", "abundance_log_effect", "method", "component",
                      "n_replications", "n_taxa", "n_signals")
setcolorder(summary, c(identity_columns, setdiff(names(summary), identity_columns)))
fwrite(summary, file.path(output, "setting_summary.csv"))
fwrite(replications, file.path(output, "replication_metrics.csv.gz"))
calibration <- merge(rbindlist(calibration_tables), tests, by = "test_id", sort = FALSE)
fwrite(calibration, file.path(output, "null_calibration.csv"))
type1_mcse <- replications[, .(raw_type1_mcse = sd(raw_type1) / sqrt(.N)), by = setting_keys]
type1 <- merge(summary[scenario == "global_null"], type1_mcse, by = setting_keys, sort = FALSE)
fwrite(type1[, c(identity_columns, "raw_type1", "raw_type1_mcse"), with = FALSE],
       file.path(output, "type1_error.csv"))
fwrite(summary[, c(identity_columns, "original_bh_power", "calibrated_bh_power"), with = FALSE],
       file.path(output, "power.csv"))
fwrite(summary[, c(identity_columns, "original_bh_fdr", "calibrated_bh_fdr"), with = FALSE],
       file.path(output, "fdr.csv"))
fwrite(summary[, c(identity_columns, "availability", "signal_availability", "null_availability"),
               with = FALSE], file.path(output, "availability.csv"))
for (name in c("design.csv", "software_versions.csv")) {
    file.copy(file.path(input, name), file.path(output, name), overwrite = TRUE)
}

core_methods <- c("SACAT", "ZINQ", "MaAsLin 3")
main <- summary[
    (scenario == "structural_only" & method %in% core_methods &
        component %in% c("structural_absence", "observed_prevalence")) |
    (scenario == "abundance_only" & method %in% core_methods &
        component %in% c("present_conditional_abundance", "detected_quantile_abundance", "detected_log_abundance")) |
    (scenario == "joint" & component %in% c("primary_omnibus", "native_joint", "general_abundance"))]
main <- main[, .(panel = "B", simulation = "metaSPARSim", setting_id, scenario,
    n_per_group = NA_integer_, covariate = NA_character_, effect_parameter = NA_real_,
    effect_index, structural_difference, abundance_log_effect, method, component,
    plotted_metric = "Power (calibrated p + BH)", plotted_value = calibrated_bh_power,
    original_bh_value = original_bh_power, include_covariate)]
stopifnot(nrow(main) == 180L, all(is.finite(main$plotted_value)),
          !anyDuplicated(main[, .(setting_id, method, component)]))
fwrite(main, file.path(output, "main_plot_data.csv"))
message("Wrote ", nrow(summary), " setting summaries and ", nrow(main), " main-figure points to ", output)
