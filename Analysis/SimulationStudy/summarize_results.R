#!/usr/bin/env Rscript

# Recompute the internal simulation tables from complete saved test families.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
study_dir <- dirname(normalizePath(sub("^--file=", "", argument[[1L]])))
data_dir <- Sys.getenv("SACAT_SIMULATION_DATA", file.path(study_dir, "data"))
table_dir <- Sys.getenv("SACAT_SIMULATION_TABLES", file.path(study_dir, "tables"))
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
base_design <- fread(file.path(data_dir, "base_design.csv"))
joint_design <- fread(file.path(data_dir, "joint_design.csv"))
added_design <- fread(file.path(data_dir, "abundance_design.csv"))
abundance_design <- rbindlist(list(
    base_design[n_per_group == 80L & scenario == "present_conditional_abundance"],
    added_design
), use.names = TRUE, fill = TRUE)
family_keys <- c("setting_id", "replication", "method", "component")
setting_keys <- c("setting_id", "method", "component")
setting_fields <- c("setting_id", "scenario", "n_per_group", "covariate",
                    "effect_parameter", "effect_index", "structural_difference",
                    "abundance_log_effect")

summarize_file <- function(filename, design, run, expected_rows) {
    data <- fread(file.path(data_dir, filename))
    stopifnot(nrow(data) == expected_rows,
              !anyDuplicated(data[, c(family_keys, "taxon"), with = FALSE]),
              !anyNA(data$available), !anyNA(data$truth_for_component),
              !anyNA(data$scenario_signal))
    design <- copy(design)
    for (field in setdiff(setting_fields, names(design))) design[, (field) := NA_real_]
    stopifnot(!anyDuplicated(design$setting_id),
              all(data$setting_id %in% design$setting_id))
    data <- merge(data, design[, ..setting_fields], by = "setting_id", sort = FALSE)
    if (!"primary_comparison" %in% names(data)) data[, primary_comparison := TRUE]
    families <- data[, .(taxa = .N), by = family_keys]
    stopifnot(all(families$taxa == 50L),
              all(families[, .N, by = setting_keys]$N == 100L),
              setequal(data$replication, 1:100))
    data[, `:=`(
        is_null = if (run == "joint") effect_index == 0L else effect_parameter == 0,
        decision_p = fifelse(available & is.finite(p_value) &
                                p_value >= 0 & p_value <= 1, p_value, 1)
    )]
    stopifnot(!any(data[is_null == TRUE]$truth_for_component))
    strata <- c("scenario", "n_per_group", "covariate", "method", "component")
    data[, calibrated_p := {
        reference <- sort(decision_p[is_null])
        null_replications <- replication[is_null]
        stopifnot(length(reference) == 5000L,
                  all(tabulate(null_replications, nbins = 100L) == 50L))
        held_out <- lapply(split(decision_p[is_null], null_replications), sort)
        by_replication <- split(seq_len(.N), replication)
        values <- numeric(.N)
        for (replication_id in names(by_replication)) {
            rows <- by_replication[[replication_id]]
            counts <- findInterval(decision_p[rows], reference) -
                findInterval(decision_p[rows], held_out[[replication_id]])
            values[rows] <- (1 + counts) / 4951
        }
        values
    }, by = strata]
    data[, raw_bh_q := p.adjust(decision_p, "BH"), by = family_keys]
    data[, calibrated_bh_q := p.adjust(calibrated_p, "BH"), by = family_keys]
    data[, `:=`(raw_reject = decision_p < 0.05,
                raw_bh_reject = raw_bh_q < 0.05,
                calibrated_raw_reject = calibrated_p < 0.05,
                calibrated_bh_reject = calibrated_bh_q < 0.05)]
    stopifnot(!any(data[available == FALSE]$raw_bh_reject),
              !any(data[available == FALSE]$calibrated_bh_reject))
    replications <- data[, .(
        taxa = .N, signals = sum(truth_for_component),
        designated_taxa = sum(scenario_signal), availability = mean(available),
        raw_type1 = mean(raw_reject[!truth_for_component]),
        raw_bh_power = mean(raw_bh_reject[truth_for_component]),
        raw_bh_fdr = sum(raw_bh_reject & !truth_for_component) / max(sum(raw_bh_reject), 1L),
        raw_bh_type1 = mean(raw_bh_reject[!truth_for_component]),
        raw_bh_designated_rejection = mean(raw_bh_reject[scenario_signal]),
        calibrated_raw_type1 = mean(calibrated_raw_reject[!truth_for_component]),
        calibrated_bh_power = mean(calibrated_bh_reject[truth_for_component]),
        calibrated_bh_fdr = sum(calibrated_bh_reject & !truth_for_component) /
            max(sum(calibrated_bh_reject), 1L),
        calibrated_bh_type1 = mean(calibrated_bh_reject[!truth_for_component]),
        calibrated_bh_designated_rejection = mean(calibrated_bh_reject[scenario_signal])
    ), by = family_keys]
    measures <- setdiff(names(replications), c(family_keys, "taxa", "signals", "designated_taxa"))
    summary <- replications[, c(list(replications = .N, taxa = unique(taxa),
                                    signals = unique(signals), designated_taxa = unique(designated_taxa)),
                               lapply(.SD, mean)), by = setting_keys, .SDcols = measures]
    mcse <- replications[, lapply(.SD, function(x) sd(x) / sqrt(.N)),
                         by = setting_keys, .SDcols = measures]
    setnames(mcse, measures, paste0(measures, "_mcse"))
    summary <- merge(summary, mcse, by = setting_keys, sort = FALSE)
    summary <- merge(summary, design[, ..setting_fields], by = "setting_id", sort = FALSE)
    primary <- unique(data[, c(setting_keys, "primary_comparison"), with = FALSE])
    stopifnot(!anyDuplicated(primary[, ..setting_keys]))
    summary <- merge(summary, primary, by = setting_keys, sort = FALSE)
    summary[, run := run]
    replications[, run := run]
    availability <- data[, .(tests = .N, available_tests = sum(available)),
                         by = .(method, component, status)][, run := run]
    list(summary = summary, replications = replications, availability = availability)
}

inputs <- list(
    list("base_n060.csv.gz", base_design, "base", 1700000L),
    list("base_n080.csv.gz", base_design, "base", 1700000L),
    list("base_n120.csv.gz", base_design, "base", 1700000L),
    list("joint.csv.gz", joint_design, "joint", 3150000L),
    list("abundance.csv.gz", abundance_design, "abundance", 160000L)
)
results <- lapply(inputs, function(input) {
    message("Summarizing ", input[[1L]])
    do.call(summarize_file, input)
})
summary <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
replications <- rbindlist(lapply(results, `[[`, "replications"), fill = TRUE)
availability <- rbindlist(lapply(results, `[[`, "availability"))[
    , .(tests = sum(tests), available_tests = sum(available_tests)),
    by = .(run, method, component, status)]
stopifnot(nrow(summary) == 1682L, nrow(replications) == 168200L,
          !anyDuplicated(summary[, c("run", setting_keys), with = FALSE]))
setorderv(summary, c("run", "setting_id", "method", "component"))
fwrite(summary, file.path(table_dir, "setting_summary.csv"), na = "")
fwrite(replications, file.path(table_dir, "replication_metrics.csv.gz"), na = "")
fwrite(availability, file.path(table_dir, "test_availability.csv"), na = "")
null_summary <- summary[(run == "joint" & effect_index == 0L) |
    (run != "joint" & effect_parameter == 0), .(
        run, setting_id, scenario, n_per_group, covariate, method, component,
        replications, type1_error = raw_type1, monte_carlo_se = raw_type1_mcse
    )]
fwrite(null_summary, file.path(table_dir, "type1_error.csv"), na = "")

structural <- summary[run == "base" & n_per_group == 80L &
    scenario %in% c("observed_prevalence_only", "structural_matched_prevalence") &
    ((method == "SACAT" & component == "structural_absence") |
     (method %in% c("ZINQ", "MaAsLin 3") & component == "observed_prevalence"))]
structural[, `:=`(
    plotted_value = calibrated_bh_designated_rejection,
    original_bh_value = raw_bh_designated_rejection,
    plotted_metric = "Designated-taxon rejection probability (calibrated p + BH)",
    effect_index = NA_integer_
)]
base_abundance <- summary[run == "base" & n_per_group == 80L &
    scenario == "present_conditional_abundance" & effect_parameter > 0 & effect_parameter < 1.4 &
    ((method == "SACAT" & component == "present_conditional_abundance") |
     (method == "MaAsLin 3" & component == "detected_log_abundance"))]
base_abundance[, effect_index := NA_integer_]
supplement <- summary[run == "abundance" & effect_parameter > 0 & effect_parameter < 1.4]
joint <- summary[run == "joint" & n_per_group == 80L & effect_index > 0L & primary_comparison]
power <- rbindlist(list(base_abundance, supplement, joint), fill = TRUE)
power[, `:=`(plotted_value = calibrated_bh_power, original_bh_value = raw_bh_power,
             plotted_metric = "Power (calibrated p + BH)")]
main <- rbindlist(list(structural, power), fill = TRUE)[, .(
    panel = "A", simulation = "Internal validation", setting_id, scenario,
    n_per_group, covariate, effect_parameter, effect_index,
    structural_difference, abundance_log_effect, method, component,
    plotted_metric, plotted_value, original_bh_value, include_covariate = NA
)]
stopifnot(nrow(main) == 204L, !anyNA(main$plotted_value),
          !anyDuplicated(main[, ..setting_keys]),
          all(main$plotted_value >= 0 & main$plotted_value <= 1))
fwrite(main, file.path(table_dir, "main_plot_data.csv"), na = "")

manuscript <- rbindlist(list(
    summary[run == "base" & method == "SACAT" & component == "structural_absence" &
        scenario %in% c("structural_only", "structural_matched_prevalence") & effect_parameter > 0,
        .(comparison = "SACAT structural FDR", settings = .N,
          raw_bh_fdr = mean(raw_bh_fdr), calibrated_bh_fdr = mean(calibrated_bh_fdr))],
    summary[run == "base" & method == "SACAT" & component == "present_conditional_abundance" &
        scenario == "present_conditional_abundance" & effect_parameter > 0,
        .(comparison = "SACAT abundance FDR", settings = .N,
          raw_bh_fdr = mean(raw_bh_fdr), calibrated_bh_fdr = mean(calibrated_bh_fdr))]
))
fwrite(manuscript, file.path(table_dir, "manuscript_summary.csv"))
message("Created complete internal simulation tables and 204 manuscript points in ", table_dir)
