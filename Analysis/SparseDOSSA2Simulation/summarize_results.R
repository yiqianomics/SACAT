#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

# Resolve the analysis module relative to this script so the file remains portable.
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("Cannot determine this script's location.", call. = FALSE)
script_file <- normalizePath(sub("^--file=", "", script_argument[1L]), mustWork = TRUE)
script_dir <- dirname(script_file)
arguments <- commandArgs(trailingOnly = TRUE)
module_relative <- if (length(arguments)) arguments[1L] else "."
module_root <- normalizePath(file.path(script_dir, module_relative), mustWork = TRUE)
simulation_file <- file.path(module_root, "simulation.R")
if (!file.exists(simulation_file)) stop("simulation.R is missing from the module root.", call. = FALSE)

assert <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
    invisible(TRUE)
}

# Source the runner; its sys.nframe() guard leaves command-line main() dormant.
load_runner <- function(path) {
    environment <- new.env(parent = globalenv())
    source(path, local = environment, chdir = FALSE)
    assert(identical(environment$RUNNER_FILE, normalizePath(path)),
           "The sourced runner did not resolve its own path.")
    environment
}

# Summarize one scalar operating characteristic across independent replications.
cluster_t <- function(data, groups, bounds = c(0, 1)) {
    output <- data[is.finite(value), {
        n <- .N
        estimate <- mean(value)
        mcse <- stats::sd(value) / sqrt(n)
        half <- stats::qt(0.975, df = n - 1L) * mcse
        list(
            n_replications = n,
            estimate = estimate,
            mcse = mcse,
            ci95_lower = max(bounds[1L], estimate - half),
            ci95_upper = min(bounds[2L], estimate + half)
        )
    }, by = groups]
    assert(nrow(output) > 0L && all(output$n_replications == 100L),
           "A replication-level summary does not contain exactly 100 values.")
    output[]
}

# Give a summarized metric a stable column prefix before joining tables.
metric_table <- function(data, metric_name, multiplicity_name, groups, prefix,
                         bounds = c(0, 1)) {
    value <- cluster_t(
        data[metric == metric_name & multiplicity == multiplicity_name],
        groups, bounds
    )
    data.table::setnames(
        value,
        c("n_replications", "estimate", "mcse", "ci95_lower", "ci95_upper"),
        paste0(prefix, c("_n", "", "_mcse", "_ci95_lower", "_ci95_upper"))
    )
    value
}

# Join several metric summaries by their shared method-family identifiers.
join_metrics <- function(parts, groups) {
    Reduce(function(left, right) merge(left, right, by = groups, all = TRUE, sort = FALSE),
           parts)
}

# Produce availability, raw power, common-BH power, and common-BH FDR columns.
operating_table <- function(data, groups) {
    join_metrics(list(
        metric_table(data, "available_fraction", "raw", groups, "availability"),
        metric_table(data, "power", "raw", groups, "raw_power"),
        metric_table(data, "power", "BH", groups, "bh_power"),
        metric_table(data, "fdp", "BH", groups, "bh_fdr")
    ), groups)
}

# Wilson intervals avoid degenerate normal intervals for all-null BH events.
wilson_summary <- function(data, groups) {
    output <- data[is.finite(value), {
        assert(all(value %in% c(0, 1)), "All-null BH FDP is not binary.")
        n <- .N
        events <- sum(value)
        estimate <- events / n
        z <- stats::qnorm(0.975)
        denominator <- 1 + z^2 / n
        center <- (estimate + z^2 / (2 * n)) / denominator
        half <- z * sqrt(estimate * (1 - estimate) / n + z^2 / (4 * n^2)) /
            denominator
        list(
            bh_fdr_fwer_n = n,
            bh_fdr_fwer_events = events,
            bh_fdr_fwer = estimate,
            bh_fdr_fwer_mcse = stats::sd(value) / sqrt(n),
            bh_fdr_fwer_ci95_lower = max(0, center - half),
            bh_fdr_fwer_ci95_upper = min(1, center + half)
        )
    }, by = groups]
    assert(nrow(output) > 0L && all(output$bh_fdr_fwer_n == 100L),
           "An all-null BH family does not contain exactly 100 replications.")
    output[]
}

engine <- load_runner(simulation_file)
raw_dir <- file.path(module_root, "results", "raw")
expected_names <- sprintf("replicate_%04d.rds", 1:100)
observed_names <- sort(list.files(raw_dir, pattern = "^replicate_[0-9]{4}\\.rds$"))
assert(identical(observed_names, expected_names),
       "The frozen analysis requires exactly replicate_0001.rds through replicate_0100.rds.")

# The frozen runner performs full design, raw-record, truth, and decision validation.
completion <- engine$summarize_results(module_root, 1:100)
summary_file <- file.path(module_root, "summary", "summary.rds")
summary_object <- readRDS(summary_file)
contract <- engine$load_contract(module_root, verify_environment = FALSE)
assert(identical(completion$completed_replications, 100L) &&
           identical(summary_object$completion$completed_replications, 100L),
       "The formal summary did not complete all 100 replications.")
assert(identical(summary_object$completion$analysis_signature,
                 contract$analysis_signature),
       "The summary signature differs from the frozen design.")

metrics <- data.table::as.data.table(summary_object$replication_metrics)
assert(setequal(unique(metrics$replication), 1:100) && all(metrics$n_taxa == 50L),
       "Replication metrics do not preserve the 100 x fixed-50-family contract.")
groups <- c("setting_id", "scenario", "depth_design", "method", "component")

# Nominal calibration and depth-imbalance robustness remain visibly distinct.
null_data <- metrics[setting_id %in% c("NULL_BALANCED", "NULL_FOURFOLD")]
null_calibration <- join_metrics(list(
    metric_table(null_data, "available_fraction", "raw", groups, "availability"),
    metric_table(null_data, "null_rejection_rate", "raw", groups, "raw_rejection"),
    wilson_summary(null_data[metric == "fdp" & multiplicity == "BH"], groups)
), groups)
null_calibration[, calibration_role := data.table::fifelse(
    setting_id == "NULL_BALANCED", "nominal global null",
    "depth-imbalance no-direct-effect robustness"
)]

# Observed-prevalence methods share an estimand; DASRA structural absence is separate.
prevalence_data <- metrics[
    setting_id %in% c("STRUCT_BALANCED", "STRUCT_FOURFOLD",
                      "JOINT_BALANCED", "JOINT_FOURFOLD") &
        ((method %in% c("ZINQ", "MaAsLin 3") & component == "observed_prevalence") |
         (method == "DASRA" & component == "structural_absence"))
]
prevalence_benchmark <- operating_table(prevalence_data, groups)
prevalence_benchmark[, `:=`(
    benchmark_role = data.table::fifelse(
        scenario == "structural_only", "isolated direct prevalence spike",
        "joint prevalence and abundance spike"
    ),
    estimand = data.table::fifelse(
        component == "observed_prevalence", "observed prevalence",
        "DASRA structural absence; reported separately"
    )
)]

# Abundance results are operational summaries for method-specific components.
abundance_components <- c(
    "present_conditional_abundance", "detected_quantile_abundance",
    "detected_log_abundance", "general_abundance"
)
abundance_data <- metrics[
    setting_id %in% c("ABUND_BALANCED", "ABUND_FOURFOLD",
                      "JOINT_BALANCED", "JOINT_FOURFOLD") &
        component %in% abundance_components
]
abundance_benchmark <- operating_table(abundance_data, groups)
abundance_benchmark[, interpretation :=
    "method-component-specific sensitivity to the native abundance perturbation"]

# Dense same-direction settings are stress tests, with 60% violating null majority.
dense_data <- metrics[
    setting_id %in% c("ABUND40_SAME_BALANCED", "ABUND60_SAME_BALANCED") &
        component %in% abundance_components
]
dense_stress <- operating_table(dense_data, groups)
dense_stress[, stress_role := data.table::fifelse(
    setting_id == "ABUND40_SAME_BALANCED", "40% same-direction abundance signals",
    "60% same-direction signals; violates DASRA strict null-majority condition"
)]
dense_stress[, interpretation :=
    "method-component-specific stress behavior; not a nominal operating point"]

audit_checks <- data.table::data.table(
    check = c(
        "frozen raw cohort", "runner validation", "frozen design signature",
        "fixed method-component family", "Monte Carlo unit", "common BH",
        "estimand separation"
    ),
    status = "PASS",
    detail = c(
        "exactly replications 1 through 100",
        "100 raw records and 1,000 frozen datasets validated by simulation.R",
        contract$analysis_signature,
        "50 taxa for every setting x method x component family",
        "replication; cluster t intervals for rates, power, and FDP",
        "fixed 50-taxon family; all-null FDR/FWER uses Wilson interval",
        "observed prevalence, DASRA structural absence, and abundance components separated"
    )
)

summary_dir <- file.path(module_root, "summary")
engine$atomic_write_csv(as.data.frame(audit_checks),
                        file.path(summary_dir, "audit_checks.csv"))
engine$atomic_write_csv(as.data.frame(null_calibration),
                        file.path(summary_dir, "null_calibration.csv"))
engine$atomic_write_csv(as.data.frame(prevalence_benchmark),
                        file.path(summary_dir, "prevalence_benchmark.csv"))
engine$atomic_write_csv(as.data.frame(abundance_benchmark),
                        file.path(summary_dir, "abundance_benchmark.csv"))
engine$atomic_write_csv(as.data.frame(dense_stress),
                        file.path(summary_dir, "dense_stress.csv"))

cat(sprintf("Reviewer summary complete: 100 replications; signature %s\n",
            contract$analysis_signature))
