#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
study_dir <- dirname(normalizePath(script))
root <- Sys.getenv("SACAT_ABUNDANCE_ROOT", file.path(study_dir, "abundance_output"))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
data_dir <- Sys.getenv("SACAT_SIMULATION_DATA", file.path(study_dir, "data"))
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
methods <- c("SACAT", "ZINQ", "MaAsLin 3")
components <- c("present_conditional_abundance", "detected_quantile_abundance",
                "detected_log_abundance")
for (directory in c("results", "datasets", "logs", "tmp", "preflight", "summary", "zinq_original")) {
    dir.create(file.path(root, directory), showWarnings = FALSE)
}

load_original <- function() {
    simulation_variables <- grep("^SACAT_SIMULATION_", names(Sys.getenv()), value = TRUE)
    Sys.unsetenv(simulation_variables)
    env <- new.env(parent = globalenv())
    source(file.path(study_dir, "run_simulation.R"), local = env)
    env
}

source_env <- load_original()
original_settings <- as.data.table(source_env$SETTINGS)
archived_settings <- fread(file.path(data_dir, "base_design.csv"))
stopifnot(isTRUE(all.equal(original_settings, archived_settings, check.attributes = FALSE)))

new_effects <- c(.30, .50, .60)
new_settings <- rbindlist(lapply(source_env$CONFIG$covariate_levels, function(covariate_value) {
    base <- original_settings[n_per_group == 80L & covariate == covariate_value &
                               scenario == "present_conditional_abundance" & effect_parameter == .40]
    stopifnot(nrow(base) == 1L)
    rbindlist(lapply(seq_along(new_effects), function(i) {
        row <- copy(base)
        effect <- new_effects[i]
        row[, `:=`(effect_parameter = effect, effect_index = 4L + i, effect_level = 4L + i,
                   base_setting_index = max(original_settings$base_setting_index) + i,
                   base_setting_id = paste0("AB_MODEL_E", source_env$format_setting_value(effect)),
                   generic_abundance = FALSE)]
        row[, setting_id := paste(design_id, base_setting_id, sep = "__")]
        row
    }))
}))
new_settings[, setting_index := max(original_settings$setting_index) + .I]
stopifnot(nrow(new_settings) == 6L, !anyDuplicated(new_settings$setting_id),
          !any(new_settings$setting_index %in% original_settings$setting_index))
jobs <- new_settings[rep(seq_len(.N), each = 100L)]
jobs[, replication := rep(seq_len(100L), nrow(new_settings))]
jobs[, job_id := .I]

run_job <- function(job, preflight = FALSE) {
    setting <- as.data.frame(job[, names(original_settings), with = FALSE])
    setting$generic_abundance <- FALSE
    replication_id <- job$replication
    name <- sprintf("%s_rep%03d", setting$setting_id, replication_id)
    destination <- file.path(root, if (preflight) "preflight" else "results", paste0(name, ".rds"))
    if (file.exists(destination)) return(invisible(destination))
    started <- Sys.time()
    design <- source_env$make_replication_design(replication_id, 80L, setting$covariate)
    template <- source_env$make_taxon_template(setting$n_signal_target, setting$template_role)
    simulation <- source_env$simulate_setting(replication_id, setting, design, template)
    stopifnot(ncol(simulation$counts) == 160L, length(simulation$evaluation_taxa) == 50L,
              all(colSums(simulation$counts) == simulation$metadata$reads),
              !any(simulation$truth$truth_structural),
              sum(simulation$truth$truth_abundance) == 10L,
              max(abs(abs(simulation$truth$oracle_present_conditional_log_relative_effect[
                  simulation$truth$truth_abundance]) - setting$effect_parameter)) < 1e-7)
    temporary_root <- file.path(root, "tmp", name)
    dir.create(temporary_root, showWarnings = FALSE)
    result <- as.data.table(source_env$analyze_setting(simulation, replication_id, temporary_root))
    stopifnot(nrow(result) == 350L, setequal(result$method, methods),
              all(result[, .N, by = .(method, component)]$N == 50L))
    errors <- result[startsWith(status, "error:")]
    if (nrow(errors)) stop(paste(unique(errors$status), collapse = "\n"))
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    saveRDS(list(setting = setting, counts = simulation$counts, metadata = simulation$metadata,
                 truth = simulation$truth),
            file.path(root, if (preflight) "preflight" else "datasets", paste0(name, "_data.rds")),
            compress = "gzip")
    saveRDS(list(results = result, elapsed_seconds = elapsed), destination, compress = "gzip")
    message(sprintf("%s | replication %d | %.1f seconds", setting$setting_id, replication_id, elapsed))
    invisible(destination)
}

prepare <- function() {
    versions <- data.table(package = c("SACAT", "ZINQ", "maaslin3", "quantreg", "logistf",
                                      "Matrix", "statmod", "data.table"))
    versions[, version := vapply(package, function(p) as.character(packageVersion(p)), character(1))]
    stopifnot(versions[package == "SACAT", version] == "0.7.1",
              versions[package == "ZINQ", version] == "2.0",
              versions[package == "maaslin3", version] == "1.5.3")
    fwrite(versions, file.path(root, "package_versions.csv"))
    fwrite(data.table(R = R.version.string, platform = R.version$platform),
           file.path(root, "runtime.csv"))
    fwrite(new_settings, file.path(root, "simulation_design.csv"))
    fwrite(jobs[, .(job_id, setting_id, replication)], file.path(root, "jobs.csv"))
    columns <- c("setting_id", "replication", "method", "component", "taxon", "p_value",
                 "q_value_bh", "available", "truth_for_component", "scenario_signal",
                 "scenario", "effect_parameter", "n_per_group", "covariate")
    original <- fread(file.path(data_dir, "base_n080.csv.gz"))
    original <- merge(original, archived_settings[, .(setting_id, scenario, effect_parameter,
                                                     n_per_group, covariate)], by = "setting_id")
    original <- original[scenario == "present_conditional_abundance" & method %in% methods]
    original[, q_value_bh := p.adjust(fifelse(available & is.finite(p_value), p_value, 1), "BH"),
             by = .(setting_id, replication, method, component)]
    original <- original[, ..columns]
    saveRDS(original[effect_parameter == 0 & component %in% components],
            file.path(root, "preflight/null_reference.rds"), compress = "gzip")
    saveRDS(original[replication == 1L & effect_parameter == .40],
            file.path(root, "preflight/original_results.rds"), compress = "gzip")
    truth <- fread(file.path(data_dir, "base_truth.csv.gz"))[
        n_per_group == 80L & scenario == "present_conditional_abundance" &
            effect_parameter == .40 & replication == 1L]
    saveRDS(truth, file.path(root, "preflight/original_truth.rds"), compress = "gzip")
    print(new_settings[, .(setting_id, effect_parameter, covariate, n_per_group, setting_index)])
}

check_preflight <- function() {
    reference <- readRDS(file.path(root, "preflight/original_results.rds"))
    truth_reference <- readRDS(file.path(root, "preflight/original_truth.rds"))
    checks <- list()
    for (covariate_value in source_env$CONFIG$covariate_levels) {
        job <- original_settings[n_per_group == 80L & covariate == covariate_value &
                                   scenario == "present_conditional_abundance" & effect_parameter == .40]
        job[, replication := 1L]
        run_job(job, preflight = TRUE)
        name <- sprintf("%s_rep001", job$setting_id)
        result <- readRDS(file.path(root, "preflight", paste0(name, ".rds")))$results
        current_truth <- readRDS(file.path(root, "preflight", paste0(name, "_data.rds")))$truth
        previous_truth <- truth_reference[setting_id == job$setting_id]
        previous_truth <- previous_truth[match(current_truth$taxon, taxon)]
        stopifnot(isTRUE(all.equal(current_truth,
            as.data.frame(previous_truth[, names(current_truth), with = FALSE]),
            check.attributes = FALSE, tolerance = 1e-12)))
        compared <- merge(result, reference[setting_id == job$setting_id],
                          by = c("setting_id", "replication", "method", "component", "taxon"),
                          suffixes = c("_local", "_original"))
        stopifnot(nrow(compared) == 350L)
        checks[[covariate_value]] <- compared[, .(
            p_max_difference = max(abs(p_value_local - p_value_original), na.rm = TRUE),
            availability_differences = sum(available_local != available_original),
            bh_decision_differences = sum((q_value_bh_local < .05) != (q_value_bh_original < .05))
        ), by = .(method, component)][, covariate := covariate_value]
    }
    checks <- rbindlist(checks)
    fwrite(checks, file.path(root, "preflight/reproduction_checks.csv"))
    print(checks)
    comparable <- checks[!(method == "ZINQ" & component == "detected_quantile_abundance")]
    stopifnot(max(comparable$p_max_difference) < 1e-5,
              sum(checks$availability_differences) == 0L,
              sum(comparable$bh_decision_differences) == 0L)
}

run_original_zinq <- function(worker, workers) {
    settings <- original_settings[n_per_group == 80L & scenario == "present_conditional_abundance"]
    tasks <- settings[rep(seq_len(.N), each = 100L)]
    tasks[, replication := rep(seq_len(100L), nrow(settings))]
    for (i in seq.int(worker, nrow(tasks), by = workers)) {
        job <- tasks[i]
        setting <- as.data.frame(job[, names(original_settings), with = FALSE])
        replication_id <- job$replication
        name <- sprintf("%s_rep%03d", setting$setting_id, replication_id)
        destination <- file.path(root, "zinq_original", paste0(name, ".rds"))
        if (file.exists(destination)) next
        started <- Sys.time()
        design <- source_env$make_replication_design(replication_id, 80L, setting$covariate)
        template <- source_env$make_taxon_template(setting$n_signal_target, setting$template_role)
        simulation <- source_env$simulate_setting(replication_id, setting, design, template)
        result <- source_env$run_zinq(simulation, replication_id)
        result <- as.data.table(source_env$apply_common_multiplicity(result))
        result[, replication := replication_id]
        result <- merge(result, as.data.table(simulation$truth), by = "taxon", sort = FALSE)
        result[, truth_for_component := fifelse(component == "observed_prevalence",
                                                truth_observed_prevalence, truth_abundance)]
        stopifnot(nrow(result) == 100L, all(result$method == "ZINQ"),
                  all(result[, .N, by = component]$N == 50L),
                  !any(startsWith(result$status, "error:")))
        elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        saveRDS(list(results = result, elapsed_seconds = elapsed), destination, compress = "gzip")
        if (job$replication %% 10L == 0L) {
            message(sprintf("ZINQ | %s | replication %d | %.1f seconds",
                            setting$setting_id, replication_id, elapsed))
        }
    }
}

if (length(arguments)) {
    mode <- arguments[[1L]]
    if (mode == "prepare") prepare()
    else if (mode == "preflight") check_preflight()
    else if (mode == "worker") {
        worker <- as.integer(arguments[[2L]])
        workers <- as.integer(arguments[[3L]])
        for (i in seq.int(worker, nrow(jobs), by = workers)) run_job(jobs[i])
    } else if (mode == "zinq") {
        run_original_zinq(as.integer(arguments[[2L]]), as.integer(arguments[[3L]]))
    } else stop("Use prepare, preflight, worker <index> <workers>, or zinq <index> <workers>.")
}
