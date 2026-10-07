#!/usr/bin/env Rscript

# Export compact analysis inputs from saved simulation results.
# Usage: Rscript prepare_inputs.R [simulation_output_directory]
# Requires data.table and R.utils.

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Install the data.table R package.", call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript.", call. = FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Supply at most one simulation output directory.", call. = FALSE)
source_dir <- if (length(args)) normalizePath(args[1L], mustWork = TRUE) else file.path(root, "output")
summary_dir <- file.path(source_dir, "summary")
results_dir <- file.path(source_dir, "results")
output <- file.path(root, "data")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
design <- fread(file.path(summary_dir, "design.csv"))
index <- fread(file.path(summary_dir, "replication_index.csv"))
stopifnot(nrow(design) == 38L, nrow(index) == 4600L,
          !anyDuplicated(index[, .(setting_id, replication)]))

columns <- c("setting_id", "replication", "taxon", "method", "component", "p_value",
             "available", "truth_for_component", "truth_any")
data <- fread(file.path(summary_dir, "taxon_results.csv.gz"), select = columns)
tests <- unique(data[, .(method, component)])
tests <- rbind(tests, data.table(method = c("MaAsLin 3", "ZINQ"), component = "native_joint"))
tests[, test_id := .I]
setcolorder(tests, c("test_id", "method", "component"))
stopifnot(nrow(data) == 5980000L, nrow(tests) == 15L)
fwrite(tests, file.path(output, "tests.csv"))
precise <- function(x) fifelse(is.na(x), NA_character_, sprintf("%.17g", x))
groups <- list(sacat = "SACAT", components = c("ZINQ", "MaAsLin 3"),
               abundance = c("ANCOM-BC2", "DESeq2", "edgeR", "LinDA", "corncob", "metagenomeSeq"))
for (name in names(groups)) {
    selected <- merge(data[method %in% groups[[name]]], tests,
                      by = c("method", "component"), sort = FALSE)
    compact <- selected[, .(
        setting_index = as.integer(sub("setting_", "", setting_id, fixed = TRUE)),
        replication, taxon_index = as.integer(sub("Taxon_", "", taxon, fixed = TRUE)),
        test_id, p_value = precise(p_value), available = as.integer(available),
        truth_for_component = as.integer(truth_for_component), truth_any = as.integer(truth_any))]
    fwrite(compact, file.path(output, paste0("raw_", name, ".csv.gz")))
}
rm(data, selected, compact)
invisible(gc())

# Native joint tests require the original component precision and zero rates.
selected_settings <- design[scenario %in% c("global_null", "joint"), setting_id]
selected_index <- index[setting_id %in% selected_settings]
result_files <- file.path(results_dir, selected_index$result_file)
if (!all(file.exists(result_files))) {
    stop("The complete 2,200 original null/joint replication files are required in results/.", call. = FALSE)
}
components <- vector("list", length(result_files))
zero_rates <- vector("list", length(result_files))
keep <- c("setting_id", "replication", "taxon", "method", "component", "p_value",
          "p_tau_025", "p_tau_050", "p_tau_075", "available", "method_detail", "truth_any")
for (i in seq_along(result_files)) {
    saved <- readRDS(result_files[i])
    result <- as.data.table(saved$results)[method %in% c("ZINQ", "MaAsLin 3"), ..keep]
    stopifnot(nrow(result) == 400L,
              all(result$setting_id == selected_index$setting_id[i]),
              all(result$replication == selected_index$replication[i]),
              !anyDuplicated(result[, .(taxon, method, component)]))
    components[[i]] <- result
    counts <- saved$simulation$counts[saved$simulation$evaluation_taxa, , drop = FALSE]
    stopifnot(nrow(counts) == 100L, ncol(counts) == 243L)
    zero_rates[[i]] <- data.table(setting_id = selected_index$setting_id[i],
        replication = selected_index$replication[i], taxon = rownames(counts),
        zero_rate = rowMeans(counts == 0L))
}
components <- merge(rbindlist(components), tests, by = c("method", "component"), sort = FALSE)
components[!is.na(method_detail) & !nzchar(trimws(method_detail)), method_detail := NA_character_]
components <- components[, .(
    setting_index = as.integer(sub("setting_", "", setting_id, fixed = TRUE)), replication,
    taxon_index = as.integer(sub("Taxon_", "", taxon, fixed = TRUE)), test_id,
    p_value, p_tau_025, p_tau_050, p_tau_075,
    available = as.integer(available), method_detail, truth_any = as.integer(truth_any))]
zero_rates <- rbindlist(zero_rates)[, .(
    setting_index = as.integer(sub("setting_", "", setting_id, fixed = TRUE)), replication,
    taxon_index = as.integer(sub("Taxon_", "", taxon, fixed = TRUE)), zero_rate)]
stopifnot(nrow(components) == 880000L, nrow(zero_rates) == 220000L)
saveRDS(as.data.frame(components), file.path(output, "native_joint_components.rds"))
saveRDS(as.data.frame(zero_rates), file.path(output, "native_joint_zero_rates.rds"))
for (name in c("design.csv", "taxon_template.csv", "replication_index.csv", "software_versions.csv")) {
    file.copy(file.path(summary_dir, name), file.path(output, name), overwrite = TRUE)
}
configuration <- readRDS(file.path(summary_dir, "configuration.rds"))
configuration$root <- NULL
fwrite(data.table(parameter = names(configuration),
    value = vapply(configuration, function(x) paste(x, collapse = ","), character(1L))),
    file.path(output, "configuration.csv"))
message("Exported inputs to ", output, ". Run reconstruct_joint_tests.R, then summarize_results.R.")
