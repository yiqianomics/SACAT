#!/usr/bin/env Rscript

# Keep complete p-value families and their truth labels for table reproduction.
suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || !args[[1L]] %in% c("base", "joint", "abundance")) {
    stop("Usage: Rscript prepare_results.R base|joint|abundance INPUT_DIRECTORY [DATA_DIRECTORY]")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
study_dir <- dirname(normalizePath(script))
mode <- args[[1L]]
input <- normalizePath(args[[2L]], mustWork = TRUE)
output <- if (length(args) >= 3L) args[[3L]] else file.path(study_dir, "data")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
columns <- c("setting_id", "replication", "taxon", "method", "component",
             "p_value", "available", "truth_for_component", "scenario_signal", "status")
write_compact <- function(data, filename, keep = columns) {
    data <- copy(data[, ..keep])
    stopifnot(!anyDuplicated(data[, .(setting_id, replication, taxon, method, component)]),
              all(data[, .N, by = .(setting_id, replication, method, component)]$N == 50L))
    # Retain CSV p-value text; use full double precision for fitted R objects.
    if (!is.character(data$p_value)) {
        data[, p_value := fifelse(is.na(p_value), "", sprintf("%.17g", p_value))]
    }
    path <- file.path(output, filename)
    fwrite(data, path, na = "")
    message(filename, ": ", nrow(data), " tests")
}
copy_input <- function(from, to) {
    stopifnot(file.copy(file.path(input, from), file.path(output, to), overwrite = TRUE))
}
versions_path <- file.path(output, "software_versions.csv")
versions <- if (file.exists(versions_path)) fread(versions_path) else
    data.table(run = character(), phase = character(), package = character(), version = character())
versions <- versions[run != mode]

if (mode == "base") {
    data <- fread(file.path(input, "combined_taxon_method_results.csv.gz"),
                  select = c(columns, "n_per_group"), colClasses = c(p_value = "character"))
    stopifnot(nrow(data) == 5100000L)
    for (n in c(60L, 80L, 120L)) {
        write_compact(data[n_per_group == n], sprintf("base_n%03d.csv.gz", n))
    }
    copy_input("simulation_design.csv", "base_design.csv")
    copy_input("combined_simulation_truth.csv.gz", "base_truth.csv.gz")
    current <- fread(file.path(input, "package_versions.csv"))[, .(
        run = "base", phase = "Original 120-setting run", package, version)]
    repair_file <- file.path(input, "profile_repair_versions.csv")
    if (file.exists(repair_file)) {
        repair <- fread(repair_file)
        current <- rbindlist(list(current, data.table(
            run = "base", phase = paste(repair$reanalyzed_datasets, "datasets reanalyzed"),
            package = c("SACAT", "R"),
            version = c(repair$reanalysis_version,
                        sub("^R version ([^ ]+).*", "\\1", repair$R_reanalysis)))))
    }
} else if (mode == "joint") {
    keep <- c(setdiff(columns, "scenario_signal"), "truth_any", "primary_comparison")
    data <- fread(file.path(input, "taxon_results.csv.gz"), select = keep,
                  colClasses = c(p_value = "character"))
    stopifnot(nrow(data) == 3150000L)
    setnames(data, "truth_any", "scenario_signal")
    write_compact(data, "joint.csv.gz", c(columns, "primary_comparison"))
    copy_input("design.csv", "joint_design.csv")
    copy_input("simulation_truth.csv.gz", "joint_truth.csv.gz")
    current <- unique(fread(file.path(input, "software_versions.csv.gz"))[, .(
        run = "joint", phase = "42-setting joint run", package, version)])
} else {
    fitted <- list.files(file.path(input, "results"), pattern = "\\.rds$", full.names = TRUE)
    rerun <- list.files(file.path(input, "zinq_original"), pattern = "\\.rds$", full.names = TRUE)
    stopifnot(length(fitted) == 600L, length(rerun) == 1000L)
    pairs <- data.table(method = c("SACAT", "ZINQ", "MaAsLin 3"),
                        component = c("present_conditional_abundance",
                                      "detected_quantile_abundance", "detected_log_abundance"))
    data <- rbindlist(lapply(c(fitted, rerun), function(path) {
        result <- as.data.table(readRDS(path)$results)
        result <- merge(result, pairs, by = c("method", "component"), sort = FALSE)
        result[, ..columns]
    }))
    stopifnot(nrow(data) == 140000L)
    data[, p_value := fifelse(is.na(p_value), "", sprintf("%.17g", p_value))]
    design <- fread(file.path(output, "base_design.csv"))
    null_settings <- design[n_per_group == 80L & scenario == "present_conditional_abundance" &
                               effect_parameter == 0, setting_id]
    nulls <- fread(file.path(output, "base_n080.csv.gz"), colClasses = c(p_value = "character"))[
        setting_id %in% null_settings & method %in% c("SACAT", "MaAsLin 3")]
    nulls <- merge(nulls, pairs, by = c("method", "component"), sort = FALSE)
    stopifnot(nrow(nulls) == 20000L)
    nulls[, status := "saved_null_reference"]
    data <- rbindlist(list(data, nulls[, ..columns]), use.names = TRUE)
    write_compact(data, "abundance.csv.gz")
    copy_input("simulation_design.csv", "abundance_design.csv")
    current <- fread(file.path(input, "package_versions.csv"))[, .(
        run = "abundance", phase = "Six added settings and ZINQ grid rerun", package, version)]
    runtime <- fread(file.path(input, "runtime.csv"))
    current <- rbindlist(list(current, data.table(
        run = "abundance", phase = "Six added settings and ZINQ grid rerun", package = "R",
        version = sub("^R version ([^ ]+).*", "\\1", runtime$R))))
}
versions <- rbindlist(list(versions, current), use.names = TRUE)
fwrite(versions, versions_path)
