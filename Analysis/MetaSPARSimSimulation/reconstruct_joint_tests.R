#!/usr/bin/env Rscript
# Reconstruct native joint tests from saved component results.
# Usage: Rscript reconstruct_joint_tests.R [DIRECTORY]

suppressPackageStartupMessages(library(data.table))
setDTthreads(2L)
stopifnot(as.character(packageVersion("ZINQ")) == "2.0",
          as.character(packageVersion("maaslin3")) == "1.5.3")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 1L) stop("Supply at most one working directory.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
root <- if (length(arguments)) normalizePath(arguments[1L], mustWork = TRUE) else
    dirname(normalizePath(script, mustWork = TRUE))
data_dir <- file.path(root, "data")
table_dir <- file.path(root, "tables")
taxon_keys <- c("setting_index", "replication", "taxon_index")
valid_p <- function(x) is.finite(x) & x >= 0 & x <= 1
normalize_error <- function(x) {
    x <- as.character(x)
    x[!is.na(x) & !nzchar(trimws(x))] <- NA_character_
    x
}

tests <- fread(file.path(data_dir, "tests.csv"))
design <- fread(file.path(data_dir, "design.csv"))
stopifnot(!anyDuplicated(tests$test_id), !anyDuplicated(design$setting_index))
selected_design <- design[scenario %chin% c("global_null", "joint")]
stopifnot(nrow(selected_design) == 14L,
          selected_design[scenario == "global_null", .N] == 2L,
          selected_design[scenario == "joint", .N] == 12L,
          all(selected_design$n_taxa == 100L),
          all(selected_design[scenario == "global_null"]$n_replications == 500L),
          all(selected_design[scenario == "joint"]$n_replications == 100L))
expected_taxon_rows <- sum(selected_design$n_taxa * selected_design$n_replications)

raw <- as.data.table(readRDS(file.path(data_dir, "native_joint_components.rds")))
required_columns <- c(taxon_keys, "test_id", "p_value", "p_tau_025", "p_tau_050",
                      "p_tau_075", "available", "method_detail", "truth_any")
stopifnot(all(required_columns %chin% names(raw)),
          nrow(raw) == 4L * expected_taxon_rows,
          all(raw$available %in% 0:1), all(raw$truth_any %in% 0:1),
          !anyDuplicated(raw[, c(taxon_keys, "test_id"), with = FALSE]))
raw[, `:=`(available = as.logical(available), truth_any = as.logical(truth_any))]
raw <- merge(raw, tests, by = "test_id", all.x = TRUE, sort = FALSE)
expected_components <- data.table(
    method = c("ZINQ", "ZINQ", "MaAsLin 3", "MaAsLin 3"),
    component = c("observed_prevalence", "detected_quantile_abundance",
                  "observed_prevalence", "detected_log_abundance")
)
observed_components <- raw[, .N, by = .(method, component)]
stopifnot(nrow(observed_components) == 4L,
          all(observed_components$N == expected_taxon_rows),
          nrow(merge(observed_components, expected_components,
                     by = c("method", "component"))) == 4L,
          setequal(raw$setting_index, selected_design$setting_index))

zero <- as.data.table(readRDS(file.path(data_dir, "native_joint_zero_rates.rds")))
stopifnot(all(c(taxon_keys, "zero_rate") %chin% names(zero)),
          nrow(zero) == expected_taxon_rows,
          !anyDuplicated(zero[, ..taxon_keys]),
          all(is.finite(zero$zero_rate) & zero$zero_rate >= 0 & zero$zero_rate <= 1),
          setequal(zero$setting_index, selected_design$setting_index))

pair_components <- function(method_name, abundance_component) {
    prevalence <- raw[method == method_name & component == "observed_prevalence",
        .(setting_index, replication, taxon_index, p_prevalence = p_value,
          prevalence_available = available,
          prevalence_error = normalize_error(method_detail), truth_any)]
    abundance <- raw[method == method_name & component == abundance_component,
        .(setting_index, replication, taxon_index, p_abundance = p_value,
          abundance_available = available,
          abundance_error = normalize_error(method_detail),
          p_tau_025, p_tau_050, p_tau_075, abundance_truth = truth_any)]
    paired <- merge(prevalence, abundance, by = taxon_keys, all = TRUE, sort = FALSE)
    stopifnot(nrow(paired) == expected_taxon_rows,
              !anyNA(paired$prevalence_available), !anyNA(paired$abundance_available),
              !anyNA(paired$truth_any), !anyNA(paired$abundance_truth),
              all(paired$truth_any == paired$abundance_truth))
    paired[, abundance_truth := NULL]
    paired
}

zinq <- pair_components("ZINQ", "detected_quantile_abundance")
zinq <- merge(zinq, zero, by = taxon_keys, all = TRUE, sort = FALSE)
stopifnot(nrow(zinq) == expected_taxon_rows, !anyNA(zinq$zero_rate),
          !anyNA(zinq$truth_any))
taus <- c(.25, .5, .75)
zinq[, p_value := vapply(seq_len(.N), function(i) {
    quantile_p <- c(p_tau_025[i], p_tau_050[i], p_tau_075[i])
    if (!prevalence_available[i] || !valid_p(p_prevalence[i]) ||
        !all(valid_p(quantile_p))) return(NA_real_)
    as.numeric(ZINQ::ZINQ_combination(
        list(pvalue.logistic = p_prevalence[i], pvalue.quantile = quantile_p,
             taus = taus, zerorate = zero_rate[i]),
        method = "Cauchy", taus = taus
    ))
}, numeric(1L))]
zinq[, `:=`(method = "ZINQ", native_fallback = "none",
             native_test = "ZINQ::ZINQ_combination (Cauchy)",
             native_package_version = as.character(packageVersion("ZINQ")))]

maaslin <- pair_components("MaAsLin 3", "detected_log_abundance")
create_combined_pval <- getFromNamespace("create_combined_pval", "maaslin3")
maaslin[, p_value := {
    native_input <- data.frame(linear = p_abundance, logistic = p_prevalence,
                               linear_error = abundance_error,
                               logistic_error = prevalence_error)
    create_combined_pval(native_input, correction = "BH")$pval_joint
}, by = .(setting_index, replication)]
maaslin[, native_fallback := fcase(
    (is.na(p_abundance) | !is.na(abundance_error)) &
        (is.na(p_prevalence) | !is.na(prevalence_error)), "both_unavailable",
    is.na(p_prevalence) | !is.na(prevalence_error), "abundance_only",
    is.na(p_abundance) | !is.na(abundance_error), "prevalence_only",
    default = "none"
)]
maaslin[, `:=`(method = "MaAsLin 3",
                native_test = "maaslin3::create_combined_pval",
                native_package_version = as.character(packageVersion("maaslin3")))]

joint <- rbindlist(list(zinq, maaslin), use.names = TRUE, fill = TRUE)
joint_tests <- tests[component == "native_joint", .(method, test_id)]
stopifnot(nrow(joint_tests) == 2L, !anyDuplicated(joint_tests$method),
          setequal(joint_tests$method, c("ZINQ", "MaAsLin 3")))
joint <- merge(joint, joint_tests, by = "method", all.x = TRUE, sort = FALSE)
joint <- merge(joint, selected_design[, .(setting_index, scenario, include_covariate,
                                        n_taxa, n_signal_target, n_replications)],
               by = "setting_index", all.x = TRUE, sort = FALSE)
joint[, `:=`(available = valid_p(p_value), truth_for_component = truth_any)]
stopifnot(nrow(joint) == 2L * expected_taxon_rows,
          !anyDuplicated(joint[, c(taxon_keys, "test_id"), with = FALSE]),
          all(joint$taxon_index >= 1L & joint$taxon_index <= joint$n_taxa),
          all(joint$replication >= 1L & joint$replication <= joint$n_replications))
families <- joint[, .(n_taxa_observed = .N, n_signals = sum(truth_for_component),
                      n_taxa_expected = unique(n_taxa),
                      n_signals_expected = unique(n_signal_target)),
                  by = .(setting_index, replication, test_id)]
stopifnot(nrow(families) == 2L * sum(selected_design$n_replications),
          all(families$n_taxa_observed == families$n_taxa_expected),
          all(families$n_signals == families$n_signals_expected))
setorderv(joint, c("setting_index", "replication", "test_id", "taxon_index"))

checks <- joint[, .(n_taxon_results = .N, n_unavailable = sum(!available),
                     n_abundance_fallback = sum(native_fallback == "abundance_only"),
                     n_prevalence_fallback = sum(native_fallback == "prevalence_only"),
                     n_both_unavailable = sum(native_fallback == "both_unavailable"),
                     minimum_raw_p = min(p_value, na.rm = TRUE),
                     maximum_raw_p = max(p_value, na.rm = TRUE)),
                 by = .(setting_index, scenario, include_covariate, method,
                        native_test, native_package_version)]
output <- joint[, .(setting_index, replication, taxon_index, test_id,
                    p_value = fifelse(is.na(p_value), "", sprintf("%.17g", p_value)),
                    available = as.integer(available),
                    truth_for_component = as.integer(truth_for_component),
                    truth_any = as.integer(truth_any))]
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(output, file.path(data_dir, "raw_native_joint.csv.gz"), quote = FALSE, na = "")
fwrite(checks, file.path(table_dir, "native_joint_checks.csv"))
message("Reconstructed ", nrow(output), " native joint tests.")
