#!/usr/bin/env Rscript

# Small-sample null calibration derived from the frozen SparseDOSSA2 design.
# The generator, panel, seeds, method wrappers, and method arguments are reused;
# only sample size, null settings, and the required DASRA version differ.

n50_runner_file <- function() {
    frames <- sys.frames()
    candidates <- vapply(frames, function(frame) {
        value <- frame$ofile
        if (is.null(value) || !length(value)) NA_character_ else as.character(value)[1L]
    }, character(1))
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
    if (length(candidates)) {
        return(normalizePath(tail(candidates, 1L), winslash = "/", mustWork = TRUE))
    }
    argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(argument)) stop("Could not determine the n=50 runner path.", call. = FALSE)
    normalizePath(sub("^--file=", "", argument[1L]), winslash = "/", mustWork = TRUE)
}

N50_RUNNER_FILE <- n50_runner_file()
N50_STAGE_DIR <- dirname(N50_RUNNER_FILE)
PARENT_STAGE_DIR <- dirname(N50_STAGE_DIR)
PARENT_RUNNER_FILE <- file.path(PARENT_STAGE_DIR, "simulation.R")
PARENT_CONTRACT_FILE <- file.path(PARENT_STAGE_DIR, "design", "contract.rds")

sys.source(PARENT_RUNNER_FILE, envir = .GlobalEnv)

EXPECTED_PARENT_RUNNER_SHA256 <-
    "02de5b1524e521f2b3703c72adc706624303dde2a49e34aaa8c5ddb0c43b4bdb"
REQUIRED_DASRA_VERSION <- "0.4.1"

assert_true(
    identical(sha256_file(PARENT_RUNNER_FILE), EXPECTED_PARENT_RUNNER_SHA256),
    "The parent SparseDOSSA2 runner differs from the frozen source."
)
assert_true(file.exists(PARENT_CONTRACT_FILE), "The frozen source design is missing.")

RUNNER_VERSION <- "2026-08-24-sparsedossa2-n50-null-v1-dasra-0.4.1"
RUNNER_FILE <- N50_RUNNER_FILE
STAGE_DIR <- N50_STAGE_DIR
METHOD_REFERENCE <- file.path(PARENT_STAGE_DIR, "frozen", "formal_hpc_method_reference.R")

STATIC_CONFIG$default_replications <- 100L
STATIC_CONFIG$n_per_group <- 50L
STATIC_CONFIG$depth_medians$balanced_low <- c(control = 1500, case = 1500)
STATIC_CONFIG$source_design_sha256 <- sha256_file(PARENT_CONTRACT_FILE)

# Reuse the panel and generator-only calibration from the frozen n=120 design.
run_generator_calibration <- function() {
    source_contract <- readRDS(PARENT_CONTRACT_FILE)
    assert_true(nrow(source_contract$panel) == STATIC_CONFIG$n_tested_taxa,
                "The frozen source panel does not contain 50 taxa.")
    list(
        panel = source_contract$panel,
        eligible = source_contract$calibration$eligible,
        all_features = source_contract$calibration$all_features,
        effect_taxon = source_contract$calibration$effect_taxon,
        truth_oracle = source_contract$calibration$truth_oracle,
        effect_summary = source_contract$effect_summary,
        chosen_abundance_effect = source_contract$chosen_abundance_effect
    )
}

# Three settings separate ordinary null calibration, low support, and depth imbalance.
parent_make_settings <- make_settings
make_settings <- function(abundance_effect) {
    source <- parent_make_settings(abundance_effect)
    standard <- source[source$setting_id == "NULL_BALANCED", , drop = FALSE]
    low <- standard
    low$setting_id <- "NULL_BALANCED_LOW"
    low$base_setting_id <- low$setting_id
    low$design_id <- "BALANCED_LOW"
    low$depth_design <- "balanced_low"
    low$design_panel <- paste(
        low$sample_size_label, low$signal_fraction_label, "Balanced low depth", sep = "; "
    )
    fourfold <- source[source$setting_id == "NULL_FOURFOLD", , drop = FALSE]
    settings <- rbind(standard, low, fourfold)
    settings$setting_index <- seq_len(nrow(settings))
    settings$generic_abundance <- FALSE
    rownames(settings) <- NULL
    settings
}

# Preserve the frozen wrappers byte-for-byte while replacing only their obsolete
# package-version preflight with the installed DASRA 0.4.1 contract.
load_method_environment <- function(runtime_root) {
    assert_true(
        identical(sha256_file(METHOD_REFERENCE), EXPECTED_METHOD_REFERENCE_SHA256),
        "The frozen formal-method reference has changed."
    )
    expressions <- parse(METHOD_REFERENCE, keep.source = FALSE)
    assert_true(length(expressions) == EXPECTED_METHOD_REFERENCE_EXPRESSIONS,
                "The frozen formal-method reference has an unexpected expression count.")
    last <- expressions[[length(expressions)]]
    assert_true(
        is.call(last) && identical(as.character(last[[1L]]), "main") && length(last) == 1L,
        "The frozen formal-method reference does not end in main()."
    )

    Sys.setenv(
        DASRA_FORMAL_ROOT = runtime_root,
        DASRA_FORMAL_N_TAXA = "50",
        DASRA_FORMAL_CONFOUNDER_GROUP_SHIFT = "0.80",
        DASRA_FORMAL_ALPHA = "0.05",
        DASRA_FORMAL_BASE_SEED = "202608190",
        DASRA_FORMAL_SAVE_DATASETS = "false",
        DASRA_FORMAL_OVERWRITE = "false",
        DASRA_FORMAL_DEPTH_MEDIAN = "8000",
        DASRA_FORMAL_DEPTH_SDLOG = "0.45",
        DASRA_FORMAL_DEPTH_MIN = "1500",
        DASRA_FORMAL_DEPTH_MAX = "40000",
        DASRA_FORMAL_DEPTH_CASE_MULTIPLIER = "0.55",
        DASRA_FORMAL_PROBABILITY_GUARD = "0.70"
    )
    Sys.unsetenv("DASRA_FORMAL_SETTING_FILTER")
    environment <- new.env(parent = globalenv())
    for (index in seq_len(length(expressions) - 1L)) {
        eval(expressions[[index]], envir = environment)
    }

    expected <- list(
        script_version = "2026-08-22-aoas-final-v4-dasra-0.4.0",
        base_seed = 202608190L,
        alpha = 0.05,
        zinq_taus = c(0.25, 0.50, 0.75),
        maaslin3_warn_prevalence = TRUE,
        corncob_robust = TRUE,
        metagenomeseq_maxit = 10L
    )
    for (name in names(expected)) {
        assert_true(
            identical(environment$CONFIG[[name]], expected[[name]]),
            sprintf("Frozen method setting %s changed.", name)
        )
    }
    required <- c("DASRA", "ZINQ", "maaslin3")
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    assert_true(!length(missing), sprintf("Missing packages: %s", paste(missing, collapse = ", ")))
    assert_true(
        identical(as.character(utils::packageVersion("DASRA")), REQUIRED_DASRA_VERSION),
        sprintf("This n=50 analysis requires DASRA %s.", REQUIRED_DASRA_VERSION)
    )
    environment$CONFIG$script_version <- RUNNER_VERSION
    environment
}

collect_package_provenance <- function() {
    method_packages <- c(
        "DASRA", "SparseDOSSA2", "ZINQ", "maaslin3", "data.table", "digest"
    )
    installed <- utils::installed.packages()
    sparse_dependencies <- unlist(tools::package_dependencies(
        packages = "SparseDOSSA2", db = installed,
        which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE
    ), use.names = FALSE)
    packages <- unique(c(method_packages, setdiff(sparse_dependencies, "R")))
    output <- do.call(rbind, lapply(packages, package_fingerprint))
    missing <- output$package[is.na(output$version)]
    if (length(missing)) stopf("Missing packages: %s", paste(missing, collapse = ", "))
    sparse <- output[output$package == "SparseDOSSA2", , drop = FALSE]
    assert_true(identical(sparse$version, EXPECTED_SPARSEDOSSA_VERSION),
                "SparseDOSSA2 must be version 0.99.2.")
    assert_true(identical(sparse$remote_sha, EXPECTED_SPARSEDOSSA_SHA),
                "SparseDOSSA2 does not match the frozen GitHub SHA.")
    dasra <- output[output$package == "DASRA", , drop = FALSE]
    assert_true(identical(dasra$version, REQUIRED_DASRA_VERSION),
                sprintf("DASRA must be version %s.", REQUIRED_DASRA_VERSION))
    output
}

if (sys.nframe() == 0L) main()
