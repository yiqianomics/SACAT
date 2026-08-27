#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
)

RUNNER_VERSION <- "2026-08-26-sparsedossa2-dasra-0.6.0-v2"
EXPECTED_METHOD_REFERENCE_SHA256 <-
    "7597ac45e2540043f1a2c56a824768ecdf7363b30bf77c134660c61f05c62d97"
EXPECTED_METHOD_REFERENCE_EXPRESSIONS <- 122L
EXPECTED_SPARSEDOSSA_VERSION <- "0.99.2"
EXPECTED_SPARSEDOSSA_SHA <-
    "26a998a6e3a5f04d6a86cce14d6d3229ca82633e"

STATIC_CONFIG <- list(
    default_replications = 100L,
    default_workers = 7L,
    n_per_group = 120L,
    n_tested_taxa = 50L,
    n_potential_signals = 30L,
    core_signal_count = 10L,
    stress_signal_counts = c(20L, 30L),
    template = "Stool",
    new_features = FALSE,
    generator_seed = 202608230L,
    panel_calibration_seed = 2026082301L,
    effect_calibration_seed = 2026082302L,
    signal_seed = 2026082303L,
    panel_calibration_samples = 20000L,
    panel_pi0_range = c(0.05, 0.90),
    panel_detection_depth = 1500L,
    panel_detection_threshold = 0.22,
    reservoir_block_fraction = 1 / 60,
    structural_native_effect = log(2),
    abundance_effect_candidates = log(c(1.5, 2, 2.5, 3)),
    mapped_effect_gate = c(0.30, 0.60),
    depth_medians = list(
        balanced = c(control = 4000, case = 4000),
        fourfold = c(control = 1500, case = 6000)
    ),
    depth_sdlog = 0.45,
    depth_bounds = c(300L, 30000L),
    other_taxon = "Other_unmodeled",
    alpha = 0.05,
    original_method_base_seed = 202608190L
)

`%||%` <- function(x, y) {
    if (is.null(x) || !length(x)) y else x
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

assert_true <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
    invisible(TRUE)
}

safe_dir_create <- function(path) {
    if (!dir.exists(path)) {
        ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (!ok && !dir.exists(path)) stopf("Could not create %s", path)
    }
    invisible(path)
}

atomic_save_rds <- function(object, path, compress = "gzip") {
    safe_dir_create(dirname(path))
    temporary <- paste0(path, ".tmp_", Sys.getpid(), "_", format(Sys.time(), "%OS6"))
    saveRDS(object, temporary, compress = compress)
    if (!file.rename(temporary, path)) {
        unlink(temporary, force = TRUE)
        stopf("Could not atomically write %s", path)
    }
    invisible(path)
}

atomic_write_csv <- function(object, path) {
    safe_dir_create(dirname(path))
    temporary <- paste0(path, ".tmp_", Sys.getpid())
    utils::write.csv(object, temporary, row.names = FALSE, na = "")
    if (!file.rename(temporary, path)) {
        unlink(temporary, force = TRUE)
        stopf("Could not atomically write %s", path)
    }
    invisible(path)
}

runner_file <- function() {
    frames <- sys.frames()
    candidates <- vapply(frames, function(frame) {
        value <- frame$ofile %||% NA_character_
        as.character(value)[1L]
    }, character(1))
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
    if (length(candidates)) {
        return(normalizePath(tail(candidates, 1L), winslash = "/", mustWork = TRUE))
    }
    argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (!length(argument)) stop("Could not determine simulation.R path.", call. = FALSE)
    normalizePath(sub("^--file=", "", argument[1L]), winslash = "/", mustWork = TRUE)
}

RUNNER_FILE <- runner_file()
STAGE_DIR <- dirname(RUNNER_FILE)
METHOD_REFERENCE <- file.path(STAGE_DIR, "frozen", "formal_hpc_method_reference.R")
find_local_r_library <- function() {
    override <- trimws(Sys.getenv("DASRA_SPARSEDOSSA_R_LIB", ""))
    ancestors <- unique(c(
        normalizePath(getwd(), winslash = "/", mustWork = TRUE),
        normalizePath(STAGE_DIR, winslash = "/", mustWork = TRUE)
    ))
    expand_ancestors <- function(path) {
        output <- character()
        repeat {
            output <- c(output, path)
            parent <- dirname(path)
            if (identical(parent, path)) break
            path <- parent
        }
        output
    }
    roots <- unique(unlist(lapply(ancestors, expand_ancestors)))
    candidates <- c(
        if (nzchar(override)) override else character(),
        file.path(roots, "Analysis", "RealDataAnalysis", "R_lib")
    )
    candidates <- candidates[dir.exists(candidates)]
    if (length(candidates)) {
        normalizePath(candidates[1L], winslash = "/", mustWork = TRUE)
    } else {
        NA_character_
    }
}

LOCAL_R_LIBRARY <- find_local_r_library()
if (!is.na(LOCAL_R_LIBRARY)) .libPaths(c(LOCAL_R_LIBRARY, .libPaths()))

sha256_file <- function(path) {
    if (!requireNamespace("digest", quietly = TRUE)) {
        stop("The digest package is required.", call. = FALSE)
    }
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

object_sha256 <- function(object) {
    digest::digest(object, algo = "sha256", serialize = TRUE)
}

sanitize_for_csv <- function(x) {
    paste(unlist(x), collapse = ";")
}

# The frozen method environment is parsed without evaluating its final main().
load_method_environment <- function(runtime_root) {
    observed_sha <- sha256_file(METHOD_REFERENCE)
    assert_true(
        identical(observed_sha, EXPECTED_METHOD_REFERENCE_SHA256),
        "The frozen formal-method reference does not match its prespecified SHA-256."
    )
    expressions <- parse(METHOD_REFERENCE, keep.source = FALSE)
    assert_true(
        length(expressions) == EXPECTED_METHOD_REFERENCE_EXPRESSIONS,
        "The frozen formal-method reference has an unexpected expression count."
    )
    last <- expressions[[length(expressions)]]
    assert_true(
        is.call(last) && identical(as.character(last[[1L]]), "main") &&
            length(last) == 1L,
        "The final expression in the method reference is not main()."
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
        script_version = "2026-08-26-aoas-final-v5-dasra-0.6.0",
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
            sprintf("Formal method setting %s differs from the frozen contract.", name)
        )
    }
    environment$check_packages(environment$replicate_packages)
    environment
}

# Package fingerprints record both installed versions and installed file content.
package_fingerprint <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
        return(data.frame(
            package = package, version = NA_character_, remote_sha = NA_character_,
            content_sha256 = NA_character_,
            stringsAsFactors = FALSE
        ))
    }
    path <- normalizePath(system.file(package = package), winslash = "/", mustWork = TRUE)
    files <- list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    files <- files[file.info(files)$isdir %in% FALSE]
    relative <- substring(files, nchar(path) + 2L)
    ordering <- order(relative)
    files <- files[ordering]
    relative <- relative[ordering]
    md5 <- unname(tools::md5sum(files))
    description <- utils::packageDescription(package)
    remote_sha <- description$RemoteSha %||% description$GithubSHA1 %||%
        description$GithubSHA %||% NA_character_
    data.frame(
        package = package,
        version = as.character(utils::packageVersion(package)),
        remote_sha = as.character(remote_sha)[1L],
        content_sha256 = object_sha256(data.frame(file = relative, md5 = md5)),
        stringsAsFactors = FALSE
    )
}

collect_package_provenance <- function() {
    method_packages <- c(
        "DASRA", "SparseDOSSA2", "ZINQ", "maaslin3", "edgeR", "DESeq2",
        "ANCOMBC", "TreeSummarizedExperiment", "SummarizedExperiment",
        "S4Vectors", "MicrobiomeStat", "corncob", "metagenomeSeq", "Biobase",
        "limma", "statmod", "data.table", "digest"
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
    assert_true(
        identical(sparse$version, EXPECTED_SPARSEDOSSA_VERSION),
        "SparseDOSSA2 must be version 0.99.2."
    )
    assert_true(
        identical(sparse$remote_sha, EXPECTED_SPARSEDOSSA_SHA),
        "SparseDOSSA2 does not match the frozen GitHub SHA."
    )
    dasra <- output[output$package == "DASRA", , drop = FALSE]
    assert_true(identical(dasra$version, "0.6.0"), "DASRA must be version 0.6.0.")
    output
}

# Ten settings isolate the global null, component-specific effects, joint effects,
# depth imbalance, and the 40%/60% same-direction reference stress cases.
make_settings <- function(abundance_effect) {
    specification <- data.frame(
        setting_id = c(
            "NULL_BALANCED", "NULL_FOURFOLD", "STRUCT_BALANCED", "STRUCT_FOURFOLD",
            "ABUND_BALANCED", "ABUND_FOURFOLD", "JOINT_BALANCED", "JOINT_FOURFOLD",
            "ABUND40_SAME_BALANCED", "ABUND60_SAME_BALANCED"
        ),
        scenario = c(
            "global_null", "global_null", "structural_only", "structural_only",
            "abundance_only", "abundance_only", "joint", "joint",
            "abundance_same_direction_40", "abundance_same_direction_60"
        ),
        depth_design = c(
            "balanced", "fourfold", "balanced", "fourfold", "balanced", "fourfold",
            "balanced", "fourfold", "balanced", "balanced"
        ),
        signal_count = c(0L, 0L, 10L, 10L, 10L, 10L, 10L, 10L, 20L, 30L),
        spike_structural = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
        spike_abundance = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
        same_direction = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
        generic_abundance = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
        stringsAsFactors = FALSE
    )
    specification$setting_index <- seq_len(nrow(specification))
    specification$base_setting_id <- specification$setting_id
    specification$design_id <- toupper(specification$depth_design)
    specification$n_per_group <- STATIC_CONFIG$n_per_group
    specification$total_sample_size <- 2L * STATIC_CONFIG$n_per_group
    specification$signal_fraction <- specification$signal_count / STATIC_CONFIG$n_tested_taxa
    specification$n_signal_target <- specification$signal_count
    specification$confounding <- ifelse(
        specification$depth_design == "fourfold", "depth_imbalanced", "unconfounded"
    )
    specification$confounded <- specification$depth_design == "fourfold"
    specification$sample_size_label <- sprintf("n/group = %d", STATIC_CONFIG$n_per_group)
    specification$signal_fraction_label <- sprintf(
        "%d%% direct signals", round(100 * specification$signal_fraction)
    )
    specification$confounding_label <- ifelse(
        specification$confounded, "Depth imbalanced", "Unconfounded"
    )
    specification$design_panel <- paste(
        specification$sample_size_label,
        specification$signal_fraction_label,
        ifelse(specification$depth_design == "balanced", "Balanced depth", "Fourfold depth"),
        sep = "; "
    )
    specification$template_role <- STATIC_CONFIG$template
    specification$study <- "SparseDOSSA2_external_robustness"
    specification$dgp <- "SparseDOSSA2_based_Stool_reservoir_semisynthetic"
    specification$effect_level <- ifelse(specification$signal_count == 0L, 0L, 1L)
    specification$effect_index <- specification$effect_level
    specification$effect_parameter <- ifelse(
        specification$spike_abundance, abundance_effect,
        ifelse(specification$spike_structural, STATIC_CONFIG$structural_native_effect, 0)
    )
    specification$effect_measure <- ifelse(
        specification$spike_structural & specification$spike_abundance,
        "native_prevalence_spike_coefficient_and_nonzero_absolute_logFC",
        ifelse(
            specification$spike_structural, "native_prevalence_spike_coefficient",
            ifelse(specification$spike_abundance, "native_nonzero_absolute_logFC", "null")
        )
    )
    specification
}

stable_detection_probability <- function(probability, depth) {
    probability <- pmin(pmax(probability, 0), 1)
    -expm1(depth * log1p(-probability))
}

# A generator-only calibration fixes the evaluable panel and native abundance effect.
run_generator_calibration <- function() {
    suppressPackageStartupMessages(requireNamespace("SparseDOSSA2", quietly = TRUE))
    set.seed(STATIC_CONFIG$panel_calibration_seed)
    generated <- SparseDOSSA2::SparseDOSSA2(
        template = STATIC_CONFIG$template,
        n_sample = STATIC_CONFIG$panel_calibration_samples,
        new_features = STATIC_CONFIG$new_features,
        spike_metadata = "none",
        median_read_depth = STATIC_CONFIG$panel_detection_depth,
        verbose = FALSE
    )
    absolute <- generated$simulated_matrices$a_null
    totals <- colSums(absolute)
    assert_true(all(is.finite(totals) & totals > 0), "Calibration produced an empty community.")
    relative <- sweep(absolute, 2L, totals, "/")
    parameters <- generated$params$feature_param
    feature <- rownames(parameters)
    q_detection <- rowMeans(stable_detection_probability(
        relative, STATIC_CONFIG$panel_detection_depth
    ))
    expected_absolute <- (1 - parameters[, "pi0"]) * exp(
        parameters[, "mu"] + 0.5 * parameters[, "sigma"]^2
    )
    all_features <- data.frame(
        source_feature = feature,
        pi0 = as.numeric(parameters[, "pi0"]),
        mu = as.numeric(parameters[, "mu"]),
        sigma = as.numeric(parameters[, "sigma"]),
        expected_absolute_abundance = as.numeric(expected_absolute),
        expected_detection_at_1500 = as.numeric(q_detection),
        stringsAsFactors = FALSE
    )
    all_features$within_pi0_range <-
        all_features$pi0 >= STATIC_CONFIG$panel_pi0_range[1L] &
        all_features$pi0 <= STATIC_CONFIG$panel_pi0_range[2L]
    all_features$passes_detection <-
        all_features$expected_detection_at_1500 >= STATIC_CONFIG$panel_detection_threshold
    eligible <- all_features[
        all_features$within_pi0_range & all_features$passes_detection, , drop = FALSE
    ]
    eligible <- eligible[order(
        -eligible$expected_absolute_abundance, eligible$source_feature
    ), , drop = FALSE]
    assert_true(nrow(eligible) == 53L, "The frozen Stool calibration must yield 53 eligible taxa.")
    assert_true(nrow(eligible) >= STATIC_CONFIG$n_tested_taxa, "Fewer than 50 taxa passed calibration.")
    eligible$selected <- seq_len(nrow(eligible)) <= STATIC_CONFIG$n_tested_taxa
    eligible$selection_rank <- seq_len(nrow(eligible))
    selected <- eligible[eligible$selected, , drop = FALSE]
    selected$taxon <- sprintf("Taxon_%03d", seq_len(nrow(selected)))
    selected$panel_rank <- seq_len(nrow(selected))
    selected <- selected[, c(
        "taxon", "source_feature", "panel_rank", "pi0", "mu", "sigma",
        "expected_absolute_abundance", "expected_detection_at_1500"
    )]
    assert_true(
        min(selected$expected_detection_at_1500) >= STATIC_CONFIG$panel_detection_threshold,
        "Selected panel violates the detection calibration threshold."
    )

    group_num <- rep(c(0, 1), each = STATIC_CONFIG$panel_calibration_samples / 2L)
    metadata <- matrix(group_num, ncol = 1L, dimnames = list(colnames(absolute), "group_num"))
    calibration_taxa <- selected$taxon
    source_index <- match(calibration_taxa, selected$taxon)
    source_features <- selected$source_feature[source_index]
    untested <- setdiff(rownames(absolute), selected$source_feature)
    reservoir_absolute <- colSums(absolute[untested, , drop = FALSE])
    reservoir_relative <- colSums(relative[untested, , drop = FALSE])
    spike_one <- get("spike_oneA_metadata", asNamespace("SparseDOSSA2"))
    rows <- list()
    candidate_summaries <- list()
    for (candidate_index in seq_along(STATIC_CONFIG$abundance_effect_candidates)) {
        candidate <- STATIC_CONFIG$abundance_effect_candidates[candidate_index]
        for (direction in c(1, -1)) {
            set.seed(STATIC_CONFIG$effect_calibration_seed)
            mapped <- numeric(length(source_features))
            for (j in seq_along(source_features)) {
                source <- source_features[j]
                w <- spike_one(
                    param = parameters[source, ], metadata = metadata,
                    col_abundance = 1L, effect_abundance = direction * candidate
                )
                baseline <- relative[source, ] +
                    STATIC_CONFIG$reservoir_block_fraction * reservoir_relative
                block_weight <- STATIC_CONFIG$reservoir_block_fraction * reservoir_absolute
                denominator <- w + block_weight
                final <- ifelse(denominator > 0, baseline * w / denominator, 0)
                control <- log(final[group_num == 0 & final > 0])
                case <- log(final[group_num == 1 & final > 0])
                baseline_control <- log(relative[source, group_num == 0 & relative[source, ] > 0])
                baseline_case <- log(relative[source, group_num == 1 & relative[source, ] > 0])
                baseline_contrast <- mean(baseline_case) - mean(baseline_control)
                raw_contrast <- mean(case) - mean(control)
                mapped[j] <- raw_contrast - baseline_contrast
                rows[[length(rows) + 1L]] <- data.frame(
                    candidate_native_effect = candidate,
                    candidate_fold_change = exp(candidate),
                    direction = ifelse(direction > 0, "positive", "negative"),
                    taxon = calibration_taxa[j],
                    source_feature = source,
                    sign = direction,
                    null_baseline_conditional_mean_log_relative_contrast = baseline_contrast,
                    raw_intervention_conditional_mean_log_relative_contrast = raw_contrast,
                    mapped_conditional_mean_log_relative_effect = mapped[j],
                    stringsAsFactors = FALSE
                )
            }
            candidate_summaries[[length(candidate_summaries) + 1L]] <- data.frame(
                candidate_native_effect = candidate,
                candidate_fold_change = exp(candidate),
                direction = ifelse(direction > 0, "positive", "negative"),
                median_absolute_mapped_effect = stats::median(abs(mapped)),
                q25_absolute_mapped_effect = unname(stats::quantile(abs(mapped), 0.25)),
                q75_absolute_mapped_effect = unname(stats::quantile(abs(mapped), 0.75)),
                min_absolute_mapped_effect = min(abs(mapped)),
                max_absolute_mapped_effect = max(abs(mapped)),
                stringsAsFactors = FALSE
            )
        }
    }
    effect_taxon <- do.call(rbind, rows)
    effect_summary <- do.call(rbind, candidate_summaries)
    effect_summary$passes_gate <-
        effect_summary$median_absolute_mapped_effect >= STATIC_CONFIG$mapped_effect_gate[1L] &
        effect_summary$median_absolute_mapped_effect <= STATIC_CONFIG$mapped_effect_gate[2L]
    candidate_pass <- vapply(
        STATIC_CONFIG$abundance_effect_candidates,
        function(candidate) all(effect_summary$passes_gate[
            abs(effect_summary$candidate_native_effect - candidate) < 1e-15
        ]),
        logical(1)
    )
    assert_true(any(candidate_pass),
                "No native abundance effect passed the oracle gate in both directions.")
    chosen_effect <- STATIC_CONFIG$abundance_effect_candidates[which(candidate_pass)[1L]]
    chosen <- effect_summary[
        abs(effect_summary$candidate_native_effect - chosen_effect) < 1e-15, , drop = FALSE
    ]

    oracle_rows <- list()
    for (direction in c(1, -1)) {
        for (j in seq_along(source_features)) {
            source <- source_features[j]
            baseline <- relative[source, ] +
                STATIC_CONFIG$reservoir_block_fraction * reservoir_relative
            block_weight <- STATIC_CONFIG$reservoir_block_fraction * reservoir_absolute
            set.seed(STATIC_CONFIG$effect_calibration_seed + j)
            abundance_w <- spike_one(
                param = parameters[source, ], metadata = metadata,
                col_abundance = 1L, effect_abundance = direction * chosen_effect
            )
            abundance_probability <- ifelse(
                abundance_w + block_weight > 0,
                baseline * abundance_w / (abundance_w + block_weight), 0
            )
            set.seed(STATIC_CONFIG$effect_calibration_seed + 1000L + j)
            joint_w <- spike_one(
                param = parameters[source, ], metadata = metadata,
                col_abundance = 1L, effect_abundance = direction * chosen_effect,
                col_prevalence = 1L,
                effect_prevalence = direction * STATIC_CONFIG$structural_native_effect
            )
            joint_probability <- ifelse(
                joint_w + block_weight > 0,
                baseline * joint_w / (joint_w + block_weight), 0
            )
            abundance_detection <- stable_detection_probability(
                abundance_probability, STATIC_CONFIG$panel_detection_depth
            )
            joint_detection <- stable_detection_probability(
                joint_probability, STATIC_CONFIG$panel_detection_depth
            )
            baseline_relative <- relative[source, ]
            baseline_detection <- stable_detection_probability(
                baseline_relative, STATIC_CONFIG$panel_detection_depth
            )
            baseline_log_contrast <-
                mean(log(baseline_relative[group_num == 1 & baseline_relative > 0])) -
                mean(log(baseline_relative[group_num == 0 & baseline_relative > 0]))
            abundance_log_contrast <-
                mean(log(abundance_probability[group_num == 1 & abundance_probability > 0])) -
                mean(log(abundance_probability[group_num == 0 & abundance_probability > 0]))
            baseline_detection_contrast <-
                mean(baseline_detection[group_num == 1]) -
                mean(baseline_detection[group_num == 0])
            abundance_detection_contrast <-
                mean(abundance_detection[group_num == 1]) -
                mean(abundance_detection[group_num == 0])
            joint_detection_contrast <-
                mean(joint_detection[group_num == 1]) -
                mean(joint_detection[group_num == 0])
            oracle_rows[[length(oracle_rows) + 1L]] <- data.frame(
                taxon = calibration_taxa[j],
                source_feature = source,
                direction = ifelse(direction > 0, "positive", "negative"),
                native_prevalence_spike_coefficient =
                    direction * STATIC_CONFIG$structural_native_effect,
                structural_logit_absence_contrast =
                    -direction * STATIC_CONFIG$structural_native_effect,
                native_nonzero_absolute_logFC = direction * chosen_effect,
                null_baseline_conditional_mean_log_relative_contrast =
                    baseline_log_contrast,
                abundance_raw_conditional_mean_log_relative_contrast =
                    abundance_log_contrast,
                abundance_mapped_conditional_mean_log_relative_contrast =
                    abundance_log_contrast - baseline_log_contrast,
                null_baseline_expected_detection_contrast_at_1500 =
                    baseline_detection_contrast,
                abundance_raw_expected_detection_contrast_at_1500 =
                    abundance_detection_contrast,
                abundance_expected_detection_contrast_at_1500 =
                    abundance_detection_contrast - baseline_detection_contrast,
                joint_raw_expected_detection_contrast_at_1500 =
                    joint_detection_contrast,
                joint_expected_detection_contrast_at_1500 =
                    joint_detection_contrast - baseline_detection_contrast,
                stringsAsFactors = FALSE
            )
        }
    }
    truth_oracle <- do.call(rbind, oracle_rows)
    oracle_sign <- ifelse(truth_oracle$direction == "positive", 1, -1)
    assert_true(all(
        oracle_sign * truth_oracle$abundance_mapped_conditional_mean_log_relative_contrast > 1e-6
    ), "A chosen abundance conditional-relative oracle contrast has the wrong direction.")
    assert_true(all(
        oracle_sign * truth_oracle$abundance_expected_detection_contrast_at_1500 > 1e-6
    ), "A chosen abundance detection oracle contrast has the wrong direction.")
    assert_true(all(
        oracle_sign * truth_oracle$joint_expected_detection_contrast_at_1500 > 1e-6
    ), "A chosen joint detection oracle contrast has the wrong direction.")

    eligible$taxon <- NA_character_
    eligible$taxon[eligible$selected] <- selected$taxon
    list(
        all_features = all_features,
        eligible = eligible,
        panel = selected,
        effect_taxon = effect_taxon,
        effect_summary = effect_summary,
        truth_oracle = truth_oracle,
        chosen_abundance_effect = chosen_effect,
        chosen_effect_summary = chosen
    )
}

design_sidecar_objects <- function(contract) {
    config_table <- data.frame(
        name = names(contract$static_config),
        value = vapply(contract$static_config, sanitize_for_csv, character(1)),
        stringsAsFactors = FALSE
    )
    list(
        "settings.csv" = contract$settings,
        "panel.csv" = contract$panel,
        "panel_calibration_eligible.csv" = contract$calibration$eligible,
        "panel_calibration_all_features.csv" = contract$calibration$all_features,
        "abundance_effect_mapping_taxon.csv" = contract$calibration$effect_taxon,
        "abundance_effect_mapping_summary.csv" = contract$effect_summary,
        "truth_oracle.csv" = contract$calibration$truth_oracle,
        "package_provenance.csv" = contract$package_provenance,
        "config.csv" = config_table
    )
}

write_design_sidecars <- function(contract, design_dir) {
    objects <- design_sidecar_objects(contract)
    for (name in names(objects)) {
        atomic_write_csv(objects[[name]], file.path(design_dir, name))
    }
    invisible(names(objects))
}

sidecar_hashes <- function(design_dir, names) {
    paths <- file.path(design_dir, names)
    if (any(!file.exists(paths))) return(NULL)
    setNames(vapply(paths, sha256_file, character(1)), names)
}

validate_design_sidecars <- function(contract, design_dir) {
    expected_names <- names(design_sidecar_objects(contract))
    observed <- sidecar_hashes(design_dir, expected_names)
    assert_true(!is.null(observed), "One or more frozen design sidecars are missing.")
    assert_true(
        identical(observed, contract$sidecar_sha256[expected_names]),
        "One or more frozen design sidecars failed SHA-256 validation."
    )
    invisible(TRUE)
}

ensure_design_sidecars <- function(contract, design_dir) {
    valid <- tryCatch({
        validate_design_sidecars(contract, design_dir)
        TRUE
    }, error = function(error) FALSE)
    if (!valid) {
        write_design_sidecars(contract, design_dir)
        validate_design_sidecars(contract, design_dir)
    }
    invisible(TRUE)
}

contract_payload <- function(contract) {
    list(
        runner_version = contract$runner_version,
        runner_sha256 = contract$runner_sha256,
        method_reference_sha256 = contract$method_reference_sha256,
        static_config = contract$static_config,
        settings = contract$settings,
        panel = contract$panel,
        calibration = contract$calibration,
        effect_summary = contract$effect_summary,
        chosen_abundance_effect = contract$chosen_abundance_effect,
        r_version = contract$r_version,
        platform = contract$platform,
        rng_kind = contract$rng_kind,
        sidecar_sha256 = contract$sidecar_sha256,
        package_fingerprints = contract$package_provenance[, c(
            "package", "version", "remote_sha", "content_sha256"
        )]
    )
}

# Design initialization is completed before any method is run.
initialize_design <- function(root) {
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    safe_dir_create(root)
    design_dir <- file.path(root, "design")
    safe_dir_create(design_dir)
    contract_path <- file.path(design_dir, "contract.rds")
    if (file.exists(contract_path)) {
        contract <- load_contract(
            root, verify_environment = TRUE, verify_sidecars = FALSE
        )
        ensure_design_sidecars(contract, design_dir)
        return(load_contract(root, verify_environment = TRUE, verify_sidecars = TRUE))
    }

    assert_true(
        identical(sha256_file(METHOD_REFERENCE), EXPECTED_METHOD_REFERENCE_SHA256),
        "Method reference SHA-256 mismatch."
    )
    calibration <- run_generator_calibration()
    packages <- collect_package_provenance()
    contract <- list(
        runner_version = RUNNER_VERSION,
        runner_sha256 = sha256_file(RUNNER_FILE),
        method_reference_sha256 = sha256_file(METHOD_REFERENCE),
        static_config = STATIC_CONFIG,
        settings = make_settings(calibration$chosen_abundance_effect),
        panel = calibration$panel,
        calibration = calibration[c(
            "eligible", "all_features", "effect_taxon", "truth_oracle"
        )],
        effect_summary = calibration$effect_summary,
        chosen_abundance_effect = calibration$chosen_abundance_effect,
        package_provenance = packages,
        r_version = R.version.string,
        platform = R.version$platform,
        rng_kind = RNGkind(),
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    write_design_sidecars(contract, design_dir)
    contract$sidecar_sha256 <- sidecar_hashes(
        design_dir, names(design_sidecar_objects(contract))
    )
    contract$analysis_signature <- object_sha256(contract_payload(contract))
    atomic_save_rds(contract, contract_path, compress = "xz")
    load_contract(root, verify_environment = TRUE, verify_sidecars = TRUE)
}

load_contract <- function(root, verify_environment = FALSE, verify_sidecars = TRUE) {
    path <- file.path(root, "design", "contract.rds")
    if (!file.exists(path)) stopf("Design contract not found at %s", path)
    contract <- readRDS(path)
    observed <- object_sha256(contract_payload(contract))
    assert_true(identical(observed, contract$analysis_signature), "Design contract signature mismatch.")
    assert_true(identical(contract$runner_version, RUNNER_VERSION), "Runner version mismatch.")
    assert_true(
        identical(contract$runner_sha256, sha256_file(RUNNER_FILE)),
        "simulation.R changed after design initialization."
    )
    assert_true(
        identical(contract$method_reference_sha256, sha256_file(METHOD_REFERENCE)),
        "Frozen method reference changed after design initialization."
    )
    assert_true(identical(contract$r_version, R.version.string), "R version changed.")
    assert_true(identical(contract$platform, R.version$platform), "R platform changed.")
    assert_true(identical(contract$rng_kind, RNGkind()), "R RNGkind changed.")
    if (isTRUE(verify_sidecars)) {
        validate_design_sidecars(contract, dirname(path))
    }
    if (isTRUE(verify_environment)) {
        current <- collect_package_provenance()
        fields <- c("package", "version", "remote_sha", "content_sha256")
        assert_true(
            identical(contract$package_provenance[, fields], current[, fields]),
            "Installed package content differs from the frozen design provenance."
        )
    }
    contract
}

# Signal taxa are sampled uniformly once per replication and nested as 10/20/30.
signal_plan <- function(replication, panel) {
    set.seed(STATIC_CONFIG$signal_seed + replication * 100003L)
    order <- sample(panel$taxon, STATIC_CONFIG$n_potential_signals, replace = FALSE)
    data.frame(
        taxon = order,
        potential_signal_rank = seq_along(order),
        alternating_sign = rep(c(1, -1), length.out = length(order)),
        stringsAsFactors = FALSE
    )
}

make_metadata <- function(replication, depth_design) {
    n_group <- STATIC_CONFIG$n_per_group
    n <- 2L * n_group
    sample <- sprintf("Sample_%03d", seq_len(n))
    group <- factor(rep(c("control", "case"), each = n_group), levels = c("control", "case"))
    group_num <- as.integer(group == "case")
    set.seed(STATIC_CONFIG$generator_seed + replication * 100003L + 31L)
    z <- stats::rnorm(n)
    medians <- STATIC_CONFIG$depth_medians[[depth_design]]
    set.seed(
        STATIC_CONFIG$generator_seed + replication * 100003L +
            ifelse(depth_design == "balanced", 41L, 43L)
    )
    depth <- round(stats::rlnorm(
        n, meanlog = log(unname(medians[as.character(group)])),
        sdlog = STATIC_CONFIG$depth_sdlog
    ))
    depth <- as.integer(pmin(pmax(depth, STATIC_CONFIG$depth_bounds[1L]), STATIC_CONFIG$depth_bounds[2L]))
    log_depth <- as.numeric(scale(log(depth)))
    metadata <- data.frame(
        group = group, group_num = group_num, z = z, reads = depth,
        log_depth = log_depth, row.names = sample, stringsAsFactors = FALSE
    )
    list(metadata = metadata, depth = depth)
}

make_spike_specification <- function(setting, plan, panel) {
    active <- head(plan, setting$signal_count)
    if (!nrow(active)) return(list(active = active, spike = "none"))
    active$sign <- if (isTRUE(setting$same_direction)) 1 else active$alternating_sign
    active$source_feature <- panel$source_feature[match(active$taxon, panel$taxon)]
    rows <- list()
    if (isTRUE(setting$spike_structural)) {
        rows[[length(rows) + 1L]] <- data.frame(
            metadata_datum = 1L,
            feature_spiked = active$source_feature,
            associated_property = "prevalence",
            effect_size = active$sign * STATIC_CONFIG$structural_native_effect,
            stringsAsFactors = FALSE
        )
    }
    if (isTRUE(setting$spike_abundance)) {
        rows[[length(rows) + 1L]] <- data.frame(
            metadata_datum = 1L,
            feature_spiked = active$source_feature,
            associated_property = "abundance",
            effect_size = active$sign * setting$effect_parameter,
            stringsAsFactors = FALSE
        )
    }
    list(active = active, spike = do.call(rbind, rows))
}

conditional_log_mean <- function(probability, selector) {
    values <- probability[selector & probability > 0]
    if (!length(values)) NA_real_ else mean(log(values))
}

# Each active taxon exchanges mass only with its allocated reservoir block.
construct_reservoir_probabilities <- function(generated, setting, panel, active, metadata) {
    absolute0 <- generated$simulated_matrices$a_null
    absolute_star <- generated$simulated_matrices$a_spiked
    total0 <- colSums(absolute0)
    assert_true(all(is.finite(total0) & total0 > 0), "SparseDOSSA2 generated an empty sample.")
    relative0 <- sweep(absolute0, 2L, total0, "/")
    source <- panel$source_feature
    untested <- setdiff(rownames(relative0), source)
    reservoir_relative <- colSums(relative0[untested, , drop = FALSE])
    reservoir_absolute <- colSums(absolute0[untested, , drop = FALSE])
    panel_probability <- relative0[source, , drop = FALSE]
    rownames(panel_probability) <- panel$taxon
    block_other <- matrix(0, nrow = nrow(active), ncol = ncol(relative0))
    if (nrow(active)) {
        for (index in seq_len(nrow(active))) {
            taxon <- active$taxon[index]
            feature <- panel$source_feature[match(taxon, panel$taxon)]
            baseline_block <- relative0[feature, ] +
                STATIC_CONFIG$reservoir_block_fraction * reservoir_relative
            w <- absolute_star[feature, ]
            w_other <- STATIC_CONFIG$reservoir_block_fraction * reservoir_absolute
            denominator <- w + w_other
            updated <- ifelse(denominator > 0, baseline_block * w / denominator, 0)
            panel_probability[taxon, ] <- updated
            block_other[index, ] <- baseline_block - updated
        }
    }
    other <- reservoir_relative *
        (1 - nrow(active) * STATIC_CONFIG$reservoir_block_fraction) +
        if (nrow(active)) colSums(block_other) else 0
    probability <- rbind(panel_probability, Other_unmodeled = other)
    preadjustment_error <- max(abs(colSums(probability) - 1))
    probability[nrow(probability), ] <- probability[nrow(probability), ] +
        (1 - colSums(probability))
    assert_true(all(is.finite(probability)), "Nonfinite final probabilities.")
    assert_true(min(probability) >= -1e-14, "Negative final probability.")
    probability[probability < 0] <- 0
    assert_true(max(abs(colSums(probability) - 1)) < 1e-12, "Final probabilities do not sum to one.")

    inactive <- setdiff(panel$taxon, active$taxon)
    inactive_error <- if (length(inactive)) {
        max(abs(panel_probability[inactive, , drop = FALSE] -
            relative0[panel$source_feature[match(inactive, panel$taxon)], , drop = FALSE]))
    } else 0
    list(
        probability = probability,
        relative0 = relative0,
        preadjustment_sum_error = preadjustment_error,
        inactive_probability_error = inactive_error,
        minimum_other_probability = min(other),
        reservoir_fraction_min = min(reservoir_relative),
        reservoir_fraction_max = max(reservoir_relative)
    )
}

draw_counts <- function(probability, depth, seed) {
    set.seed(seed)
    counts <- vapply(seq_along(depth), function(index) {
        as.integer(stats::rmultinom(1L, size = depth[index], prob = probability[, index]))
    }, integer(nrow(probability)))
    dimnames(counts) <- dimnames(probability)
    storage.mode(counts) <- "integer"
    counts
}

make_truth <- function(setting, panel, plan, active, probability, relative0, metadata, depth) {
    active_match <- match(panel$taxon, active$taxon)
    is_active <- !is.na(active_match)
    sign <- rep(0, nrow(panel))
    sign[is_active] <- active$sign[active_match[is_active]]
    group_control <- metadata$group == "control"
    group_case <- metadata$group == "case"
    expected_detection <- stable_detection_probability(
        probability[panel$taxon, , drop = FALSE],
        matrix(depth, nrow = nrow(panel), ncol = length(depth), byrow = TRUE)
    )
    truth <- panel
    truth$potential_signal_rank <- plan$potential_signal_rank[match(panel$taxon, plan$taxon)]
    truth$direct_structural_spike <- is_active & isTRUE(setting$spike_structural)
    truth$direct_abundance_spike <- is_active & isTRUE(setting$spike_abundance)
    truth$direct_any_spike <- truth$direct_structural_spike | truth$direct_abundance_spike
    truth$effect_sign <- sign
    truth$native_prevalence_spike_coefficient <- ifelse(
        truth$direct_structural_spike, sign * STATIC_CONFIG$structural_native_effect, 0
    )
    truth$structural_logit_absence_contrast <-
        -truth$native_prevalence_spike_coefficient
    truth$native_nonzero_absolute_logFC <- ifelse(
        truth$direct_abundance_spike, sign * setting$effect_parameter, 0
    )
    truth$truth_structural <- truth$direct_structural_spike
    truth$truth_observed_prevalence <- truth$direct_any_spike
    truth$truth_abundance <- truth$direct_abundance_spike
    truth$generator_baseline_mean_probability_control <- rowMeans(
        relative0[panel$source_feature, group_control, drop = FALSE]
    )
    truth$generator_baseline_mean_probability_case <- rowMeans(
        relative0[panel$source_feature, group_case, drop = FALSE]
    )
    truth$generator_final_mean_probability_control <- rowMeans(
        probability[panel$taxon, group_control, drop = FALSE]
    )
    truth$generator_final_mean_probability_case <- rowMeans(
        probability[panel$taxon, group_case, drop = FALSE]
    )
    truth$generator_expected_detection_control <- rowMeans(
        expected_detection[, group_control, drop = FALSE]
    )
    truth$generator_expected_detection_case <- rowMeans(
        expected_detection[, group_case, drop = FALSE]
    )
    truth$generator_mapped_conditional_mean_log_relative_effect <- vapply(
        panel$taxon, function(taxon) {
            conditional_log_mean(probability[taxon, ], group_case) -
                conditional_log_mean(probability[taxon, ], group_control)
        }, numeric(1)
    )
    context <- c(
        "setting_id", "base_setting_id", "design_id", "n_per_group",
        "total_sample_size", "signal_fraction", "n_signal_target", "confounding",
        "confounded", "sample_size_label", "signal_fraction_label",
        "confounding_label", "design_panel", "template_role", "study", "dgp",
        "scenario", "effect_level", "effect_index", "effect_parameter", "effect_measure"
    )
    for (column in rev(context)) truth <- cbind(setNames(list(setting[[column]][1L]), column), truth)
    rownames(truth) <- NULL
    truth
}

# One simulation object is frozen before any analysis method is called.
simulate_setting <- function(replication, setting, contract) {
    panel <- contract$panel
    plan <- signal_plan(replication, panel)
    metadata_object <- make_metadata(replication, setting$depth_design)
    metadata <- metadata_object$metadata
    depth <- metadata_object$depth
    spike <- make_spike_specification(setting, plan, panel)
    sparse_metadata <- as.matrix(metadata[, c("group_num", "z")])
    set.seed(STATIC_CONFIG$generator_seed + replication * 100003L + 11L)
    generated <- SparseDOSSA2::SparseDOSSA2(
        template = STATIC_CONFIG$template,
        n_sample = nrow(metadata),
        new_features = STATIC_CONFIG$new_features,
        spike_metadata = spike$spike,
        metadata_matrix = sparse_metadata,
        median_read_depth = 4000,
        verbose = FALSE
    )
    reservoir <- construct_reservoir_probabilities(
        generated, setting, panel, spike$active, metadata
    )
    colnames(reservoir$probability) <- rownames(metadata)
    counts <- draw_counts(
        reservoir$probability, depth,
        STATIC_CONFIG$generator_seed + replication * 100003L +
            setting$setting_index * 1009L + 197L
    )
    assert_true(identical(as.integer(colSums(counts)), depth), "Count totals differ from reads.")
    truth <- make_truth(
        setting, panel, plan, spike$active, reservoir$probability,
        reservoir$relative0, metadata, depth
    )
    positive_counts <- data.frame(
        taxon = panel$taxon,
        positive_control = rowSums(counts[panel$taxon, metadata$group == "control", drop = FALSE] > 0),
        positive_case = rowSums(counts[panel$taxon, metadata$group == "case", drop = FALSE] > 0),
        stringsAsFactors = FALSE
    )
    generator_contract <- data.frame(
        replication = replication,
        setting_id = setting$setting_id,
        max_probability_sum_error_before_adjustment = reservoir$preadjustment_sum_error,
        max_inactive_tested_probability_error = reservoir$inactive_probability_error,
        minimum_other_probability = reservoir$minimum_other_probability,
        reservoir_fraction_min = reservoir$reservoir_fraction_min,
        reservoir_fraction_max = reservoir$reservoir_fraction_max,
        minimum_positive_control = min(positive_counts$positive_control),
        minimum_positive_case = min(positive_counts$positive_case),
        stringsAsFactors = FALSE
    )
    assert_true(generator_contract$max_probability_sum_error_before_adjustment < 1e-12,
                "Reservoir probability sum invariant failed.")
    assert_true(identical(generator_contract$max_inactive_tested_probability_error, 0),
                "An inactive tested probability changed.")
    assert_true(generator_contract$minimum_other_probability >= -1e-14,
                "Reservoir construction produced a negative Other probability.")
    list(
        setting = setting,
        counts = counts,
        biological_counts = t(counts[panel$taxon, , drop = FALSE]),
        metadata = metadata,
        depth = depth,
        structural_absence = t(reservoir$probability[panel$taxon, , drop = FALSE] == 0),
        latent_x = t(log(pmax(reservoir$probability[panel$taxon, , drop = FALSE], 1e-300))),
        latent_probability = t(reservoir$probability[panel$taxon, , drop = FALSE]),
        observed_detection = t(counts[panel$taxon, , drop = FALSE] > 0),
        present_but_undetected = t(
            reservoir$probability[panel$taxon, , drop = FALSE] > 0 &
                counts[panel$taxon, , drop = FALSE] == 0
        ),
        truth = truth,
        evaluation_taxa = panel$taxon,
        other_taxon = STATIC_CONFIG$other_taxon,
        dgp = "SparseDOSSA2_based_Stool_reservoir_semisynthetic",
        signal_plan = plan,
        active_signals = spike$active,
        positive_counts = positive_counts,
        generator_contract = generator_contract,
        analysis_signature = contract$analysis_signature
    )
}

validate_simulation <- function(simulation, contract, replication = NULL, setting = NULL) {
    panel <- contract$panel$taxon
    assert_true(identical(simulation$evaluation_taxa, panel), "Evaluation panel order changed.")
    assert_true(identical(rownames(simulation$counts), c(panel, STATIC_CONFIG$other_taxon)),
                "Count row contract failed.")
    assert_true(identical(colnames(simulation$counts), rownames(simulation$metadata)),
                "Count/metadata sample order differs.")
    assert_true(identical(levels(simulation$metadata$group), c("control", "case")),
                "Group reference direction changed.")
    assert_true(identical(as.integer(simulation$metadata$group_num),
                          as.integer(simulation$metadata$group == "case")),
                "group_num coding changed.")
    assert_true(identical(as.integer(simulation$metadata$reads), simulation$depth),
                "reads does not equal full multinomial depth.")
    assert_true(all(is.finite(simulation$metadata$z)), "Nonfinite z covariate.")
    assert_true(all(is.finite(simulation$metadata$log_depth)), "Nonfinite log_depth covariate.")
    assert_true(all(simulation$counts >= 0 & simulation$counts == round(simulation$counts)),
                "Counts are not nonnegative integers.")
    assert_true(identical(as.integer(colSums(simulation$counts)), simulation$depth),
                "Count library totals changed.")
    assert_true(identical(simulation$analysis_signature, contract$analysis_signature),
                "Simulation signature mismatch.")
    assert_true(nrow(simulation$metadata) == 2L * STATIC_CONFIG$n_per_group,
                "Simulation sample size changed.")
    assert_true(nrow(simulation$truth) == STATIC_CONFIG$n_tested_taxa &&
                    identical(simulation$truth$taxon, panel),
                "Simulation truth panel changed.")
    assert_true(!anyDuplicated(simulation$truth$taxon), "Simulation truth has duplicate taxa.")
    if (!is.null(replication)) {
        assert_true(
            identical(as.integer(simulation$generator_contract$replication), as.integer(replication)),
            "Simulation replication identifier changed."
        )
    }
    if (!is.null(setting)) {
        assert_true(identical(simulation$setting$setting_id, setting$setting_id),
                    "Simulation setting identifier changed.")
        assert_true(identical(simulation$generator_contract$setting_id, setting$setting_id),
                    "Generator setting identifier changed.")
    }
    invisible(TRUE)
}

expected_method_pairs <- function(generic_abundance) {
    output <- data.frame(
        method = c("DASRA", "DASRA", "ZINQ", "ZINQ", "MaAsLin 3", "MaAsLin 3"),
        component = c(
            "structural_absence", "present_conditional_abundance",
            "observed_prevalence", "detected_quantile_abundance",
            "observed_prevalence", "detected_log_abundance"
        ),
        stringsAsFactors = FALSE
    )
    if (isTRUE(generic_abundance)) {
        output <- rbind(output, data.frame(
            method = c("edgeR", "DESeq2", "ANCOM-BC2", "LinDA", "corncob", "metagenomeSeq"),
            component = "general_abundance",
            stringsAsFactors = FALSE
        ))
    }
    output
}

validate_method_results <- function(results, simulation) {
    expected_rows <- if (isTRUE(simulation$setting$generic_abundance)) 600L else 300L
    assert_true(nrow(results) == expected_rows, "Unexpected number of method-result rows.")
    required <- c(
        "taxon", "method", "component", "p_value", "available", "status",
        "q_value_bh", "q_value_by", "reject_raw", "reject_bh", "reject_by",
        "replication", "setting_id", "truth_structural",
        "truth_observed_prevalence", "truth_abundance", "truth_for_component"
    )
    assert_true(all(required %in% names(results)), "Method-result schema is incomplete.")
    expected_pairs <- expected_method_pairs(simulation$setting$generic_abundance)
    observed_pairs <- unique(results[, c("method", "component")])
    observed_pairs <- observed_pairs[order(observed_pairs$method, observed_pairs$component), ]
    expected_pairs <- expected_pairs[order(expected_pairs$method, expected_pairs$component), ]
    rownames(observed_pairs) <- rownames(expected_pairs) <- NULL
    assert_true(identical(observed_pairs, expected_pairs), "Method/component family changed.")
    for (index in seq_len(nrow(expected_pairs))) {
        keep <- results$method == expected_pairs$method[index] &
            results$component == expected_pairs$component[index]
        taxa <- results$taxon[keep]
        assert_true(length(taxa) == STATIC_CONFIG$n_tested_taxa,
                    "A method/component does not have exactly 50 rows.")
        assert_true(!anyDuplicated(taxa) && setequal(taxa, simulation$evaluation_taxa),
                    "A method/component taxon family changed.")
    }
    key <- paste(results$method, results$component, results$taxon, sep = "\r")
    assert_true(!anyDuplicated(key), "Method-result keys are not unique.")
    assert_true(all(is.finite(results$q_value_bh) & results$q_value_bh >= 0 & results$q_value_bh <= 1),
                "Invalid common BH value.")
    assert_true(all(is.finite(results$q_value_by) & results$q_value_by >= 0 & results$q_value_by <= 1),
                "Invalid common BY value.")
    assert_true(!anyNA(results[, c("reject_raw", "reject_bh", "reject_by")]),
                "A rejection flag is missing.")
    valid_p <- is.finite(results$p_value) & results$p_value >= 0 & results$p_value <= 1
    assert_true(identical(
        as.logical(results$reject_raw),
        as.logical(results$available & valid_p & results$p_value < STATIC_CONFIG$alpha)
    ), "Raw rejection flags changed.")
    assert_true(identical(
        as.logical(results$reject_bh), as.logical(results$q_value_bh < STATIC_CONFIG$alpha)
    ), "BH rejection flags changed.")
    assert_true(identical(
        as.logical(results$reject_by), as.logical(results$q_value_by < STATIC_CONFIG$alpha)
    ), "BY rejection flags changed.")
    groups <- interaction(results$method, results$component, drop = TRUE, lex.order = TRUE)
    for (indices in split(seq_len(nrow(results)), groups)) {
        p <- results$p_value[indices]
        p_for_adjustment <- ifelse(is.finite(p) & p >= 0 & p <= 1, p, 1)
        assert_true(max(abs(results$q_value_bh[indices] - stats::p.adjust(p_for_adjustment, "BH"))) < 1e-14,
                    "Common BH calculation changed.")
        assert_true(max(abs(results$q_value_by[indices] - stats::p.adjust(p_for_adjustment, "BY"))) < 1e-14,
                    "Common BY calculation changed.")
    }
    expected_truth <- ifelse(
        results$component == "structural_absence", results$truth_structural,
        ifelse(results$component == "observed_prevalence",
               results$truth_observed_prevalence, results$truth_abundance)
    )
    assert_true(identical(as.logical(results$truth_for_component), as.logical(expected_truth)),
                "Component truth mapping changed.")
    assert_true(!anyNA(results[, c(
        "taxon", "method", "component", "available", "replication",
        "setting_id", "truth_structural", "truth_observed_prevalence",
        "truth_abundance", "truth_for_component"
    )]), "Required method-result fields contain missing values.")
    context <- c(
        "setting_id", "base_setting_id", "design_id", "n_per_group",
        "total_sample_size", "signal_fraction", "n_signal_target", "confounding",
        "confounded", "sample_size_label", "signal_fraction_label",
        "confounding_label", "design_panel", "template_role", "study", "dgp",
        "scenario", "effect_level", "effect_index", "effect_parameter", "effect_measure"
    )
    assert_true(all(context %in% names(results)), "Method-result context is incomplete.")
    assert_true(!anyNA(results[, context]), "Method-result context contains missing values.")
    for (column in context) {
        assert_true(length(unique(results[[column]])) == 1L,
                    sprintf("Method-result context %s is not constant.", column))
        assert_true(isTRUE(all.equal(
            results[[column]][1L], simulation$setting[[column]][1L],
            check.attributes = FALSE
        )), sprintf("Method-result context %s differs from its setting.", column))
    }
    truth <- simulation$truth[match(results$taxon, simulation$truth$taxon), , drop = FALSE]
    for (column in c("truth_structural", "truth_observed_prevalence", "truth_abundance")) {
        assert_true(identical(as.logical(results[[column]]), as.logical(truth[[column]])),
                    sprintf("Method-result %s differs from frozen truth.", column))
    }
    unavailable <- is.na(results$available) | !results$available
    assert_true(!any(unavailable & (results$reject_raw | results$reject_bh | results$reject_by)),
                "An unavailable result was marked as rejected.")
    invisible(TRUE)
}

.WORKER_STATE <- new.env(parent = emptyenv())

initialize_worker <- function(root) {
    contract <- load_contract(root, verify_environment = FALSE)
    method_environment <- load_method_environment(file.path(root, "tmp", paste0("worker_", Sys.getpid())))
    .WORKER_STATE$root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    .WORKER_STATE$contract <- contract
    .WORKER_STATE$method_environment <- method_environment
    invisible(TRUE)
}

worker_state <- function(root) {
    normalized <- normalizePath(root, winslash = "/", mustWork = TRUE)
    if (is.null(.WORKER_STATE$root) || !identical(.WORKER_STATE$root, normalized)) {
        initialize_worker(normalized)
    }
    .WORKER_STATE
}

setting_checkpoint_dir <- function(root, replication, setting_id) {
    file.path(root, "checkpoints", sprintf("replicate_%04d", replication), setting_id)
}

valid_simulation_checkpoint <- function(value, contract, replication, setting) {
    tryCatch({
        assert_true(is.list(value), "Simulation checkpoint is not a list.")
        assert_true(identical(value$analysis_signature, contract$analysis_signature),
                    "Simulation checkpoint signature mismatch.")
        assert_true(identical(value$replication, as.integer(replication)),
                    "Simulation checkpoint replication mismatch.")
        assert_true(identical(value$setting_id, setting$setting_id),
                    "Simulation checkpoint setting mismatch.")
        validate_simulation(value$simulation, contract, replication, setting)
        TRUE
    }, error = function(error) FALSE)
}

valid_result_checkpoint <- function(value, simulation, contract, replication, setting) {
    tryCatch({
        assert_true(is.list(value) && isTRUE(value$success),
                    "Result checkpoint is not successful.")
        assert_true(identical(value$analysis_signature, contract$analysis_signature),
                    "Result checkpoint signature mismatch.")
        assert_true(identical(value$replication, as.integer(replication)),
                    "Result checkpoint replication mismatch.")
        assert_true(identical(value$setting_id, setting$setting_id),
                    "Result checkpoint setting mismatch.")
        validate_method_results(value$method_results, simulation)
        assert_true(identical(value$truth, simulation$truth),
                    "Result checkpoint truth differs from the simulation checkpoint.")
        assert_true(identical(value$dataset$counts, simulation$counts) &&
                        identical(value$dataset$metadata, simulation$metadata) &&
                        identical(value$dataset$depth, simulation$depth),
                    "Result checkpoint dataset differs from its simulation checkpoint.")
        TRUE
    }, error = function(error) FALSE)
}

run_one_setting <- function(replication, setting, root, state) {
    checkpoint_dir <- setting_checkpoint_dir(root, replication, setting$setting_id)
    result_path <- file.path(checkpoint_dir, "result.rds")
    simulation_path <- file.path(checkpoint_dir, "simulation.rds")
    simulation_record <- if (file.exists(simulation_path)) {
        tryCatch(readRDS(simulation_path), error = function(e) NULL)
    } else NULL
    if (!valid_simulation_checkpoint(
        simulation_record, state$contract, replication, setting
    )) {
        simulation <- simulate_setting(replication, setting, state$contract)
        validate_simulation(simulation, state$contract, replication, setting)
        simulation_record <- list(
            analysis_signature = state$contract$analysis_signature,
            replication = as.integer(replication),
            setting_id = setting$setting_id,
            simulation = simulation
        )
        atomic_save_rds(simulation_record, simulation_path, compress = "gzip")
    }
    simulation <- simulation_record$simulation
    validate_simulation(simulation, state$contract, replication, setting)
    if (file.exists(result_path)) {
        existing <- tryCatch(readRDS(result_path), error = function(e) NULL)
        if (valid_result_checkpoint(
            existing, simulation, state$contract, replication, setting
        )) return(existing)
    }
    result <- tryCatch({
        method_results <- state$method_environment$analyze_setting(
            simulation = simulation,
            replication_id = replication,
            temporary_root = file.path(root, "tmp", paste0("worker_", Sys.getpid()),
                                       sprintf("replicate_%04d", replication))
        )
        validate_method_results(method_results, simulation)
        list(
            success = TRUE,
            analysis_signature = state$contract$analysis_signature,
            replication = as.integer(replication),
            setting_id = setting$setting_id,
            method_results = method_results,
            truth = simulation$truth,
            dataset = list(
                counts = simulation$counts,
                metadata = simulation$metadata,
                depth = simulation$depth,
                signal_plan = simulation$signal_plan,
                active_signals = simulation$active_signals,
                positive_counts = simulation$positive_counts,
                generator_contract = simulation$generator_contract,
                setting = simulation$setting,
                evaluation_taxa = simulation$evaluation_taxa,
                truth = simulation$truth
            ),
            error = NA_character_
        )
    }, error = function(error) {
        list(
            success = FALSE,
            analysis_signature = state$contract$analysis_signature,
            replication = as.integer(replication),
            setting_id = setting$setting_id,
            method_results = NULL, truth = simulation$truth, dataset = NULL,
            error = conditionMessage(error)
        )
    })
    atomic_save_rds(result, result_path, compress = "gzip")
    result
}

replication_path <- function(root, replication) {
    file.path(root, "results", "raw", sprintf("replicate_%04d.rds", replication))
}

validate_frozen_dataset <- function(dataset, contract, setting, replication) {
    required <- c(
        "counts", "metadata", "depth", "signal_plan", "active_signals",
        "positive_counts", "generator_contract", "setting", "evaluation_taxa", "truth"
    )
    assert_true(is.list(dataset) && all(required %in% names(dataset)),
                "Frozen dataset schema is incomplete.")
    assert_true(identical(dataset$setting$setting_id, setting$setting_id),
                "Frozen dataset setting mismatch.")
    assert_true(identical(dataset$evaluation_taxa, contract$panel$taxon),
                "Frozen dataset panel mismatch.")
    assert_true(identical(rownames(dataset$counts),
                          c(contract$panel$taxon, STATIC_CONFIG$other_taxon)),
                "Frozen dataset count rows changed.")
    assert_true(identical(colnames(dataset$counts), rownames(dataset$metadata)),
                "Frozen dataset sample order changed.")
    assert_true(identical(as.integer(colSums(dataset$counts)), dataset$depth) &&
                    identical(as.integer(dataset$metadata$reads), dataset$depth),
                "Frozen dataset library totals changed.")
    assert_true(identical(levels(dataset$metadata$group), c("control", "case")),
                "Frozen dataset group reference changed.")
    assert_true(nrow(dataset$truth) == STATIC_CONFIG$n_tested_taxa &&
                    identical(dataset$truth$taxon, contract$panel$taxon),
                "Frozen dataset truth panel changed.")
    assert_true(identical(as.integer(dataset$generator_contract$replication),
                          as.integer(replication)) &&
                    identical(dataset$generator_contract$setting_id, setting$setting_id),
                "Frozen generator identifiers changed.")
    invisible(TRUE)
}

validate_replication_record <- function(value, contract, replication) {
    settings <- contract$settings
    assert_true(is.list(value) && isTRUE(value$success),
                "Replication aggregate is not successful.")
    assert_true(identical(value$analysis_signature, contract$analysis_signature),
                "Replication aggregate signature mismatch.")
    assert_true(identical(value$replication, as.integer(replication)),
                "Replication aggregate identifier mismatch.")
    assert_true(identical(value$setting_ids, settings$setting_id),
                "Replication aggregate setting order changed.")
    assert_true(identical(names(value$datasets), settings$setting_id),
                "Replication dataset names changed.")
    expected_rows <- sum(ifelse(settings$generic_abundance, 600L, 300L))
    assert_true(nrow(value$method_results) == expected_rows,
                "Replication aggregate method-result row count changed.")
    assert_true(nrow(value$truth) == nrow(settings) * STATIC_CONFIG$n_tested_taxa,
                "Replication aggregate truth row count changed.")
    result_key <- paste(
        value$method_results$setting_id, value$method_results$method,
        value$method_results$component, value$method_results$taxon, sep = "\r"
    )
    truth_key <- paste(value$truth$setting_id, value$truth$taxon, sep = "\r")
    assert_true(!anyDuplicated(result_key), "Replication method-result keys are duplicated.")
    assert_true(!anyDuplicated(truth_key), "Replication truth keys are duplicated.")
    for (index in seq_len(nrow(settings))) {
        setting <- settings[index, , drop = FALSE]
        dataset <- value$datasets[[setting$setting_id]]
        validate_frozen_dataset(dataset, contract, setting, replication)
        results <- value$method_results[
            value$method_results$setting_id == setting$setting_id, , drop = FALSE
        ]
        simulation_stub <- list(
            setting = setting,
            evaluation_taxa = contract$panel$taxon,
            truth = dataset$truth
        )
        validate_method_results(results, simulation_stub)
        truth <- value$truth[value$truth$setting_id == setting$setting_id, , drop = FALSE]
        dataset_truth <- dataset$truth
        rownames(truth) <- rownames(dataset_truth) <- NULL
        assert_true(identical(truth, dataset_truth),
                    "Replication truth differs from its frozen dataset.")
        assert_true(all(results$replication == replication),
                    "Method-result replication identifier changed.")
    }
    invisible(TRUE)
}

valid_replication_file <- function(path, contract, replication) {
    if (!file.exists(path)) return(FALSE)
    value <- tryCatch(readRDS(path), error = function(e) NULL)
    tryCatch({
        validate_replication_record(value, contract, replication)
        TRUE
    }, error = function(error) FALSE)
}

run_replication <- function(replication, root) {
    started <- proc.time()[["elapsed"]]
    state <- worker_state(root)
    settings <- state$contract$settings
    path <- replication_path(root, replication)
    if (valid_replication_file(path, state$contract, replication)) {
        return(data.frame(replication = replication, success = TRUE, resumed = TRUE,
                          elapsed_seconds = 0, error = NA_character_))
    }
    records <- lapply(seq_len(nrow(settings)), function(index) {
        run_one_setting(replication, settings[index, , drop = FALSE], root, state)
    })
    failed <- !vapply(records, function(record) isTRUE(record$success), logical(1))
    if (any(failed)) {
        errors <- paste(vapply(records[failed], `[[`, character(1), "error"), collapse = " | ")
        return(data.frame(replication = replication, success = FALSE, resumed = FALSE,
                          elapsed_seconds = proc.time()[["elapsed"]] - started, error = errors))
    }
    output <- list(
        success = TRUE,
        analysis_signature = state$contract$analysis_signature,
        replication = as.integer(replication),
        setting_ids = settings$setting_id,
        method_results = do.call(rbind, lapply(records, `[[`, "method_results")),
        truth = do.call(rbind, lapply(records, `[[`, "truth")),
        datasets = setNames(lapply(records, `[[`, "dataset"), settings$setting_id),
        completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    atomic_save_rds(output, path, compress = "xz")
    assert_true(
        valid_replication_file(path, state$contract, replication),
        "Replication aggregate failed read-back validation."
    )
    checkpoint_root <- file.path(root, "checkpoints", sprintf("replicate_%04d", replication))
    unlink(checkpoint_root, recursive = TRUE, force = TRUE)
    data.frame(replication = replication, success = TRUE, resumed = FALSE,
               elapsed_seconds = proc.time()[["elapsed"]] - started, error = NA_character_)
}

completed_replications <- function(root, contract, replications) {
    replications[vapply(replications, function(replication) {
        valid_replication_file(replication_path(root, replication), contract, replication)
    }, logical(1))]
}

# Replication-level summaries retain unavailable taxa as non-rejections.
summarize_results <- function(root, replications) {
    contract <- load_contract(root, verify_environment = FALSE)
    complete <- completed_replications(root, contract, replications)
    assert_true(identical(complete, as.integer(replications)),
                "Summary requires every requested replication to be complete.")
    decisions <- c(raw = "reject_raw", BH = "reject_bh", BY = "reject_by")
    replication_parts <- vector("list", length(replications))
    availability_parts <- vector("list", length(replications))
    generator_parts <- vector("list", length(replications))
    support_parts <- vector("list", length(replications))
    structural_parts <- vector("list", length(replications))
    runtime_parts <- vector("list", length(replications))
    status_parts <- vector("list", length(replications))
    for (replication_index in seq_along(replications)) {
        record <- readRDS(replication_path(root, replications[replication_index]))
        results <- data.table::as.data.table(record$method_results)
        results[, available_for_summary := !is.na(available) & available]
        results[, truth_for_component := !is.na(truth_for_component) & truth_for_component]
        replication_rows <- lapply(names(decisions), function(multiplicity) {
            decision_column <- decisions[[multiplicity]]
            working <- data.table::copy(results)
            working[, primary_decision :=
                available_for_summary & !is.na(get(decision_column)) &
                    as.logical(get(decision_column))]
            working[, {
                null <- !truth_for_component
                signal <- truth_for_component
                rejected <- primary_decision
                list(
                    n_taxa = .N,
                    n_available = sum(available_for_summary),
                    available_fraction = mean(available_for_summary),
                    rejection_rate = mean(rejected),
                    null_rejection_rate = if (any(null)) mean(rejected[null]) else NA_real_,
                    power = if (any(signal)) mean(rejected[signal]) else NA_real_,
                    fdp = if (sum(rejected) == 0) 0 else sum(rejected & null) / sum(rejected),
                    available_conditional_raw_null_rate = if (
                        multiplicity == "raw" && any(null & available_for_summary)
                    ) mean(rejected[null & available_for_summary]) else NA_real_
                )
            }, by = .(
                replication, setting_id, scenario, depth_design = design_id,
                method, component
            )][, multiplicity := multiplicity]
        })
        replication_parts[[replication_index]] <- data.table::rbindlist(
            replication_rows, fill = TRUE
        )
        availability_parts[[replication_index]] <- results[, .(
            n = .N,
            n_available = sum(available_for_summary),
            n_finite_p = sum(is.finite(p_value) & p_value >= 0 & p_value <= 1)
        ), by = .(setting_id, method, component, status)]
        structural_parts[[replication_index]] <- results[
            method == "DASRA" & component == "structural_absence",
            .(
                returned = sum(available_for_summary),
                regular_score = sum(available_for_summary & is.finite(statistic)),
                conservative_nonregular = sum(available_for_summary & !is.finite(statistic)),
                unavailable = sum(!available_for_summary)
            ),
            by = .(replication, setting_id)
        ]
        runtime_parts[[replication_index]] <- results[, {
            total <- runtime_seconds_total_setting[
                is.finite(runtime_seconds_total_setting)
            ]
            direct <- runtime_seconds[is.finite(runtime_seconds)]
            list(runtime_seconds = if (length(total)) total[1L] else if (length(direct)) direct[1L] else NA_real_)
        }, by = .(replication, setting_id, method)]
        status_parts[[replication_index]] <- results[, .(
            n = .N,
            warning_nonempty = sum(!is.na(warning) & nzchar(trimws(warning)))
        ), by = .(
            setting_id, method, component,
            status = ifelse(is.na(status) | !nzchar(trimws(status)), "<missing>", status)
        )]
        generator_parts[[replication_index]] <- data.table::rbindlist(
            lapply(record$datasets, `[[`, "generator_contract"), fill = TRUE
        )
        support_parts[[replication_index]] <- data.table::rbindlist(
            lapply(names(record$datasets), function(setting_id) {
                value <- data.table::as.data.table(
                    record$datasets[[setting_id]]$positive_counts
                )
                value[, `:=`(replication = record$replication, setting_id = setting_id)]
                value
            }), fill = TRUE
        )
        rm(record, results, replication_rows)
    }
    replication_metrics_wide <- data.table::rbindlist(replication_parts, fill = TRUE)
    generator <- data.table::rbindlist(generator_parts, fill = TRUE)
    positive_support <- data.table::rbindlist(support_parts, fill = TRUE)
    structural_replication <- data.table::rbindlist(structural_parts, fill = TRUE)
    runtime_replication <- data.table::rbindlist(runtime_parts, fill = TRUE)
    replication_metrics <- data.table::melt(
        replication_metrics_wide,
        id.vars = c(
            "replication", "setting_id", "scenario", "depth_design",
            "method", "component", "multiplicity", "n_taxa", "n_available"
        ),
        measure.vars = c(
            "available_fraction", "rejection_rate", "null_rejection_rate",
            "power", "fdp", "available_conditional_raw_null_rate"
        ),
        variable.name = "metric", value.name = "value"
    )
    performance <- replication_metrics[is.finite(value), {
        estimate <- mean(value)
        n <- .N
        mcse <- if (n > 1L) stats::sd(value) / sqrt(n) else NA_real_
        list(
            n_replications = n,
            estimate = estimate,
            mcse = mcse,
            ci95_lower = if (is.finite(mcse)) max(0, estimate - stats::qnorm(0.975) * mcse) else NA_real_,
            ci95_upper = if (is.finite(mcse)) min(1, estimate + stats::qnorm(0.975) * mcse) else NA_real_
        )
    }, by = .(setting_id, scenario, depth_design, method, component, multiplicity, metric)]

    availability <- data.table::rbindlist(availability_parts, fill = TRUE)[, .(
        n = sum(n),
        available_fraction = sum(n_available) / sum(n),
        finite_p_fraction = sum(n_finite_p) / sum(n)
    ), by = .(setting_id, method, component, status)]
    support_summary <- positive_support[, .(
        n_replications = .N,
        mean_positive_control = mean(positive_control),
        mean_positive_case = mean(positive_case),
        probability_control_below_20 = mean(positive_control < 20),
        probability_case_below_20 = mean(positive_case < 20),
        probability_either_group_below_20 = mean(
            positive_control < 20 | positive_case < 20
        )
    ), by = .(setting_id, taxon)]
    structural_contract <- structural_replication[, .(
        n_replications = .N,
        mean_returned = mean(returned),
        mean_regular_score = mean(regular_score),
        mean_conservative_nonregular = mean(conservative_nonregular),
        mean_unavailable = mean(unavailable),
        fraction_regular_score = sum(regular_score) / (50 * .N),
        fraction_conservative_nonregular = sum(conservative_nonregular) / (50 * .N),
        fraction_unavailable = sum(unavailable) / (50 * .N)
    ), by = setting_id]
    runtime_summary <- runtime_replication[is.finite(runtime_seconds), .(
        n_replications = .N,
        median_runtime_seconds = stats::median(runtime_seconds),
        q25_runtime_seconds = unname(stats::quantile(runtime_seconds, 0.25)),
        q75_runtime_seconds = unname(stats::quantile(runtime_seconds, 0.75))
    ), by = .(setting_id, method)]
    status_summary <- data.table::rbindlist(status_parts, fill = TRUE)[, .(
        n = sum(n), warning_nonempty = sum(warning_nonempty)
    ), by = .(setting_id, method, component, status)][, `:=`(
        proportion = n / sum(n),
        warning_nonempty_proportion = warning_nonempty / n
    ), by = .(setting_id, method, component)]
    completion <- data.frame(
        expected_replications = length(replications),
        completed_replications = length(complete),
        analysis_signature = contract$analysis_signature,
        summarized_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        stringsAsFactors = FALSE
    )
    summary_dir <- file.path(root, "summary")
    atomic_write_csv(as.data.frame(performance), file.path(summary_dir, "method_performance.csv"))
    atomic_write_csv(as.data.frame(replication_metrics), file.path(summary_dir, "replication_metrics.csv"))
    atomic_write_csv(as.data.frame(availability), file.path(summary_dir, "availability.csv"))
    atomic_write_csv(as.data.frame(generator), file.path(summary_dir, "generator_contracts.csv"))
    atomic_write_csv(as.data.frame(support_summary), file.path(summary_dir, "positive_support.csv"))
    atomic_write_csv(as.data.frame(structural_contract),
                     file.path(summary_dir, "dasra_structural_contract.csv"))
    atomic_write_csv(as.data.frame(runtime_summary), file.path(summary_dir, "runtime_summary.csv"))
    atomic_write_csv(as.data.frame(status_summary), file.path(summary_dir, "status_summary.csv"))
    atomic_write_csv(completion, file.path(summary_dir, "completion.csv"))
    atomic_save_rds(
        list(
            performance = performance,
            replication_metrics = replication_metrics,
            availability = availability,
            generator_contracts = generator,
            positive_support = support_summary,
            dasra_structural_contract = structural_contract,
            runtime_summary = runtime_summary,
            status_summary = status_summary,
            completion = completion
        ),
        file.path(summary_dir, "summary.rds"), compress = "xz"
    )
    readback <- readRDS(file.path(summary_dir, "summary.rds"))
    assert_true(identical(readback$completion$analysis_signature, contract$analysis_signature),
                "Summary read-back signature mismatch.")
    invisible(completion)
}

cleanup_successful_run <- function(root) {
    for (name in c("logs", "checkpoints", "tmp")) {
        path <- file.path(root, name)
        if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
    }
    invisible(TRUE)
}

run_parallel <- function(root, replications, workers) {
    contract <- initialize_design(root)
    safe_dir_create(file.path(root, "results", "raw"))
    safe_dir_create(file.path(root, "logs"))
    complete <- completed_replications(root, contract, replications)
    pending <- setdiff(replications, complete)
    if (!length(pending)) {
        summarize_results(root, replications)
        return(invisible(data.frame()))
    }
    workers <- min(as.integer(workers), length(pending))
    assert_true(workers >= 1L, "workers must be positive.")
    if (workers == 1L) {
        initialize_worker(root)
        status <- do.call(rbind, lapply(pending, run_replication, root = root))
    } else {
        cluster <- parallel::makePSOCKcluster(
            workers, outfile = file.path(root, "logs", "workers.log")
        )
        on.exit(parallel::stopCluster(cluster), add = TRUE)
        parallel::clusterCall(cluster, function(script, analysis_root, library_path) {
            if (!is.na(library_path) && dir.exists(library_path)) {
                .libPaths(c(library_path, .libPaths()))
            }
            source(script, local = .GlobalEnv)
            initialize_worker(analysis_root)
            TRUE
        }, RUNNER_FILE, root, LOCAL_R_LIBRARY)
        status <- do.call(rbind, parallel::parLapplyLB(
            cluster, pending, function(replication, analysis_root) {
                run_replication(replication, analysis_root)
            }, analysis_root = root
        ))
    }
    atomic_write_csv(status, file.path(root, "logs", "last_run_status.csv"))
    if (any(!status$success)) {
        stopf("%d replication(s) failed; checkpoints were retained for resume.", sum(!status$success))
    }
    complete <- completed_replications(root, contract, replications)
    if (identical(complete, as.integer(replications))) summarize_results(root, replications)
    invisible(status)
}

parse_cli <- function(arguments) {
    mode <- if (length(arguments) && !startsWith(arguments[1L], "--")) arguments[1L] else "design"
    if (length(arguments) && identical(arguments[1L], mode)) arguments <- arguments[-1L]
    values <- list(
        mode = mode,
        root = if (mode == "smoke") {
            file.path(tempdir(), "dasra_sparsedossa2_smoke")
        } else {
            STAGE_DIR
        },
        workers = if (mode == "smoke") 1L else STATIC_CONFIG$default_workers,
        replications = if (mode == "smoke") 1L else STATIC_CONFIG$default_replications
    )
    index <- 1L
    while (index <= length(arguments)) {
        argument <- arguments[index]
        if (grepl("^--[^=]+=", argument)) {
            pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
            name <- pieces[1L]
            value <- paste(pieces[-1L], collapse = "=")
        } else if (startsWith(argument, "--") && index < length(arguments)) {
            name <- sub("^--", "", argument)
            value <- arguments[index + 1L]
            index <- index + 1L
        } else {
            stopf("Unknown command-line argument: %s", argument)
        }
        if (name == "root") values$root <- value
        else if (name == "workers") values$workers <- as.integer(value)
        else if (name == "replications") values$replications <- as.integer(value)
        else stopf("Unknown option --%s", name)
        index <- index + 1L
    }
    values$root <- normalizePath(values$root, winslash = "/", mustWork = FALSE)
    assert_true(is.finite(values$workers) && values$workers >= 1L, "Invalid workers value.")
    assert_true(is.finite(values$replications) && values$replications >= 1L, "Invalid replications value.")
    values
}

main <- function() {
    arguments <- parse_cli(commandArgs(trailingOnly = TRUE))
    replications <- seq_len(arguments$replications)
    if (arguments$mode == "design") {
        contract <- initialize_design(arguments$root)
        print(contract$settings[, c(
            "setting_id", "scenario", "depth_design", "signal_count",
            "generic_abundance", "effect_parameter", "effect_measure"
        )], row.names = FALSE)
        cat(sprintf("\nAnalysis signature: %s\n", contract$analysis_signature))
    } else if (arguments$mode %in% c("run", "resume", "smoke")) {
        if (arguments$mode %in% c("run", "resume")) {
            assert_true(arguments$workers == STATIC_CONFIG$default_workers,
                        "Formal run/resume requires exactly seven PSOCK workers.")
            assert_true(arguments$replications == STATIC_CONFIG$default_replications,
                        "Formal run/resume requires the frozen 100 replications.")
        }
        run_parallel(arguments$root, replications, arguments$workers)
        cleanup_successful_run(arguments$root)
    } else if (arguments$mode == "summarize") {
        summarize_results(arguments$root, replications)
        cleanup_successful_run(arguments$root)
    } else {
        stopf("Unknown mode: %s", arguments$mode)
    }
}

if (sys.nframe() == 0L) main()
