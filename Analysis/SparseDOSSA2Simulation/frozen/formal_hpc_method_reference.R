#!/usr/bin/env Rscript

# Final AOAS-oriented HPC simulation study for DASRA (version 5)
#
# The prespecified design crosses three sample sizes, two signal densities and
# confounded/unconfounded covariate structures with an exact model-validation
# study, component-specificity checks, a depth-imbalance negative control and a
# correlated-community robustness study. One Slurm array task runs one complete
# independent replication over the full setting grid.
#
# Usage
#   Rscript dasra_formal_hpc_simulation.R replicate <replication_id>
#   Rscript dasra_formal_hpc_simulation.R summarize
#   Rscript dasra_formal_hpc_simulation.R design
#   Rscript dasra_formal_hpc_simulation.R preflight
#
# Raw taxon-level p-values, adjusted p-values, diagnostics, simulation truth,
# latent states and count data are saved for reproducibility.

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
)

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0L) y else x
}

stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)
messagef <- function(fmt, ...) message(sprintf(fmt, ...))

clamp <- function(x, lower, upper) pmin(pmax(x, lower), upper)
expit <- function(x) stats::plogis(x)
logit <- function(p) stats::qlogis(clamp(p, 1e-12, 1 - 1e-12))

safe_dir_create <- function(path) {
    if (!dir.exists(path)) {
        ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (!ok && !dir.exists(path)) stopf("Could not create directory: %s", path)
    }
    invisible(path)
}

env_int <- function(name, default, minimum = NULL) {
    value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
    if (length(value) != 1L || is.na(value)) stopf("%s must be an integer.", name)
    if (!is.null(minimum) && value < minimum) {
        stopf("%s must be at least %d.", name, minimum)
    }
    value
}

env_num <- function(name, default, minimum = NULL) {
    value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
    if (length(value) != 1L || !is.finite(value)) stopf("%s must be numeric.", name)
    if (!is.null(minimum) && value < minimum) {
        stopf("%s must be at least %g.", name, minimum)
    }
    value
}

env_flag <- function(name, default = FALSE) {
    raw <- tolower(Sys.getenv(name, if (isTRUE(default)) "true" else "false"))
    if (!raw %in% c("true", "false", "1", "0", "yes", "no")) {
        stopf("%s must be TRUE or FALSE.", name)
    }
    raw %in% c("true", "1", "yes")
}

collapse_text <- function(x) {
    x <- unique(trimws(as.character(x)))
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) NA_character_ else paste(x, collapse = " | ")
}

safe_capture <- function(expr) {
    warnings <- character()
    started <- proc.time()[["elapsed"]]
    value <- tryCatch(
        withCallingHandlers(
            expr,
            warning = function(w) {
                warnings <<- c(warnings, conditionMessage(w))
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) structure(
            list(message = conditionMessage(e)), class = "dasra_sim_error"
        )
    )
    elapsed <- proc.time()[["elapsed"]] - started
    list(value = value, warnings = unique(warnings), elapsed = elapsed)
}

is_sim_error <- function(x) inherits(x, "dasra_sim_error")

atomic_save_rds <- function(object, path, compress = "xz") {
    safe_dir_create(dirname(path))
    temporary <- paste0(
        path, ".tmp_", Sys.getpid(), "_", sample.int(.Machine$integer.max, 1L)
    )
    saveRDS(object, temporary, compress = compress)
    if (!file.rename(temporary, path)) {
        unlink(temporary, force = TRUE)
        stopf("Could not atomically write %s", path)
    }
    invisible(path)
}

package_version_or_na <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package))
}

cauchy_combine <- function(p_values, weights = NULL) {
    p_values <- as.numeric(p_values)
    keep <- is.finite(p_values) & p_values >= 0 & p_values <= 1
    if (!any(keep)) return(NA_real_)
    p_values <- clamp(p_values[keep], 1e-15, 1 - 1e-15)
    if (is.null(weights)) {
        weights <- rep(1 / length(p_values), length(p_values))
    } else {
        weights <- as.numeric(weights)[keep]
        if (length(weights) != length(p_values) || any(!is.finite(weights)) ||
            any(weights < 0) || sum(weights) <= 0) {
            return(NA_real_)
        }
        weights <- weights / sum(weights)
    }
    statistic <- sum(weights * tan((0.5 - p_values) * pi))
    clamp(1 - stats::pcauchy(statistic), 0, 1)
}

first_existing_column <- function(data, candidates) {
    candidates <- candidates[candidates %in% names(data)]
    if (!length(candidates)) return(NULL)
    candidates[[1L]]
}

as_numeric_column <- function(data, candidates) {
    column <- first_existing_column(data, candidates)
    if (is.null(column)) return(rep(NA_real_, nrow(data)))
    suppressWarnings(as.numeric(data[[column]]))
}

CONFIG <- list(
    script_version = "2026-08-26-aoas-final-v5-dasra-0.6.0",
    root = Sys.getenv(
        "DASRA_FORMAL_ROOT",
        "/home/zhang.16383/DORAM/dasra_formal_simulation"
    ),
    n_taxa = env_int("DASRA_FORMAL_N_TAXA", 50L, 30L),
    sample_sizes_per_group = c(60L, 80L, 120L),
    signal_fractions = c(0.20, 0.40),
    confounding_levels = c("unconfounded", "confounded"),
    confounder_group_shift = env_num(
        "DASRA_FORMAL_CONFOUNDER_GROUP_SHIFT", 0.80, 0
    ),
    primary_n_per_group = 120L,
    primary_signal_fraction = 0.20,
    alpha = env_num("DASRA_FORMAL_ALPHA", 0.05, 0),
    base_seed = env_int("DASRA_FORMAL_BASE_SEED", 202608190L, 1L),
    save_datasets = env_flag("DASRA_FORMAL_SAVE_DATASETS", TRUE),
    overwrite = env_flag("DASRA_FORMAL_OVERWRITE", FALSE),
    depth_median = env_num("DASRA_FORMAL_DEPTH_MEDIAN", 8000, 500),
    depth_sdlog = env_num("DASRA_FORMAL_DEPTH_SDLOG", 0.45, 0),
    depth_min = env_int("DASRA_FORMAL_DEPTH_MIN", 1500L, 100L),
    depth_max = env_int("DASRA_FORMAL_DEPTH_MAX", 40000L, 1000L),
    depth_case_multiplier = env_num(
        "DASRA_FORMAL_DEPTH_CASE_MULTIPLIER", 0.55, 0.05
    ),
    probability_guard = env_num(
        "DASRA_FORMAL_PROBABILITY_GUARD", 0.70, 0.1
    ),
    zinq_taus = c(0.25, 0.50, 0.75),
    gh_order_truth = 81L,
    observed_prevalence_differences = c(
        0, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32
    ),
    matched_structural_differences = c(
        0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
        0.50
    ),
    abundance_effects = c(
        0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 1.05,
        1.20, 1.40
    ),
    specificity_structural_differences = c(
        0, 0.10, 0.20, 0.30, 0.40, 0.45, 0.50
    ),
    joint_abundance_effects = c(0.40, 0.80),
    joint_structural_differences = c(0.20, 0.35),
    joint_correlation = 0.40,
    joint_other_weight = 1,
    maaslin3_warn_prevalence = TRUE,
    corncob_robust = TRUE,
    metagenomeseq_maxit = 10L
)

if (CONFIG$alpha <= 0 || CONFIG$alpha >= 1) {
    stop("DASRA_FORMAL_ALPHA must lie strictly between zero and one.",
         call. = FALSE)
}
if (any(CONFIG$sample_sizes_per_group < 20L) ||
    any(CONFIG$sample_sizes_per_group != as.integer(CONFIG$sample_sizes_per_group))) {
    stop("All per-group sample sizes must be integers of at least 20.",
         call. = FALSE)
}
if (any(CONFIG$signal_fractions <= 0) ||
    any(CONFIG$signal_fractions >= 0.5)) {
    stop("Signal fractions must lie strictly between zero and one half.",
         call. = FALSE)
}
CONFIG$signal_counts <- as.integer(round(
    CONFIG$n_taxa * CONFIG$signal_fractions
))
if (any(abs(CONFIG$signal_counts / CONFIG$n_taxa -
            CONFIG$signal_fractions) > 1e-12)) {
    stop(
        "Each signal fraction must correspond to an integer number of taxa.",
        call. = FALSE
    )
}
if (any(CONFIG$signal_counts %% 2L != 0L)) {
    stop(
        "Each signal count must be even for balanced positive/negative effects.",
        call. = FALSE
    )
}
if (max(CONFIG$signal_counts) >= CONFIG$n_taxa / 2) {
    stop(
        "The largest signal set must leave a strict majority of null taxa.",
        call. = FALSE
    )
}
if (!all(CONFIG$confounding_levels %in% c("unconfounded", "confounded"))) {
    stop("Unknown confounding level.", call. = FALSE)
}
if (CONFIG$probability_guard >= 0.95) {
    stop("DASRA_FORMAL_PROBABILITY_GUARD must be below 0.95.",
         call. = FALSE)
}

CONFIG$results_dir <- file.path(CONFIG$root, "results")
CONFIG$datasets_dir <- file.path(CONFIG$root, "datasets")
CONFIG$summary_dir <- file.path(CONFIG$root, "summary")
CONFIG$tmp_dir <- file.path(CONFIG$root, "tmp")
CONFIG$logs_dir <- file.path(CONFIG$root, "logs")

replicate_packages <- c(
    "DASRA", "ZINQ", "maaslin3", "edgeR", "DESeq2", "ANCOMBC",
    "TreeSummarizedExperiment", "SummarizedExperiment", "S4Vectors",
    "MicrobiomeStat", "corncob", "metagenomeSeq", "Biobase", "limma",
    "statmod", "data.table"
)
summary_packages <- c("data.table", "ggplot2", "scales")

check_packages <- function(packages) {
    missing <- packages[
        !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing)) {
        stopf(
            "Missing required packages: %s. Run the supplied r451 precheck first.",
            paste(missing, collapse = ", ")
        )
    }

    if ("DASRA" %in% packages) {
        installed_dasra <- as.character(utils::packageVersion("DASRA"))
        if (!identical(installed_dasra, "0.6.0")) {
            stopf(
                "This simulation requires DASRA 0.6.0, but DASRA %s was found in %s.",
                installed_dasra,
                normalizePath(
                    system.file(package = "DASRA"),
                    winslash = "/",
                    mustWork = FALSE
                )
            )
        }
    }

    invisible(TRUE)
}

format_setting_value <- function(x, digits = 3L) {
    gsub("\\.", "p", formatC(x, format = "f", digits = digits))
}

make_base_setting_grid <- function() {
    rows <- list()
    add <- function(base_setting_id, study, dgp, scenario, index, effect,
                    effect_measure, generic_abundance, template_role,
                    description) {
        rows[[length(rows) + 1L]] <<- data.frame(
            base_setting_id = base_setting_id,
            study = study,
            dgp = dgp,
            scenario = scenario,
            effect_level = as.integer(index),
            effect_index = as.integer(index),
            effect_parameter = as.numeric(effect),
            effect_measure = effect_measure,
            generic_abundance = isTRUE(generic_abundance),
            template_role = template_role,
            description = description,
            stringsAsFactors = FALSE
        )
    }

    for (k in seq_along(CONFIG$observed_prevalence_differences)) {
        value <- CONFIG$observed_prevalence_differences[k]
        add(
            paste0("SA_PREV_D", format_setting_value(value)),
            "structural_estimand", "taxonwise",
            "observed_prevalence_only", k - 1L, value,
            "standardized_expected_observed_prevalence_difference", FALSE,
            "detection",
            paste(
                "A conditional-present abundance shift is calibrated to the",
                "target standardized observed-prevalence difference while",
                "structural absence remains unchanged."
            )
        )
    }
    for (k in seq_along(CONFIG$matched_structural_differences)) {
        value <- CONFIG$matched_structural_differences[k]
        add(
            paste0("SA_MATCHED_D", format_setting_value(value)),
            "structural_estimand", "taxonwise",
            "structural_matched_prevalence", k - 1L, value,
            "standardized_structural_absence_probability_difference", FALSE,
            "structural",
            paste(
                "The standardized structural-absence probability changes by",
                "the target amount and conditional-present abundance is",
                "calibrated to match standardized expected observed prevalence."
            )
        )
    }
    for (k in seq_along(CONFIG$abundance_effects)) {
        value <- CONFIG$abundance_effects[k]
        add(
            paste0("AB_MODEL_E", format_setting_value(value)),
            "abundance_comparison", "taxonwise",
            "present_conditional_abundance", k - 1L, value,
            "present_conditional_mean_log_relative_abundance_difference", TRUE,
            "abundance",
            paste(
                "Balanced positive and negative present-conditional abundance",
                "effects among commonly detected taxa, with structural absence",
                "unchanged."
            )
        )
    }
    for (k in seq_along(CONFIG$specificity_structural_differences)) {
        value <- CONFIG$specificity_structural_differences[k]
        add(
            paste0("COMP_STRUCT_D", format_setting_value(value)),
            "component_specificity", "taxonwise",
            "structural_only", k - 1L, value,
            "standardized_structural_absence_probability_difference", FALSE,
            "structural",
            paste(
                "Structural absence changes by the target amount while the",
                "conditional-present abundance distribution remains unchanged."
            )
        )
    }

    add(
        "CAL_DEPTH_NULL", "calibration", "taxonwise",
        "depth_imbalanced_null", 0L, 0,
        "null", TRUE, "neutral",
        paste(
            "No biological group effect; the two groups have different",
            "library-size distributions."
        )
    )
    add(
        "ROBUST_JOINT_NULL", "joint_robustness", "joint_lognormal",
        "joint_global_null", 0L, 0,
        "null", TRUE, "neutral",
        paste(
            "Correlated lognormal absolute abundances, structural states,",
            "closure and multinomial sequencing under the global null."
        )
    )
    for (k in seq_along(CONFIG$joint_abundance_effects)) {
        value <- CONFIG$joint_abundance_effects[k]
        add(
            paste0("ROBUST_JOINT_ABUND_E", format_setting_value(value)),
            "joint_robustness", "joint_lognormal",
            "joint_abundance", k, value,
            "direct_log_absolute_abundance_effect", TRUE, "abundance",
            paste(
                "Balanced direct log-absolute-abundance effects under a",
                "correlated lognormal community followed by closure and",
                "multinomial sequencing."
            )
        )
    }
    for (k in seq_along(CONFIG$joint_structural_differences)) {
        value <- CONFIG$joint_structural_differences[k]
        add(
            paste0("ROBUST_JOINT_STRUCT_D", format_setting_value(value)),
            "joint_robustness", "joint_lognormal",
            "joint_structural", k, value,
            "standardized_structural_absence_probability_difference", FALSE,
            "structural",
            paste(
                "Structural-absence effects under a correlated lognormal",
                "community followed by closure and multinomial sequencing."
            )
        )
    }

    out <- do.call(rbind, rows)
    out$base_setting_index <- seq_len(nrow(out))
    rownames(out) <- NULL
    out
}

make_design_grid <- function() {
    design <- expand.grid(
        n_per_group = CONFIG$sample_sizes_per_group,
        signal_fraction = CONFIG$signal_fractions,
        confounding = CONFIG$confounding_levels,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    design$n_per_group <- as.integer(design$n_per_group)
    design$total_sample_size <- 2L * design$n_per_group
    design$n_signal_target <- as.integer(round(
        CONFIG$n_taxa * design$signal_fraction
    ))
    design$confounded <- design$confounding == "confounded"
    design$design_id <- sprintf(
        "N%03d_S%02d_%s",
        design$n_per_group,
        as.integer(round(100 * design$signal_fraction)),
        ifelse(design$confounded, "C", "U")
    )
    design$sample_size_label <- sprintf(
        "n/group = %d", design$n_per_group
    )
    design$signal_fraction_label <- sprintf(
        "%d%% signal taxa", as.integer(round(100 * design$signal_fraction))
    )
    design$confounding_label <- ifelse(
        design$confounded, "Confounded", "Unconfounded"
    )
    design$design_panel <- paste(
        design$sample_size_label,
        design$signal_fraction_label,
        design$confounding_label,
        sep = "; "
    )
    design$design_index <- seq_len(nrow(design))
    design
}

make_setting_grid <- function() {
    base <- make_base_setting_grid()
    design <- make_design_grid()
    rows <- vector("list", nrow(design))
    for (i in seq_len(nrow(design))) {
        block <- cbind(
            design[rep(i, nrow(base)), , drop = FALSE],
            base
        )
        block$setting_id <- paste(
            block$design_id, block$base_setting_id, sep = "__"
        )
        block$template_id <- paste(
            block$template_role,
            sprintf("S%02d", block$n_signal_target),
            sep = "_"
        )
        block$description <- paste0(
            block$description, " Design stratum: ", block$design_panel, "."
        )
        rows[[i]] <- block
    }
    out <- do.call(rbind, rows)
    out$setting_index <- seq_len(nrow(out))
    rownames(out) <- NULL

    filter_text <- trimws(Sys.getenv("DASRA_FORMAL_SETTING_FILTER", ""))
    if (nzchar(filter_text)) {
        requested <- trimws(strsplit(filter_text, ",", fixed = TRUE)[[1L]])
        keep <- out$setting_id %in% requested |
            out$base_setting_id %in% requested |
            out$design_id %in% requested
        out <- out[keep, , drop = FALSE]
        if (!nrow(out)) {
            stop("DASRA_FORMAL_SETTING_FILTER selected no settings.",
                 call. = FALSE)
        }
    }
    out
}

BASE_SETTINGS <- make_base_setting_grid()
DESIGN_STRATA <- make_design_grid()
SETTINGS <- make_setting_grid()
CONFIG$n_base_settings <- nrow(BASE_SETTINGS)
CONFIG$n_design_strata <- nrow(DESIGN_STRATA)
CONFIG$n_total_settings <- nrow(SETTINGS)

package_versions <- function() {
    packages <- unique(c(replicate_packages, summary_packages))
    data.frame(
        package = packages,
        version = vapply(packages, package_version_or_na, character(1)),
        remote_sha = vapply(packages, function(package) {
            if (!requireNamespace(package, quietly = TRUE)) {
                return(NA_character_)
            }
            description <- utils::packageDescription(package)
            value <- description$RemoteSha %||%
                description$GithubSHA1 %||%
                description$GithubSHA %||%
                NA_character_
            as.character(value)[1L]
        }, character(1)),
        library_path = vapply(packages, function(package) {
            if (!requireNamespace(package, quietly = TRUE)) {
                return(NA_character_)
            }
            normalizePath(
                system.file(package = package),
                winslash = "/", mustWork = FALSE
            )
        }, character(1)),
        stringsAsFactors = FALSE
    )
}
method_configuration <- function() {
    data.frame(
        method = c(
            "DASRA", "ZINQ", "MaAsLin 3", "edgeR", "DESeq2",
            "ANCOM-BC2", "LinDA", "corncob", "metagenomeSeq"
        ),
        configuration = c(
            paste(
                "component='all'; formula=~group+z in both confounded and",
                "unconfounded strata; original total library size supplied;",
                "package formation diagnostics retained"
            ),
            paste(
                "official Firth prevalence and quantile rank-score components;",
                "taus=0.25,0.50,0.75; z and standardized log depth adjusted;",
                "quantile-only weighted Cauchy combination"
            ),
            paste(
                "TSS normalization; LOG transform; augmentation;",
                "standardization; median-comparison abundance correction;",
                "prevalence warning retained; z and log depth adjusted"
            ),
            paste(
                "design=~group+z; filterByExpr; TMM normalization; robust",
                "dispersion estimation and quasi-likelihood test"
            ),
            paste(
                "design=~z+group; positive-count size factors; Wald test;",
                "Cook's handling and independent filtering retained"
            ),
            paste(
                "fix_formula='group+z'; ANCOM-BC2 structural-zero detection;",
                "negative lower bound and pseudocount sensitivity retained;",
                "no manual prevalence filtering"
            ),
            paste(
                "formula=~group+z; count input; default winsorization and",
                "adaptive zero handling; compositional-bias correction;",
                "no manual prevalence filtering"
            ),
            paste(
                "beta-binomial mean model ~group+z and dispersion model",
                "~group+z; robust Wald inference"
            ),
            paste(
                "CSS normalization; model=~group+z in the zero-inflated",
                "Gaussian fitZig procedure; MRcoefs extraction and moderated",
                "inference"
            )
        ),
        stringsAsFactors = FALSE
    )
}

make_taxon_template <- local({
    base_template <- NULL
    cache <- new.env(parent = emptyenv())

    initialize_base <- function() {
        if (!is.null(base_template)) return(base_template)
        set.seed(CONFIG$base_seed + 1701L)
        J <- CONFIG$n_taxa
        taxa <- sprintf("Taxon_%03d", seq_len(J))
        probability <- exp(stats::runif(J, log(1.0e-4), log(2.8e-3)))
        base_template <<- data.frame(
            taxon = taxa,
            taxon_index = seq_len(J),
            baseline_eta = logit(probability),
            baseline_probability = probability,
            baseline_rho = stats::runif(J, 0.15, 0.45),
            sigma = stats::runif(J, 0.50, 0.90),
            eta_z = stats::runif(J, -0.16, 0.16),
            rho_z = stats::runif(J, -0.22, 0.22),
            joint_sd = stats::runif(J, 0.45, 0.75),
            stringsAsFactors = FALSE
        )
        base_template
    }

    function(n_signal, role = c(
        "neutral", "detection", "abundance", "structural"
    )) {
        role <- match.arg(role)
        n_signal <- as.integer(n_signal)
        maximum_signal <- max(CONFIG$signal_counts)
        if (length(n_signal) != 1L || is.na(n_signal) ||
            n_signal < 0L || n_signal > maximum_signal) {
            stop("Invalid signal count for the taxon template.", call. = FALSE)
        }
        key <- paste(role, n_signal, sep = "::")
        if (exists(key, envir = cache, inherits = FALSE)) {
            return(get(key, envir = cache, inherits = FALSE))
        }

        template <- initialize_base()
        candidate_index <- seq_len(maximum_signal)
        if (role == "detection") {
            probability <- exp(seq(
                log(4.0e-5), log(8.0e-5), length.out = maximum_signal
            ))
            template$baseline_probability[candidate_index] <- probability
            template$baseline_eta[candidate_index] <- logit(probability)
            template$baseline_rho[candidate_index] <- seq(
                0.18, 0.28, length.out = maximum_signal
            )
            template$sigma[candidate_index] <- seq(
                0.60, 0.82, length.out = maximum_signal
            )
        } else if (role == "abundance") {
            probability <- exp(seq(
                log(5.0e-4), log(2.0e-3), length.out = maximum_signal
            ))
            template$baseline_probability[candidate_index] <- probability
            template$baseline_eta[candidate_index] <- logit(probability)
            template$baseline_rho[candidate_index] <- seq(
                0.15, 0.35, length.out = maximum_signal
            )
            template$sigma[candidate_index] <- seq(
                0.50, 0.78, length.out = maximum_signal
            )
        } else if (role == "structural") {
            probability <- exp(seq(
                log(4.0e-4), log(1.5e-3), length.out = maximum_signal
            ))
            template$baseline_probability[candidate_index] <- probability
            template$baseline_eta[candidate_index] <- logit(probability)
            template$baseline_rho[candidate_index] <- seq(
                0.55, 0.68, length.out = maximum_signal
            )
            template$sigma[candidate_index] <- seq(
                0.52, 0.82, length.out = maximum_signal
            )
        }

        signal_mask <- rep(FALSE, CONFIG$n_taxa)
        if (n_signal > 0L) signal_mask[seq_len(n_signal)] <- TRUE
        template$detection_signal <- role == "detection" & signal_mask
        template$abundance_signal <- role == "abundance" & signal_mask
        template$structural_signal <- role == "structural" & signal_mask
        template$candidate_signal_pool <- seq_len(CONFIG$n_taxa) %in%
            candidate_index
        template$template_role <- role
        template$n_signal_target <- n_signal
        assign(key, template, envir = cache)
        template
    }
})

normal_gh_rule <- local({
    cache <- new.env(parent = emptyenv())
    function(order = CONFIG$gh_order_truth) {
        key <- as.character(order)
        if (exists(key, envir = cache, inherits = FALSE)) {
            return(get(key, envir = cache, inherits = FALSE))
        }
        rule <- statmod::gauss.quad(order, kind = "hermite")
        value <- list(
            z = sqrt(2) * as.numeric(rule$nodes),
            w = as.numeric(rule$weights) / sqrt(pi)
        )
        assign(key, value, envir = cache)
        value
    }
})

expected_nondetection <- function(eta, sigma, depth,
                                  gh = normal_gh_rule()) {
    eta <- as.numeric(eta)
    depth <- as.numeric(depth)
    if (length(eta) == 1L) eta <- rep(eta, length(depth))
    if (length(eta) != length(depth)) {
        stop("eta and depth must have compatible lengths.", call. = FALSE)
    }
    vapply(seq_along(depth), function(i) {
        probability <- expit(eta[i] + sigma * gh$z)
        sum(gh$w * exp(depth[i] * log1p(-probability)))
    }, numeric(1))
}

expected_detection <- function(eta, sigma, depth,
                               gh = normal_gh_rule()) {
    clamp(1 - expected_nondetection(eta, sigma, depth, gh), 0, 1)
}

expected_log_probability <- function(eta, sigma,
                                     gh = normal_gh_rule()) {
    eta <- as.numeric(eta)
    vapply(eta, function(mu) {
        sum(gh$w * log(expit(mu + sigma * gh$z)))
    }, numeric(1))
}

solve_scalar_root <- function(objective, lower, upper,
                              tolerance = 1e-10) {
    f_lower <- objective(lower)
    f_upper <- objective(upper)
    if (!is.finite(f_lower) || !is.finite(f_upper) ||
        f_lower * f_upper > 0) {
        return(NA_real_)
    }
    stats::uniroot(
        objective, interval = c(lower, upper), tol = tolerance
    )$root
}

solve_observed_prevalence_shift <- function(
        target_difference, eta_reference, eta_case_base,
        rho_reference, rho_case, sigma, depth_reference, depth_case) {
    reference_prevalence <- mean(
        (1 - rho_reference) *
            expected_detection(eta_reference, sigma, depth_reference)
    )
    objective <- function(shift) {
        case_prevalence <- mean(
            (1 - rho_case) *
                expected_detection(eta_case_base + shift, sigma, depth_case)
        )
        case_prevalence - reference_prevalence - target_difference
    }
    solve_scalar_root(objective, lower = -8, upper = 8)
}

solve_structural_shift <- function(
        target_difference, baseline_rho, rho_z, standardized_z) {
    rho_zero <- expit(logit(baseline_rho) + rho_z * standardized_z)
    objective <- function(delta) {
        rho_one <- expit(
            logit(baseline_rho) + rho_z * standardized_z + delta
        )
        mean(rho_zero - rho_one) - target_difference
    }
    solve_scalar_root(objective, lower = -12, upper = 0)
}

solve_present_log_effect_shift <- function(
        target_effect, eta_base, sigma) {
    objective <- function(shift) {
        mean(
            expected_log_probability(eta_base + shift, sigma) -
                expected_log_probability(eta_base, sigma)
        ) - target_effect
    }
    solve_scalar_root(objective, lower = -8, upper = 8)
}

make_replication_design <- function(
        replication_id, n_per_group, confounding) {
    n_per_group <- as.integer(n_per_group)
    confounding <- match.arg(
        as.character(confounding), CONFIG$confounding_levels
    )
    set.seed(
        CONFIG$base_seed + replication_id * 100003L +
            n_per_group * 1009L + 11L
    )
    n0 <- n_per_group
    n <- 2L * n0
    group_numeric <- rep(c(0, 1), each = n0)
    group <- factor(
        ifelse(group_numeric == 0, "control", "case"),
        levels = c("control", "case")
    )

    base_noise <- stats::rnorm(n)
    group_centered <- group_numeric - mean(group_numeric)
    base_noise <- stats::residuals(stats::lm(base_noise ~ group_centered))
    z_raw <- base_noise
    if (confounding == "confounded") {
        z_raw <- z_raw + CONFIG$confounder_group_shift * group_centered
    }
    z <- as.numeric(scale(z_raw))

    reference_depth <- as.integer(round(exp(stats::rnorm(
        n0, log(CONFIG$depth_median), CONFIG$depth_sdlog
    ))))
    reference_depth <- as.integer(clamp(
        reference_depth, CONFIG$depth_min, CONFIG$depth_max
    ))
    balanced_case_depth <- sample(reference_depth, n0, replace = FALSE)
    imbalanced_case_depth <- as.integer(round(exp(stats::rnorm(
        n0,
        log(CONFIG$depth_median * CONFIG$depth_case_multiplier),
        CONFIG$depth_sdlog
    ))))
    imbalanced_case_depth <- as.integer(clamp(
        imbalanced_case_depth, CONFIG$depth_min, CONFIG$depth_max
    ))

    sample_names <- sprintf("Sample_%03d", seq_len(n))
    group_z_correlation <- suppressWarnings(stats::cor(group_numeric, z))
    z_mean_reference <- mean(z[group_numeric == 0])
    z_mean_comparison <- mean(z[group_numeric == 1])
    make_metadata <- function(depth) {
        data.frame(
            group = group,
            group_num = group_numeric,
            z = z,
            reads = as.integer(depth),
            log_depth = as.numeric(scale(log(depth))),
            confounding = confounding,
            confounded = confounding == "confounded",
            n_per_group = n_per_group,
            total_sample_size = n,
            group_z_correlation = group_z_correlation,
            z_mean_reference = z_mean_reference,
            z_mean_comparison = z_mean_comparison,
            z_mean_difference = z_mean_comparison - z_mean_reference,
            row.names = sample_names,
            check.names = FALSE
        )
    }

    list(
        replication = replication_id,
        n = n,
        n_per_group = n_per_group,
        confounding = confounding,
        confounded = confounding == "confounded",
        group_numeric = group_numeric,
        group = group,
        z = z,
        group_z_correlation = group_z_correlation,
        z_mean_reference = z_mean_reference,
        z_mean_comparison = z_mean_comparison,
        sample_names = sample_names,
        balanced_depth = c(reference_depth, balanced_case_depth),
        imbalanced_depth = c(reference_depth, imbalanced_case_depth),
        balanced_metadata = make_metadata(
            c(reference_depth, balanced_case_depth)
        ),
        imbalanced_metadata = make_metadata(
            c(reference_depth, imbalanced_case_depth)
        )
    )
}

scenario_seed_offset <- function(scenario) {
    scenarios <- c(
        "observed_prevalence_only",
        "structural_matched_prevalence",
        "present_conditional_abundance",
        "structural_only",
        "depth_imbalanced_null",
        "joint_global_null",
        "joint_abundance",
        "joint_structural"
    )
    index <- match(scenario, scenarios)
    if (is.na(index)) stopf("Unknown scenario: %s", scenario)
    index * 1000003L
}

make_effect_vectors <- function(setting, metadata, template) {
    J <- nrow(template)
    delta <- numeric(J)
    zeta <- numeric(J)
    absolute_effect <- numeric(J)
    signal <- rep(FALSE, J)
    target <- setting$effect_parameter
    scenario <- setting$scenario

    all_z <- metadata$z
    all_depth <- metadata$reads

    if (scenario == "observed_prevalence_only") {
        signal <- template$detection_signal
        if (target > 0) {
            for (j in which(signal)) {
                eta_reference <- template$baseline_eta[j] +
                    template$eta_z[j] * all_z
                eta_case_base <- eta_reference
                rho_reference <- expit(
                    logit(template$baseline_rho[j]) +
                        template$rho_z[j] * all_z
                )
                rho_case <- rho_reference
                zeta[j] <- solve_observed_prevalence_shift(
                    target_difference = target,
                    eta_reference = eta_reference,
                    eta_case_base = eta_case_base,
                    rho_reference = rho_reference,
                    rho_case = rho_case,
                    sigma = template$sigma[j],
                    depth_reference = all_depth,
                    depth_case = all_depth
                )
            }
        }
    } else if (scenario == "structural_matched_prevalence") {
        signal <- template$structural_signal
        if (target > 0) {
            for (j in which(signal)) {
                delta[j] <- solve_structural_shift(
                    target_difference = target,
                    baseline_rho = template$baseline_rho[j],
                    rho_z = template$rho_z[j],
                    standardized_z = all_z
                )
                eta_reference <- template$baseline_eta[j] +
                    template$eta_z[j] * all_z
                eta_case_base <- eta_reference
                rho_reference <- expit(
                    logit(template$baseline_rho[j]) +
                        template$rho_z[j] * all_z
                )
                rho_case <- expit(
                    logit(template$baseline_rho[j]) +
                        template$rho_z[j] * all_z + delta[j]
                )
                zeta[j] <- solve_observed_prevalence_shift(
                    target_difference = 0,
                    eta_reference = eta_reference,
                    eta_case_base = eta_case_base,
                    rho_reference = rho_reference,
                    rho_case = rho_case,
                    sigma = template$sigma[j],
                    depth_reference = all_depth,
                    depth_case = all_depth
                )
            }
        }
    } else if (scenario == "present_conditional_abundance") {
        signal <- template$abundance_signal
        directions <- rep(c(1, -1), length.out = sum(signal))
        if (target > 0) {
            signal_index <- which(signal)
            for (k in seq_along(signal_index)) {
                j <- signal_index[k]
                eta_base <- template$baseline_eta[j] +
                    template$eta_z[j] * all_z
                zeta[j] <- solve_present_log_effect_shift(
                    target_effect = directions[k] * target,
                    eta_base = eta_base,
                    sigma = template$sigma[j]
                )
            }
        }
    } else if (scenario == "structural_only") {
        signal <- template$structural_signal
        if (target > 0) {
            for (j in which(signal)) {
                delta[j] <- solve_structural_shift(
                    target_difference = target,
                    baseline_rho = template$baseline_rho[j],
                    rho_z = template$rho_z[j],
                    standardized_z = all_z
                )
            }
        }
    } else if (scenario == "depth_imbalanced_null" ||
               scenario == "joint_global_null") {
        signal[] <- FALSE
    } else if (scenario == "joint_abundance") {
        signal <- template$abundance_signal
        directions <- rep(c(1, -1), length.out = sum(signal))
        absolute_effect[signal] <- target * directions
    } else if (scenario == "joint_structural") {
        signal <- template$structural_signal
        for (j in which(signal)) {
            delta[j] <- solve_structural_shift(
                target_difference = target,
                baseline_rho = template$baseline_rho[j],
                rho_z = template$rho_z[j],
                standardized_z = all_z
            )
        }
    } else {
        stopf("Unknown scenario: %s", scenario)
    }

    if (any(signal & !is.finite(delta)) ||
        any(signal & !is.finite(zeta)) ||
        any(signal & !is.finite(absolute_effect))) {
        stopf("Effect calibration failed in setting %s.", setting$setting_id)
    }

    list(
        delta = delta,
        zeta = zeta,
        absolute_effect = absolute_effect,
        signal = signal
    )
}

expected_detected_log_probability <- function(
        eta, sigma, depth, rho, gh = normal_gh_rule()) {
    eta <- as.numeric(eta)
    depth <- as.numeric(depth)
    rho <- as.numeric(rho)
    numerator <- numeric(length(depth))
    denominator <- numeric(length(depth))
    for (i in seq_along(depth)) {
        probability <- expit(eta[i] + sigma * gh$z)
        detection <- 1 - exp(depth[i] * log1p(-probability))
        numerator[i] <- (1 - rho[i]) * sum(
            gh$w * log(probability) * detection
        )
        denominator[i] <- (1 - rho[i]) * sum(gh$w * detection)
    }
    if (sum(denominator) <= 0) return(NA_real_)
    sum(numerator) / sum(denominator)
}

draw_multinomial_counts <- function(probability, depth, other_probability,
                                    seed) {
    n <- nrow(probability)
    J <- ncol(probability)
    set.seed(seed)
    count_by_sample <- vapply(seq_len(n), function(i) {
        as.numeric(stats::rmultinom(
            1L,
            size = depth[i],
            prob = c(probability[i, ], other_probability[i])
        )[, 1L])
    }, numeric(J + 1L))
    count_by_sample
}

make_truth_rows <- function(
        setting, template, metadata, effects, rho, eta,
        structural_absence, latent_probability, observed_detection,
        present_but_undetected, depth, dgp) {
    J <- nrow(template)
    reference <- metadata$group_num == 0
    comparison <- metadata$group_num == 1
    all_z <- metadata$z
    all_depth <- metadata$reads
    rows <- vector("list", J)

    for (j in seq_len(J)) {
        eta_zero_all <- template$baseline_eta[j] +
            template$eta_z[j] * all_z
        eta_one_all <- eta_zero_all + effects$zeta[j]
        rho_zero_all <- expit(
            logit(template$baseline_rho[j]) +
                template$rho_z[j] * all_z
        )
        rho_one_all <- expit(
            logit(template$baseline_rho[j]) +
                template$rho_z[j] * all_z + effects$delta[j]
        )
        detection_zero_all <- expected_detection(
            eta_zero_all, template$sigma[j], all_depth
        )
        detection_one_all <- expected_detection(
            eta_one_all, template$sigma[j], all_depth
        )

        expected_prevalence_reference <- mean(
            (1 - rho[reference, j]) * expected_detection(
                eta[reference, j], template$sigma[j], depth[reference]
            )
        )
        expected_prevalence_comparison <- mean(
            (1 - rho[comparison, j]) * expected_detection(
                eta[comparison, j], template$sigma[j], depth[comparison]
            )
        )
        standardized_prevalence_difference <- mean(
            (1 - rho_one_all) * detection_one_all -
                (1 - rho_zero_all) * detection_zero_all
        )

        if (identical(dgp, "taxonwise")) {
            present_effect <- mean(
                expected_log_probability(
                    eta_one_all, template$sigma[j]
                ) -
                    expected_log_probability(
                        eta_zero_all, template$sigma[j]
                    )
            )
            detected_reference <- expected_detected_log_probability(
                eta = eta_zero_all,
                sigma = template$sigma[j],
                depth = all_depth,
                rho = rho_zero_all
            )
            detected_comparison <- expected_detected_log_probability(
                eta = eta_one_all,
                sigma = template$sigma[j],
                depth = all_depth,
                rho = rho_one_all
            )
            detected_effect <- detected_comparison - detected_reference
        } else {
            present_reference <- reference & !structural_absence[, j] &
                latent_probability[, j] > 0
            present_comparison <- comparison & !structural_absence[, j] &
                latent_probability[, j] > 0
            detected_reference_index <- reference & observed_detection[, j]
            detected_comparison_index <- comparison & observed_detection[, j]
            present_effect <- if (
                any(present_reference) && any(present_comparison)
            ) {
                mean(log(latent_probability[present_comparison, j])) -
                    mean(log(latent_probability[present_reference, j]))
            } else {
                NA_real_
            }
            detected_effect <- if (
                any(detected_reference_index) &&
                any(detected_comparison_index)
            ) {
                mean(log(latent_probability[detected_comparison_index, j])) -
                    mean(log(latent_probability[detected_reference_index, j]))
            } else {
                NA_real_
            }
            expected_prevalence_reference <- mean(
                observed_detection[reference, j]
            )
            expected_prevalence_comparison <- mean(
                observed_detection[comparison, j]
            )
            standardized_prevalence_difference <-
                expected_prevalence_comparison -
                expected_prevalence_reference
        }

        truth_observed <- effects$signal[j] &&
            setting$effect_parameter > 0 &&
            setting$scenario != "structural_matched_prevalence"

        rows[[j]] <- data.frame(
            setting_id = setting$setting_id,
            base_setting_id = setting$base_setting_id,
            design_id = setting$design_id,
            n_per_group = setting$n_per_group,
            total_sample_size = setting$total_sample_size,
            signal_fraction = setting$signal_fraction,
            n_signal_target = setting$n_signal_target,
            confounding = setting$confounding,
            confounded = setting$confounded,
            sample_size_label = setting$sample_size_label,
            signal_fraction_label = setting$signal_fraction_label,
            confounding_label = setting$confounding_label,
            design_panel = setting$design_panel,
            template_role = setting$template_role,
            group_z_correlation = metadata$group_z_correlation[1L],
            z_mean_difference = metadata$z_mean_difference[1L],
            study = setting$study,
            dgp = setting$dgp,
            scenario = setting$scenario,
            effect_level = setting$effect_level,
            effect_index = setting$effect_index,
            effect_parameter = setting$effect_parameter,
            effect_measure = setting$effect_measure,
            taxon = template$taxon[j],
            taxon_index = j,
            evaluation_taxon = TRUE,
            scenario_signal = effects$signal[j],
            truth_structural = abs(effects$delta[j]) > 1e-12,
            truth_abundance = abs(effects$zeta[j]) > 1e-12 ||
                abs(effects$absolute_effect[j]) > 1e-12,
            truth_observed_prevalence = truth_observed,
            structural_log_odds_effect = effects$delta[j],
            conditional_logit_abundance_effect = effects$zeta[j],
            direct_log_absolute_abundance_effect =
                effects$absolute_effect[j],
            standardized_structural_probability_difference = mean(
                rho_one_all - rho_zero_all
            ),
            standardized_expected_prevalence_difference =
                standardized_prevalence_difference,
            expected_prevalence_reference =
                expected_prevalence_reference,
            expected_prevalence_comparison =
                expected_prevalence_comparison,
            expected_prevalence_difference =
                expected_prevalence_comparison -
                expected_prevalence_reference,
            oracle_present_conditional_log_relative_effect =
                present_effect,
            oracle_detected_log_relative_effect = detected_effect,
            actual_structural_absence_reference = mean(
                structural_absence[reference, j]
            ),
            actual_structural_absence_comparison = mean(
                structural_absence[comparison, j]
            ),
            actual_observed_prevalence_reference = mean(
                observed_detection[reference, j]
            ),
            actual_observed_prevalence_comparison = mean(
                observed_detection[comparison, j]
            ),
            actual_nondetection_given_presence_reference = mean(
                present_but_undetected[reference, j]
            ) / max(
                mean(!structural_absence[reference, j]),
                .Machine$double.eps
            ),
            actual_nondetection_given_presence_comparison = mean(
                present_but_undetected[comparison, j]
            ) / max(
                mean(!structural_absence[comparison, j]),
                .Machine$double.eps
            ),
            baseline_relative_abundance =
                template$baseline_probability[j],
            baseline_structural_absence =
                template$baseline_rho[j],
            sigma = template$sigma[j],
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, rows)
}

simulate_taxonwise_setting <- function(
        replication_id, setting, design, template) {
    metadata <- if (setting$scenario == "depth_imbalanced_null") {
        design$imbalanced_metadata
    } else {
        design$balanced_metadata
    }
    depth <- metadata$reads
    n <- nrow(metadata)
    J <- nrow(template)
    effects <- make_effect_vectors(setting, metadata, template)

    eta <- outer(rep(1, n), template$baseline_eta) +
        outer(metadata$z, template$eta_z) +
        outer(metadata$group_num, effects$zeta)
    rho_linear <- outer(rep(1, n), logit(template$baseline_rho)) +
        outer(metadata$z, template$rho_z) +
        outer(metadata$group_num, effects$delta)
    rho <- expit(rho_linear)

    set.seed(
        CONFIG$base_seed + replication_id * 100003L +
            scenario_seed_offset(setting$scenario)
    )
    structural_uniform <- matrix(
        stats::runif(n * J), nrow = n, ncol = J
    )
    latent_normal <- matrix(
        stats::rnorm(n * J), nrow = n, ncol = J
    )
    structural_absence <- structural_uniform < rho
    latent_x <- eta + sweep(
        latent_normal, 2L, template$sigma, "*"
    )
    latent_probability <- expit(latent_x)
    latent_probability[structural_absence] <- 0

    total_probability <- rowSums(latent_probability)
    if (any(!is.finite(total_probability)) ||
        any(total_probability >= CONFIG$probability_guard)) {
        stopf(
            paste(
                "Generated taxon probabilities exceeded the guard in %s;",
                "maximum sum was %.6f."
            ),
            setting$setting_id, max(total_probability)
        )
    }
    other_probability <- 1 - total_probability
    counts_full <- draw_multinomial_counts(
        probability = latent_probability,
        depth = depth,
        other_probability = other_probability,
        seed = CONFIG$base_seed + replication_id * 100003L +
            setting$setting_index * 1009L + 97L
    )
    rownames(counts_full) <- c(template$taxon, "Other_unmodeled")
    colnames(counts_full) <- rownames(metadata)

    biological_counts <- t(
        counts_full[template$taxon, , drop = FALSE]
    )
    observed_detection <- biological_counts > 0
    present_but_undetected <- !structural_absence &
        !observed_detection

    truth <- make_truth_rows(
        setting = setting,
        template = template,
        metadata = metadata,
        effects = effects,
        rho = rho,
        eta = eta,
        structural_absence = structural_absence,
        latent_probability = latent_probability,
        observed_detection = observed_detection,
        present_but_undetected = present_but_undetected,
        depth = depth,
        dgp = "taxonwise"
    )
    rownames(truth) <- NULL

    list(
        setting = setting,
        counts = counts_full,
        biological_counts = biological_counts,
        metadata = metadata,
        depth = depth,
        structural_absence = structural_absence,
        latent_x = latent_x,
        latent_probability = latent_probability,
        observed_detection = observed_detection,
        present_but_undetected = present_but_undetected,
        truth = truth,
        evaluation_taxa = template$taxon,
        other_taxon = "Other_unmodeled",
        dgp = "taxonwise"
    )
}

simulate_joint_setting <- function(
        replication_id, setting, design, template) {
    metadata <- design$balanced_metadata
    depth <- metadata$reads
    n <- nrow(metadata)
    J <- nrow(template)
    effects <- make_effect_vectors(setting, metadata, template)

    rho_linear <- outer(rep(1, n), logit(template$baseline_rho)) +
        outer(metadata$z, template$rho_z) +
        outer(metadata$group_num, effects$delta)
    rho <- expit(rho_linear)

    set.seed(
        CONFIG$base_seed + replication_id * 100003L +
            scenario_seed_offset(setting$scenario)
    )
    structural_uniform <- matrix(
        stats::runif(n * J), nrow = n, ncol = J
    )
    structural_absence <- structural_uniform < rho

    index <- seq_len(J)
    correlation <- CONFIG$joint_correlation^abs(
        outer(index, index, "-")
    )
    chol_correlation <- chol(correlation)
    independent_normal <- matrix(
        stats::rnorm(n * J), nrow = n, ncol = J
    )
    correlated_normal <- independent_normal %*% chol_correlation
    correlated_normal <- sweep(
        correlated_normal, 2L, template$joint_sd, "*"
    )

    log_absolute <- outer(
        rep(1, n), log(template$baseline_probability)
    ) +
        outer(metadata$z, template$eta_z) +
        outer(metadata$group_num, effects$absolute_effect) +
        correlated_normal
    absolute_weight <- exp(log_absolute)
    absolute_weight[structural_absence] <- 0

    other_weight <- rep(CONFIG$joint_other_weight, n)
    denominator <- rowSums(absolute_weight) + other_weight
    latent_probability <- absolute_weight / denominator
    other_probability <- other_weight / denominator
    latent_x <- logit(clamp(
        latent_probability, 1e-12, 1 - 1e-12
    ))
    latent_x[structural_absence] <- NA_real_

    counts_full <- draw_multinomial_counts(
        probability = latent_probability,
        depth = depth,
        other_probability = other_probability,
        seed = CONFIG$base_seed + replication_id * 100003L +
            setting$setting_index * 1009L + 197L
    )
    rownames(counts_full) <- c(template$taxon, "Other_unmodeled")
    colnames(counts_full) <- rownames(metadata)
    biological_counts <- t(
        counts_full[template$taxon, , drop = FALSE]
    )
    observed_detection <- biological_counts > 0
    present_but_undetected <- !structural_absence &
        !observed_detection

    eta_working <- outer(
        rep(1, n), template$baseline_eta
    ) + outer(metadata$z, template$eta_z)

    truth <- make_truth_rows(
        setting = setting,
        template = template,
        metadata = metadata,
        effects = effects,
        rho = rho,
        eta = eta_working,
        structural_absence = structural_absence,
        latent_probability = latent_probability,
        observed_detection = observed_detection,
        present_but_undetected = present_but_undetected,
        depth = depth,
        dgp = "joint_lognormal"
    )
    rownames(truth) <- NULL

    list(
        setting = setting,
        counts = counts_full,
        biological_counts = biological_counts,
        metadata = metadata,
        depth = depth,
        structural_absence = structural_absence,
        latent_x = latent_x,
        latent_probability = latent_probability,
        observed_detection = observed_detection,
        present_but_undetected = present_but_undetected,
        truth = truth,
        evaluation_taxa = template$taxon,
        other_taxon = "Other_unmodeled",
        dgp = "joint_lognormal"
    )
}

simulate_setting <- function(replication_id, setting, design, template) {
    if (identical(setting$dgp, "taxonwise")) {
        simulate_taxonwise_setting(
            replication_id, setting, design, template
        )
    } else if (identical(setting$dgp, "joint_lognormal")) {
        simulate_joint_setting(
            replication_id, setting, design, template
        )
    } else {
        stopf("Unknown data-generating mechanism: %s", setting$dgp)
    }
}

method_result_template <- function(taxa, method, component, status,
                                   elapsed = NA_real_) {
    data.frame(
        taxon = taxa,
        method = method,
        component = component,
        p_value = NA_real_,
        p_value_unfiltered = NA_real_,
        package_q_value = NA_real_,
        estimate = NA_real_,
        standard_error = NA_real_,
        statistic = NA_real_,
        available = FALSE,
        status = status,
        warning = NA_character_,
        method_detail = NA_character_,
        abundance_induced_flag = NA,
        passed_sensitivity = NA,
        runtime_seconds = elapsed,
        runtime_seconds_total_setting = NA_real_,
        p_tau_025 = NA_real_,
        p_tau_050 = NA_real_,
        p_tau_075 = NA_real_,
        stringsAsFactors = FALSE
    )
}

valid_p <- function(x) is.finite(x) & x >= 0 & x <= 1

run_dasra <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture(
        DASRA::dasra(
            counts = simulation$counts[taxa, , drop = FALSE],
            metadata = simulation$metadata,
            formula = ~ group + z,
            group = "group",
            library_size = "reads",
            taxa_are_rows = TRUE,
            reference = "control",
            p_adjust_method = "BH",
            component = "all",
            full_output = FALSE
        )
    )
    if (is_sim_error(captured$value)) {
        status <- paste0("error: ", captured$value$message)
        return(rbind(
            method_result_template(
                taxa, "DASRA", "structural_absence", status, captured$elapsed
            ),
            method_result_template(
                taxa, "DASRA", "present_conditional_abundance", status,
                captured$elapsed
            )
        ))
    }

    fit <- captured$value
    result <- fit$results
    diagnostics <- fit$diagnostics
    result <- result[match(taxa, result$taxon), , drop = FALSE]
    diagnostics <- diagnostics[match(taxa, diagnostics$taxon), , drop = FALSE]

    structural <- method_result_template(
        taxa, "DASRA", "structural_absence", "unavailable", captured$elapsed
    )
    structural$p_value <- suppressWarnings(as.numeric(
        result$p_structural_absence
    ))
    structural$package_q_value <- suppressWarnings(as.numeric(
        result$p_adj_structural_absence
    ))
    structural$statistic <- suppressWarnings(as.numeric(
        result$z_structural_absence
    ))
    structural$available <- valid_p(structural$p_value) &
        as.logical(diagnostics$formed_structural_absence)
    structural$status <- ifelse(
        structural$available &
            as.logical(diagnostics$regular_structural_absence),
        "ok",
        as.character(diagnostics$reason_structural_absence)
    )

    abundance <- method_result_template(
        taxa, "DASRA", "present_conditional_abundance", "unavailable",
        captured$elapsed
    )
    abundance$p_value <- suppressWarnings(as.numeric(
        result$p_relative_abundance
    ))
    abundance$package_q_value <- suppressWarnings(as.numeric(
        result$p_adj_relative_abundance
    ))
    abundance$estimate <- suppressWarnings(as.numeric(
        result$estimate_relative_abundance
    ))
    abundance$standard_error <- suppressWarnings(as.numeric(
        result$se_relative_abundance
    ))
    abundance$statistic <- suppressWarnings(as.numeric(
        result$z_relative_abundance
    ))
    abundance$available <- valid_p(abundance$p_value) &
        as.logical(diagnostics$formed_relative_abundance)
    abundance$status <- ifelse(
        abundance$available,
        "ok",
        as.character(diagnostics$reason_relative_abundance)
    )

    structural$warning <- as.character(
        diagnostics$warning_structural_absence
    )
    abundance$warning <- as.character(
        diagnostics$warning_relative_abundance
    )
    rbind(structural, abundance)
}

run_zinq <- function(simulation, replication_id) {
    taxa <- simulation$evaluation_taxa
    count <- simulation$counts[taxa, , drop = FALSE]
    relative_abundance <- sweep(count, 2L, simulation$depth, "/")
    prevalence_rows <- vector("list", length(taxa))
    abundance_rows <- vector("list", length(taxa))
    tau_weights <- pmin(CONFIG$zinq_taus, 1 - CONFIG$zinq_taus)
    tau_weights <- tau_weights / sum(tau_weights)

    started_all <- proc.time()[["elapsed"]]
    for (j in seq_along(taxa)) {
        data <- data.frame(
            y = as.numeric(relative_abundance[j, ]),
            group_num = simulation$metadata$group_num,
            z = simulation$metadata$z,
            log_depth = simulation$metadata$log_depth
        )
        set.seed(
            CONFIG$base_seed + replication_id * 100003L +
                simulation$setting$setting_index * 1009L + j
        )
        captured <- safe_capture(
            ZINQ::ZINQ_tests(
                formula.logistic = y ~ group_num + z + log_depth,
                formula.quantile = y ~ group_num + z + log_depth,
                C = "group_num",
                y_CorD = "C",
                data = data,
                taus = CONFIG$zinq_taus,
                seed = CONFIG$base_seed + replication_id * 1009L + j
            )
        )

        prevalence <- method_result_template(
            taxa[j], "ZINQ", "observed_prevalence", "unavailable",
            captured$elapsed
        )
        abundance <- method_result_template(
            taxa[j], "ZINQ", "detected_quantile_abundance", "unavailable",
            captured$elapsed
        )
        warning_text <- collapse_text(captured$warnings)
        prevalence$warning <- warning_text
        abundance$warning <- warning_text

        if (is_sim_error(captured$value)) {
            status <- paste0("error: ", captured$value$message)
            prevalence$status <- status
            abundance$status <- status
        } else {
            fit <- captured$value
            prevalence$p_value <- suppressWarnings(as.numeric(
                fit$pvalue.logistic
            )[1L])
            prevalence$available <- valid_p(prevalence$p_value)
            prevalence$status <- ifelse(
                prevalence$available, "ok", "nonfinite_prevalence_p_value"
            )

            quantile_p <- suppressWarnings(as.numeric(fit$pvalue.quantile))
            if (length(quantile_p) == length(CONFIG$zinq_taus)) {
                abundance$p_tau_025 <- quantile_p[1L]
                abundance$p_tau_050 <- quantile_p[2L]
                abundance$p_tau_075 <- quantile_p[3L]
                if (all(valid_p(quantile_p))) {
                    abundance$p_value <- cauchy_combine(
                        quantile_p, tau_weights
                    )
                }
            }
            abundance$available <- valid_p(abundance$p_value)
            abundance$status <- ifelse(
                abundance$available, "ok", "quantile_component_unavailable"
            )
        }
        prevalence_rows[[j]] <- prevalence
        abundance_rows[[j]] <- abundance
    }
    elapsed_all <- proc.time()[["elapsed"]] - started_all
    prevalence <- do.call(rbind, prevalence_rows)
    abundance <- do.call(rbind, abundance_rows)
    prevalence$runtime_seconds_total_setting <- elapsed_all
    abundance$runtime_seconds_total_setting <- elapsed_all
    rbind(prevalence, abundance)
}

maaslin_result_frame <- function(fit, component, output_dir) {
    candidate <- if (component == "abundance") {
        fit$fit_data_abundance %||% fit$abundance
    } else {
        fit$fit_data_prevalence %||% fit$prevalence
    }
    if (is.data.frame(candidate)) return(candidate)
    if (is.list(candidate) && is.data.frame(candidate$results)) {
        return(candidate$results)
    }
    if (is.list(fit) && is.data.frame(fit$results)) return(fit$results)

    files <- list.files(
        output_dir,
        pattern = "(all_results|results).*\\.(tsv|txt|csv)$",
        recursive = TRUE, full.names = TRUE, ignore.case = TRUE
    )
    for (file in files) {
        data <- tryCatch(
            if (grepl("\\.csv$", file, ignore.case = TRUE)) {
                utils::read.csv(file, check.names = FALSE)
            } else {
                utils::read.delim(file, check.names = FALSE)
            },
            error = function(e) NULL
        )
        if (is.data.frame(data) && nrow(data)) return(data)
    }
    NULL
}

select_maaslin_row <- function(data, taxon, component) {
    if (is.null(data) || !is.data.frame(data) || !nrow(data)) return(NULL)
    feature_column <- first_existing_column(
        data, c("feature", "taxon", "name_feature")
    )
    candidates <- data
    if (!is.null(feature_column)) {
        candidates <- candidates[
            as.character(candidates[[feature_column]]) == taxon,
            , drop = FALSE
        ]
    }
    if (!nrow(candidates)) return(NULL)

    if ("model" %in% names(candidates)) {
        model <- tolower(trimws(as.character(candidates$model)))
        accepted <- if (component == "abundance") {
            c("linear", "abundance")
        } else {
            c("logistic", "prevalence")
        }
        candidates <- candidates[model %in% accepted, , drop = FALSE]
    }
    if (!nrow(candidates)) return(NULL)

    if ("metadata" %in% names(candidates)) {
        keep <- as.character(candidates$metadata) == "group"
        if (any(keep, na.rm = TRUE)) candidates <- candidates[keep, , drop = FALSE]
    }
    if (nrow(candidates) > 1L && "value" %in% names(candidates)) {
        keep <- tolower(as.character(candidates$value)) == "case"
        if (any(keep, na.rm = TRUE)) candidates <- candidates[keep, , drop = FALSE]
    }
    if (nrow(candidates) > 1L && "name" %in% names(candidates)) {
        normalized <- gsub("[^a-z0-9]", "", tolower(as.character(candidates$name)))
        keep <- grepl("group.*case|case.*group", normalized)
        if (any(keep, na.rm = TRUE)) candidates <- candidates[keep, , drop = FALSE]
    }
    if (!nrow(candidates)) NULL else candidates[1L, , drop = FALSE]
}

extract_maaslin_component <- function(data, taxa, component, elapsed, warnings) {
    component_label <- if (component == "abundance") {
        "detected_log_abundance"
    } else {
        "observed_prevalence"
    }
    output <- method_result_template(
        taxa, "MaAsLin 3", component_label, "row_unavailable", elapsed
    )
    output$warning <- collapse_text(warnings)
    for (j in seq_along(taxa)) {
        row <- select_maaslin_row(data, taxa[j], component)
        if (is.null(row)) next
        output$p_value[j] <- as_numeric_column(
            row, c("pval_individual", "pval", "p_value", "pvalue")
        )[1L]
        output$package_q_value[j] <- as_numeric_column(
            row, c("qval_individual", "qval", "q_value", "qvalue")
        )[1L]
        output$estimate[j] <- as_numeric_column(
            row, c("coef", "coefficient", "estimate")
        )[1L]
        output$standard_error[j] <- as_numeric_column(
            row, c("stderr", "std_error", "se")
        )[1L]
        output$statistic[j] <- as_numeric_column(
            row, c("stat", "test_statistic", "z")
        )[1L]
        output$available[j] <- valid_p(output$p_value[j])
        error_column <- first_existing_column(row, c("error", "warning"))
        error_text <- if (is.null(error_column)) NA_character_ else {
            collapse_text(row[[error_column]])
        }
        output$method_detail[j] <- error_text
        if (component == "prevalence") {
            output$abundance_induced_flag[j] <- !is.na(error_text) &&
                grepl("abundance|induced", error_text, ignore.case = TRUE)
        }
        output$status[j] <- if (output$available[j]) {
            "ok"
        } else if (!is.na(error_text)) {
            error_text
        } else {
            "nonfinite_p_value"
        }
    }
    output
}

run_maaslin3 <- function(simulation, replication_id, temporary_root) {
    taxa <- simulation$evaluation_taxa
    output_dir <- file.path(
        temporary_root,
        paste0("maaslin3_", simulation$setting$setting_id)
    )
    unlink(output_dir, recursive = TRUE, force = TRUE)
    safe_dir_create(output_dir)
    on.exit({
        try(maaslin3::maaslin_log_reset(), silent = TRUE)
        unlink(output_dir, recursive = TRUE, force = TRUE)
    }, add = TRUE)

    input_data <- as.data.frame(t(simulation$counts))
    metadata <- simulation$metadata[, c("group", "z", "log_depth"), drop = FALSE]
    arguments <- list(
        input_data = input_data,
        input_metadata = metadata,
        output = output_dir,
        formula = "~ group + z + log_depth",
        min_abundance = 0,
        min_prevalence = 0,
        max_prevalence = 1.01,
        zero_threshold = 0,
        min_variance = 0,
        normalization = "TSS",
        transform = "LOG",
        correction = "BH",
        standardize = TRUE,
        median_comparison_abundance = TRUE,
        median_comparison_prevalence = FALSE,
        warn_prevalence = isTRUE(CONFIG$maaslin3_warn_prevalence),
        augment = TRUE,
        evaluate_only = NULL,
        plot_summary_plot = FALSE,
        plot_associations = FALSE,
        max_pngs = 1,
        cores = 1,
        save_models = FALSE,
        save_plots_rds = FALSE,
        verbosity = "ERROR",
        reference = "group,control"
    )
    available <- names(formals(maaslin3::maaslin3))
    arguments <- arguments[names(arguments) %in% available]

    set.seed(
        CONFIG$base_seed + replication_id * 100003L +
            simulation$setting$setting_index * 4001L
    )
    try(maaslin3::maaslin_log_reset(), silent = TRUE)
    captured <- safe_capture({
        invisible(capture.output(
            fit_value <- suppressMessages(do.call(maaslin3::maaslin3, arguments))
        ))
        fit_value
    })
    if (is_sim_error(captured$value)) {
        status <- paste0("error: ", captured$value$message)
        return(rbind(
            method_result_template(
                taxa, "MaAsLin 3", "observed_prevalence", status,
                captured$elapsed
            ),
            method_result_template(
                taxa, "MaAsLin 3", "detected_log_abundance", status,
                captured$elapsed
            )
        ))
    }

    abundance_data <- maaslin_result_frame(
        captured$value, "abundance", output_dir
    )
    prevalence_data <- maaslin_result_frame(
        captured$value, "prevalence", output_dir
    )
    rbind(
        extract_maaslin_component(
            prevalence_data, taxa, "prevalence", captured$elapsed,
            captured$warnings
        ),
        extract_maaslin_component(
            abundance_data, taxa, "abundance", captured$elapsed,
            captured$warnings
        )
    )
}

run_edger <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        design <- stats::model.matrix(
            ~ group + z, data = simulation$metadata
        )
        coefficient <- grep("^group", colnames(design))
        if (length(coefficient) != 1L) {
            stop("Could not identify the edgeR group coefficient.")
        }
        dge <- edgeR::DGEList(counts = simulation$counts)
        keep <- edgeR::filterByExpr(dge, design = design)
        if (!any(keep)) stop("edgeR filterByExpr retained no taxa.")
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        dge <- edgeR::estimateDisp(dge, design = design, robust = TRUE)
        fit <- edgeR::glmQLFit(dge, design = design, robust = TRUE)
        test <- edgeR::glmQLFTest(fit, coef = coefficient)
        table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
        list(table = table, keep = keep)
    })
    output <- method_result_template(
        taxa, "edgeR", "general_abundance", "filtered_or_unavailable",
        captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }
    table <- captured$value$table
    matched <- match(taxa, rownames(table))
    present <- !is.na(matched)
    output$p_value[present] <- suppressWarnings(as.numeric(
        table$PValue[matched[present]]
    ))
    output$package_q_value[present] <- suppressWarnings(as.numeric(
        table$FDR[matched[present]]
    ))
    output$estimate[present] <- suppressWarnings(as.numeric(
        table$logFC[matched[present]]
    ))
    output$statistic[present] <- suppressWarnings(as.numeric(
        table$F[matched[present]]
    ))
    output$available <- valid_p(output$p_value)
    output$status <- ifelse(output$available, "ok", "filterByExpr_excluded")
    output
}

run_deseq2 <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        col_data <- simulation$metadata[, c("group", "z"), drop = FALSE]
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = round(simulation$counts),
            colData = col_data,
            design = ~ z + group
        )
        dds <- DESeq2::DESeq(
            dds,
            test = "Wald",
            fitType = "parametric",
            sfType = "poscounts",
            betaPrior = FALSE,
            quiet = TRUE
        )
        coefficient <- grep(
            "group.*case.*control|group_case_vs_control",
            DESeq2::resultsNames(dds), value = TRUE, ignore.case = TRUE
        )
        if (length(coefficient) != 1L) {
            stop("Could not identify the DESeq2 group coefficient.")
        }
        result <- DESeq2::results(
            dds,
            name = coefficient,
            alpha = CONFIG$alpha,
            independentFiltering = TRUE,
            cooksCutoff = TRUE
        )
        as.data.frame(result)
    })
    output <- method_result_template(
        taxa, "DESeq2", "general_abundance", "unavailable", captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }
    table <- captured$value
    matched <- match(taxa, rownames(table))
    present <- !is.na(matched)
    output$p_value[present] <- suppressWarnings(as.numeric(
        table$pvalue[matched[present]]
    ))
    output$package_q_value[present] <- suppressWarnings(as.numeric(
        table$padj[matched[present]]
    ))
    output$estimate[present] <- suppressWarnings(as.numeric(
        table$log2FoldChange[matched[present]]
    ))
    output$standard_error[present] <- suppressWarnings(as.numeric(
        table$lfcSE[matched[present]]
    ))
    output$statistic[present] <- suppressWarnings(as.numeric(
        table$stat[matched[present]]
    ))
    output$available <- valid_p(output$p_value)
    output$status <- ifelse(
        output$available, "ok", "independent_filter_or_nonfinite_p"
    )
    output
}

find_ancom_column <- function(names_vector, prefix, term_pattern) {
    normalized <- tolower(names_vector)
    prefix_pattern <- paste0("^", tolower(prefix), "[_\\.]")
    candidates <- which(
        grepl(prefix_pattern, normalized) & grepl(term_pattern, normalized)
    )
    if (!length(candidates)) return(NA_integer_)
    candidates[[1L]]
}

run_ancombc2 <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        arguments <- list(
            data = simulation$counts,
            taxa_are_rows = TRUE,
            assay_name = "counts",
            tax_level = NULL,
            aggregate_data = NULL,
            meta_data = simulation$metadata[, c("group", "z"), drop = FALSE],
            fix_formula = "group + z",
            rand_formula = NULL,
            p_adj_method = "BH",
            pseudo = 0,
            pseudo_sens = TRUE,
            prv_cut = 0,
            lib_cut = 0,
            s0_perc = 0.05,
            group = "group",
            struc_zero = TRUE,
            neg_lb = TRUE,
            alpha = CONFIG$alpha,
            n_cl = 1,
            verbose = FALSE,
            global = FALSE,
            pairwise = FALSE,
            dunnet = FALSE,
            trend = FALSE
        )
        available_arguments <- names(formals(ANCOMBC::ancombc2))
        arguments <- arguments[names(arguments) %in% available_arguments]
        fit <- suppressMessages(do.call(ANCOMBC::ancombc2, arguments))
        result <- fit$res
        if (!is.data.frame(result)) result <- as.data.frame(result)
        result
    })
    output <- method_result_template(
        taxa, "ANCOM-BC2", "general_abundance", "unavailable",
        captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }

    table <- captured$value
    taxon_column <- first_existing_column(
        table, c("taxon", "feature", "otu")
    )
    table_taxa <- if (is.null(taxon_column)) {
        rownames(table)
    } else {
        as.character(table[[taxon_column]])
    }
    matched <- match(taxa, table_taxa)
    present <- !is.na(matched)
    term_pattern <- "group.*case|case.*group"
    p_col <- find_ancom_column(names(table), "p", term_pattern)
    q_col <- find_ancom_column(names(table), "q", term_pattern)
    lfc_col <- find_ancom_column(names(table), "lfc", term_pattern)
    se_col <- find_ancom_column(names(table), "se", term_pattern)
    w_col <- find_ancom_column(names(table), "w", term_pattern)
    sensitivity_col <- find_ancom_column(
        names(table), "passed_ss", term_pattern
    )

    if (is.na(p_col)) {
        output$status <- "group_p_value_column_unavailable"
        return(output)
    }

    raw_p <- rep(NA_real_, length(taxa))
    raw_p[present] <- suppressWarnings(as.numeric(
        table[matched[present], p_col]
    ))
    output$p_value_unfiltered <- raw_p
    output$p_value <- raw_p

    if (!is.na(sensitivity_col)) {
        passed <- rep(NA, length(taxa))
        passed[present] <- as.logical(
            table[matched[present], sensitivity_col]
        )
        output$passed_sensitivity <- passed
        failed <- !is.na(passed) & !passed & valid_p(raw_p)
        output$p_value[failed] <- 1
        output$method_detail[failed] <-
            "pseudocount_sensitivity_failed; p-value set to one"
    } else {
        output$method_detail[present] <-
            "pseudocount_sensitivity_column_not_returned"
    }

    if (!is.na(q_col)) {
        output$package_q_value[present] <- suppressWarnings(as.numeric(
            table[matched[present], q_col]
        ))
    }
    if (!is.na(lfc_col)) {
        output$estimate[present] <- suppressWarnings(as.numeric(
            table[matched[present], lfc_col]
        ))
    }
    if (!is.na(se_col)) {
        output$standard_error[present] <- suppressWarnings(as.numeric(
            table[matched[present], se_col]
        ))
    }
    if (!is.na(w_col)) {
        output$statistic[present] <- suppressWarnings(as.numeric(
            table[matched[present], w_col]
        ))
    }

    output$available <- valid_p(output$p_value)
    output$status <- ifelse(
        output$available,
        ifelse(
            !is.na(output$passed_sensitivity) &
                !output$passed_sensitivity,
            "ok_after_sensitivity_exclusion",
            "ok"
        ),
        "nonfinite_p_value"
    )
    output
}
run_linda <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    linda_formals <- names(formals(MicrobiomeStat::linda))

    arguments <- if ("feature.dat" %in% linda_formals) {
        list(
            feature.dat = simulation$counts,
            meta.dat = simulation$metadata[, c("group", "z"), drop = FALSE],
            formula = "~ group + z",
            feature.dat.type = "count",
            prev.filter = 0,
            mean.abund.filter = 0,
            max.abund.filter = 0,
            adaptive = TRUE,
            p.adj.method = "BH",
            alpha = CONFIG$alpha,
            n.cores = 1,
            verbose = FALSE
        )
    } else {
        list(
            otu.tab = simulation$counts,
            meta = simulation$metadata[, c("group", "z"), drop = FALSE],
            formula = "~ group + z",
            type = "count",
            adaptive = TRUE,
            p.adj.method = "BH",
            alpha = CONFIG$alpha,
            prev.cut = 0,
            lib.cut = 0,
            n.cores = 1
        )
    }
    arguments <- arguments[names(arguments) %in% linda_formals]

    captured <- safe_capture(
        suppressMessages(do.call(MicrobiomeStat::linda, arguments))
    )
    output <- method_result_template(
        taxa, "LinDA", "general_abundance", "unavailable",
        captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }

    fit <- captured$value
    if (is.null(fit$output) || !is.list(fit$output)) {
        output$status <- "output_list_unavailable"
        return(output)
    }
    output_names <- names(fit$output)
    coefficient_name <- output_names[
        grepl("group.*case|case.*group", output_names, ignore.case = TRUE)
    ]
    if (length(coefficient_name) > 1L) {
        normalized <- gsub(
            "[^a-z0-9]", "", tolower(coefficient_name)
        )
        exact <- which(normalized %in% c("groupcase", "groupcasevscontrol"))
        if (length(exact)) coefficient_name <- coefficient_name[exact[1L]]
    }
    if (length(coefficient_name) != 1L) {
        output$status <- "group_output_unavailable"
        return(output)
    }

    table <- as.data.frame(fit$output[[coefficient_name]])
    matched <- match(taxa, rownames(table))
    present <- !is.na(matched)
    selected <- table[matched[present], , drop = FALSE]
    output$p_value[present] <- as_numeric_column(
        selected, c("pvalue", "p_value", "p.value")
    )
    output$package_q_value[present] <- as_numeric_column(
        selected, c("padj", "qvalue", "q_value", "adj.p")
    )
    output$estimate[present] <- as_numeric_column(
        selected, c("log2FoldChange", "coef", "estimate")
    )
    output$standard_error[present] <- as_numeric_column(
        selected, c("lfcSE", "se", "stderr", "standard_error")
    )
    output$statistic[present] <- as_numeric_column(
        selected, c("stat", "t", "z")
    )
    output$available <- valid_p(output$p_value)
    output$status <- ifelse(
        output$available, "ok", "nonfinite_p_value"
    )
    output
}
run_corncob <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    rows <- vector("list", length(taxa))
    started <- proc.time()[["elapsed"]]
    for (j in seq_along(taxa)) {
        data <- data.frame(
            W = as.numeric(simulation$counts[taxa[j], ]),
            M = as.numeric(simulation$depth),
            group = simulation$metadata$group,
            z = simulation$metadata$z
        )
        captured <- safe_capture(
            corncob::bbdml(
                formula = cbind(W, M - W) ~ group + z,
                phi.formula = ~ group + z,
                data = data,
                method = "trust",
                robust = isTRUE(CONFIG$corncob_robust)
            )
        )
        row <- method_result_template(
            taxa[j], "corncob", "general_abundance", "unavailable",
            captured$elapsed
        )
        row$warning <- collapse_text(captured$warnings)
        if (is_sim_error(captured$value)) {
            row$status <- paste0("error: ", captured$value$message)
            rows[[j]] <- row
            next
        }
        fit <- captured$value
        summary_fit <- tryCatch(summary(fit), error = function(e) NULL)
        coefficient_table <- if (is.null(summary_fit)) NULL else summary_fit$coefficients
        if (is.null(coefficient_table) || !is.matrix(coefficient_table)) {
            row$status <- "coefficient_table_unavailable"
            rows[[j]] <- row
            next
        }
        abundance_table <- coefficient_table[
            seq_len(min(fit$np.mu, nrow(coefficient_table))), , drop = FALSE
        ]
        coefficient_row <- grep(
            "group.*case|case.*group", rownames(abundance_table),
            ignore.case = TRUE
        )
        if (length(coefficient_row) != 1L) {
            row$status <- "group_coefficient_unavailable"
            rows[[j]] <- row
            next
        }
        columns <- colnames(abundance_table)
        estimate_col <- grep("estimate", columns, ignore.case = TRUE)[1L]
        se_col <- grep("std.*error|standard.*error", columns, ignore.case = TRUE)[1L]
        stat_col <- grep("t value|z value|stat", columns, ignore.case = TRUE)[1L]
        p_col <- grep("pr\\(|p.value|p value", columns, ignore.case = TRUE)[1L]
        if (is.na(p_col)) {
            row$status <- "group_p_value_unavailable"
            rows[[j]] <- row
            next
        }
        row$p_value <- suppressWarnings(as.numeric(
            abundance_table[coefficient_row, p_col]
        ))
        if (!is.na(estimate_col)) {
            row$estimate <- suppressWarnings(as.numeric(
                abundance_table[coefficient_row, estimate_col]
            ))
        }
        if (!is.na(se_col)) {
            row$standard_error <- suppressWarnings(as.numeric(
                abundance_table[coefficient_row, se_col]
            ))
        }
        if (!is.na(stat_col)) {
            row$statistic <- suppressWarnings(as.numeric(
                abundance_table[coefficient_row, stat_col]
            ))
        }
        row$available <- valid_p(row$p_value)
        row$status <- ifelse(row$available, "ok", "nonfinite_p_value")
        rows[[j]] <- row
    }
    output <- do.call(rbind, rows)
    output$runtime_seconds_total_setting <- proc.time()[["elapsed"]] - started
    output
}

run_metagenomeseq <- function(simulation) {
    taxa <- simulation$evaluation_taxa
    captured <- safe_capture({
        pheno <- Biobase::AnnotatedDataFrame(
            simulation$metadata[, c("group", "z"), drop = FALSE]
        )
        object <- metagenomeSeq::newMRexperiment(
            counts = simulation$counts,
            phenoData = pheno
        )
        percentile <- metagenomeSeq::cumNormStatFast(
            object, pFlag = FALSE
        )
        object <- metagenomeSeq::cumNorm(object, p = percentile)
        design <- stats::model.matrix(
            ~ group + z, data = simulation$metadata
        )
        coefficient <- grep("^group", colnames(design))
        if (length(coefficient) != 1L) {
            stop("Could not identify the metagenomeSeq group coefficient.")
        }
        coefficient_name <- colnames(design)[coefficient]
        control <- metagenomeSeq::zigControl(
            maxit = CONFIG$metagenomeseq_maxit,
            verbose = FALSE
        )
        fit <- metagenomeSeq::fitZig(
            obj = object, mod = design, control = control
        )

        table <- tryCatch(
            metagenomeSeq::MRcoefs(
                fit,
                by = coefficient_name,
                coef = coefficient_name,
                number = Inf,
                group = 4,
                adjustMethod = "BH",
                eff = 0,
                counts = 0
            ),
            error = function(e) NULL
        )
        if (is.null(table)) {
            fit_slot <- methods::slot(fit, "fit")
            table <- limma::topTable(
                fit_slot,
                coef = coefficient,
                number = Inf,
                sort.by = "none",
                adjust.method = "BH"
            )
        }
        list(
            table = as.data.frame(table),
            coefficient_name = coefficient_name
        )
    })
    output <- method_result_template(
        taxa, "metagenomeSeq", "general_abundance", "unavailable",
        captured$elapsed
    )
    output$warning <- collapse_text(captured$warnings)
    if (is_sim_error(captured$value)) {
        output$status <- paste0("error: ", captured$value$message)
        return(output)
    }

    table <- captured$value$table
    coefficient_name <- captured$value$coefficient_name
    taxon_column <- first_existing_column(
        table, c("taxon", "feature", "otu")
    )
    table_taxa <- if (is.null(taxon_column)) {
        rownames(table)
    } else {
        as.character(table[[taxon_column]])
    }
    matched <- match(taxa, table_taxa)
    present <- !is.na(matched)
    selected <- table[matched[present], , drop = FALSE]

    output$p_value[present] <- as_numeric_column(
        selected, c("pvalues", "P.Value", "pvalue", "p_value")
    )
    output$package_q_value[present] <- as_numeric_column(
        selected, c("adjPvalues", "adj.P.Val", "padj", "q_value")
    )
    estimate_candidates <- c(
        "logFC", "coef", "coefficient", "estimate"
    )
    if (!is.null(first_existing_column(selected, coefficient_name))) {
        estimate_candidates <- c(coefficient_name, estimate_candidates)
    }
    output$estimate[present] <- as_numeric_column(
        selected, estimate_candidates
    )
    output$standard_error[present] <- as_numeric_column(
        selected, c("se", "SE", "stderr", "standard_error")
    )
    output$statistic[present] <- as_numeric_column(
        selected, c("t", "stat", "z")
    )
    output$available <- valid_p(output$p_value)
    output$status <- ifelse(
        output$available, "ok", "nonfinite_p_value"
    )
    output
}
run_generic_abundance_methods <- function(simulation) {
    rbind(
        run_edger(simulation),
        run_deseq2(simulation),
        run_ancombc2(simulation),
        run_linda(simulation),
        run_corncob(simulation),
        run_metagenomeseq(simulation)
    )
}

apply_common_multiplicity <- function(results) {
    data <- data.table::as.data.table(results)
    data[, p_value_for_multiplicity := ifelse(
        valid_p(p_value), p_value, 1
    )]
    data[, q_value_bh := stats::p.adjust(
        p_value_for_multiplicity, method = "BH"
    ), by = .(method, component)]
    data[, q_value_by := stats::p.adjust(
        p_value_for_multiplicity, method = "BY"
    ), by = .(method, component)]
    data[, q_value := q_value_bh]
    data[, reject_raw := available & valid_p(p_value) &
        p_value < CONFIG$alpha]
    data[, reject_bh := q_value_bh < CONFIG$alpha]
    data[, reject_by := q_value_by < CONFIG$alpha]
    as.data.frame(data)
}
analyze_setting <- function(simulation, replication_id, temporary_root) {
    setting <- simulation$setting
    method_seed <- CONFIG$base_seed + replication_id * 100003L +
        setting$setting_index * 10007L
    set.seed(method_seed + 1L)
    dasra_result <- run_dasra(simulation)
    set.seed(method_seed + 2L)
    zinq_result <- run_zinq(simulation, replication_id)
    set.seed(method_seed + 3L)
    maaslin_result <- run_maaslin3(
        simulation, replication_id, temporary_root
    )
    method_results <- list(dasra_result, zinq_result, maaslin_result)
    if (isTRUE(setting$generic_abundance)) {
        set.seed(method_seed + 101L)
        method_results[[length(method_results) + 1L]] <-
            run_generic_abundance_methods(simulation)
    }
    result <- do.call(rbind, method_results)
    rownames(result) <- NULL
    result <- apply_common_multiplicity(result)
    result$replication <- replication_id
    context_columns <- c(
        "setting_id", "base_setting_id", "design_id", "n_per_group",
        "total_sample_size", "signal_fraction", "n_signal_target",
        "confounding", "confounded", "sample_size_label",
        "signal_fraction_label", "confounding_label", "design_panel",
        "template_role", "study", "dgp", "scenario", "effect_level",
        "effect_index", "effect_parameter", "effect_measure"
    )
    for (column in context_columns) {
        result[[column]] <- setting[[column]][1L]
    }

    truth <- simulation$truth
    result <- merge(
        result,
        truth,
        by = c(context_columns, "taxon"),
        all.x = TRUE,
        sort = FALSE
    )
    result$truth_for_component <- ifelse(
        result$component == "structural_absence",
        result$truth_structural,
        ifelse(
            result$component == "observed_prevalence",
            result$truth_observed_prevalence,
            result$truth_abundance
        )
    )
    result
}

setting_checkpoint_path <- function(replication_id, setting_id) {
    file.path(
        CONFIG$results_dir,
        "setting_checkpoints",
        sprintf("replicate_%04d", replication_id),
        paste0(setting_id, ".rds")
    )
}

run_one_setting <- function(replication_id, setting, design, template,
                            temporary_root) {
    checkpoint <- setting_checkpoint_path(replication_id, setting$setting_id)
    if (file.exists(checkpoint) && !isTRUE(CONFIG$overwrite)) {
        existing <- tryCatch(readRDS(checkpoint), error = function(e) NULL)
        if (!is.null(existing) && isTRUE(existing$success) &&
            identical(existing$script_version, CONFIG$script_version)) {
            return(existing)
        }
    }

    result <- tryCatch({
        simulation <- simulate_setting(
            replication_id = replication_id,
            setting = setting,
            design = design,
            template = template
        )
        method_results <- analyze_setting(
            simulation = simulation,
            replication_id = replication_id,
            temporary_root = temporary_root
        )
        dataset <- if (isTRUE(CONFIG$save_datasets)) {
            list(
                setting = simulation$setting,
                counts = simulation$counts,
                metadata = simulation$metadata,
                structural_absence = simulation$structural_absence,
                latent_x = simulation$latent_x,
                latent_probability = simulation$latent_probability,
                observed_detection = simulation$observed_detection,
                present_but_undetected = simulation$present_but_undetected,
                truth = simulation$truth
            )
        } else {
            NULL
        }
        list(
            success = TRUE,
            script_version = CONFIG$script_version,
            replication = replication_id,
            setting_id = setting$setting_id,
            method_results = method_results,
            truth = simulation$truth,
            dataset = dataset,
            error = NA_character_
        )
    }, error = function(e) {
        list(
            success = FALSE,
            script_version = CONFIG$script_version,
            replication = replication_id,
            setting_id = setting$setting_id,
            method_results = NULL,
            truth = NULL,
            dataset = NULL,
            error = conditionMessage(e)
        )
    })
    atomic_save_rds(result, checkpoint, compress = "gzip")
    result
}

run_replication <- function(replication_id) {
    check_packages(replicate_packages)
    for (directory in c(
        CONFIG$root, CONFIG$results_dir, CONFIG$datasets_dir,
        CONFIG$summary_dir, CONFIG$tmp_dir, CONFIG$logs_dir
    )) {
        safe_dir_create(directory)
    }

    final_path <- file.path(
        CONFIG$results_dir, sprintf("replicate_%04d.rds", replication_id)
    )
    if (file.exists(final_path) && !isTRUE(CONFIG$overwrite)) {
        existing <- tryCatch(readRDS(final_path), error = function(e) NULL)
        if (!is.null(existing) && isTRUE(existing$success) &&
            identical(existing$script_version, CONFIG$script_version)) {
            messagef("Replication %d already completed: %s", replication_id, final_path)
            return(invisible(existing))
        }
    }

    messagef(
        paste(
            "Starting formal DASRA replication %d with %d settings",
            "(%d base settings x %d design strata)."
        ),
        replication_id, nrow(SETTINGS),
        CONFIG$n_base_settings, CONFIG$n_design_strata
    )
    started <- proc.time()[["elapsed"]]

    design_specs <- unique(SETTINGS[, c(
        "design_id", "n_per_group", "confounding"
    ), drop = FALSE])
    designs <- setNames(
        lapply(seq_len(nrow(design_specs)), function(i) {
            make_replication_design(
                replication_id = replication_id,
                n_per_group = design_specs$n_per_group[i],
                confounding = design_specs$confounding[i]
            )
        }),
        design_specs$design_id
    )

    template_specs <- unique(SETTINGS[, c(
        "template_id", "n_signal_target", "template_role"
    ), drop = FALSE])
    templates <- setNames(
        lapply(seq_len(nrow(template_specs)), function(i) {
            make_taxon_template(
                n_signal = template_specs$n_signal_target[i],
                role = template_specs$template_role[i]
            )
        }),
        template_specs$template_id
    )

    temporary_root <- file.path(
        CONFIG$tmp_dir, sprintf("replicate_%04d", replication_id)
    )
    unlink(temporary_root, recursive = TRUE, force = TRUE)
    safe_dir_create(temporary_root)
    on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)

    setting_results <- vector("list", nrow(SETTINGS))
    for (i in seq_len(nrow(SETTINGS))) {
        setting <- SETTINGS[i, , drop = FALSE]
        messagef(
            "Replication %d: setting %d/%d, %s",
            replication_id, i, nrow(SETTINGS), setting$setting_id
        )
        setting_results[[i]] <- run_one_setting(
            replication_id = replication_id,
            setting = setting,
            design = designs[[setting$design_id]],
            template = templates[[setting$template_id]],
            temporary_root = temporary_root
        )
    }

    successful <- vapply(setting_results, function(x) isTRUE(x$success), logical(1))
    method_results <- data.table::rbindlist(
        lapply(setting_results[successful], `[[`, "method_results"),
        fill = TRUE, use.names = TRUE
    )
    truth <- data.table::rbindlist(
        lapply(setting_results[successful], `[[`, "truth"),
        fill = TRUE, use.names = TRUE
    )
    errors <- data.frame(
        setting_id = vapply(setting_results, `[[`, character(1), "setting_id"),
        success = successful,
        error = vapply(
            setting_results,
            function(x) if (is.null(x$error)) NA_character_ else x$error,
            character(1)
        ),
        stringsAsFactors = FALSE
    )
    errors <- merge(
        errors,
        SETTINGS[, c(
            "setting_id", "base_setting_id", "design_id", "n_per_group",
            "signal_fraction", "n_signal_target", "confounding",
            "scenario", "effect_parameter"
        )],
        by = "setting_id", all.x = TRUE, sort = FALSE
    )

    dataset_path <- NA_character_
    if (isTRUE(CONFIG$save_datasets)) {
        dataset_path <- file.path(
            CONFIG$datasets_dir, sprintf("replicate_%04d.rds", replication_id)
        )
        datasets <- setNames(
            lapply(setting_results[successful], `[[`, "dataset"),
            vapply(
                setting_results[successful], `[[`, character(1), "setting_id"
            )
        )
        atomic_save_rds(
            list(
                script_version = CONFIG$script_version,
                replication = replication_id,
                designs = designs,
                taxon_templates = templates,
                settings = SETTINGS,
                datasets = datasets
            ),
            dataset_path,
            compress = "gzip"
        )
    }

    elapsed <- proc.time()[["elapsed"]] - started
    final <- list(
        success = all(successful),
        script_version = CONFIG$script_version,
        replication = replication_id,
        elapsed_seconds = elapsed,
        configuration = CONFIG,
        settings = SETTINGS,
        design_strata = DESIGN_STRATA,
        package_versions = package_versions(),
        method_configuration = method_configuration(),
        session_info = utils::capture.output(sessionInfo()),
        method_results = as.data.frame(method_results),
        truth = as.data.frame(truth),
        setting_status = errors,
        dataset_path = dataset_path
    )
    atomic_save_rds(final, final_path, compress = "gzip")

    messagef(
        "Replication %d completed in %.1f minutes; %d/%d settings succeeded.",
        replication_id, elapsed / 60, sum(successful), length(successful)
    )
    if (!all(successful)) {
        failed <- errors$setting_id[!successful]
        stopf(
            "Replication %d finished with failed settings: %s",
            replication_id, paste(failed, collapse = ", ")
        )
    }
    invisible(final)
}

method_display_name <- function(method, component) {
    key <- paste(method, component, sep = "::")
    labels <- c(
        "DASRA::structural_absence" = "DASRA structural absence",
        "DASRA::present_conditional_abundance" =
            "DASRA present-conditional abundance",
        "ZINQ::observed_prevalence" = "ZINQ Firth prevalence",
        "ZINQ::detected_quantile_abundance" = "ZINQ quantile abundance",
        "MaAsLin 3::observed_prevalence" = "MaAsLin 3 prevalence",
        "MaAsLin 3::detected_log_abundance" = "MaAsLin 3 abundance",
        "edgeR::general_abundance" = "edgeR",
        "DESeq2::general_abundance" = "DESeq2",
        "ANCOM-BC2::general_abundance" = "ANCOM-BC2",
        "LinDA::general_abundance" = "LinDA",
        "corncob::general_abundance" = "corncob",
        "metagenomeSeq::general_abundance" = "metagenomeSeq"
    )
    answer <- unname(labels[key])
    missing <- is.na(answer)
    answer[missing] <- paste(method[missing], component[missing], sep = ": ")
    answer
}

abundance_component <- function(component) {
    component %in% c(
        "present_conditional_abundance",
        "detected_quantile_abundance",
        "detected_log_abundance",
        "general_abundance"
    )
}

structural_component <- function(component) {
    component %in% c("structural_absence", "observed_prevalence")
}

make_replicate_metrics <- function(results) {
    data <- data.table::as.data.table(results)
    data[, method_label := method_display_name(method, component)]
    data[, .(
        n_taxa = .N,
        n_available = sum(available, na.rm = TRUE),
        availability = mean(available, na.rm = TRUE),
        n_signal = sum(truth_for_component, na.rm = TRUE),
        n_null = sum(!truth_for_component, na.rm = TRUE),
        signal_availability = if (sum(truth_for_component) > 0) {
            mean(available[truth_for_component], na.rm = TRUE)
        } else {
            NA_real_
        },
        null_availability = if (sum(!truth_for_component) > 0) {
            mean(available[!truth_for_component], na.rm = TRUE)
        } else {
            NA_real_
        },
        raw_rejections_signal = sum(
            truth_for_component & reject_raw, na.rm = TRUE
        ),
        raw_rejections_null = sum(
            !truth_for_component & reject_raw, na.rm = TRUE
        ),
        raw_signal_rejection_rate = if (sum(truth_for_component) > 0) {
            sum(truth_for_component & reject_raw, na.rm = TRUE) /
                sum(truth_for_component)
        } else {
            NA_real_
        },
        type1_error = if (sum(!truth_for_component) > 0) {
            sum(!truth_for_component & reject_raw, na.rm = TRUE) /
                sum(!truth_for_component)
        } else {
            NA_real_
        },
        discoveries = sum(reject_bh, na.rm = TRUE),
        true_discoveries = sum(
            truth_for_component & reject_bh, na.rm = TRUE
        ),
        false_discoveries = sum(
            !truth_for_component & reject_bh, na.rm = TRUE
        ),
        power = if (sum(truth_for_component) > 0) {
            sum(truth_for_component & reject_bh, na.rm = TRUE) /
                sum(truth_for_component)
        } else {
            NA_real_
        },
        fdp = if (sum(reject_bh, na.rm = TRUE) > 0) {
            sum(!truth_for_component & reject_bh, na.rm = TRUE) /
                sum(reject_bh, na.rm = TRUE)
        } else {
            0
        },
        discoveries_by = sum(reject_by, na.rm = TRUE),
        power_by = if (sum(truth_for_component) > 0) {
            sum(truth_for_component & reject_by, na.rm = TRUE) /
                sum(truth_for_component)
        } else {
            NA_real_
        },
        fdp_by = if (sum(reject_by, na.rm = TRUE) > 0) {
            sum(!truth_for_component & reject_by, na.rm = TRUE) /
                sum(reject_by, na.rm = TRUE)
        } else {
            0
        },
        designated_signal_raw_rejection = if (any(scenario_signal)) {
            mean(reject_raw[scenario_signal], na.rm = TRUE)
        } else {
            NA_real_
        },
        designated_signal_bh_rejection = if (any(scenario_signal)) {
            mean(reject_bh[scenario_signal], na.rm = TRUE)
        } else {
            NA_real_
        },
        designated_signal_by_rejection = if (any(scenario_signal)) {
            mean(reject_by[scenario_signal], na.rm = TRUE)
        } else {
            NA_real_
        }
    ), by = .(
        replication, setting_id, base_setting_id, design_id,
        n_per_group, total_sample_size, signal_fraction, n_signal_target,
        confounding, confounded, sample_size_label,
        signal_fraction_label, confounding_label, design_panel,
        template_role, study, dgp, scenario,
        effect_level, effect_index, effect_parameter, effect_measure,
        method, component, method_label
    )]
}

aggregate_metric <- function(metrics, value_columns) {
    data <- data.table::as.data.table(metrics)
    groups <- c(
        "setting_id", "base_setting_id", "design_id",
        "n_per_group", "total_sample_size", "signal_fraction",
        "n_signal_target", "confounding", "confounded",
        "sample_size_label", "signal_fraction_label",
        "confounding_label", "design_panel", "template_role",
        "study", "dgp", "scenario", "effect_level", "effect_index",
        "effect_parameter", "effect_measure", "method", "component",
        "method_label"
    )
    long <- data.table::melt(
        data,
        id.vars = c("replication", groups),
        measure.vars = value_columns,
        variable.name = "metric",
        value.name = "value"
    )
    long[, .(
        n_replications = sum(is.finite(value)),
        mean = if (any(is.finite(value))) {
            mean(value, na.rm = TRUE)
        } else {
            NA_real_
        },
        sd = if (sum(is.finite(value)) > 1L) {
            stats::sd(value, na.rm = TRUE)
        } else {
            NA_real_
        },
        mcse = if (sum(is.finite(value)) > 1L) {
            stats::sd(value, na.rm = TRUE) /
                sqrt(sum(is.finite(value)))
        } else {
            NA_real_
        }
    ), by = c(groups, "metric")][
        , `:=`(
            ci_lower = pmax(0, mean - 1.96 * mcse),
            ci_upper = pmin(1, mean + 1.96 * mcse)
        )
    ]
}

METHOD_COLORS <- c(
    "DASRA structural absence" = "#0072B2",
    "DASRA present-conditional abundance" = "#0072B2",
    "ZINQ Firth prevalence" = "#D55E00",
    "ZINQ quantile abundance" = "#D55E00",
    "MaAsLin 3 prevalence" = "#009E73",
    "MaAsLin 3 abundance" = "#009E73",
    "edgeR" = "#CC79A7",
    "DESeq2" = "#E69F00",
    "ANCOM-BC2" = "#56B4E9",
    "LinDA" = "#7A7A00",
    "corncob" = "#000000",
    "metagenomeSeq" = "#7F7F7F"
)

METHOD_SHAPES <- c(
    "DASRA structural absence" = 16,
    "DASRA present-conditional abundance" = 16,
    "ZINQ Firth prevalence" = 17,
    "ZINQ quantile abundance" = 17,
    "MaAsLin 3 prevalence" = 15,
    "MaAsLin 3 abundance" = 15,
    "edgeR" = 18,
    "DESeq2" = 0,
    "ANCOM-BC2" = 1,
    "LinDA" = 2,
    "corncob" = 5,
    "metagenomeSeq" = 6
)

METHOD_LINETYPES <- c(
    "DASRA structural absence" = "solid",
    "DASRA present-conditional abundance" = "solid",
    "ZINQ Firth prevalence" = "longdash",
    "ZINQ quantile abundance" = "longdash",
    "MaAsLin 3 prevalence" = "dotdash",
    "MaAsLin 3 abundance" = "dotdash",
    "edgeR" = "twodash",
    "DESeq2" = "dashed",
    "ANCOM-BC2" = "dotted",
    "LinDA" = "solid",
    "corncob" = "dashed",
    "metagenomeSeq" = "dotted"
)

method_scales <- function() {
    list(
        ggplot2::scale_color_manual(
            values = METHOD_COLORS, drop = TRUE
        ),
        ggplot2::scale_shape_manual(
            values = METHOD_SHAPES, drop = TRUE
        ),
        ggplot2::scale_linetype_manual(
            values = METHOD_LINETYPES, drop = TRUE
        )
    )
}

publication_theme <- function() {
    ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
            legend.position = "top",
            legend.title = ggplot2::element_blank(),
            legend.key.width = grid::unit(1.6, "lines"),
            panel.spacing = grid::unit(1.2, "lines"),
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 11),
            axis.text.x = ggplot2::element_text(
                angle = 0, hjust = 0.5, vjust = 0.5
            ),
            plot.margin = ggplot2::margin(7, 8, 7, 8)
        )
}

save_publication_plot <- function(plot, stem, width = 8, height = 5) {
    pdf_path <- file.path(CONFIG$summary_dir, paste0(stem, ".pdf"))
    png_path <- file.path(CONFIG$summary_dir, paste0(stem, ".png"))
    pdf_device <- if (capabilities("cairo")) {
        grDevices::cairo_pdf
    } else {
        "pdf"
    }
    ggplot2::ggsave(
        pdf_path, plot = plot, width = width, height = height,
        units = "in", device = pdf_device
    )
    ggplot2::ggsave(
        png_path, plot = plot, width = width, height = height,
        units = "in", dpi = 600, bg = "white"
    )
}

summarize_curve <- function(data, value, groups) {
    value_name <- deparse(substitute(value))
    data[, .(
        mean = mean(get(value_name), na.rm = TRUE),
        mcse = stats::sd(get(value_name), na.rm = TRUE) /
            sqrt(sum(is.finite(get(value_name))))
    ), by = groups]
}

is_primary_sample_signal <- function(data) {
    data$n_per_group == CONFIG$primary_n_per_group &
        abs(data$signal_fraction - CONFIG$primary_signal_fraction) < 1e-12
}

primary_design_subtitle <- function() {
    sprintf(
        "%d samples per group; %d%% signal taxa",
        CONFIG$primary_n_per_group,
        as.integer(round(100 * CONFIG$primary_signal_fraction))
    )
}

order_design_panel <- function(data) {
    levels <- DESIGN_STRATA[order(
        DESIGN_STRATA$n_per_group,
        DESIGN_STRATA$signal_fraction,
        DESIGN_STRATA$confounded
    ), "design_panel"]
    factor(data$design_panel, levels = unique(levels))
}

summarize_rejection_curve <- function(data, grouping) {
    per_rep <- data[, .(
        rejection_rate = mean(reject_raw, na.rm = TRUE)
    ), by = c("replication", grouping)]
    per_rep[, .(
        mean = mean(rejection_rate, na.rm = TRUE),
        mcse = stats::sd(rejection_rate, na.rm = TRUE) /
            sqrt(sum(is.finite(rejection_rate)))
    ), by = grouping]
}

plot_structural_estimand <- function(results) {
    data <- data.table::as.data.table(results)
    data <- data[
        scenario %in% c(
            "observed_prevalence_only",
            "structural_matched_prevalence"
        ) & scenario_signal & structural_component(component)
    ]
    data[, method_label := method_display_name(method, component)]
    data[, scenario_label := ifelse(
        scenario == "observed_prevalence_only",
        paste(
            "Observed-prevalence difference",
            "with no structural-absence difference",
            sep = "\n"
        ),
        paste(
            "Structural-absence difference",
            "with matched observed prevalence",
            sep = "\n"
        )
    )]

    primary <- data[is_primary_sample_signal(data)]
    summary <- summarize_rejection_curve(
        primary,
        c("confounding_label", "scenario_label", "effect_parameter",
          "method_label")
    )
    plot <- ggplot2::ggplot(
        summary,
        ggplot2::aes(
            x = effect_parameter, y = mean,
            group = method_label, color = method_label,
            linetype = method_label, shape = method_label
        )
    ) +
        ggplot2::geom_hline(
            yintercept = CONFIG$alpha,
            linetype = "dotted", linewidth = 0.5
        ) +
        ggplot2::geom_line(linewidth = 0.75) +
        ggplot2::geom_point(size = 2.25) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = pmax(0, mean - 1.96 * mcse),
                ymax = pmin(1, mean + 1.96 * mcse)
            ),
            width = 0.008, linewidth = 0.42
        ) +
        ggplot2::facet_grid(
            confounding_label ~ scenario_label, scales = "free_x"
        ) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::scale_x_continuous(
            breaks = scales::breaks_pretty(n = 6)
        ) +
        ggplot2::labs(
            x = "Target standardized probability difference",
            y = "Empirical rejection rate",
            subtitle = primary_design_subtitle()
        ) +
        method_scales() +
        publication_theme()
    save_publication_plot(
        plot, "fig_structural_estimand_separation",
        width = 10.8, height = 7.8
    )

    sensitivity <- summarize_rejection_curve(
        data,
        c("design_panel", "scenario_label", "effect_parameter",
          "method_label")
    )
    sensitivity[, design_panel := order_design_panel(sensitivity)]
    sensitivity_plot <- ggplot2::ggplot(
        sensitivity,
        ggplot2::aes(
            x = effect_parameter, y = mean,
            group = method_label, color = method_label,
            linetype = method_label, shape = method_label
        )
    ) +
        ggplot2::geom_hline(
            yintercept = CONFIG$alpha,
            linetype = "dotted", linewidth = 0.4
        ) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::geom_point(size = 1.6) +
        ggplot2::facet_grid(
            design_panel ~ scenario_label, scales = "free_x"
        ) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(
            x = "Target standardized probability difference",
            y = "Empirical rejection rate"
        ) +
        method_scales() +
        publication_theme() +
        ggplot2::theme(strip.text.y = ggplot2::element_text(size = 8))
    save_publication_plot(
        sensitivity_plot, "fig_structural_estimand_separation_sensitivity",
        width = 12.5, height = 20.5
    )
}

plot_structural_design <- function(truth) {
    data <- data.table::as.data.table(truth)
    data <- data[
        scenario %in% c(
            "observed_prevalence_only",
            "structural_matched_prevalence"
        ) & scenario_signal
    ]
    data[, scenario_label := ifelse(
        scenario == "observed_prevalence_only",
        paste(
            "Observed-prevalence difference",
            "with no structural-absence difference",
            sep = "\n"
        ),
        paste(
            "Structural-absence difference",
            "with matched observed prevalence",
            sep = "\n"
        )
    )]
    per_rep <- data[, .(
        achieved_observed_prevalence_difference = mean(
            abs(standardized_expected_prevalence_difference), na.rm = TRUE
        ),
        achieved_structural_probability_difference = mean(
            abs(standardized_structural_probability_difference), na.rm = TRUE
        )
    ), by = .(
        replication, design_panel, n_per_group, signal_fraction,
        confounding_label, scenario_label, effect_parameter
    )]
    long <- data.table::melt(
        per_rep,
        id.vars = c(
            "replication", "design_panel", "n_per_group", "signal_fraction",
            "confounding_label", "scenario_label", "effect_parameter"
        ),
        measure.vars = c(
            "achieved_observed_prevalence_difference",
            "achieved_structural_probability_difference"
        ),
        variable.name = "quantity",
        value.name = "value"
    )
    long[, quantity_label := ifelse(
        quantity == "achieved_observed_prevalence_difference",
        "Standardized observed-prevalence difference",
        "Structural-absence probability difference"
    )]
    primary <- long[
        n_per_group == CONFIG$primary_n_per_group &
            abs(signal_fraction - CONFIG$primary_signal_fraction) < 1e-12
    ]
    summary <- primary[, .(
        mean = mean(value, na.rm = TRUE),
        mcse = stats::sd(value, na.rm = TRUE) /
            sqrt(sum(is.finite(value)))
    ), by = .(
        confounding_label, scenario_label, effect_parameter, quantity_label
    )]
    plot <- ggplot2::ggplot(
        summary,
        ggplot2::aes(
            x = effect_parameter, y = mean,
            group = quantity_label,
            linetype = quantity_label, shape = quantity_label
        )
    ) +
        ggplot2::geom_abline(
            slope = 1, intercept = 0,
            linetype = "dotted", linewidth = 0.45
        ) +
        ggplot2::geom_line(linewidth = 0.75) +
        ggplot2::geom_point(size = 2.25) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = pmax(0, mean - 1.96 * mcse),
                ymax = mean + 1.96 * mcse
            ),
            width = 0.008, linewidth = 0.42
        ) +
        ggplot2::facet_grid(
            confounding_label ~ scenario_label, scales = "free"
        ) +
        ggplot2::labs(
            x = "Target standardized probability difference",
            y = "Achieved absolute probability difference",
            subtitle = primary_design_subtitle()
        ) +
        publication_theme()
    save_publication_plot(
        plot, "fig_structural_design_verification",
        width = 10.8, height = 7.4
    )
}

plot_component_specificity <- function(results) {
    data <- data.table::as.data.table(results)
    structural_under_abundance <- data[
        scenario == "present_conditional_abundance" &
            scenario_signal & method == "DASRA" &
            component == "structural_absence"
    ]
    abundance_under_structural <- data[
        scenario == "structural_only" &
            scenario_signal & method == "DASRA" &
            component == "present_conditional_abundance"
    ]
    structural_under_abundance[, panel := paste(
        "Structural component under",
        "present-conditional abundance effects", sep = "\n"
    )]
    abundance_under_structural[, panel := paste(
        "Abundance component under",
        "structural-absence effects", sep = "\n"
    )]
    combined <- data.table::rbindlist(
        list(structural_under_abundance, abundance_under_structural),
        fill = TRUE
    )
    primary <- combined[is_primary_sample_signal(combined)]
    summary <- summarize_rejection_curve(
        primary,
        c("confounding_label", "panel", "effect_parameter")
    )
    upper <- max(0.15, summary$mean + 2 * summary$mcse, na.rm = TRUE)
    plot <- ggplot2::ggplot(
        summary,
        ggplot2::aes(x = effect_parameter, y = mean, group = 1)
    ) +
        ggplot2::geom_hline(
            yintercept = CONFIG$alpha,
            linetype = "dotted", linewidth = 0.5
        ) +
        ggplot2::geom_line(
            linewidth = 0.75,
            color = METHOD_COLORS["DASRA present-conditional abundance"]
        ) +
        ggplot2::geom_point(
            size = 2.25, shape = 16,
            color = METHOD_COLORS["DASRA present-conditional abundance"]
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = pmax(0, mean - 1.96 * mcse),
                ymax = pmin(1, mean + 1.96 * mcse)
            ),
            width = 0.015, linewidth = 0.42,
            color = METHOD_COLORS["DASRA present-conditional abundance"]
        ) +
        ggplot2::facet_grid(
            confounding_label ~ panel, scales = "free_x"
        ) +
        ggplot2::coord_cartesian(ylim = c(0, upper)) +
        ggplot2::labs(
            x = "Target effect size",
            y = "Non-target component Type I error (raw p < 0.05)",
            subtitle = primary_design_subtitle()
        ) +
        publication_theme()
    save_publication_plot(
        plot, "fig_component_specificity",
        width = 10.4, height = 7.0
    )

    sensitivity <- summarize_rejection_curve(
        combined,
        c("design_panel", "panel", "effect_parameter")
    )
    sensitivity[, design_panel := order_design_panel(sensitivity)]
    sensitivity_plot <- ggplot2::ggplot(
        sensitivity,
        ggplot2::aes(x = effect_parameter, y = mean, group = 1)
    ) +
        ggplot2::geom_hline(
            yintercept = CONFIG$alpha,
            linetype = "dotted", linewidth = 0.4
        ) +
        ggplot2::geom_line(
            linewidth = 0.58,
            color = METHOD_COLORS["DASRA present-conditional abundance"]
        ) +
        ggplot2::geom_point(
            size = 1.5,
            color = METHOD_COLORS["DASRA present-conditional abundance"]
        ) +
        ggplot2::facet_grid(design_panel ~ panel, scales = "free_x") +
        ggplot2::coord_cartesian(ylim = c(0, max(0.15, upper))) +
        ggplot2::labs(
            x = "Target effect size",
            y = "Non-target component Type I error (raw p < 0.05)"
        ) +
        publication_theme() +
        ggplot2::theme(strip.text.y = ggplot2::element_text(size = 8))
    save_publication_plot(
        sensitivity_plot, "fig_component_specificity_sensitivity",
        width = 12.0, height = 18.5
    )
}

summarize_metric_curve <- function(data, metric_name, grouping) {
    data[, .(
        mean = mean(get(metric_name), na.rm = TRUE),
        mcse = stats::sd(get(metric_name), na.rm = TRUE) /
            sqrt(sum(is.finite(get(metric_name))))
    ), by = grouping]
}

plot_abundance_metric <- function(metrics, metric_name, y_label, stem,
                                  reference_line = FALSE) {
    data <- data.table::as.data.table(metrics)
    data <- data[
        scenario == "present_conditional_abundance" &
            effect_parameter > 0 & abundance_component(component)
    ]
    primary <- data[is_primary_sample_signal(data)]
    summary <- summarize_metric_curve(
        primary, metric_name,
        c("confounding_label", "effect_parameter", "method_label")
    )
    plot <- ggplot2::ggplot(
        summary,
        ggplot2::aes(
            x = effect_parameter, y = mean,
            group = method_label, color = method_label,
            linetype = method_label, shape = method_label
        )
    )
    if (isTRUE(reference_line)) {
        plot <- plot + ggplot2::geom_hline(
            yintercept = CONFIG$alpha,
            linetype = "dotted", linewidth = 0.5
        )
    }
    plot <- plot +
        ggplot2::geom_line(linewidth = 0.68) +
        ggplot2::geom_point(size = 2.0) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = pmax(0, mean - 1.96 * mcse),
                ymax = pmin(1, mean + 1.96 * mcse)
            ),
            width = 0.015, linewidth = 0.38
        ) +
        ggplot2::facet_wrap(~ confounding_label, nrow = 1) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::scale_x_continuous(
            breaks = scales::breaks_pretty(n = 8)
        ) +
        ggplot2::labs(
            x = paste(
                "Absolute present-conditional mean",
                "log-relative-abundance difference"
            ),
            y = y_label,
            subtitle = primary_design_subtitle()
        ) +
        method_scales() + publication_theme() +
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
            shape = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
            linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
        )
    save_publication_plot(plot, stem, width = 12.0, height = 7.1)

    sensitivity <- summarize_metric_curve(
        data, metric_name,
        c("design_panel", "effect_parameter", "method_label")
    )
    sensitivity[, design_panel := order_design_panel(sensitivity)]
    sensitivity_plot <- ggplot2::ggplot(
        sensitivity,
        ggplot2::aes(
            x = effect_parameter, y = mean,
            group = method_label, color = method_label,
            linetype = method_label, shape = method_label
        )
    )
    if (isTRUE(reference_line)) {
        sensitivity_plot <- sensitivity_plot + ggplot2::geom_hline(
            yintercept = CONFIG$alpha,
            linetype = "dotted", linewidth = 0.4
        )
    }
    sensitivity_plot <- sensitivity_plot +
        ggplot2::geom_line(linewidth = 0.55) +
        ggplot2::geom_point(size = 1.35) +
        ggplot2::facet_wrap(~ design_panel, ncol = 2) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(
            x = paste(
                "Absolute present-conditional mean",
                "log-relative-abundance difference"
            ),
            y = y_label
        ) +
        method_scales() + publication_theme() +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 8)) +
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
            shape = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
            linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
        )
    save_publication_plot(
        sensitivity_plot, paste0(stem, "_sensitivity"),
        width = 13.0, height = 16.5
    )
}

plot_abundance_power <- function(metrics) {
    plot_abundance_metric(
        metrics, "power", "Power at BH FDR 0.05",
        "fig_abundance_power", FALSE
    )
}

plot_abundance_fdr <- function(metrics) {
    plot_abundance_metric(
        metrics, "fdp", "FDR relative to directly perturbed taxa",
        "fig_abundance_fdr", TRUE
    )
}

plot_availability <- function(metrics) {
    plot_abundance_metric(
        metrics, "signal_availability",
        "Proportion of signal-taxon tests available",
        "fig_method_availability", FALSE
    )
}

plot_global_null_qq <- function(results) {
    data <- data.table::as.data.table(results)
    data <- data[
        scenario == "present_conditional_abundance" &
            effect_parameter == 0 & abundance_component(component) &
            valid_p(p_value) & is_primary_sample_signal(data)
    ]
    data[, method_label := method_display_name(method, component)]
    qq <- data[, {
        observed <- sort(p_value)
        expected <- stats::ppoints(length(observed))
        .(expected = expected, observed = observed)
    }, by = .(confounding_label, method_label)]

    plot <- ggplot2::ggplot(
        qq, ggplot2::aes(x = expected, y = observed)
    ) +
        ggplot2::geom_abline(linetype = "dotted", linewidth = 0.5) +
        ggplot2::geom_point(size = 0.56, alpha = 0.50) +
        ggplot2::facet_grid(confounding_label ~ method_label) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(
            x = "Expected uniform quantile",
            y = "Observed abundance p-value",
            subtitle = primary_design_subtitle()
        ) +
        publication_theme() +
        ggplot2::theme(
            strip.text.x = ggplot2::element_text(size = 8),
            strip.text.y = ggplot2::element_text(size = 9)
        )
    save_publication_plot(
        plot, "fig_global_null_qq",
        width = 17.0, height = 6.6
    )
}

plot_type1_error <- function(metrics) {
    data <- data.table::as.data.table(metrics)
    data <- data[
        (
            scenario == "present_conditional_abundance" &
                effect_parameter == 0
        ) |
            scenario %in% c("depth_imbalanced_null", "joint_global_null")
    ]
    data <- data[
        abs(signal_fraction - CONFIG$primary_signal_fraction) < 1e-12
    ]
    data[, null_scenario_label := factor(
        scenario,
        levels = c(
            "present_conditional_abundance",
            "depth_imbalanced_null",
            "joint_global_null"
        ),
        labels = c(
            "Taxonwise global null",
            "Depth-imbalanced global null",
            "Correlated-community global null"
        )
    )]
    data[, component_family := ifelse(
        structural_component(component),
        "Structural / prevalence tests",
        "Abundance tests"
    )]
    summary <- data[, .(
        mean = mean(type1_error, na.rm = TRUE),
        mcse = stats::sd(type1_error, na.rm = TRUE) /
            sqrt(sum(is.finite(type1_error)))
    ), by = .(
        n_per_group, confounding_label, null_scenario_label,
        component_family, method_label
    )]
    upper <- max(0.15, summary$mean + 2 * summary$mcse, na.rm = TRUE)
    plot <- ggplot2::ggplot(
        summary,
        ggplot2::aes(
            x = factor(n_per_group), y = mean,
            group = method_label, color = method_label,
            linetype = method_label, shape = method_label
        )
    ) +
        ggplot2::geom_hline(
            yintercept = CONFIG$alpha,
            linetype = "dotted", linewidth = 0.5
        ) +
        ggplot2::geom_line(linewidth = 0.65) +
        ggplot2::geom_point(size = 1.9) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = pmax(0, mean - 1.96 * mcse),
                ymax = pmin(1, mean + 1.96 * mcse)
            ),
            width = 0.12, linewidth = 0.38
        ) +
        ggplot2::facet_grid(
            component_family + confounding_label ~ null_scenario_label
        ) +
        ggplot2::coord_cartesian(ylim = c(0, upper)) +
        ggplot2::labs(
            x = "Samples per group",
            y = "Empirical Type I error (raw p < 0.05)",
            subtitle = "Dotted line: nominal alpha = 0.05"
        ) +
        method_scales() + publication_theme() +
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
            shape = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
            linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
        )
    save_publication_plot(
        plot, "fig_type1_error",
        width = 15.0, height = 12.0
    )
}

plot_joint_community_robustness <- function(metrics) {
    data <- data.table::as.data.table(metrics)
    data <- data[is_primary_sample_signal(data)]

    abundance <- data[
        scenario == "joint_abundance" & abundance_component(component)
    ]
    abundance_summary <- abundance[, .(
        power = mean(power, na.rm = TRUE),
        power_mcse = stats::sd(power, na.rm = TRUE) /
            sqrt(sum(is.finite(power))),
        fdr = mean(fdp, na.rm = TRUE),
        fdr_mcse = stats::sd(fdp, na.rm = TRUE) /
            sqrt(sum(is.finite(fdp)))
    ), by = .(confounding_label, effect_parameter, method_label)]

    if (nrow(abundance_summary)) {
        power_plot <- ggplot2::ggplot(
            abundance_summary,
            ggplot2::aes(
                x = effect_parameter, y = power,
                group = method_label, color = method_label,
                linetype = method_label, shape = method_label
            )
        ) +
            ggplot2::geom_line(linewidth = 0.65) +
            ggplot2::geom_point(size = 1.9) +
            ggplot2::geom_errorbar(
                ggplot2::aes(
                    ymin = pmax(0, power - 1.96 * power_mcse),
                    ymax = pmin(1, power + 1.96 * power_mcse)
                ),
                width = 0.02, linewidth = 0.38
            ) +
            ggplot2::facet_wrap(~ confounding_label, nrow = 1) +
            ggplot2::coord_cartesian(ylim = c(0, 1)) +
            ggplot2::labs(
                x = "Direct log-absolute-abundance effect",
                y = "Power at BH FDR 0.05",
                subtitle = primary_design_subtitle()
            ) +
            method_scales() + publication_theme() +
            ggplot2::guides(
                color = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
                shape = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
                linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
            )
        save_publication_plot(
            power_plot, "fig_joint_abundance_power",
            width = 12.0, height = 7.0
        )

        fdr_plot <- ggplot2::ggplot(
            abundance_summary,
            ggplot2::aes(
                x = effect_parameter, y = fdr,
                group = method_label, color = method_label,
                linetype = method_label, shape = method_label
            )
        ) +
            ggplot2::geom_hline(
                yintercept = CONFIG$alpha,
                linetype = "dotted", linewidth = 0.5
            ) +
            ggplot2::geom_line(linewidth = 0.65) +
            ggplot2::geom_point(size = 1.9) +
            ggplot2::geom_errorbar(
                ggplot2::aes(
                    ymin = pmax(0, fdr - 1.96 * fdr_mcse),
                    ymax = pmin(1, fdr + 1.96 * fdr_mcse)
                ),
                width = 0.02, linewidth = 0.38
            ) +
            ggplot2::facet_wrap(~ confounding_label, nrow = 1) +
            ggplot2::labs(
                x = "Direct log-absolute-abundance effect",
                y = "FDR relative to directly perturbed taxa",
                subtitle = primary_design_subtitle()
            ) +
            method_scales() + publication_theme() +
            ggplot2::guides(
                color = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
                shape = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
                linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
            )
        save_publication_plot(
            fdr_plot, "fig_joint_abundance_fdr",
            width = 12.0, height = 7.0
        )
    }

    structural <- data[
        scenario == "joint_structural" & structural_component(component)
    ]
    structural_summary <- structural[, .(
        rejection_rate = mean(
            designated_signal_raw_rejection, na.rm = TRUE
        ),
        mcse = stats::sd(
            designated_signal_raw_rejection, na.rm = TRUE
        ) / sqrt(sum(is.finite(designated_signal_raw_rejection)))
    ), by = .(confounding_label, effect_parameter, method_label)]

    if (nrow(structural_summary)) {
        structural_plot <- ggplot2::ggplot(
            structural_summary,
            ggplot2::aes(
                x = effect_parameter, y = rejection_rate,
                group = method_label, color = method_label,
                linetype = method_label, shape = method_label
            )
        ) +
            ggplot2::geom_hline(
                yintercept = CONFIG$alpha,
                linetype = "dotted", linewidth = 0.5
            ) +
            ggplot2::geom_line(linewidth = 0.7) +
            ggplot2::geom_point(size = 2.1) +
            ggplot2::geom_errorbar(
                ggplot2::aes(
                    ymin = pmax(0, rejection_rate - 1.96 * mcse),
                    ymax = pmin(1, rejection_rate + 1.96 * mcse)
                ),
                width = 0.015, linewidth = 0.4
            ) +
            ggplot2::facet_wrap(~ confounding_label, nrow = 1) +
            ggplot2::coord_cartesian(ylim = c(0, 1)) +
            ggplot2::labs(
                x = "Structural-absence probability difference",
                y = "Empirical rejection rate",
                subtitle = primary_design_subtitle()
            ) +
            method_scales() + publication_theme()
        save_publication_plot(
            structural_plot, "fig_joint_structural_robustness",
            width = 11.0, height = 6.2
        )
    }
}

plot_confounding_design <- function(truth) {
    data <- data.table::as.data.table(truth)
    design_values <- unique(data[, .(
        replication, n_per_group, confounding_label,
        group_z_correlation, z_mean_difference
    )])
    summary <- design_values[, .(
        mean_correlation = mean(group_z_correlation, na.rm = TRUE),
        correlation_mcse = stats::sd(group_z_correlation, na.rm = TRUE) /
            sqrt(sum(is.finite(group_z_correlation))),
        mean_standardized_difference = mean(z_mean_difference, na.rm = TRUE),
        difference_mcse = stats::sd(z_mean_difference, na.rm = TRUE) /
            sqrt(sum(is.finite(z_mean_difference)))
    ), by = .(n_per_group, confounding_label)]

    long <- data.table::rbindlist(list(
        summary[, .(
            n_per_group,
            confounding_label,
            quantity_label =
                "Correlation between group and adjustment covariate",
            mean = mean_correlation,
            mcse = correlation_mcse
        )],
        summary[, .(
            n_per_group,
            confounding_label,
            quantity_label = "Standardized covariate mean difference",
            mean = mean_standardized_difference,
            mcse = difference_mcse
        )]
    ))
    plot <- ggplot2::ggplot(
        long,
        ggplot2::aes(
            x = factor(n_per_group), y = mean,
            group = confounding_label, color = confounding_label,
            shape = confounding_label
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0, linetype = "dotted", linewidth = 0.45
        ) +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::geom_point(size = 2.2) +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = mean - 1.96 * mcse,
                ymax = mean + 1.96 * mcse
            ),
            width = 0.08, linewidth = 0.4
        ) +
        ggplot2::facet_wrap(~ quantity_label, scales = "free_y", nrow = 1) +
        ggplot2::labs(
            x = "Samples per group",
            y = NULL,
            color = NULL,
            shape = NULL
        ) +
        publication_theme()
    save_publication_plot(
        plot, "fig_confounding_design_verification",
        width = 10.5, height = 4.8
    )
    summary
}

summarize_completed_replications <- function() {
    check_packages(summary_packages)
    safe_dir_create(CONFIG$summary_dir)
    previous_outputs <- list.files(
        CONFIG$summary_dir,
        pattern = paste0(
            "^(all_|base_|confounding_|design_|replicate_|setting_|method_|",
            "truth_|simulation_|package_|session_|fig_)"
        ),
        full.names = TRUE
    )
    if (length(previous_outputs)) {
        unlink(previous_outputs, recursive = TRUE, force = TRUE)
    }
    files <- list.files(
        CONFIG$results_dir,
        pattern = "^replicate_[0-9]+\\.rds$",
        full.names = TRUE
    )
    if (!length(files)) {
        stopf("No completed replicate files were found in %s", CONFIG$results_dir)
    }
    raw_results_path <- file.path(
        CONFIG$summary_dir, "all_taxon_method_results.csv.gz"
    )
    raw_truth_path <- file.path(
        CONFIG$summary_dir, "all_simulation_truth.csv.gz"
    )
    metrics_parts <- list()
    failure_parts <- list()
    plot_result_parts <- list()
    truth_parts <- list()
    setting_status_parts <- list()
    replication_ids <- integer()
    first_metadata <- NULL
    n_completed <- 0L

    for (file in files) {
        object <- tryCatch(readRDS(file), error = function(e) NULL)
        if (is.null(object) || !isTRUE(object$success) ||
            !identical(object$script_version, CONFIG$script_version)) {
            next
        }
        n_completed <- n_completed + 1L
        replication_ids[n_completed] <- object$replication
        if (is.null(first_metadata)) {
            first_metadata <- list(
                package_versions = object$package_versions,
                method_configuration = object$method_configuration,
                session_info = object$session_info
            )
        }

        result <- data.table::as.data.table(object$method_results)
        result[, method_label := method_display_name(method, component)]
        data.table::fwrite(
            result, raw_results_path,
            append = n_completed > 1L,
            col.names = n_completed == 1L,
            compress = "gzip"
        )
        metrics_parts[[n_completed]] <- make_replicate_metrics(result)
        failure_parts[[n_completed]] <- result[, .N, by = .(
            method_label, component, design_id, n_per_group,
            signal_fraction, n_signal_target, confounding, dgp, scenario,
            effect_parameter, status, available
        )]
        plot_result_parts[[n_completed]] <- result[
            (
                scenario %in% c(
                    "observed_prevalence_only",
                    "structural_matched_prevalence"
                ) & scenario_signal & structural_component(component)
            ) |
                (
                    scenario == "present_conditional_abundance" &
                        scenario_signal & method == "DASRA" &
                        component == "structural_absence"
                ) |
                (
                    scenario == "structural_only" & scenario_signal &
                        method == "DASRA" &
                        component == "present_conditional_abundance"
                ) |
                (
                    scenario == "present_conditional_abundance" &
                        effect_parameter == 0 &
                        abundance_component(component) & valid_p(p_value) &
                        n_per_group == CONFIG$primary_n_per_group &
                        abs(signal_fraction -
                            CONFIG$primary_signal_fraction) < 1e-12
                )
        ]

        truth_value <- data.table::as.data.table(object$truth)
        truth_value[, replication := object$replication]
        data.table::fwrite(
            truth_value, raw_truth_path,
            append = n_completed > 1L,
            col.names = n_completed == 1L,
            compress = "gzip"
        )
        truth_parts[[n_completed]] <- truth_value

        status_value <- data.table::as.data.table(object$setting_status)
        status_value[, replication := object$replication]
        setting_status_parts[[n_completed]] <- status_value
        rm(object, result, truth_value, status_value)
        gc(verbose = FALSE)
    }
    if (!n_completed) {
        stop("No replicate files match the current script version.", call. = FALSE)
    }

    metrics <- data.table::rbindlist(
        metrics_parts, fill = TRUE, use.names = TRUE
    )
    plot_results <- data.table::rbindlist(
        plot_result_parts, fill = TRUE, use.names = TRUE
    )
    truth <- data.table::rbindlist(
        truth_parts, fill = TRUE, use.names = TRUE
    )
    setting_status <- data.table::rbindlist(
        setting_status_parts, fill = TRUE, use.names = TRUE
    )
    failure_summary <- data.table::rbindlist(
        failure_parts, fill = TRUE, use.names = TRUE
    )[, .(N = sum(N)), by = .(
        method_label, component, design_id, n_per_group, signal_fraction,
        n_signal_target, confounding, dgp, scenario,
        effect_parameter, status, available
    )][order(method_label, design_id, scenario, effect_parameter, -N)]
    rm(metrics_parts, failure_parts, plot_result_parts,
       truth_parts, setting_status_parts)
    gc(verbose = FALSE)

    metric_summary <- aggregate_metric(
        metrics,
        c(
            "availability", "signal_availability", "null_availability",
            "raw_signal_rejection_rate", "type1_error",
            "power", "fdp", "power_by", "fdp_by",
            "designated_signal_raw_rejection",
            "designated_signal_bh_rejection",
            "designated_signal_by_rejection"
        )
    )
    truth_summary <- truth[, .(
        n_replications = data.table::uniqueN(replication),
        mean_structural_probability_difference = mean(
            abs(standardized_structural_probability_difference), na.rm = TRUE
        ),
        mean_expected_prevalence_difference = mean(
            abs(expected_prevalence_difference), na.rm = TRUE
        ),
        mean_standardized_expected_prevalence_difference = mean(
            abs(standardized_expected_prevalence_difference), na.rm = TRUE
        ),
        mean_group_z_correlation = mean(
            group_z_correlation, na.rm = TRUE
        ),
        mean_z_mean_difference = mean(
            z_mean_difference, na.rm = TRUE
        ),
        mean_present_conditional_effect = mean(
            abs(oracle_present_conditional_log_relative_effect), na.rm = TRUE
        ),
        mean_detected_effect = mean(
            abs(oracle_detected_log_relative_effect), na.rm = TRUE
        ),
        mean_nondetection_reference = mean(
            actual_nondetection_given_presence_reference, na.rm = TRUE
        ),
        mean_nondetection_comparison = mean(
            actual_nondetection_given_presence_comparison, na.rm = TRUE
        )
    ), by = .(
        setting_id, base_setting_id, design_id,
        n_per_group, total_sample_size, signal_fraction, n_signal_target,
        confounding, confounded, sample_size_label,
        signal_fraction_label, confounding_label, design_panel,
        template_role, study, dgp, scenario,
        effect_level, effect_index, effect_parameter,
        effect_measure, scenario_signal
    )]

    confounding_summary <- unique(truth[, .(
        replication, n_per_group, confounding, confounding_label,
        group_z_correlation, z_mean_difference
    )])[, .(
        n_replications = data.table::uniqueN(replication),
        mean_group_z_correlation = mean(group_z_correlation, na.rm = TRUE),
        sd_group_z_correlation = stats::sd(group_z_correlation, na.rm = TRUE),
        mean_z_mean_difference = mean(z_mean_difference, na.rm = TRUE),
        sd_z_mean_difference = stats::sd(z_mean_difference, na.rm = TRUE)
    ), by = .(n_per_group, confounding, confounding_label)]

    data.table::fwrite(
        metrics,
        file.path(CONFIG$summary_dir, "replicate_metrics.csv")
    )
    data.table::fwrite(
        metric_summary,
        file.path(CONFIG$summary_dir, "setting_summary.csv")
    )
    data.table::fwrite(
        failure_summary,
        file.path(CONFIG$summary_dir, "method_failure_summary.csv")
    )
    data.table::fwrite(
        truth_summary,
        file.path(CONFIG$summary_dir, "truth_summary.csv")
    )
    data.table::fwrite(
        confounding_summary,
        file.path(CONFIG$summary_dir, "confounding_design_summary.csv")
    )
    data.table::fwrite(
        setting_status,
        file.path(CONFIG$summary_dir, "setting_completion_status.csv")
    )
    data.table::fwrite(
        SETTINGS,
        file.path(CONFIG$summary_dir, "simulation_design.csv")
    )
    data.table::fwrite(
        BASE_SETTINGS,
        file.path(CONFIG$summary_dir, "base_simulation_design.csv")
    )
    data.table::fwrite(
        DESIGN_STRATA,
        file.path(CONFIG$summary_dir, "design_strata.csv")
    )
    data.table::fwrite(
        first_metadata$package_versions,
        file.path(CONFIG$summary_dir, "package_versions.csv")
    )
    data.table::fwrite(
        first_metadata$method_configuration,
        file.path(CONFIG$summary_dir, "method_configuration.csv")
    )

    writeLines(
        c(
            sprintf("Script version: %s", CONFIG$script_version),
            sprintf("Completed replications: %d", n_completed),
            sprintf("Base settings: %d", CONFIG$n_base_settings),
            sprintf("Design strata: %d", CONFIG$n_design_strata),
            sprintf("Total settings per replication: %d", CONFIG$n_total_settings),
            sprintf(
                "Replication IDs: %s",
                paste(sort(replication_ids), collapse = ", ")
            ),
            "",
            first_metadata$session_info
        ),
        file.path(CONFIG$summary_dir, "session_information.txt")
    )

    plot_structural_estimand(plot_results)
    plot_structural_design(truth)
    plot_component_specificity(plot_results)
    plot_abundance_power(metrics)
    plot_abundance_fdr(metrics)
    plot_availability(metrics)
    plot_global_null_qq(plot_results)
    plot_type1_error(metrics)
    plot_joint_community_robustness(metrics)
    plot_confounding_design(truth)

    messagef(
        "Combined %d completed replications into %s",
        n_completed, CONFIG$summary_dir
    )
    invisible(list(
        results_file = raw_results_path,
        truth = truth,
        metrics = metrics,
        summary = metric_summary
    ))
}

run_preflight <- function() {
    check_packages(replicate_packages)

    endpoint_settings <- SETTINGS[
        (SETTINGS$scenario == "observed_prevalence_only" &
            SETTINGS$effect_parameter == max(
                CONFIG$observed_prevalence_differences
            )) |
        (SETTINGS$scenario == "structural_matched_prevalence" &
            SETTINGS$effect_parameter == max(
                CONFIG$matched_structural_differences
            )) |
        (SETTINGS$scenario == "present_conditional_abundance" &
            SETTINGS$effect_parameter == max(CONFIG$abundance_effects)) |
        (SETTINGS$scenario == "structural_only" &
            SETTINGS$effect_parameter == max(
                CONFIG$specificity_structural_differences
            )),
        , drop = FALSE
    ]
    calibration_rows <- vector("list", nrow(endpoint_settings))
    for (i in seq_len(nrow(endpoint_settings))) {
        setting <- endpoint_settings[i, , drop = FALSE]
        design <- make_replication_design(
            replication_id = 9991L,
            n_per_group = setting$n_per_group,
            confounding = setting$confounding
        )
        template <- make_taxon_template(
            n_signal = setting$n_signal_target,
            role = setting$template_role
        )
        simulation <- simulate_setting(
            replication_id = 9991L,
            setting = setting,
            design = design,
            template = template
        )
        signal_truth <- simulation$truth[
            simulation$truth$scenario_signal, , drop = FALSE
        ]
        if (!nrow(signal_truth)) {
            stopf("Preflight found no designated signals in %s.", setting$setting_id)
        }
        target <- setting$effect_parameter
        achieved <- switch(
            setting$scenario,
            observed_prevalence_only = max(abs(
                abs(signal_truth$standardized_expected_prevalence_difference) -
                    target
            )),
            structural_matched_prevalence = max(c(
                abs(
                    abs(signal_truth$standardized_structural_probability_difference) -
                        target
                ),
                abs(signal_truth$standardized_expected_prevalence_difference)
            )),
            present_conditional_abundance = max(abs(
                abs(signal_truth$oracle_present_conditional_log_relative_effect) -
                    target
            )),
            structural_only = max(abs(
                abs(signal_truth$standardized_structural_probability_difference) -
                    target
            )),
            NA_real_
        )
        calibration_rows[[i]] <- data.frame(
            setting_id = setting$setting_id,
            scenario = setting$scenario,
            n_per_group = setting$n_per_group,
            signal_fraction = setting$signal_fraction,
            confounding = setting$confounding,
            target = target,
            maximum_calibration_error = achieved,
            stringsAsFactors = FALSE
        )
    }
    calibration_status <- do.call(rbind, calibration_rows)
    print(calibration_status, row.names = FALSE)
    if (any(!is.finite(calibration_status$maximum_calibration_error)) ||
        any(calibration_status$maximum_calibration_error > 5e-05)) {
        stop(
            "At least one endpoint setting failed the simulation-calibration precheck.",
            call. = FALSE
        )
    }

    selected <- SETTINGS[
        SETTINGS$scenario == "present_conditional_abundance" &
            SETTINGS$effect_parameter == 0 &
            SETTINGS$n_per_group == min(CONFIG$sample_sizes_per_group) &
            abs(SETTINGS$signal_fraction - max(CONFIG$signal_fractions)) < 1e-12 &
            SETTINGS$confounding == "confounded",
        , drop = FALSE
    ]
    if (nrow(selected) != 1L) {
        stop("Could not identify the hardest abundance global-null setting.",
             call. = FALSE)
    }
    design <- make_replication_design(
        replication_id = 9999L,
        n_per_group = selected$n_per_group,
        confounding = selected$confounding
    )
    template <- make_taxon_template(
        n_signal = selected$n_signal_target,
        role = selected$template_role
    )
    temporary_root <- tempfile("dasra_formal_preflight_")
    safe_dir_create(temporary_root)
    on.exit(
        unlink(temporary_root, recursive = TRUE, force = TRUE),
        add = TRUE
    )
    simulation <- simulate_setting(
        replication_id = 9999L,
        setting = selected,
        design = design,
        template = template
    )
    result <- analyze_setting(
        simulation = simulation,
        replication_id = 9999L,
        temporary_root = temporary_root
    )
    status <- data.table::as.data.table(result)[, .(
        n_taxa = .N,
        n_available = sum(available, na.rm = TRUE),
        availability = mean(available, na.rm = TRUE),
        representative_status = collapse_text(status)
    ), by = .(method, component)]
    print(selected[, c(
        "setting_id", "n_per_group", "signal_fraction", "confounding"
    )], row.names = FALSE)
    print(status)
    failed_methods <- status[
        n_available == 0,
        paste(method, component, sep = "::")
    ]
    if (length(failed_methods)) {
        stopf(
            "Preflight failed for: %s",
            paste(failed_methods, collapse = ", ")
        )
    }
    message(
        "All endpoint calibrations and method/component interfaces passed preflight."
    )
    invisible(list(
        calibration = calibration_status,
        method_status = status
    ))
}

print_design <- function() {
    cat(sprintf(
        "Base settings: %d\nDesign strata: %d\nTotal settings: %d\n\n",
        nrow(BASE_SETTINGS), nrow(DESIGN_STRATA), nrow(SETTINGS)
    ))
    cat("Design strata:\n")
    print(DESIGN_STRATA, row.names = FALSE)
    cat("\nBase settings:\n")
    print(BASE_SETTINGS, row.names = FALSE)
    cat("\nFull setting grid:\n")
    print(SETTINGS, row.names = FALSE)
    cat("\nConfiguration:\n")
    print(CONFIG)
}

main <- function() {
    arguments <- commandArgs(trailingOnly = TRUE)
    if (!length(arguments)) {
        stop(
            paste(
                "Usage: Rscript dasra_formal_hpc_simulation.R",
                "replicate <id> | summarize | design | preflight"
            ),
            call. = FALSE
        )
    }
    mode <- tolower(arguments[[1L]])
    if (mode == "replicate") {
        if (length(arguments) < 2L) {
            stop("replicate mode requires a replication ID.", call. = FALSE)
        }
        replication_id <- suppressWarnings(as.integer(arguments[[2L]]))
        if (length(replication_id) != 1L || is.na(replication_id) ||
            replication_id < 1L) {
            stop("The replication ID must be a positive integer.", call. = FALSE)
        }
        run_replication(replication_id)
    } else if (mode == "summarize") {
        summarize_completed_replications()
    } else if (mode == "design") {
        print_design()
    } else if (mode == "preflight") {
        run_preflight()
    } else {
        stopf("Unknown mode: %s", mode)
    }
}

main()
