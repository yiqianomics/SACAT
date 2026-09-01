# Run the Rauer defined mock-community benchmark.

options(stringsAsFactors = FALSE, warn = 1)

locate_dataset_directory <- function() {
    file_argument <- grep(
        "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
    )
    if (length(file_argument) == 1L) {
        script_path <- sub("^--file=", "", file_argument)
        return(dirname(normalizePath(script_path)))
    }
    candidates <- unique(c(
        getwd(),
        file.path(getwd(), "Analysis", "RealDataAnalysis", "rauer_mock")
    ))
    candidates <- normalizePath(candidates, mustWork = FALSE)
    matches <- candidates[file.exists(file.path(candidates, "analysis.R"))]
    if (length(matches) != 1L) {
        stop("Could not identify the rauer_mock directory.", call. = FALSE)
    }
    matches[[1L]]
}

arguments <- commandArgs(trailingOnly = TRUE)
design_argument <- grep("^--design=", arguments, value = TRUE)
if (length(design_argument) > 1L) {
    stop("Specify at most one --design argument.", call. = FALSE)
}
unknown_arguments <- setdiff(arguments, design_argument)
if (length(unknown_arguments)) {
    stop(sprintf("Unknown argument: %s", unknown_arguments[[1L]]),
         call. = FALSE)
}
requested_design <- if (length(design_argument)) {
    sub("^--design=", "", design_argument)
} else {
    "all"
}
if (!requested_design %in% c("all", "nonnull", "null")) {
    stop("--design must be all, nonnull, or null.", call. = FALSE)
}

workers <- suppressWarnings(as.integer(Sys.getenv("DASRA_WORKERS", "1")))
if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("DASRA_WORKERS must be a positive integer.", call. = FALSE)
}

required_packages <- c("DASRA", "maaslin3", "ZINQ")
missing_packages <- required_packages[!vapply(
    required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
    stop(
        sprintf(
            "Required R packages are unavailable: %s",
            paste(missing_packages, collapse = ", ")
        ),
        call. = FALSE
    )
}

dataset_directory <- locate_dataset_directory()
input_file <- file.path(
    dataset_directory, "processed", "rauer_mock_dasra_input.rds"
)
if (!file.exists(input_file)) {
    stop("Run prepare_data.R before the benchmark analysis.", call. = FALSE)
}
work_directory <- file.path(dataset_directory, "work")
table_directory <- file.path(dataset_directory, "table")
dir.create(work_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

prepared <- readRDS(input_file)
required_input_names <- c(
    "counts", "metadata", "taxon_metadata", "preprocessing"
)
if (!identical(sort(names(prepared)), sort(required_input_names))) {
    stop("The prepared input has an unexpected structure.", call. = FALSE)
}

taxa <- as.character(prepared$taxon_metadata$taxon)
d6300_taxa <- prepared$taxon_metadata$taxon[
    prepared$taxon_metadata$source_standard == "D6300"
]
d6321_taxa <- prepared$taxon_metadata$taxon[
    prepared$taxon_metadata$source_standard == "D6321"
]
target_counts <- t(prepared$counts[taxa, , drop = FALSE])
metadata <- prepared$metadata[colnames(prepared$counts), , drop = FALSE]
if (nrow(target_counts) != 32L || ncol(target_counts) != 11L ||
    length(d6300_taxa) != 8L || length(d6321_taxa) != 3L ||
    !identical(rownames(target_counts), metadata$sample_id)) {
    stop("The prepared count and metadata tables failed validation.",
         call. = FALSE)
}

other_counts <- metadata$library_size - rowSums(target_counts)
if (any(other_counts < 0L)) {
    stop("Reference-taxon counts exceed the recorded library size.",
         call. = FALSE)
}
full_counts <- cbind(
    target_counts,
    Other_unmodeled = as.integer(other_counts)
)
storage.mode(full_counts) <- "integer"

metadata$source_standard <- factor(
    metadata$source_standard, levels = c("D6300", "D6321")
)
metadata$input_stratum <- factor(
    metadata$input_stratum, levels = c("higher", "lower")
)
metadata$extraction_kit <- factor(metadata$extraction_kit)
metadata$lysis <- factor(metadata$lysis)
metadata$extraction_buffer <- factor(metadata$extraction_buffer)

d6300_metadata <- metadata[metadata$source_standard == "D6300", , drop = FALSE]
cell_metadata <- d6300_metadata[
    match(1:16, d6300_metadata$pair_cell),
    c("input_stratum", "extraction_kit", "lysis", "extraction_buffer"),
    drop = FALSE
]
if (anyNA(cell_metadata)) {
    stop("The matched experimental cells are incomplete.", call. = FALSE)
}

nonnull_candidates <- utils::combn(16L, 4L)
nonnull_balanced <- apply(nonnull_candidates, 2L, function(cells) {
    chosen <- cell_metadata[cells, , drop = FALSE]
    all(vapply(chosen, function(values) {
        counts <- table(values)
        length(counts) == 2L && all(counts == 2L)
    }, logical(1)))
})
nonnull_candidates <- nonnull_candidates[, nonnull_balanced, drop = FALSE]
if (ncol(nonnull_candidates) != 52L) {
    stop("The expected 52 nonnull allocations were not recovered.",
         call. = FALSE)
}
set.seed(20260831L)
nonnull_allocations <- nonnull_candidates[
    , sample.int(ncol(nonnull_candidates)), drop = FALSE
]

null_candidates <- utils::combn(16L, 8L)
null_balanced <- apply(null_candidates, 2L, function(cells) {
    chosen <- cell_metadata[cells, , drop = FALSE]
    all(vapply(chosen, function(values) {
        counts <- table(values)
        length(counts) == 2L && all(counts == 4L)
    }, logical(1)))
})
null_candidates <- null_candidates[, null_balanced, drop = FALSE]
null_allocations <- null_candidates[
    , apply(null_candidates, 2L, function(cells) 1L %in% cells),
    drop = FALSE
]
if (ncol(null_allocations) != 111L) {
    stop("The expected 111 complete-null allocations were not recovered.",
         call. = FALSE)
}

depths <- c(
    native = NA_real_, `1000` = 1000, `750` = 750, `500` = 500,
    `250` = 250, `100` = 100, `50` = 50
)
depth_seed_ids <- c(
    native = 1L, `1000` = 2L, `750` = 7L, `500` = 3L,
    `250` = 4L, `100` = 5L, `50` = 6L
)
if (min(rowSums(full_counts)) < max(depths, na.rm = TRUE)) {
    stop("At least one library is too small for the requested depths.",
         call. = FALSE)
}

rarefy_row <- function(counts, target_depth) {
    counts <- as.integer(counts)
    if (sum(counts) == target_depth) return(counts)
    draws_left <- as.integer(target_depth)
    total_left <- sum(counts)
    output <- integer(length(counts))
    for (index in seq_len(length(counts) - 1L)) {
        output[[index]] <- stats::rhyper(
            1L,
            counts[[index]],
            total_left - counts[[index]],
            draws_left
        )
        draws_left <- draws_left - output[[index]]
        total_left <- total_left - counts[[index]]
    }
    output[[length(counts)]] <- draws_left
    output
}

rarefy_table <- function(counts, target_depth) {
    if (is.na(target_depth)) return(counts)
    rarefied <- t(apply(
        counts, 1L, rarefy_row, target_depth = target_depth
    ))
    dimnames(rarefied) <- dimnames(counts)
    storage.mode(rarefied) <- "integer"
    if (any(rowSums(rarefied) != target_depth)) {
        stop("Downsampling did not produce the requested depth.",
             call. = FALSE)
    }
    rarefied
}

adjust_family <- function(p_value, available) {
    p_for_adjustment <- ifelse(
        available & is.finite(p_value), p_value, 1
    )
    stats::p.adjust(p_for_adjustment, method = "BH")
}

extract_maaslin_component <- function(results, component) {
    results <- results[
        results$metadata == "group" & results$value == "comparison" &
            results$feature %in% taxa,
        ,
        drop = FALSE
    ]
    p_value <- estimate <- stats::setNames(rep(NA_real_, length(taxa)), taxa)
    available <- stats::setNames(rep(FALSE, length(taxa)), taxa)
    reason <- stats::setNames(
        rep("model did not return a finite p-value", length(taxa)), taxa
    )
    for (taxon in taxa) {
        row <- results[results$feature == taxon, , drop = FALSE]
        if (nrow(row) != 1L) next
        p_value[[taxon]] <- row$pval_individual[[1L]]
        estimate[[taxon]] <- row$coef[[1L]]
        available[[taxon]] <- is.finite(p_value[[taxon]])
        if (available[[taxon]]) reason[[taxon]] <- "available"
    }
    data.frame(
        taxon = taxa,
        method = paste("MaAsLin3", component),
        component = component,
        available = unname(available),
        reason = unname(reason),
        p_value = unname(p_value),
        estimate = unname(estimate),
        stringsAsFactors = FALSE
    )
}

run_maaslin <- function(counts, model_metadata, output_directory) {
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    on.exit({
        try(maaslin3::maaslin_log_reset(), silent = TRUE)
        unlink(output_directory, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    use_depth <- stats::sd(model_metadata$library_size) > 0
    model_formula <- if (use_depth) {
        ~ extraction_kit + lysis + extraction_buffer + input_stratum +
            log_library_size_z + group
    } else {
        ~ extraction_kit + lysis + extraction_buffer + input_stratum + group
    }
    fit <- maaslin3::maaslin3(
        input_data = as.data.frame(counts, check.names = FALSE),
        input_metadata = model_metadata,
        output = output_directory,
        formula = model_formula,
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
        warn_prevalence = TRUE,
        augment = TRUE,
        cores = 1,
        max_significance = 1,
        plot_summary_plot = FALSE,
        plot_associations = FALSE,
        save_models = FALSE,
        save_plots_rds = FALSE,
        verbosity = "ERROR",
        reference = "group,reference"
    )
    maaslin3::maaslin_log_reset()
    rbind(
        extract_maaslin_component(
            fit$fit_data_prevalence$results, "prevalence"
        ),
        extract_maaslin_component(
            fit$fit_data_abundance$results, "abundance"
        )
    )
}

quantile_cauchy_p <- function(p_values, taus) {
    if (length(p_values) != length(taus) || any(!is.finite(p_values)) ||
        any(p_values < 0 | p_values > 1)) return(NA_real_)
    p_values <- pmin(pmax(p_values, 1e-15), 1 - 1e-15)
    weights <- ifelse(taus <= 0.5, taus, 1 - taus)
    weights <- weights / sum(weights)
    statistic <- sum(weights * tan((0.5 - p_values) * pi))
    min(max(1 - stats::pcauchy(statistic), 0), 1)
}

run_zinq <- function(counts, model_metadata, seed) {
    relative_abundance <- sweep(
        counts[, taxa, drop = FALSE],
        1L,
        model_metadata$library_size,
        "/"
    )
    group_binary <- as.integer(model_metadata$group == "comparison")
    use_depth <- stats::sd(model_metadata$library_size) > 0
    covariate_formula <- if (use_depth) {
        ~ extraction_kit + lysis + extraction_buffer + input_stratum +
            log_library_size_z
    } else {
        ~ extraction_kit + lysis + extraction_buffer + input_stratum
    }
    covariates <- stats::model.matrix(
        covariate_formula, model_metadata
    )[, -1L, drop = FALSE]
    colnames(covariates) <- paste0("covariate", seq_len(ncol(covariates)))
    model_formula <- stats::as.formula(paste(
        "y ~", paste(c("group", colnames(covariates)), collapse = " + ")
    ))
    taus <- c(0.25, 0.50, 0.75)

    do.call(rbind, lapply(seq_along(taxa), function(index) {
        taxon <- taxa[[index]]
        model_data <- data.frame(
            y = as.numeric(relative_abundance[, taxon]),
            group = group_binary,
            covariates,
            check.names = FALSE
        )
        fit <- tryCatch(
            ZINQ::ZINQ_tests(
                formula.logistic = model_formula,
                formula.quantile = model_formula,
                C = "group",
                y_CorD = "C",
                data = model_data,
                taus = taus,
                seed = seed + index
            ),
            error = function(error) error
        )
        if (inherits(fit, "error")) {
            return(data.frame(
                taxon = taxon,
                method = c("ZINQ prevalence", "ZINQ abundance"),
                component = c("prevalence", "abundance"),
                available = FALSE,
                reason = "model did not return a finite p-value",
                p_value = NA_real_,
                estimate = NA_real_,
                stringsAsFactors = FALSE
            ))
        }
        p_values <- c(
            as.numeric(fit$pvalue.logistic)[1L],
            quantile_cauchy_p(as.numeric(fit$pvalue.quantile), taus)
        )
        data.frame(
            taxon = taxon,
            method = c("ZINQ prevalence", "ZINQ abundance"),
            component = c("prevalence", "abundance"),
            available = is.finite(p_values),
            reason = ifelse(
                is.finite(p_values),
                "available",
                "model did not return a finite p-value"
            ),
            p_value = p_values,
            estimate = NA_real_,
            stringsAsFactors = FALSE
        )
    }))
}

make_group_metadata <- function(design, allocation) {
    output <- metadata
    output$group <- ifelse(
        (output$source_standard == "D6300" &
            output$pair_cell %in% allocation) |
            (output$source_standard == "D6321" &
                !(output$pair_cell %in% allocation)),
        "comparison",
        "reference"
    )
    output$group <- factor(
        output$group, levels = c("reference", "comparison")
    )
    if (!all(table(output$group) == 16L)) {
        stop("An allocation did not produce two groups of 16 libraries.",
             call. = FALSE)
    }
    if (design == "nonnull") {
        if (!identical(
            as.integer(table(output$group, output$source_standard)),
            c(12L, 4L, 4L, 12L)
        )) {
            stop("A nonnull allocation has an unexpected source balance.",
                 call. = FALSE)
        }
    } else if (!all(table(output$group, output$source_standard) == 8L)) {
        stop("A null allocation has an unexpected source balance.",
             call. = FALSE)
    }
    if (!all(table(
        output$group,
        output$input_stratum,
        output$extraction_kit,
        output$lysis,
        output$extraction_buffer
    ) == 1L)) {
        stop("An allocation is not balanced across experimental factors.",
             call. = FALSE)
    }
    if (design == "null") {
        for (factor_name in c(
            "input_stratum", "extraction_kit", "lysis",
            "extraction_buffer"
        )) {
            if (!all(table(
                output$group,
                output$source_standard,
                output[[factor_name]]
            ) == 4L)) {
                stop(
                    "A null allocation is not balanced within each mock community.",
                    call. = FALSE
                )
            }
        }
    }
    output
}

run_one_depth <- function(design, allocation_id, allocation,
                          depth_name, target_depth, model_metadata) {
    seed <- 630000L + allocation_id * 100L +
        depth_seed_ids[[depth_name]]
    set.seed(seed)
    counts <- rarefy_table(full_counts, target_depth)
    library_size <- rowSums(counts)
    model_metadata$library_size <- library_size
    model_metadata$log_library_size_z <- if (stats::sd(library_size) > 0) {
        as.numeric(scale(log(library_size)))
    } else {
        0
    }

    dasra_fit <- DASRA::dasra(
        counts = counts[, taxa, drop = FALSE],
        metadata = model_metadata,
        formula = ~ group + extraction_kit + lysis + extraction_buffer +
            input_stratum,
        group = "group",
        library_size = library_size,
        taxa_are_rows = FALSE,
        reference = "reference",
        component = "all",
        full_output = FALSE,
        structural_conditional_present_starts = "full",
        workers = 1L,
        verbose = FALSE
    )
    dasra_rows <- rbind(
        data.frame(
            taxon = taxa,
            method = "DASRA structural absence",
            component = "structural absence",
            available = dasra_fit$diagnostics$formed_structural_absence,
            reason = ifelse(
                dasra_fit$diagnostics$formed_structural_absence,
                "available",
                dasra_fit$diagnostics$reason_structural_absence
            ),
            p_value = dasra_fit$results$p_structural_absence,
            estimate = dasra_fit$results$z_structural_absence,
            stringsAsFactors = FALSE
        ),
        data.frame(
            taxon = taxa,
            method = "DASRA present-conditional abundance",
            component = "abundance",
            available = dasra_fit$diagnostics$formed_relative_abundance,
            reason = ifelse(
                dasra_fit$diagnostics$formed_relative_abundance,
                "available",
                dasra_fit$diagnostics$reason_relative_abundance
            ),
            p_value = dasra_fit$results$p_relative_abundance,
            estimate = dasra_fit$results$estimate_relative_abundance,
            stringsAsFactors = FALSE
        )
    )

    maaslin_rows <- tryCatch(
        run_maaslin(
            counts,
            model_metadata,
            tempfile(
                sprintf(
                    "maaslin-%s-%03d-%s-",
                    design, allocation_id, depth_name
                ),
                tmpdir = work_directory
            )
        ),
        error = function(error) {
            expand.grid(
                taxon = taxa,
                method = c(
                    "MaAsLin3 prevalence", "MaAsLin3 abundance"
                ),
                KEEP.OUT.ATTRS = FALSE,
                stringsAsFactors = FALSE
            ) |>
                transform(
                    component = ifelse(
                        method == "MaAsLin3 prevalence",
                        "prevalence",
                        "abundance"
                    ),
                    available = FALSE,
                    reason = "model did not return a finite p-value",
                    p_value = NA_real_,
                    estimate = NA_real_
                )
        }
    )
    zinq_rows <- run_zinq(counts, model_metadata, seed + 10000L)
    rows <- rbind(dasra_rows, maaslin_rows, zinq_rows)
    rows$q_value <- ave(
        seq_len(nrow(rows)),
        rows$method,
        FUN = function(indices) {
            adjust_family(rows$p_value[indices], rows$available[indices])
        }
    )
    rows$significant <- rows$available & rows$q_value <= 0.05
    rows$allocation_id <- allocation_id
    rows$design <- design
    rows$depth <- depth_name
    rows$target_library_size <- if (is.na(target_depth)) {
        NA_real_
    } else {
        target_depth
    }
    rows$median_library_size <- stats::median(library_size)
    rows$target_set <- ifelse(
        rows$taxon %in% d6300_taxa,
        "Structural targets",
        "Abundance targets"
    )
    rows$d6300_comparison_cells <- paste(allocation, collapse = ";")
    rows
}

method_order <- c(
    "DASRA structural absence",
    "DASRA present-conditional abundance",
    "MaAsLin3 prevalence",
    "MaAsLin3 abundance",
    "ZINQ prevalence",
    "ZINQ abundance"
)
checkpoint_row_columns <- c(
    "taxon", "method", "component", "available", "reason", "p_value",
    "estimate", "q_value", "significant", "allocation_id", "design",
    "depth", "target_library_size", "median_library_size", "target_set",
    "d6300_comparison_cells"
)

run_allocation <- function(design, allocation_id, allocations) {
    design_directory <- file.path(work_directory, design)
    dir.create(design_directory, recursive = TRUE, showWarnings = FALSE)
    output_file <- file.path(
        design_directory, sprintf("allocation_%03d.rds", allocation_id)
    )
    allocation <- allocations[, allocation_id]
    analysis_id <- paste0("rauer_mock_component_diagnostics_", design)
    if (file.exists(output_file)) {
        existing <- readRDS(output_file)
        expected_rows <- length(depths) * length(method_order) * length(taxa)
        valid_existing <- is.list(existing) &&
            identical(existing$analysis_id, analysis_id) &&
            identical(as.integer(existing$allocation),
                      as.integer(allocation)) &&
            identical(existing$depths, depths) &&
            is.data.frame(existing$rows) &&
            identical(names(existing$rows), checkpoint_row_columns) &&
            nrow(existing$rows) == expected_rows &&
            !anyDuplicated(existing$rows[c("depth", "method", "taxon")])
        if (valid_existing) {
            return(existing)
        }
    }

    model_metadata <- make_group_metadata(design, allocation)
    rows <- do.call(rbind, lapply(names(depths), function(depth_name) {
        run_one_depth(
            design,
            allocation_id,
            allocation,
            depth_name,
            depths[[depth_name]],
            model_metadata
        )
    }))
    output <- list(
        analysis_id = analysis_id,
        allocation_id = allocation_id,
        allocation = allocation,
        depths = depths,
        rows = rows
    )
    saveRDS(output, output_file, compress = "xz")
    message(sprintf("Completed %s allocation %d", design, allocation_id))
    output
}

run_design <- function(design, allocations) {
    allocation_results <- parallel::mclapply(
        seq_len(ncol(allocations)),
        run_allocation,
        design = design,
        allocations = allocations,
        mc.cores = min(workers, ncol(allocations)),
        mc.preschedule = FALSE,
        mc.set.seed = FALSE
    )
    results <- do.call(rbind, lapply(allocation_results, `[[`, "rows"))
    rownames(results) <- NULL
    expected_rows <- ncol(allocations) * length(depths) *
        length(method_order) * length(taxa)
    if (nrow(results) != expected_rows || anyDuplicated(results[
        c("allocation_id", "depth", "method", "taxon")
    ])) {
        stop(sprintf("The %s result grid is incomplete.", design),
             call. = FALSE)
    }
    results$depth <- factor(results$depth, levels = names(depths))
    results$method <- factor(results$method, levels = method_order)
    results <- results[order(
        results$allocation_id,
        results$depth,
        results$method,
        match(results$taxon, taxa)
    ), ]
    results$depth <- as.character(results$depth)
    results$method <- as.character(results$method)
    output_columns <- c(
        "design", "allocation_id", "depth", "target_library_size",
        "median_library_size", "taxon", "target_set", "method",
        "component", "available", "reason", "p_value", "estimate",
        "q_value", "significant", "d6300_comparison_cells"
    )
    results <- results[, output_columns]
    utils::write.csv(
        results,
        file.path(
            table_directory,
            sprintf("rauer_mock_%s_results.csv", design)
        ),
        row.names = FALSE,
        quote = TRUE,
        na = ""
    )
    results
}

if (requested_design %in% c("all", "nonnull")) {
    invisible(run_design("nonnull", nonnull_allocations))
}
if (requested_design %in% c("all", "null")) {
    invisible(run_design("null", null_allocations))
}
