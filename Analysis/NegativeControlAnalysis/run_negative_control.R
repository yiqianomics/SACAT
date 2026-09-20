#!/usr/bin/env Rscript

# Native counts and full-table library sizes stay fixed; only labels change.
# Source this file to use the reader, single-replicate runner, or CSV exporter.

negative_control_datasets <- c(
    "BaxterE_2016", "cdi_schubert", "gems_pediatric_diarrhea",
    "korean_hypertension", "qiita_1939_pediatric_crohn", "NogueraJulianHIV",
    "ArtPrize_2015_Forehead", "HMPV35Throat", "Qiita11993Colombia",
    "Qiita13631ASD", "NHANESOral_2009_2010", "NHANESOral_2011_2012",
    "Goodrich_2014", "Zeevi_2015_16S", "Singh_2015_EDD", "AmericanGut_2018",
    "HillBurns_2017", "Heinzel_2021_TREND", "Wallen_2020_Birmingham", "iMSMS_2022"
)

read_negative_control_input <- function(dataset_directory) {
    input <- readRDS(file.path(dataset_directory, "analysis_input.rds"))
    stopifnot(
        is.list(input), length(input$dataset_id) == 1L,
        input$dataset_id %in% negative_control_datasets,
        identical(input$dataset_id, basename(normalizePath(dataset_directory))),
        length(input$label) == 1L, nzchar(input$label),
        length(input$dataset_index) == 1L, is.finite(input$dataset_index),
        input$dataset_index > 0, input$dataset_index == as.integer(input$dataset_index),
        is.matrix(input$counts), is.numeric(input$counts),
        identical(dim(input$counts), c(30L, 200L)),
        length(input$evaluation_taxa) == 30L,
        !anyNA(input$evaluation_taxa), all(nzchar(input$evaluation_taxa)),
        !anyDuplicated(input$evaluation_taxa),
        identical(rownames(input$counts), input$evaluation_taxa),
        identical(colnames(input$counts), rownames(input$metadata)),
        !anyDuplicated(colnames(input$counts)),
        !anyNA(input$counts), all(is.finite(input$counts)),
        all(input$counts >= 0), all(input$counts == round(input$counts)),
        length(input$metadata$library_size) == 200L,
        all(is.finite(input$metadata$library_size)),
        all(input$metadata$library_size > 0),
        all(input$metadata$library_size == round(input$metadata$library_size)),
        all(colSums(input$counts) <= input$metadata$library_size)
    )
    # The complete genus table was used to select this panel before randomization.
    stopifnot(
        all(rowMeans(input$counts > 0) >= 0.10),
        all(rowMeans(sweep(input$counts, 2L, input$metadata$library_size, "/")) >= 1e-5)
    )
    input
}

negative_control_groups <- function(input, replicate) {
    stopifnot(length(replicate) == 1L, is.finite(replicate),
              replicate == as.integer(replicate), replicate >= 1L, replicate <= 100L)
    seed <- 202609190L + input$dataset_index * 100000L + replicate * 1009L
    set.seed(seed)
    list(seed = seed, group = factor(sample(rep(c("A", "B"), each = 100L)),
                                    levels = c("A", "B")))
}

validate_negative_control_result <- function(result, input, replicate) {
    grouping <- negative_control_groups(input, replicate)
    stopifnot(
        identical(result$dataset, input$dataset_id),
        identical(result$package_version, "0.7.1"),
        result$replicate == replicate, result$group_seed == grouping$seed,
        identical(result$taxa, input$evaluation_taxa),
        nrow(result$results) == 30L, nrow(result$diagnostics) == 30L,
        setequal(result$results$taxon, input$evaluation_taxa),
        setequal(result$diagnostics$taxon, input$evaluation_taxa),
        identical(result$assignments$sample_id, colnames(input$counts)),
        identical(as.character(result$assignments$group), as.character(grouping$group)),
        identical(result$assignments$library_size, input$metadata$library_size),
        identical(result$settings$group, "group"),
        identical(result$settings$component, "all"),
        identical(result$settings$p_adjust_method, "BH"),
        identical(all.vars(result$settings$formula), "group")
    )
    invisible(result)
}

run_negative_control_replicate <- function(input, replicate, result_directory) {
    path <- file.path(result_directory, "replicates",
                      sprintf("replicate_%03d.rds", replicate))
    if (file.exists(path)) {
        result <- readRDS(path)
        validate_negative_control_result(result, input, replicate)
        return(result)
    }
    if (!requireNamespace("SACAT", quietly = TRUE) ||
        as.character(utils::packageVersion("SACAT")) != "0.7.1") {
        stop("Install SACAT 0.7.1 before running the negative controls.", call. = FALSE)
    }
    grouping <- negative_control_groups(input, replicate)
    metadata <- input$metadata
    metadata$group <- grouping$group
    started <- proc.time()[["elapsed"]]
    fit <- SACAT::sacat(
        counts = input$counts, metadata = metadata,
        formula = ~ group, group = "group", library_size = "library_size",
        taxa_are_rows = TRUE, reference = "A", p_adjust_method = "BH",
        component = "all", full_output = FALSE,
        structural_conditional_present_starts = 1L
    )
    depth_a <- metadata$library_size[metadata$group == "A"]
    depth_b <- metadata$library_size[metadata$group == "B"]
    result <- list(
        dataset = input$dataset_id, label = input$label, replicate = replicate,
        taxa = input$evaluation_taxa, results = fit$results,
        diagnostics = fit$diagnostics, settings = fit$settings,
        package_version = "0.7.1", group_seed = grouping$seed,
        assignments = data.frame(sample_id = rownames(metadata),
            group = metadata$group, library_size = metadata$library_size),
        depth = data.frame(dataset = input$dataset_id, replicate = replicate,
            n_A = length(depth_a), n_B = length(depth_b),
            median_A = median(depth_a), median_B = median(depth_b),
            mean_A = mean(depth_a), mean_B = mean(depth_b)),
        elapsed_seconds = proc.time()[["elapsed"]] - started
    )
    validate_negative_control_result(result, input, replicate)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(result, paste0(path, ".tmp"))
    stopifnot(file.rename(paste0(path, ".tmp"), path))
    result
}

negative_control_setting_text <- function(value) {
    if (inherits(value, "formula")) return(paste(deparse(value), collapse = " "))
    if (!is.null(names(value))) value <- paste(names(value), value, sep = "=")
    paste(value, collapse = "; ")
}

# 'results' may be saved result objects or paths to the original per-replicate fits.
# Exporting requires no SACAT installation and never fits a model.
export_negative_control_results <- function(input, results, result_directory,
                                            expected_replicates = 100L) {
    if (is.character(results)) results <- lapply(results, readRDS)
    stopifnot(is.list(results), length(results) == expected_replicates,
              expected_replicates >= 1L, expected_replicates <= 100L)
    replicate_ids <- vapply(results, function(x) as.integer(x$replicate), integer(1))
    stopifnot(!anyDuplicated(replicate_ids),
              setequal(replicate_ids, seq_len(expected_replicates)))
    results <- results[order(replicate_ids)]
    components <- c(SA = "structural_absence", RA = "relative_abundance",
                    omnibus = "omnibus")
    taxon_tables <- diagnostic_tables <- depth_tables <- setting_tables <-
        vector("list", length(results))

    for (i in seq_along(results)) {
        result <- results[[i]]
        validate_negative_control_result(result, input, i)
        p_index <- match(input$evaluation_taxa, result$results$taxon)
        d_index <- match(input$evaluation_taxa, result$diagnostics$taxon)
        diagnostics <- result$diagnostics[d_index, , drop = FALSE]
        taxon_tables[[i]] <- do.call(rbind, lapply(names(components), function(component) {
            suffix <- components[[component]]
            p <- result$results[[paste0("p_", suffix)]][p_index]
            available <- diagnostics$retained &
                diagnostics[[paste0("formed_", suffix)]] &
                is.finite(p) & p >= 0 & p <= 1
            stopifnot(length(p) == 30L, !anyNA(available))
            # Keep every selected genus in the denominator; unavailable tests do not reject.
            analysis_p <- ifelse(available, p, 1)
            data.frame(dataset = input$dataset_id, label = input$label,
                replicate = i, taxon = input$evaluation_taxa, component = component,
                raw_p = p, available = available, analysis_p = analysis_p,
                bh_adjusted_p = stats::p.adjust(analysis_p, "BH"),
                rejected = analysis_p < 0.05, stringsAsFactors = FALSE)
        }))
        diagnostic_tables[[i]] <- cbind(data.frame(
            dataset = input$dataset_id, label = input$label, replicate = i), diagnostics)
        depth_tables[[i]] <- result$depth
        depth_tables[[i]]$group_seed <- result$group_seed
        settings <- c(list(package_version = result$package_version,
            planned_randomizations = 100L,
            reported_randomizations = expected_replicates,
            full_output = FALSE, type1_error_threshold = 0.05,
            type1_error_uses = "raw component and omnibus p-values; unavailable tests count as no rejection",
            library_size_definition = input$source$library_size),
            result$settings)
        setting_tables[[i]] <- data.frame(dataset = input$dataset_id, replicate = i,
            setting = names(settings),
            value = vapply(settings, negative_control_setting_text, character(1)),
            stringsAsFactors = FALSE)
    }
    taxon <- do.call(rbind, taxon_tables)
    replicate <- do.call(rbind, lapply(seq_along(results), function(i) {
        do.call(rbind, lapply(names(components), function(component) {
            rows <- taxon[taxon$replicate == i & taxon$component == component, ]
            data.frame(dataset = input$dataset_id, label = input$label,
                replicate = i, component = component, taxa = nrow(rows),
                rejected_taxa = sum(rows$rejected), type1_error = mean(rows$rejected),
                availability = mean(rows$available), stringsAsFactors = FALSE)
        }))
    }))
    tables <- list(taxon_pvalues = taxon, replicate_metrics = replicate,
        diagnostics = do.call(rbind, diagnostic_tables),
        depth_diagnostics = do.call(rbind, depth_tables),
        analysis_settings = do.call(rbind, setting_tables))
    # Constant settings appear once; retain replicate IDs only for settings that vary.
    settings <- tables$analysis_settings
    tables$analysis_settings <- do.call(rbind, lapply(unique(settings$setting), function(name) {
        rows <- settings[settings$setting == name, , drop = FALSE]
        if (length(unique(rows$value)) == 1L) {
            rows <- rows[1L, , drop = FALSE]
            rows$replicate <- NA_integer_
        }
        rows
    }))
    dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
    for (name in names(tables)) {
        utils::write.csv(tables[[name]], file.path(result_directory, paste0(name, ".csv")),
                         row.names = FALSE, na = "")
    }
    invisible(tables)
}

run_negative_control <- function(dataset_directory, workers = 7L, replicates = 100L) {
    stopifnot(length(workers) == 1L, is.finite(workers), workers >= 1L,
              workers == as.integer(workers), length(replicates) == 1L,
              is.finite(replicates), replicates >= 1L, replicates <= 100L,
              replicates == as.integer(replicates))
    input <- read_negative_control_input(dataset_directory)
    result_directory <- file.path(dataset_directory, "results")
    # Each worker fits one randomization at a time.
    thread_variables <- c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
                          "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
    previous <- Sys.getenv(thread_variables, unset = NA_character_)
    do.call(Sys.setenv, as.list(setNames(rep("1", length(thread_variables)), thread_variables)))
    on.exit({
        Sys.unsetenv(thread_variables[is.na(previous)])
        do.call(Sys.setenv, as.list(previous[!is.na(previous)]))
    }, add = TRUE)
    if (.Platform$OS.type == "windows" && workers > 1L) {
        message("Forked workers are unavailable on Windows; running sequentially.")
        workers <- 1L
    }
    run <- function(i) run_negative_control_replicate(input, i, result_directory)
    results <- if (workers == 1L) {
        lapply(seq_len(replicates), run)
    } else {
        parallel::mclapply(seq_len(replicates), run,
            mc.cores = min(workers, replicates), mc.preschedule = FALSE, mc.set.seed = FALSE)
    }
    failed <- vapply(results, function(x) inherits(x, "try-error") || is.null(x), logical(1))
    if (any(failed)) {
        stop("Randomizations failed; completed caches were retained. Re-run to resume.", call. = FALSE)
    }
    tables <- export_negative_control_results(input, results, result_directory, replicates)
    message(sprintf("%s: %d randomizations exported.", input$label, replicates))
    invisible(tables)
}

if (sys.nframe() == 0L) {
    arguments <- commandArgs(trailingOnly = TRUE)
    if (!length(arguments) || arguments[[1L]] %in% c("-h", "--help")) {
        cat("Usage: Rscript run_negative_control.R DATASET [WORKERS]\n",
            "       Rscript run_negative_control.R all [WORKERS]\n",
            "Defaults: 7 workers, 100 randomizations per dataset.\n", sep = "")
    } else {
        script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
        analysis_directory <- dirname(normalizePath(script_file))
        selected <- if (arguments[[1L]] == "all") negative_control_datasets else arguments[[1L]]
        stopifnot(all(selected %in% negative_control_datasets), length(arguments) <= 2L)
        workers <- if (length(arguments) == 2L) as.integer(arguments[[2L]]) else 7L
        for (dataset in selected) {
            run_negative_control(file.path(analysis_directory, dataset), workers = workers)
        }
    }
}
