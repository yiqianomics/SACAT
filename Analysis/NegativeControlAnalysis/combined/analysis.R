#!/usr/bin/env Rscript

# Summarize SACAT raw-p rejection rates with fixed native counts and library sizes.
run_negative_control_summary <- function(root, make_figure = TRUE) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package data.table is required.", call. = FALSE)
    }
    root <- normalizePath(root, mustWork = TRUE)
    output_directory <- file.path(root, "combined")
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    components <- c("SA", "RA", "omnibus")
    alpha <- 0.05
    replicates <- 100L
    taxa_per_replicate <- 30L

    input_paths <- file.path(list.dirs(root, recursive = FALSE),
                             "analysis_input.rds")
    input_paths <- input_paths[file.exists(input_paths)]
    if (length(input_paths) != 20L) {
        stop("Expected 20 dataset directories with native-count inputs.",
             call. = FALSE)
    }
    inputs <- lapply(input_paths, readRDS)
    catalog <- data.table::rbindlist(lapply(seq_along(inputs), function(i) {
        x <- inputs[[i]]
        counts <- x$counts
        depth <- x$metadata$library_size
        stopifnot(identical(x$dataset_id, basename(dirname(input_paths[[i]]))),
                  length(x$label) == 1L, !is.na(x$label), nzchar(x$label),
                  length(x$dataset_index) == 1L, is.finite(x$dataset_index),
                  identical(dim(counts), c(30L, 200L)), is.numeric(counts),
                  identical(rownames(counts), x$evaluation_taxa),
                  !anyDuplicated(x$evaluation_taxa),
                  !anyDuplicated(colnames(counts)),
                  identical(colnames(counts), rownames(x$metadata)),
                  all(is.finite(counts)), all(counts >= 0),
                  all(counts == round(counts)), length(depth) == 200L,
                  all(is.finite(depth)), all(depth > 0),
                  all(depth == round(depth)), all(colSums(counts) <= depth))
        data.table::data.table(
            dataset = x$dataset_id, label = x$label,
            dataset_index = as.integer(x$dataset_index),
            samples = ncol(counts), genera = nrow(counts),
            planned_randomizations = replicates,
            minimum_library_size = min(depth),
            median_library_size = median(depth),
            mean_library_size = mean(depth),
            maximum_library_size = max(depth),
            sample_selection_seed = if (is.null(x$sample_selection_seed)) {
                NA_integer_
            } else x$sample_selection_seed
        )
    }))
    stopifnot(!anyDuplicated(catalog$dataset),
              !anyDuplicated(catalog$label),
              !anyDuplicated(catalog$dataset_index))
    names(inputs) <- vapply(inputs, function(x) x$dataset_id, character(1))
    data.table::setorder(catalog, dataset_index)

    tables <- lapply(catalog$dataset, function(dataset) {
        input <- inputs[[dataset]]
        path <- file.path(root, dataset, "results", "taxon_pvalues.csv")
        if (!file.exists(path)) stop("Missing results: ", path, call. = FALSE)
        result <- data.table::fread(path, showProgress = FALSE)
        required <- c("dataset", "label", "replicate", "taxon", "component",
                      "raw_p", "available", "analysis_p", "rejected")
        stopifnot(all(required %in% names(result)),
                  nrow(result) == replicates * taxa_per_replicate * 3L,
                  identical(unique(result$dataset), dataset),
                  identical(unique(result$label), input$label),
                  setequal(unique(result$replicate), seq_len(replicates)),
                  setequal(unique(result$component), components),
                  setequal(unique(result$taxon), input$evaluation_taxa),
                  !anyNA(result[, c("replicate", "taxon", "component"),
                                with = FALSE]),
                  !anyDuplicated(result[, c("replicate", "taxon", "component"),
                                        with = FALSE]),
                  is.logical(result$available), !anyNA(result$available))
        valid_p <- is.finite(result$raw_p) &
            result$raw_p >= 0 & result$raw_p <= 1
        stopifnot(all(valid_p[result$available]))
        analysis_p <- ifelse(result$available, result$raw_p, 1)
        rejected <- analysis_p < alpha
        stopifnot(all(is.finite(result$analysis_p)),
                  all(abs(result$analysis_p - analysis_p) < 1e-12),
                  !anyNA(result$rejected),
                  all(result$rejected == rejected))
        # Derive decisions from raw p-values, retaining every panel member.
        data.table::set(result, j = "analysis_p", value = analysis_p)
        data.table::set(result, j = "rejected", value = rejected)
        result[, dataset_index := input$dataset_index]
        result[, method := "SACAT"]
        result
    })
    taxon <- data.table::rbindlist(tables, use.names = TRUE, fill = TRUE)
    group_columns <- c("dataset", "label", "dataset_index", "method", "component")
    replicate_summary <- taxon[, list(
        taxa = .N, rejected_taxa = sum(rejected),
        type1_error = sum(rejected) / taxa_per_replicate,
        available_taxa = sum(available), unavailable_taxa = sum(!available),
        availability = mean(available)
    ), by = c(group_columns, "replicate")]
    stopifnot(nrow(replicate_summary) == 20L * replicates * 3L,
              all(replicate_summary$taxa == taxa_per_replicate))
    dataset_summary <- replicate_summary[, list(
        replications = .N, taxa_per_replication = taxa_per_replicate,
        type1_error = mean(type1_error),
        monte_carlo_sd = stats::sd(type1_error),
        monte_carlo_se = stats::sd(type1_error) / sqrt(.N),
        availability = mean(availability)
    ), by = group_columns]
    stopifnot(all(dataset_summary$replications == replicates))
    dataset_summary[, mc_lower_95 := pmax(0, type1_error - 1.96 * monte_carlo_se)]
    dataset_summary[, mc_upper_95 := pmin(1, type1_error + 1.96 * monte_carlo_se)]
    overall_summary <- dataset_summary[, list(
        datasets = .N, replications_per_dataset = replicates,
        taxa_per_replication = taxa_per_replicate,
        type1_error = mean(type1_error),
        minimum_type1_error = min(type1_error),
        maximum_type1_error = max(type1_error),
        availability = mean(availability),
        minimum_availability = min(availability),
        maximum_availability = max(availability)
    ), by = c("method", "component")]
    availability_summary <- replicate_summary[, list(
        replications = .N, tests = sum(taxa),
        available_tests = sum(available_taxa),
        unavailable_tests = sum(unavailable_taxa),
        availability = sum(available_taxa) / sum(taxa)
    ), by = group_columns]

    order_components <- function(table) {
        order(table$dataset_index, match(table$component, components),
              if ("replicate" %in% names(table)) table$replicate else rep(0L, nrow(table)))
    }
    dataset_summary <- dataset_summary[order_components(dataset_summary)]
    replicate_summary <- replicate_summary[order_components(replicate_summary)]
    availability_summary <- availability_summary[order_components(availability_summary)]
    overall_summary <- overall_summary[match(components, component)]
    outputs <- list(
        dataset_type1_error = dataset_summary,
        replicate_type1_error = replicate_summary,
        overall_type1_error = overall_summary,
        availability_summary = availability_summary,
        dataset_design = catalog
    )
    for (name in names(outputs)) {
        data.table::fwrite(outputs[[name]],
                          file.path(output_directory, paste0(name, ".csv")))
    }

    figure_paths <- NULL
    if (isTRUE(make_figure)) {
        if (!requireNamespace("ggplot2", quietly = TRUE) ||
            !requireNamespace("scales", quietly = TRUE)) {
            stop("Packages ggplot2 and scales are required for the figure.",
                 call. = FALSE)
        }
        plot_data <- as.data.frame(dataset_summary)
        plot_data$label <- factor(plot_data$label, levels = rev(catalog$label))
        plot_data$component <- factor(plot_data$component, levels = components)
        strip_labels <- c(SA = "Structural absence (SA)",
                          RA = "Relative abundance (RA)", omnibus = "Omnibus")
        x_max <- max(0.08, ceiling(max(plot_data$mc_upper_95) * 50) / 50)
        plot <- ggplot2::ggplot(plot_data,
            ggplot2::aes(x = type1_error, y = label, colour = component)) +
            ggplot2::geom_vline(xintercept = alpha, linetype = "dashed",
                               colour = "grey45", linewidth = 0.5) +
            ggplot2::geom_segment(ggplot2::aes(x = mc_lower_95,
                xend = mc_upper_95, yend = label), linewidth = 0.65) +
            ggplot2::geom_point(size = 2.2) +
            ggplot2::facet_grid(cols = ggplot2::vars(component),
                labeller = ggplot2::labeller(component = strip_labels)) +
            ggplot2::scale_colour_manual(
                values = c(SA = "#0072B2", RA = "#009E73", omnibus = "#D55E00"),
                guide = "none") +
            ggplot2::scale_x_continuous(limits = c(0, x_max),
                breaks = seq(0, x_max, by = 0.02),
                labels = scales::label_percent(accuracy = 1),
                expand = ggplot2::expansion(mult = c(0.025, 0.035))) +
            ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.7)) +
            ggplot2::labs(
                title = "Native-count negative controls",
                subtitle = "SACAT | 20 cohorts | 100 label randomizations per cohort | 30 genera per replication",
                x = "Raw-p Type I error", y = NULL,
                caption = paste(
                    "Points: mean rejection rate at raw p < 0.05; bars: 95% Monte Carlo intervals.",
                    "Unavailable tests are retained as p = 1; dashed line: nominal 5%.",
                    sep = "\n")) +
            ggplot2::theme_bw(base_size = 10, base_family = "sans") +
            ggplot2::theme(
                panel.grid.major.y = ggplot2::element_line(colour = "grey92"),
                panel.grid.major.x = ggplot2::element_line(colour = "grey94"),
                panel.grid.minor = ggplot2::element_blank(),
                panel.border = ggplot2::element_rect(colour = "grey75"),
                panel.spacing.x = grid::unit(7, "mm"),
                strip.background = ggplot2::element_rect(fill = "grey96",
                                                        colour = "grey75"),
                strip.text = ggplot2::element_text(face = "bold", size = 10),
                axis.text.y = ggplot2::element_text(size = 9, colour = "grey15"),
                axis.text.x = ggplot2::element_text(colour = "grey20"),
                axis.ticks.y = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(face = "bold", size = 14),
                plot.subtitle = ggplot2::element_text(size = 10,
                    margin = ggplot2::margin(b = 10)),
                plot.caption = ggplot2::element_text(hjust = 0, size = 9,
                    margin = ggplot2::margin(t = 10)),
                plot.margin = ggplot2::margin(12, 14, 10, 10))
        pdf_path <- file.path(output_directory, "negative_control_type1_error.pdf")
        ggplot2::ggsave(pdf_path, plot, device = grDevices::cairo_pdf,
                       width = 325, height = 210, units = "mm")
        figure_paths <- c(pdf = pdf_path)
    }
    cat("Summarized 20 cohorts, 2,000 replications, and 180,000 raw-p tests.\n")
    if (length(figure_paths)) cat("PDF:", figure_paths[["pdf"]], "\n")
    invisible(list(tables = outputs, figures = figure_paths))
}

if (sys.nframe() == 0L) {
    arguments <- commandArgs(trailingOnly = TRUE)
    file_argument <- grep("^--file=", commandArgs(), value = TRUE)
    default_root <- dirname(dirname(normalizePath(
        sub("^--file=", "", file_argument[[1L]]), mustWork = TRUE)))
    run_negative_control_summary(if (length(arguments)) arguments[[1L]] else default_root)
}
