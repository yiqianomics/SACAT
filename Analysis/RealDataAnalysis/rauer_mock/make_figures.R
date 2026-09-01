# Create the Rauer mock-community figure and summary tables.

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
    matches <- candidates[file.exists(file.path(candidates, "make_figures.R"))]
    if (length(matches) != 1L) {
        stop("Could not identify the rauer_mock directory.", call. = FALSE)
    }
    matches[[1L]]
}

required_packages <- c("ggplot2", "patchwork", "scales")
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
table_directory <- file.path(dataset_directory, "table")
figure_directory <- file.path(dataset_directory, "figs")
input_file <- file.path(
    dataset_directory, "processed", "rauer_mock_dasra_input.rds"
)
nonnull_file <- file.path(table_directory, "rauer_mock_nonnull_results.csv")
null_file <- file.path(table_directory, "rauer_mock_null_results.csv")
required_files <- c(input_file, nonnull_file, null_file)
if (any(!file.exists(required_files))) {
    stop("Prepared input and both completed result tables are required.",
         call. = FALSE)
}
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

prepared <- readRDS(input_file)
nonnull_results <- utils::read.csv(
    nonnull_file, check.names = FALSE, stringsAsFactors = FALSE
)
null_results <- utils::read.csv(
    null_file, check.names = FALSE, stringsAsFactors = FALSE
)

depth_levels <- c("native", "1000", "750", "500", "250", "100", "50")
depth_positions <- stats::setNames(seq_along(depth_levels) - 1L, depth_levels)
depth_labels <- c(
    native = "Native", `1000` = "1,000", `750` = "750", `500` = "500",
    `250` = "250", `100` = "100", `50` = "50"
)
method_order <- c(
    "DASRA structural absence",
    "DASRA present-conditional abundance",
    "MaAsLin3 prevalence",
    "MaAsLin3 abundance",
    "ZINQ prevalence",
    "ZINQ abundance"
)
taxa <- as.character(prepared$taxon_metadata$taxon)

validate_results <- function(results, design, allocations) {
    required_columns <- c(
        "design", "allocation_id", "depth", "taxon", "target_set",
        "method", "available", "reason", "p_value", "q_value",
        "significant"
    )
    if (length(setdiff(required_columns, names(results))) ||
        nrow(results) != allocations * length(depth_levels) *
            length(method_order) * length(taxa) ||
        !identical(unique(results$design), design) ||
        !setequal(results$depth, depth_levels) ||
        !setequal(results$method, method_order) ||
        !setequal(results$taxon, taxa) ||
        anyDuplicated(results[
            c("allocation_id", "depth", "method", "taxon")
        ])) {
        stop(sprintf("The %s result table is incomplete.", design),
             call. = FALSE)
    }
    block_indices <- split(
        seq_len(nrow(results)),
        interaction(
            results$allocation_id,
            results$depth,
            results$method,
            drop = TRUE
        )
    )
    adjustment_matches <- vapply(block_indices, function(indices) {
        expected <- stats::p.adjust(
            ifelse(
                results$available[indices] &
                    is.finite(results$p_value[indices]),
                results$p_value[indices],
                1
            ),
            method = "BH"
        )
        max(abs(expected - results$q_value[indices])) < 1e-12
    }, logical(1))
    if (!all(adjustment_matches)) {
        stop(sprintf("BH adjustment is inconsistent in the %s results.",
                     design), call. = FALSE)
    }
    expected_significant <- results$available & results$q_value <= 0.05
    if (!identical(as.logical(results$significant), expected_significant)) {
        stop(sprintf("Discovery indicators are inconsistent in the %s results.",
                     design), call. = FALSE)
    }
    if (anyNA(results$reason) || any(!nzchar(results$reason)) ||
        any(results$available != (results$reason == "available"))) {
        stop(sprintf("Availability reasons are inconsistent in the %s results.",
                     design), call. = FALSE)
    }
}

validate_results(nonnull_results, "nonnull", 52L)
validate_results(null_results, "null", 111L)

method_colours <- c(
    "DASRA" = "#0072B2",
    "MaAsLin3" = "#D55E00",
    "ZINQ" = "#009E73"
)
method_shapes <- c("DASRA" = 16, "MaAsLin3" = 17, "ZINQ" = 15)

base_theme <- ggplot2::theme_classic(
    base_family = "Helvetica", base_size = 9
) +
    ggplot2::theme(
        plot.title = ggplot2::element_text(size = 10, hjust = 0),
        axis.title = ggplot2::element_text(size = 9),
        axis.text = ggplot2::element_text(size = 8, colour = "#222222"),
        axis.line = ggplot2::element_line(
            linewidth = 0.45, colour = "#333333"
        ),
        axis.ticks = ggplot2::element_line(
            linewidth = 0.4, colour = "#333333"
        ),
        axis.ticks.length = grid::unit(2.5, "pt"),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.width = grid::unit(17, "pt"),
        legend.spacing.x = grid::unit(4, "pt"),
        plot.margin = ggplot2::margin(7, 8, 5, 7)
    )

add_horizontal_guides <- function(breaks) {
    ggplot2::geom_hline(
        yintercept = breaks,
        colour = "#E5E5E5",
        linewidth = 0.35
    )
}

percent_scale <- function() {
    ggplot2::scale_y_continuous(
        limits = c(0, 103),
        breaks = c(0, 25, 50, 75, 100),
        labels = c(0, 25, 50, 75, 100),
        expand = ggplot2::expansion(mult = c(0, 0.01))
    )
}

# Native library sizes

library_data <- prepared$metadata
library_data$source <- ifelse(
    library_data$source_standard == "D6300", "Even mock", "Spike-in mock"
)
library_data$input_load <- ifelse(
    library_data$input_cells == "100000",
    "10⁵",
    ifelse(library_data$input_cells == "10000", "10⁴", "6×10³")
)
library_data$source <- factor(
    library_data$source, levels = c("Even mock", "Spike-in mock")
)
library_data$input_load <- factor(
    library_data$input_load, levels = c("10⁵", "10⁴", "6×10³")
)
if (nrow(library_data) != 32L ||
    !identical(as.numeric(stats::median(library_data$library_size)), 21398) ||
    !identical(as.numeric(range(library_data$library_size)), c(1176, 37683))) {
    stop("Native library sizes do not match the prepared benchmark.",
         call. = FALSE)
}

source_colours <- c("Even mock" = "#0072B2", "Spike-in mock" = "#6F4C9B")
library_panel <- ggplot2::ggplot(
    library_data,
    ggplot2::aes(x = input_load, y = library_size, colour = source)
) +
    ggplot2::geom_boxplot(
        width = 0.42,
        outlier.shape = NA,
        fill = "white",
        linewidth = 0.55
    ) +
    ggplot2::geom_jitter(
        position = ggplot2::position_jitter(
            width = 0.095,
            height = 0,
            seed = 20260831L
        ),
        size = 1.65,
        alpha = 0.88,
        stroke = 0
    ) +
    ggplot2::facet_grid(
        cols = ggplot2::vars(source), scales = "free_x", space = "free_x"
    ) +
    ggplot2::scale_colour_manual(values = source_colours, guide = "none") +
    ggplot2::scale_y_log10(
        breaks = c(1000, 3000, 10000, 30000),
        labels = c("1,000", "3,000", "10,000", "30,000"),
        limits = c(1000, 42000),
        expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::labs(
        title = "A   Native library sizes",
        x = "Input cells",
        y = "Reads per sample\n(log scale)"
    ) +
    base_theme +
    ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
            size = 8.5,
            face = "bold",
            colour = "#222222",
            margin = ggplot2::margin(b = 3)
        ),
        panel.grid.major.y = ggplot2::element_line(
            colour = "#E5E5E5", linewidth = 0.35
        ),
        axis.line.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 7.7),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
        plot.title = ggplot2::element_text(size = 9.5, face = "bold"),
        plot.title.position = "plot",
        plot.margin = ggplot2::margin(3, 5, 3, 4)
    )

# Complete-null familywise Type I error

null_discovery <- stats::aggregate(
    significant ~ allocation_id + depth + method,
    null_results,
    max
)
names(null_discovery)[names(null_discovery) == "significant"] <-
    "any_discovery"
null_availability <- stats::aggregate(
    available ~ allocation_id + depth + method,
    null_results,
    mean
)
null_blocks <- merge(
    null_discovery,
    null_availability,
    by = c("allocation_id", "depth", "method")
)
null_operating <- stats::aggregate(
    cbind(any_discovery, available) ~ depth + method,
    null_blocks,
    mean
)
null_definitions <- data.frame(
    method = method_order,
    method_component = c(
        "DASRA structural",
        "DASRA abundance",
        "MaAsLin3 prevalence",
        "MaAsLin3 abundance",
        "ZINQ prevalence",
        "ZINQ abundance"
    ),
    family = c(
        "Absence-oriented", "Abundance-oriented",
        "Absence-oriented", "Abundance-oriented",
        "Absence-oriented", "Abundance-oriented"
    ),
    stringsAsFactors = FALSE
)
null_operating <- merge(null_operating, null_definitions, by = "method")
null_operating$fwer_percent <- 100 * null_operating$any_discovery
null_operating$availability_percent <- 100 * null_operating$available
null_operating$value_label <- ifelse(
    null_operating$method == "DASRA present-conditional abundance",
    sprintf(
        "%.1f (%d%%)",
        null_operating$fwer_percent,
        round(null_operating$availability_percent)
    ),
    sprintf("%.1f", null_operating$fwer_percent)
)
null_operating$text_colour <- ifelse(
    null_operating$fwer_percent >= 8, "white", "black"
)
null_operating$depth <- factor(
    null_operating$depth, levels = depth_levels
)
null_operating$method_component <- factor(
    null_operating$method_component,
    levels = c(
        "ZINQ abundance",
        "MaAsLin3 abundance",
        "DASRA abundance",
        "ZINQ prevalence",
        "MaAsLin3 prevalence",
        "DASRA structural"
    )
)

type1_panel <- ggplot2::ggplot(
    null_operating,
    ggplot2::aes(x = depth, y = method_component, fill = fwer_percent)
) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
    ggplot2::geom_tile(
        data = null_operating[null_operating$fwer_percent > 5, , drop = FALSE],
        fill = NA,
        colour = "#A33A2B",
        linewidth = 0.70
    ) +
    ggplot2::geom_text(
        ggplot2::aes(label = value_label, colour = text_colour),
        family = "Helvetica",
        size = 2.15,
        show.legend = FALSE
    ) +
    ggplot2::geom_hline(
        yintercept = 3.5, colour = "#666666", linewidth = 0.45
    ) +
    ggplot2::geom_vline(
        xintercept = 1.5, colour = "#BEBEBE", linewidth = 0.4
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_discrete(
        labels = unname(depth_labels), expand = c(0, 0)
    ) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_fill_gradientn(
        colours = c(
            "#F7FBFF", "#D5E8F2", "#9ECAE1",
            "#FEE0B6", "#EF8A62", "#B2182B"
        ),
        values = c(0, 0.25, 0.4999, 0.5001, 0.75, 1),
        limits = c(0, 10),
        oob = scales::squish
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::labs(
        title = "B   Complete-null familywise Type I error (%)",
        x = NULL,
        y = NULL
    ) +
    ggplot2::coord_fixed(ratio = 0.72, clip = "off") +
    ggplot2::theme_minimal(base_family = "Helvetica", base_size = 9) +
    ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(size = 7.2, colour = "#222222"),
        axis.text.y = ggplot2::element_text(size = 6.8, colour = "#222222"),
        axis.ticks = ggplot2::element_blank(),
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 9.5, face = "bold"),
        plot.title.position = "plot",
        plot.margin = ggplot2::margin(3, 3, 1, 3)
    )

# Target recovery and availability

recovery_definitions <- data.frame(
    target_set = rep(c("Structural targets", "Abundance targets"), each = 3L),
    method = c(
        "DASRA structural absence",
        "MaAsLin3 prevalence",
        "ZINQ prevalence",
        "DASRA present-conditional abundance",
        "MaAsLin3 abundance",
        "ZINQ abundance"
    ),
    display_method = rep(c("DASRA", "MaAsLin3", "ZINQ"), 2L),
    target = rep(c("Structural", "Abundance"), each = 3L),
    stringsAsFactors = FALSE
)
recovery_rows <- merge(
    nonnull_results,
    recovery_definitions,
    by = c("target_set", "method")
)
recovery_allocation <- stats::aggregate(
    cbind(significant, available) ~
        allocation_id + depth + target + display_method,
    recovery_rows,
    mean
)
recovery_summary <- stats::aggregate(
    cbind(significant, available) ~ depth + target + display_method,
    recovery_allocation,
    mean
)
recovery_summary$discovery_percent <- 100 * recovery_summary$significant
recovery_summary$availability_percent <- 100 * recovery_summary$available
recovery_summary$depth <- factor(
    recovery_summary$depth, levels = depth_levels
)
recovery_summary$display_method <- factor(
    recovery_summary$display_method, levels = names(method_colours)
)
recovery_summary$depth_position <- unname(
    depth_positions[as.character(recovery_summary$depth)]
)

make_recovery_panel <- function(target_name, title, show_y_title) {
    panel_data <- recovery_summary[
        recovery_summary$target == target_name, , drop = FALSE
    ]
    availability_data <- panel_data[
        as.character(panel_data$display_method) == "DASRA",
        ,
        drop = FALSE
    ]
    fixed_depth_data <- panel_data[
        as.character(panel_data$depth) != "native", , drop = FALSE
    ]
    native_transition_data <- panel_data[
        as.character(panel_data$depth) %in% c("native", "1000"),
        ,
        drop = FALSE
    ]

    ggplot2::ggplot(
        panel_data,
        ggplot2::aes(
            x = depth_position,
            y = discovery_percent,
            colour = display_method,
            shape = display_method,
            group = display_method
        )
    ) +
        ggplot2::annotate(
            "rect",
            xmin = -0.50,
            xmax = 0.50,
            ymin = -Inf,
            ymax = Inf,
            fill = "#F5F5F5",
            colour = NA
        ) +
        add_horizontal_guides(c(0, 25, 50, 75, 100)) +
        ggplot2::geom_vline(
            xintercept = 0.50, colour = "#C8C8C8", linewidth = 0.35
        ) +
        ggplot2::geom_line(
            data = availability_data,
            ggplot2::aes(
                x = depth_position,
                y = availability_percent,
                group = 1,
                linetype = "DASRA availability"
            ),
            inherit.aes = FALSE,
            colour = "#777777",
            linewidth = 0.70
        ) +
        ggplot2::geom_point(
            data = availability_data,
            ggplot2::aes(x = depth_position, y = availability_percent),
            inherit.aes = FALSE,
            shape = 21,
            fill = "white",
            colour = "#777777",
            size = 1.45,
            stroke = 0.55,
            show.legend = FALSE
        ) +
        ggplot2::geom_line(
            data = native_transition_data,
            linewidth = 0.65,
            linetype = "22",
            show.legend = FALSE
        ) +
        ggplot2::geom_line(data = fixed_depth_data, linewidth = 0.75) +
        ggplot2::geom_point(size = 2.15, stroke = 0.8) +
        ggplot2::scale_x_continuous(
            breaks = unname(depth_positions),
            labels = unname(depth_labels),
            limits = c(-0.50, 6.50),
            expand = ggplot2::expansion(mult = c(0, 0))
        ) +
        percent_scale() +
        ggplot2::scale_colour_manual(values = method_colours, drop = FALSE) +
        ggplot2::scale_shape_manual(values = method_shapes, drop = FALSE) +
        ggplot2::scale_linetype_manual(
            values = c("DASRA availability" = "22"), name = NULL
        ) +
        ggplot2::guides(
            colour = ggplot2::guide_legend(
                order = 1,
                override.aes = list(
                    shape = unname(method_shapes), linewidth = 0.75
                )
            ),
            shape = "none",
            linetype = ggplot2::guide_legend(
                order = 2,
                override.aes = list(colour = "#777777", linewidth = 0.70)
            )
        ) +
        ggplot2::labs(
            title = title,
            x = "Reads per sample",
            y = if (show_y_title) "Target recovery (%)" else NULL,
            colour = NULL,
            shape = NULL,
            linetype = NULL
        ) +
        base_theme +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 9.5, face = "bold"),
            plot.title.position = "plot",
            axis.title.y = ggplot2::element_text(
                margin = ggplot2::margin(r = 3)
            ),
            axis.title.x = ggplot2::element_text(
                margin = ggplot2::margin(t = 4)
            ),
            axis.text.x = ggplot2::element_text(size = 7.2),
            legend.position = "bottom",
            legend.margin = ggplot2::margin(t = 0, b = 0),
            legend.box.margin = ggplot2::margin(t = -2),
            plot.margin = ggplot2::margin(3, 5, 3, 4)
        )
}

structural_panel <- make_recovery_panel(
    "Structural", "C   Structural targets", TRUE
)
abundance_panel <- make_recovery_panel(
    "Abundance", "D   Abundance targets", FALSE
)

# Summary tables

component_output_names <- c(
    "DASRA structural absence" = "DASRA structural absence",
    "DASRA present-conditional abundance" = "DASRA abundance",
    "MaAsLin3 prevalence" = "MaAsLin3 prevalence",
    "MaAsLin3 abundance" = "MaAsLin3 abundance",
    "ZINQ prevalence" = "ZINQ prevalence",
    "ZINQ abundance" = "ZINQ abundance"
)
depth_method_name <- function(method) {
    if (method == "DASRA present-conditional abundance") {
        return("DASRA present-conditional abundance")
    }
    component_output_names[[method]]
}
native_rows <- nonnull_results[
    nonnull_results$depth == "native", , drop = FALSE
]
native_frequency <- stats::aggregate(
    significant ~ target_set + taxon + method,
    native_rows,
    mean
)
native_frequency$discovery_percent <- 100 * native_frequency$significant

make_native_row <- function(target_set, taxon) {
    output <- data.frame(
        target_set = target_set,
        taxon = taxon,
        stringsAsFactors = FALSE
    )
    for (method in names(component_output_names)) {
        if (taxon == "Mean across taxa") {
            values <- native_frequency$discovery_percent[
                native_frequency$target_set == target_set &
                    native_frequency$method == method
            ]
        } else {
            values <- native_frequency$discovery_percent[
                native_frequency$target_set == target_set &
                    native_frequency$taxon == taxon &
                    native_frequency$method == method
            ]
        }
        if (!length(values)) {
            stop("The native-depth summary is incomplete.", call. = FALSE)
        }
        output[[component_output_names[[method]]]] <- mean(values)
    }
    output
}

structural_taxa <- sort(prepared$taxon_metadata$taxon[
    prepared$taxon_metadata$target_set == "structural absence"
])
abundance_taxa <- sort(prepared$taxon_metadata$taxon[
    prepared$taxon_metadata$target_set == "present-conditional abundance"
])
native_taxon_summary <- do.call(rbind, c(
    lapply(
        c(structural_taxa, "Mean across taxa"),
        make_native_row,
        target_set = "Structural targets"
    ),
    lapply(
        c(abundance_taxa, "Mean across taxa"),
        make_native_row,
        target_set = "Abundance targets"
    )
))

calculate_fdr <- function(target_set, method, display_method) {
    rows <- nonnull_results[nonnull_results$method == method, , drop = FALSE]
    blocks <- split(
        rows,
        interaction(rows$allocation_id, rows$depth, drop = TRUE)
    )
    allocation_fdp <- do.call(rbind, lapply(blocks, function(block) {
        discoveries <- sum(block$significant)
        false_discoveries <- sum(
            block$significant & block$target_set != target_set
        )
        data.frame(
            allocation_id = block$allocation_id[[1L]],
            depth = block$depth[[1L]],
            fdp = if (discoveries) false_discoveries / discoveries else 0,
            stringsAsFactors = FALSE
        )
    }))
    output <- stats::aggregate(fdp ~ depth, allocation_fdp, mean)
    output$target_set <- target_set
    output$display_method <- display_method
    output
}

fdr_definitions <- data.frame(
    target_set = rep(c("Structural targets", "Abundance targets"), each = 3L),
    method = recovery_definitions$method,
    display_method = recovery_definitions$display_method,
    stringsAsFactors = FALSE
)
fdr_summary <- do.call(rbind, lapply(seq_len(nrow(fdr_definitions)), function(i) {
    calculate_fdr(
        fdr_definitions$target_set[[i]],
        fdr_definitions$method[[i]],
        fdr_definitions$display_method[[i]]
    )
}))
fdr_summary$fdp_percent <- 100 * fdr_summary$fdp

nonnull_availability <- stats::aggregate(
    available ~ depth + target_set + method,
    nonnull_results[nonnull_results$method %in% c(
        "DASRA structural absence",
        "DASRA present-conditional abundance"
    ), ],
    mean
)
nonnull_availability$availability_percent <-
    100 * nonnull_availability$available

depth_values <- function(data, value, keep) {
    selected <- data[keep, , drop = FALSE]
    selected <- selected[
        match(depth_levels, as.character(selected$depth)),
        ,
        drop = FALSE
    ]
    if (nrow(selected) != length(depth_levels) ||
        !identical(as.character(selected$depth), depth_levels)) {
        stop("A depth-summary series is incomplete.", call. = FALSE)
    }
    as.numeric(selected[[value]])
}

make_depth_row <- function(section, target_set, method_component, values) {
    output <- data.frame(
        section = section,
        target_set = target_set,
        method_component = method_component,
        stringsAsFactors = FALSE
    )
    display_columns <- c("Native", "1,000", "750", "500", "250", "100", "50")
    for (index in seq_along(values)) {
        output[[display_columns[[index]]]] <- values[[index]]
    }
    output
}

depth_rows <- list()
append_depth_row <- function(section, target_set, method_component, values) {
    depth_rows[[length(depth_rows) + 1L]] <<- make_depth_row(
        section, target_set, method_component, values
    )
}

for (target_name in c("Structural", "Abundance")) {
    target_set <- paste(target_name, "targets")
    definitions <- recovery_definitions[
        recovery_definitions$target == target_name, , drop = FALSE
    ]
    for (index in seq_len(nrow(definitions))) {
        append_depth_row(
            "Empirical power",
            target_set,
            depth_method_name(definitions$method[[index]]),
            depth_values(
                recovery_summary,
                "discovery_percent",
                recovery_summary$target == target_name &
                    as.character(recovery_summary$display_method) ==
                        definitions$display_method[[index]]
            )
        )
    }
}

for (target_set in c("Structural targets", "Abundance targets")) {
    definitions <- fdr_definitions[
        fdr_definitions$target_set == target_set, , drop = FALSE
    ]
    for (index in seq_len(nrow(definitions))) {
        append_depth_row(
            "Mechanism-specific empirical FDR",
            sub(" targets$", " family", target_set),
            depth_method_name(definitions$method[[index]]),
            depth_values(
                fdr_summary,
                "fdp_percent",
                fdr_summary$target_set == target_set &
                    fdr_summary$display_method ==
                        definitions$display_method[[index]]
            )
        )
    }
}

for (index in seq_len(nrow(null_definitions))) {
    definition <- null_definitions[index, ]
    append_depth_row(
        "Complete-null familywise Type I error",
        if (definition$family == "Absence-oriented") {
            "Structural family"
        } else {
            "Abundance family"
        },
        depth_method_name(definition$method),
        depth_values(
            null_operating,
            "fwer_percent",
            null_operating$method == definition$method
        )
    )
}

availability_order <- data.frame(
    method = c(
        "DASRA structural absence",
        "DASRA structural absence",
        "DASRA present-conditional abundance",
        "DASRA present-conditional abundance"
    ),
    target_set = c(
        "Structural targets",
        "Abundance targets",
        "Abundance targets",
        "Structural targets"
    ),
    component = c(
        "Structural-absence component",
        "Structural-absence component",
        "Abundance component",
        "Abundance component"
    ),
    stringsAsFactors = FALSE
)
for (index in seq_len(nrow(availability_order))) {
    definition <- availability_order[index, ]
    append_depth_row(
        "DASRA availability: nonnull allocations",
        definition$component,
        if (definition$target_set == "Structural targets") {
            "Structural target set (8 taxa)"
        } else {
            "Abundance target set (3 taxa)"
        },
        depth_values(
            nonnull_availability,
            "availability_percent",
            nonnull_availability$method == definition$method &
                nonnull_availability$target_set == definition$target_set
        )
    )
}

for (method in c(
    "DASRA structural absence", "DASRA present-conditional abundance"
)) {
    append_depth_row(
        "DASRA availability: complete-null allocations",
        "All 11 taxa",
        if (method == "DASRA structural absence") {
            "Structural-absence component"
        } else {
            "Abundance component"
        },
        depth_values(
            null_operating,
            "availability_percent",
            null_operating$method == method
        )
    )
}

depth_summary <- do.call(rbind, depth_rows)
native_numeric <- vapply(native_taxon_summary, is.numeric, logical(1))
native_taxon_summary[native_numeric] <- lapply(
    native_taxon_summary[native_numeric], round, digits = 2
)
depth_numeric <- vapply(depth_summary, is.numeric, logical(1))
depth_summary[depth_numeric] <- lapply(
    depth_summary[depth_numeric], round, digits = 2
)
utils::write.csv(
    native_taxon_summary,
    file.path(table_directory, "rauer_mock_native_taxon_summary.csv"),
    row.names = FALSE,
    quote = TRUE
)
utils::write.csv(
    depth_summary,
    file.path(table_directory, "rauer_mock_depth_summary.csv"),
    row.names = FALSE,
    quote = TRUE
)

top_row <- patchwork::free(library_panel, side = "l") + type1_panel +
    patchwork::plot_layout(ncol = 2L, widths = c(0.74, 1.26))
bottom_row <- structural_panel + abundance_panel +
    patchwork::plot_layout(ncol = 2L, guides = "collect")
bottom_row <- bottom_row & ggplot2::theme(legend.position = "bottom")
main_figure <- top_row / bottom_row +
    patchwork::plot_layout(nrow = 2L, heights = c(0.96, 1.04))

ggplot2::ggsave(
    file.path(figure_directory, "rauer_mock_main.pdf"),
    main_figure,
    width = 7.2,
    height = 5.8,
    units = "in",
    device = grDevices::cairo_pdf
)
