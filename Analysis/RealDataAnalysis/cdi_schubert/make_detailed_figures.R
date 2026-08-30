#!/usr/bin/env Rscript

# Create the detailed CDI component and method-intersection figures.

required_packages <- c("ggplot2", "patchwork", "scales", "tidyr")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
    stop(
        "Install the required packages: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument)) {
    sub("^--file=", "", script_argument[[1L]])
} else {
    "Analysis/RealDataAnalysis/cdi_schubert/make_detailed_figures.R"
}
analysis_directory <- normalizePath(dirname(script_path), mustWork = TRUE)
table_directory <- file.path(analysis_directory, "table")
figure_directory <- file.path(analysis_directory, "figs")
input_path <- file.path(
    analysis_directory, "processed", "cdi_schubert_dasra_input.rds"
)

result_path <- file.path(
    table_directory, "schubert_cdi_method_results_all_taxa.csv"
)
required_files <- c(result_path, input_path)
if (!all(file.exists(required_files))) {
    stop(
        "Prepare the CDI input and complete the real-data analysis first.",
        call. = FALSE
    )
}

results <- utils::read.csv(result_path, stringsAsFactors = FALSE)
input <- readRDS(input_path)

analysis_method_order <- c(
    "DASRA structural absence",
    "DASRA present-conditional abundance",
    "DASRA combined",
    "MaAsLin3 prevalence",
    "MaAsLin3 abundance",
    "MaAsLin3 combined",
    "ZINQ prevalence",
    "ZINQ abundance",
    "ZINQ combined",
    "ANCOM-BC2",
    "LinDA",
    "corncob",
    "edgeR",
    "DESeq2",
    "metagenomeSeq"
)
display_method_sources <- c(
    "DASRA structural absence" = "DASRA structural absence",
    "DASRA present-conditional abundance" =
        "DASRA present-conditional abundance",
    "MaAsLin3 combined" = "MaAsLin3",
    "ZINQ combined" = "ZINQ",
    "ANCOM-BC2" = "ANCOM-BC2",
    "LinDA" = "LinDA",
    "corncob" = "corncob",
    "edgeR" = "edgeR",
    "DESeq2" = "DESeq2",
    "metagenomeSeq" = "metagenomeSeq"
)
display_method_order <- unname(display_method_sources)
component_taxa <- c("Coprococcus", "Enterococcus", "Faecalibacterium")
intersection_taxa <- component_taxa
focus_colors <- c(
    "Coprococcus" = "#0072B2",
    "Enterococcus" = "#D55E00",
    "Faecalibacterium" = "#009E73"
)
group_colors <- c("Healthy" = "#4C78A8", "CDI" = "#E45756")

if (!identical(unique(results$method), analysis_method_order) ||
    !all(intersection_taxa %in% results$taxon) ||
    !all(component_taxa %in% rownames(input$counts))) {
    stop("The completed CDI results do not match the expected analysis.",
         call. = FALSE)
}

wilson_interval <- function(successes, total, level = 0.95) {
    z <- stats::qnorm(1 - (1 - level) / 2)
    proportion <- successes / total
    denominator <- 1 + z ^ 2 / total
    center <- (proportion + z ^ 2 / (2 * total)) / denominator
    half_width <- z * sqrt(
        proportion * (1 - proportion) / total + z ^ 2 / (4 * total ^ 2)
    ) / denominator
    c(lower = center - half_width, upper = center + half_width)
}

bh_z_gate <- function(p_value, q_value, alpha = 0.05) {
    tested <- is.finite(p_value)
    rejected <- tested & is.finite(q_value) & q_value <= alpha
    m <- sum(tested)
    k <- sum(rejected)
    if (!m || !k) return(NA_real_)
    critical_p <- alpha * k / m
    stats::qnorm(critical_p / 2, lower.tail = FALSE)
}

base_theme <- ggplot2::theme_classic(base_size = 9.5) +
    ggplot2::theme(
        axis.title = ggplot2::element_text(size = 9),
        axis.text = ggplot2::element_text(size = 8, color = "grey15"),
        plot.title = ggplot2::element_text(size = 9, face = "italic"),
        plot.subtitle = ggplot2::element_text(size = 7.5, color = "grey25"),
        plot.tag = ggplot2::element_text(size = 11, face = "bold"),
        plot.margin = ggplot2::margin(5, 5, 5, 5)
    )

# Method intersections ----------------------------------------------------

display_results <- results[
    results$method %in% names(display_method_sources),
    c("taxon", "method", "significant")
]
display_results$display_method <- unname(
    display_method_sources[display_results$method]
)
display_signature <- tidyr::pivot_wider(
    display_results[, c("taxon", "display_method", "significant")],
    names_from = "display_method",
    values_from = "significant",
    values_fill = FALSE
)
display_signature <- display_signature[
    , c("taxon", display_method_order), drop = FALSE
]
significant_matrix <- as.matrix(
    display_signature[, display_method_order, drop = FALSE]
)
storage.mode(significant_matrix) <- "logical"
has_discovery <- rowSums(significant_matrix) > 0L
display_patterns <- apply(
    significant_matrix[has_discovery, , drop = FALSE],
    1L,
    function(row) paste(display_method_order[row], collapse = " | ")
)
display_membership <- data.frame(
    taxon = display_signature$taxon[has_discovery],
    methods = unname(display_patterns),
    stringsAsFactors = FALSE
)
pattern_table <- as.data.frame(
    table(display_membership$methods), stringsAsFactors = FALSE
)
names(pattern_table) <- c("methods", "n")
pattern_table <- pattern_table[
    order(-pattern_table$n, pattern_table$methods), , drop = FALSE
]
pattern_table$intersection <- paste0("I", seq_len(nrow(pattern_table)))
display_membership$intersection <- pattern_table$intersection[
    match(display_membership$methods, pattern_table$methods)
]
displayed_intersections <- pattern_table[
    , c("intersection", "n", "methods")
]
intersection_levels <- displayed_intersections$intersection
method_intersections_export <- displayed_intersections
method_membership_export <- display_membership[
    order(
        match(display_membership$intersection, intersection_levels),
        display_membership$taxon,
        method = "radix"
    ),
    c("intersection", "taxon", "methods"),
    drop = FALSE
]
displayed_intersections$intersection <- factor(
    displayed_intersections$intersection, levels = intersection_levels
)

positions <- data.frame(
    method = display_method_order,
    y = rev(seq_along(display_method_order)),
    stringsAsFactors = FALSE
)
all_points <- merge(
    expand.grid(
        intersection = intersection_levels,
        method = display_method_order,
        stringsAsFactors = FALSE
    ),
    positions,
    by = "method",
    sort = FALSE
)
all_points$intersection <- factor(
    all_points$intersection, levels = intersection_levels
)

active_points <- do.call(rbind, lapply(
    seq_len(nrow(displayed_intersections)),
    function(index) {
        methods <- strsplit(
            displayed_intersections$methods[[index]], " | ", fixed = TRUE
        )[[1L]]
        data.frame(
            intersection = displayed_intersections$intersection[[index]],
            method = methods,
            stringsAsFactors = FALSE
        )
    }
))
active_points <- merge(active_points, positions, by = "method", sort = FALSE)
active_points$intersection <- factor(
    active_points$intersection, levels = intersection_levels
)
segment_data <- stats::aggregate(
    y ~ intersection, active_points,
    function(values) c(minimum = min(values), maximum = max(values))
)
segment_data <- data.frame(
    intersection = segment_data$intersection,
    minimum = segment_data$y[, "minimum"],
    maximum = segment_data$y[, "maximum"]
)

row_background <- data.frame(
    y = positions$y,
    fill = rep(c("grey98", "grey94"), length.out = nrow(positions))
)

focus_membership <- display_membership[
    display_membership$taxon %in% intersection_taxa,
]
focus_membership <- merge(
    focus_membership,
    displayed_intersections[, c("intersection", "n")],
    by = "intersection"
)
if (nrow(focus_membership) != length(intersection_taxa)) {
    stop("One or more highlighted taxa are absent from the displayed intersections.",
         call. = FALSE)
}
focus_membership$intersection <- factor(
    focus_membership$intersection, levels = intersection_levels
)

intersection_index <- stats::setNames(
    seq_along(intersection_levels), intersection_levels
)
displayed_intersections$intersection_index <- unname(
    intersection_index[as.character(displayed_intersections$intersection)]
)
all_points$intersection_index <- unname(
    intersection_index[as.character(all_points$intersection)]
)
active_points$intersection_index <- unname(
    intersection_index[as.character(active_points$intersection)]
)
segment_data$intersection_index <- unname(
    intersection_index[as.character(segment_data$intersection)]
)
focus_membership$intersection_index <- unname(
    intersection_index[as.character(focus_membership$intersection)]
)

annotated_intersections <- c("I2", "I4", "I5")
annotated_members <- display_membership[
    display_membership$intersection %in% annotated_intersections,
    c("intersection", "taxon")
]
annotated_members <- annotated_members[order(
    match(annotated_members$intersection, annotated_intersections),
    annotated_members$taxon,
    method = "radix"
), ]
expected_members <- list(
    I2 = c(
        "Enterococcus", "Lactobacillus", "Prevotella", "Streptococcus",
        "Unclassified_Enterobacteriaceae", "Veillonella"
    ),
    I4 = c(
        "Butyricicoccus", "Coprococcus", "Parabacteroides", "Ruminococcus2"
    ),
    I5 = c(
        "Faecalibacterium", "Lachnospiracea_incertae_sedis",
        "Unclassified_Clostridiales"
    )
)
observed_members <- split(
    annotated_members$taxon,
    factor(annotated_members$intersection, levels = annotated_intersections)
)
if (!identical(observed_members, expected_members)) {
    stop("The annotated intersections do not match the completed CDI results.",
         call. = FALSE)
}

member_y <- list(
    I2 = c(5.70, 4.75, 3.80, 2.85, 1.90, 0.95),
    I4 = c(5.70, 4.55, 3.40, 2.25),
    I5 = c(5.70, 4.30, 2.65)
)
annotated_members$y <- unlist(
    lapply(annotated_intersections, function(intersection) {
        member_y[[intersection]]
    }),
    use.names = FALSE
)
annotated_members$label <- gsub(
    "_", " ", annotated_members$taxon, fixed = TRUE
)
annotated_members$label[
    annotated_members$taxon == "Lachnospiracea_incertae_sedis"
] <- "Lachnospiracea\nincertae sedis"
annotated_members$highlighted <- annotated_members$taxon %in% intersection_taxa
annotated_members$italicize_label <- !grepl(
    "^Unclassified_|incertae_sedis$|^Ruminococcus2$",
    annotated_members$taxon
)
ordinary_italic_labels <-
    !annotated_members$highlighted & annotated_members$italicize_label
ordinary_plain_labels <-
    !annotated_members$highlighted & !annotated_members$italicize_label
key_columns <- c(I2 = 0.02, I4 = 0.36, I5 = 0.68)
annotated_members$text_x <- unname(
    key_columns[annotated_members$intersection]
) + 0.035
annotated_members$marker_x <- unname(
    key_columns[annotated_members$intersection]
)
annotation_headers <- data.frame(
    intersection = annotated_intersections,
    label = sprintf(
        "%s (n = %d)",
        annotated_intersections,
        lengths(expected_members)
    ),
    x = unname(key_columns),
    y = rep(6.90, length(annotated_intersections)),
    taxon = c("Enterococcus", "Coprococcus", "Faecalibacterium"),
    stringsAsFactors = FALSE
)
arrow_data <- merge(
    focus_membership[, c(
        "taxon", "intersection", "intersection_index", "n"
    )],
    data.frame(
        intersection = annotated_intersections,
        label_x = c(2.55, 4.55, 5.65),
        label_y = c(9.2, 7.3, 5.8),
        stringsAsFactors = FALSE
    ),
    by = "intersection",
    sort = FALSE
)
arrow_data <- arrow_data[
    match(annotated_intersections, arrow_data$intersection),
]
arrow_data$arrow_x <- arrow_data$intersection_index + 0.23
arrow_data$arrow_y <- arrow_data$n + 0.55

intersection_bar <- ggplot2::ggplot(
    displayed_intersections,
    ggplot2::aes(x = intersection_index, y = n)
) +
    ggplot2::geom_col(width = 0.66, fill = "grey25") +
    ggplot2::geom_text(
        data = displayed_intersections[
            displayed_intersections$n >= 2 &
                !as.character(displayed_intersections$intersection) %in%
                    annotated_intersections,
        ],
        ggplot2::aes(
            x = intersection_index, y = n, label = n
        ),
        vjust = -0.25,
        size = 2.6, color = "grey15"
    ) +
    ggplot2::geom_col(
        data = focus_membership,
        ggplot2::aes(x = intersection_index, y = n, color = taxon),
        fill = NA, width = 0.66, linewidth = 0.8,
        inherit.aes = FALSE
    ) +
    ggplot2::geom_curve(
        data = arrow_data,
        ggplot2::aes(
            x = label_x, y = label_y - 0.35,
            xend = arrow_x, yend = arrow_y,
            color = taxon
        ),
        curvature = 0.12, linewidth = 0.55,
        arrow = grid::arrow(
            length = grid::unit(1.7, "mm"), type = "closed"
        ),
        inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
        data = arrow_data,
        ggplot2::aes(
            x = label_x, y = label_y, label = intersection, color = taxon
        ),
        size = 2.8, fontface = "bold", inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(values = focus_colors) +
    ggplot2::scale_x_continuous(
        limits = c(0.5, length(intersection_levels) + 0.5),
        breaks = seq_along(intersection_levels),
        labels = intersection_levels,
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
        limits = c(0, max(displayed_intersections$n) * 1.18),
        expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(x = NULL, y = "Intersection size") +
    base_theme +
    ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(
            color = "grey90", linewidth = 0.3
        ),
        plot.margin = ggplot2::margin(5, 5, 2, 3)
    )

intersection_matrix <- ggplot2::ggplot() +
    ggplot2::geom_rect(
        data = row_background,
        ggplot2::aes(
            xmin = -Inf, xmax = Inf,
            ymin = y - 0.5, ymax = y + 0.5, fill = fill
        ),
        color = NA
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_point(
        data = all_points,
        ggplot2::aes(x = intersection_index, y = y),
        color = "grey82", size = 2.05
    ) +
    ggplot2::geom_segment(
        data = segment_data,
        ggplot2::aes(
            x = intersection_index, xend = intersection_index,
            y = minimum, yend = maximum
        ),
        color = "grey15", linewidth = 0.4
    ) +
    ggplot2::geom_point(
        data = active_points,
        ggplot2::aes(x = intersection_index, y = y),
        color = "grey15", size = 2.45
    ) +
    ggplot2::geom_rect(
        data = focus_membership,
        ggplot2::aes(
            xmin = intersection_index - 0.42,
            xmax = intersection_index + 0.42,
            ymin = 0.55, ymax = length(display_method_order) + 0.45,
            color = taxon
        ),
        fill = NA, linewidth = 0.45, inherit.aes = FALSE
    ) +
    ggplot2::scale_color_manual(values = focus_colors) +
    ggplot2::scale_x_continuous(
        limits = c(0.5, length(intersection_levels) + 0.5),
        breaks = seq_along(intersection_levels),
        labels = intersection_levels,
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
        limits = c(0.5, length(display_method_order) + 0.5),
        breaks = positions$y,
        labels = positions$method,
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
        x = "Intersection pattern",
        y = NULL
    ) +
    base_theme +
    ggplot2::theme(
        legend.position = "none",
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(
            size = 7.8, hjust = 1, margin = ggplot2::margin(r = 4)
        ),
        axis.text.x = ggplot2::element_text(
            size = 6.6, angle = 90, hjust = 1, vjust = 0.5
        ),
        plot.margin = ggplot2::margin(4, 5, 8, 5)
    )

intersection_key <- ggplot2::ggplot() +
    ggplot2::geom_text(
        data = annotation_headers,
        ggplot2::aes(
            x = x, y = y, label = label, color = taxon
        ),
        hjust = 0, vjust = 0.5, size = 3.05, fontface = "bold"
    ) +
    ggplot2::geom_text(
        data = annotated_members[ordinary_italic_labels, ],
        ggplot2::aes(x = text_x, y = y, label = label),
        hjust = 0, size = 2.55, lineheight = 0.9,
        fontface = "italic", color = "grey25"
    ) +
    ggplot2::geom_text(
        data = annotated_members[ordinary_plain_labels, ],
        ggplot2::aes(x = text_x, y = y, label = label),
        hjust = 0, size = 2.55, lineheight = 0.9, color = "grey25"
    ) +
    ggplot2::geom_point(
        data = annotated_members[annotated_members$highlighted, ],
        ggplot2::aes(x = marker_x, y = y, color = taxon),
        shape = 1, size = 3.0, stroke = 0.85
    ) +
    ggplot2::geom_text(
        data = annotated_members[annotated_members$highlighted, ],
        ggplot2::aes(
            x = text_x, y = y, label = label, color = taxon
        ),
        hjust = 0, size = 2.7, fontface = "bold.italic"
    ) +
    ggplot2::scale_color_manual(values = focus_colors) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 7.6)) +
    ggplot2::theme_void() +
    ggplot2::theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(5, 5, 3, 9)
    )

intersection_figure <- patchwork::wrap_plots(
    intersection_key,
    intersection_bar,
    intersection_matrix,
    ncol = 1,
    heights = c(0.72, 0.65, 1.55)
)

grDevices::pdf(
    file.path(figure_directory, "schubert_cdi_method_intersections.pdf"),
    width = 225 / 25.4, height = 105 / 25.4, onefile = TRUE,
    family = "Helvetica", useDingbats = FALSE, version = "1.5"
)
print(intersection_figure)
grDevices::dev.off()

# DASRA component evidence ------------------------------------------------

component_rows <- results[results$method %in% c(
    "DASRA structural absence",
    "DASRA present-conditional abundance",
    "DASRA combined"
), c("taxon", "method", "available", "p_value", "q_value",
     "significant", "estimate", "statistic")]
component_wide <- tidyr::pivot_wider(
    component_rows,
    names_from = method,
    values_from = c(available, p_value, q_value, significant,
                    estimate, statistic),
    names_sep = "__"
)
names(component_wide) <- gsub(
    "DASRA structural absence", "structural", names(component_wide),
    fixed = TRUE
)
names(component_wide) <- gsub(
    "DASRA present-conditional abundance", "abundance",
    names(component_wide), fixed = TRUE
)
names(component_wide) <- gsub(
    "DASRA combined", "combined", names(component_wide), fixed = TRUE
)

structural_gate <- bh_z_gate(
    component_wide$p_value__structural,
    component_wide$q_value__structural
)
abundance_gate <- bh_z_gate(
    component_wide$p_value__abundance,
    component_wide$q_value__abundance
)
component_evidence <- rbind(
    data.frame(
        taxon = component_wide$taxon,
        component = "Structural absence",
        statistic = component_wide$statistic__structural,
        adjusted_p = component_wide$q_value__structural,
        gate = structural_gate,
        stringsAsFactors = FALSE
    ),
    data.frame(
        taxon = component_wide$taxon,
        component = "Present-conditional abundance",
        statistic = component_wide$statistic__abundance,
        adjusted_p = component_wide$q_value__abundance,
        gate = abundance_gate,
        stringsAsFactors = FALSE
    )
)
component_evidence <- component_evidence[
    component_evidence$taxon %in% component_taxa &
        is.finite(component_evidence$statistic),
]
component_evidence$taxon <- factor(
    component_evidence$taxon, levels = rev(component_taxa)
)
component_evidence$component <- factor(
    component_evidence$component,
    levels = c("Structural absence", "Present-conditional abundance")
)
component_evidence$selected <- component_evidence$adjusted_p <= 0.05
component_gate_data <- unique(
    component_evidence[, c("component", "gate")]
)
component_gate_data <- rbind(
    transform(component_gate_data, boundary = -gate),
    transform(component_gate_data, boundary = gate)
)
component_facet_labels <- c(
    "Structural absence" = "Structural\nabsence",
    "Present-conditional abundance" = "Present-conditional\nabundance"
)

# Taxon-level observed summaries -----------------------------------------

group_label <- c("H" = "Healthy", "CDI" = "CDI")
taxon_observations <- do.call(rbind, lapply(component_taxa, function(taxon) {
    data.frame(
        taxon = taxon,
        group = unname(group_label[as.character(input$metadata$group)]),
        count = as.numeric(input$counts[taxon, ]),
        relative_abundance = 100 * as.numeric(input$counts[taxon, ]) /
            input$metadata$library_size,
        stringsAsFactors = FALSE
    )
}))
taxon_observations$group <- factor(
    taxon_observations$group, levels = c("Healthy", "CDI")
)

detection_summary <- do.call(rbind, lapply(
    split(taxon_observations, list(
        taxon_observations$taxon, taxon_observations$group
    ), drop = TRUE),
    function(data) {
        successes <- sum(data$count > 0)
        total <- nrow(data)
        interval <- wilson_interval(successes, total)
        data.frame(
            taxon = data$taxon[[1L]], group = data$group[[1L]],
            positive = successes, total = total,
            proportion = successes / total,
            lower = interval[["lower"]], upper = interval[["upper"]],
            stringsAsFactors = FALSE
        )
    }
))
detection_summary$group <- factor(
    detection_summary$group, levels = c("Healthy", "CDI")
)

focus_results <- component_wide[
    match(component_taxa, component_wide$taxon),
]
names(focus_results)[names(focus_results) == "taxon"] <- "feature"
focus_intersections <- display_membership[
    match(component_taxa, display_membership$taxon),
    c("taxon", "intersection", "methods")
]
focus_summary <- merge(
    detection_summary,
    focus_results[, c(
        "feature", "q_value__structural", "q_value__abundance",
        "q_value__combined", "statistic__structural",
        "statistic__abundance", "estimate__abundance"
    )],
    by.x = "taxon", by.y = "feature", all.x = TRUE, sort = FALSE
)
focus_summary <- merge(
    focus_summary, focus_intersections,
    by = "taxon", all.x = TRUE, sort = FALSE
)
focus_summary <- focus_summary[
    order(match(focus_summary$taxon, component_taxa), focus_summary$group),
]
focus_summary <- focus_summary[, c(
    "taxon", "group", "positive", "total", "proportion", "lower", "upper",
    "q_value__structural", "q_value__abundance", "q_value__combined",
    "statistic__structural", "statistic__abundance", "estimate__abundance",
    "intersection", "methods"
)]
names(focus_summary) <- c(
    "taxon", "group", "observed_positive_samples", "total_samples",
    "observed_detection_proportion", "detection_ci_lower",
    "detection_ci_upper", "structural_bh_adjusted_p",
    "abundance_bh_adjusted_p", "combined_bh_adjusted_p",
    "structural_z_statistic", "abundance_z_statistic",
    "present_conditional_abundance_estimate", "primary_result_intersection",
    "primary_results"
)
utils::write.csv(
    focus_summary,
    file.path(table_directory, "schubert_cdi_component_examples.csv"),
    row.names = FALSE, na = ""
)
utils::write.csv(
    method_intersections_export,
    file.path(table_directory, "schubert_cdi_method_intersections.csv"),
    row.names = FALSE, na = ""
)
utils::write.csv(
    method_membership_export,
    file.path(
        table_directory,
        "schubert_cdi_method_intersection_membership.csv"
    ),
    row.names = FALSE, na = ""
)

detection_summary$taxon <- factor(
    detection_summary$taxon, levels = rev(component_taxa)
)
taxon_observations$taxon <- factor(
    taxon_observations$taxon, levels = rev(component_taxa)
)
component_plot <- ggplot2::ggplot(
    component_evidence,
    ggplot2::aes(x = statistic, y = taxon)
) +
    ggplot2::geom_vline(
        xintercept = 0, color = "grey65", linewidth = 0.35
    ) +
    ggplot2::geom_vline(
        data = component_gate_data,
        ggplot2::aes(xintercept = boundary),
        color = "grey55", linetype = "dashed", linewidth = 0.45
    ) +
    ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = statistic, yend = taxon),
        color = "grey70", linewidth = 0.45
    ) +
    ggplot2::geom_point(
        ggplot2::aes(fill = selected),
        shape = 21, color = "grey15", size = 3.0, stroke = 0.7
    ) +
    ggplot2::facet_wrap(
        ggplot2::vars(component), nrow = 1,
        labeller = ggplot2::as_labeller(component_facet_labels)
    ) +
    ggplot2::scale_fill_manual(
        values = c("FALSE" = "white", "TRUE" = "grey15"),
        labels = c(
            "FALSE" = "BH-adjusted p > 0.05",
            "TRUE" = "BH-adjusted p <= 0.05"
        )
    ) +
    ggplot2::scale_y_discrete(
        labels = function(values) {
            parse(text = paste0("italic(", values, ")"))
        }
    ) +
    ggplot2::scale_x_continuous(
        limits = c(-6.4, 6.4), breaks = c(-6, -3, 0, 3, 6),
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
        x = "Signed component\nZ-statistic", y = NULL,
        fill = NULL, tag = "A"
    ) +
    base_theme +
    ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
            size = 7.6, face = "bold", lineheight = 0.92
        ),
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 6.8),
        legend.key.width = grid::unit(8, "pt"),
        legend.spacing.x = grid::unit(2, "pt"),
        panel.spacing.x = grid::unit(7, "pt"),
        panel.grid.major.y = ggplot2::element_line(
            color = "grey92", linewidth = 0.3
        ),
        plot.margin = ggplot2::margin(7, 4, 4, 7)
    )

detection_summary$taxon_index <- as.numeric(detection_summary$taxon)
detection_summary$y <- detection_summary$taxon_index + ifelse(
    detection_summary$group == "Healthy", 0.13, -0.13
)

detection_plot <- ggplot2::ggplot(detection_summary) +
    ggplot2::geom_segment(
        ggplot2::aes(
            x = lower, xend = upper, y = y, yend = y, color = group
        ),
        linewidth = 0.7
    ) +
    ggplot2::geom_segment(
        ggplot2::aes(
            x = lower, xend = lower, y = y - 0.045, yend = y + 0.045,
            color = group
        ),
        linewidth = 0.55
    ) +
    ggplot2::geom_segment(
        ggplot2::aes(
            x = upper, xend = upper, y = y - 0.045, yend = y + 0.045,
            color = group
        ),
        linewidth = 0.55
    ) +
    ggplot2::geom_point(
        ggplot2::aes(x = proportion, y = y, color = group), size = 2.7
    ) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_y_continuous(
        breaks = seq_along(component_taxa),
        labels = parse(text = paste0("italic(", rev(component_taxa), ")")),
        limits = c(0.55, length(component_taxa) + 0.45),
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_x_continuous(
        limits = c(0, 1.04), breaks = c(0, 0.5, 1),
        labels = scales::label_percent(accuracy = 1),
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
        x = "Observed detection", y = NULL,
        color = NULL, tag = "B"
    ) +
    base_theme +
    ggplot2::theme(
        legend.position = "bottom",
        legend.text = ggplot2::element_text(size = 7.0),
        legend.key.width = grid::unit(9, "pt"),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(
            color = "grey92", linewidth = 0.3
        ),
        plot.margin = ggplot2::margin(7, 4, 4, 4)
    )

positive_observations <- taxon_observations[
    taxon_observations$relative_abundance > 0,
]
positive_observations$taxon_index <- as.numeric(positive_observations$taxon)
positive_observations$y <- positive_observations$taxon_index + ifelse(
    positive_observations$group == "Healthy", 0.13, -0.13
)

abundance_plot <- ggplot2::ggplot(
    positive_observations,
    ggplot2::aes(
        x = relative_abundance, y = y,
        color = group, fill = group,
        group = interaction(taxon, group)
    )
) +
    ggplot2::geom_boxplot(
        orientation = "y", width = 0.20, outlier.shape = NA,
        alpha = 0.17, linewidth = 0.55
    ) +
    ggplot2::geom_point(
        position = ggplot2::position_jitter(
            width = 0, height = 0.045, seed = 1L
        ),
        size = 0.65, alpha = 0.30, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::guides(color = "none", fill = "none") +
    ggplot2::scale_y_continuous(
        breaks = seq_along(component_taxa),
        labels = parse(text = paste0("italic(", rev(component_taxa), ")")),
        limits = c(0.55, length(component_taxa) + 0.45),
        expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_x_log10(
        limits = c(0.01, 100),
        breaks = c(0.01, 0.1, 1, 10, 100),
        labels = c("0.01", "0.1", "1", "10", "100")
    ) +
    ggplot2::labs(
        x = "Relative abundance among\nobserved-positive samples (%)",
        y = NULL, tag = "C"
    ) +
    base_theme +
    ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(
            color = "grey92", linewidth = 0.3
        ),
        plot.margin = ggplot2::margin(7, 7, 4, 4)
    )

component_figure <- patchwork::wrap_plots(
    component_plot, detection_plot, abundance_plot,
    nrow = 1, widths = c(1.55, 0.82, 1.18)
)

grDevices::pdf(
    file.path(figure_directory, "schubert_cdi_component_examples.pdf"),
    width = 180 / 25.4, height = 90 / 25.4, onefile = TRUE,
    family = "Helvetica", useDingbats = FALSE, version = "1.5"
)
print(component_figure)
grDevices::dev.off()
