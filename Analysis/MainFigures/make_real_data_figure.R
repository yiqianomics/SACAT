#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(grid)
})

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
directory <- dirname(normalizePath(script))
data_directory <- file.path(directory, "data", "real_data")
read_table <- function(name) {
    read.csv(file.path(data_directory, name), stringsAsFactors = FALSE)
}
recovery <- read_table("rauer_mock_recovery.csv")
results <- read_table("schubert_cdi_figure_method_results.csv")
samples <- read_table("schubert_cdi_figure_sample_values.csv")
fitted <- read_table("schubert_cdi_figure_fitted_abundance.csv")
taxa <- c("Coprococcus", "Enterococcus", "Faecalibacterium")
scale <- 183 / 210
font_families <- if (requireNamespace("systemfonts", quietly = TRUE)) {
    systemfonts::system_fonts()$family
} else character()
figure_font <- if ("Arial" %in% font_families) "Arial" else if (
    "Liberation Sans" %in% font_families
) "Liberation Sans" else "sans"

# Discoveries retain the original 88-taxon BH families.
for (method in unique(results$method)) {
    family <- results[results$method == method, ]
    stopifnot(nrow(family) == 88L, !anyDuplicated(family$taxon),
              max(abs(p.adjust(family$p_for_bh, "BH") - family$q_value)) < 1e-12,
              all(family$significant == (family$available & family$q_value <= .05)))
}
stopifnot(nrow(recovery) == 42L,
          max(abs(recovery$recovery_percent - 100 * recovery$n_discoveries /
                  (recovery$n_allocations * recovery$n_targets))) < 1e-10,
          all(samples$gamma >= 0 & samples$gamma <= 1),
          all(samples$gamma[samples$count > 0] == 0))
recovery$method[recovery$method == "MaAsLin3"] <- "MaAsLin 3"

components <- do.call(rbind, lapply(c(
    SA = "SACAT structural absence",
    RA = "SACAT present-conditional abundance"
), function(method) {
    family <- results[results$method == method, ]
    rows <- family[family$taxon %in% taxa, ]
    data.frame(taxon = rows$taxon,
               component = if (method == "SACAT structural absence") "SA" else "RA",
               statistic = rows$statistic, adjusted_p = rows$q_value,
               gate = qnorm(.05 * sum(family$significant) / nrow(family) / 2,
                            lower.tail = FALSE))
}))

# Zero-source shares average posterior weights over all samples in each group.
allocations <- do.call(rbind, lapply(taxa, function(taxon) {
    do.call(rbind, lapply(c("Healthy", "CDI"), function(group) {
        rows <- samples[samples$taxon == taxon & samples$group == group, ]
        data.frame(taxon = taxon, group = group,
                   observed_detection = mean(rows$count > 0),
                   estimated_sa = mean(rows$gamma),
                   estimated_nd = mean((rows$count == 0) - rows$gamma))
    }))
}))
stopifnot(max(abs(rowSums(allocations[, 3:5]) - 1)) < 1e-12)
observations <- samples[samples$count > 0, c("taxon", "group", "count", "library_size")]
observations$relative_abundance <- 100 * observations$count / observations$library_size
names(fitted)[names(fitted) == "abundance_percent"] <- "abundance"

# Shared styling

combined_text_theme <- function(theme, scale = 1) {
    for (name in names(theme)) {
        element <- theme[[name]]
        if (inherits(element, "element_text")) {
            element$family <- figure_font
            element$size <- if (name == "plot.tag") {
                8 / (183 / 210) * scale
            } else {
                7 / (183 / 210) * scale
            }
            theme[[name]] <- element
        }
    }
    theme
}

combined_tag_theme <- function(scale = 1) {
    ggplot2::theme(
        plot.tag = ggplot2::element_text(
            family = figure_font, size = 8 / (183 / 210) * scale,
            face = "bold", hjust = 0, vjust = 1,
            margin = ggplot2::margin(0, 0, 0, 0)
        ),
        plot.tag.position = c(0.008, 0.995),
        plot.tag.location = "plot"
    )
}

# Spike-in recovery

make_mock_panels <- function(mock, scale = 1) {
    depth_levels <- c("native", "1000", "750", "500", "250", "100", "50")
    method_colors <- c("SACAT" = "#0072B2", "MaAsLin 3" = "#D55E00",
                       "ZINQ" = "#009E73")
    method_shapes <- c("SACAT" = 16, "MaAsLin 3" = 17, "ZINQ" = 15)
    mock$depth_position <- match(as.character(mock$depth), depth_levels) - 1L
    mock$method <- factor(mock$method, levels = names(method_colors))
    legend_theme <- combined_text_theme(ggplot2::theme(
        text = ggplot2::element_text(),
        legend.text = ggplot2::element_text(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.width = grid::unit(8 * scale, "pt"),
        legend.key.height = grid::unit(8 * scale, "pt"),
        legend.spacing.x = grid::unit(6 * scale, "pt"),
        legend.spacing.y = grid::unit(0, "pt"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(11 * scale, "pt"),
        legend.box.margin = ggplot2::margin(-2 * scale, 0, 5 * scale, 0)
    ), scale)
    panel_theme <- ggplot2::theme_classic(
        base_family = figure_font, base_size = 9 * scale
    ) + ggplot2::theme(
        plot.title.position = "panel",
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        strip.text = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(),
        axis.text = ggplot2::element_text(color = "#222222"),
        axis.title.x = ggplot2::element_text(
            margin = ggplot2::margin(t = 4 * scale)
        ),
        axis.title.y = ggplot2::element_text(
            margin = ggplot2::margin(r = 3 * scale)
        ),
        axis.line = ggplot2::element_line(
            linewidth = 0.45 * scale, color = "#333333"
        ),
        axis.ticks = ggplot2::element_line(
            linewidth = 0.4 * scale, color = "#333333"
        ),
        axis.ticks.length = grid::unit(2.5 * scale, "pt")
    ) + legend_theme + combined_tag_theme(scale)
    panel_theme <- combined_text_theme(panel_theme, scale)

    recovery_panel <- function(target, tag, show_y_title) {
        data <- mock[mock$target == target, , drop = FALSE]
        fixed_depth <- data[data$depth != "native", , drop = FALSE]
        native_transition <- data[
            data$depth %in% c("native", "1000"), , drop = FALSE
        ]
        ggplot2::ggplot(
            data, ggplot2::aes(
                x = depth_position, y = recovery_percent,
                color = method, shape = method, group = method
            )
        ) +
            ggplot2::annotate(
                "rect", xmin = -0.5, xmax = 0.5, ymin = -Inf,
                ymax = Inf, fill = "#F5F5F5", color = NA
            ) +
            ggplot2::geom_hline(
                yintercept = c(0, 25, 50, 75, 100),
                color = "#E5E5E5", linewidth = 0.35 * scale
            ) +
            ggplot2::geom_vline(
                xintercept = 0.5, color = "#C8C8C8", linewidth = 0.35 * scale
            ) +
            ggplot2::geom_line(
                data = native_transition, linewidth = 0.65 * scale,
                linetype = "22", show.legend = FALSE
            ) +
            ggplot2::geom_line(data = fixed_depth, linewidth = 0.75 * scale) +
            ggplot2::geom_point(size = 2.15 * scale, stroke = 0.8 * scale) +
            ggplot2::scale_x_continuous(
                breaks = 0:6,
                labels = c("Native", "1k", "750", "500", "250", "100", "50"),
                limits = c(-0.5, 6.5), expand = c(0, 0)
            ) +
            ggplot2::scale_y_continuous(
                limits = c(0, 103), breaks = c(0, 25, 50, 75, 100),
                labels = c(0, 25, 50, 75, 100),
                expand = ggplot2::expansion(mult = c(0, 0.01))
            ) +
            ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
            ggplot2::scale_shape_manual(values = method_shapes, drop = FALSE) +
            ggplot2::guides(
                color = ggplot2::guide_legend(
                    order = 1, override.aes = list(
                        shape = unname(method_shapes), linewidth = 0.75 * scale
                    )
                ),
                shape = "none"
            ) +
            ggplot2::labs(
                title = paste(target, "targets"), tag = tag,
                x = "Reads per sample",
                y = if (show_y_title) "Target recovery (%)" else NULL,
                color = NULL, shape = NULL
            ) + panel_theme + ggplot2::theme(
                plot.margin = if (show_y_title) {
                    ggplot2::margin(4, 7, 3, 4) * scale
                } else {
                    ggplot2::margin(4, 4, 3, 7) * scale
                }
            )
    }

    panels <- patchwork::wrap_plots(
        recovery_panel("Structural", "A", TRUE),
        recovery_panel("Abundance", "B", FALSE),
        nrow = 1, widths = c(1, 1),
        guides = "collect", axis_titles = "collect_x"
    )
    panels <- panels & legend_theme
    patchwork::free(panels, side = "l")
}

# CDI component statistics

make_component_panel <- function(components, scale = 1) {
    taxa <- c("Coprococcus", "Enterococcus", "Faecalibacterium")
    components$taxon <- factor(components$taxon, levels = rev(taxa))
    components$component <- factor(components$component, levels = c("SA", "RA"))
    components$selected <- components$adjusted_p <= 0.05
    gates <- unique(components[, c("component", "gate")])
    gates <- gates[is.finite(gates$gate), , drop = FALSE]
    boundaries <- rbind(
        transform(gates, boundary = -gate), transform(gates, boundary = gate)
    )
    panel_theme <- ggplot2::theme_classic(base_size = 9.5 * scale) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(),
            axis.text = ggplot2::element_text(color = "grey15"),
            axis.text.y = ggplot2::element_text(
                face = "italic", angle = 32, hjust = 1, vjust = 0.5,
                margin = ggplot2::margin(r = 2 * scale)
            ),
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(color = "grey25"),
            strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(face = "bold", lineheight = 0.92),
            legend.position = "bottom",
            legend.key.width = grid::unit(8 * scale, "pt"),
            legend.key.size = grid::unit(1.2 * scale, "lines"),
            legend.spacing.x = grid::unit(2 * scale, "pt"),
            legend.margin = ggplot2::margin(0, 0, 0, 0),
            legend.box.margin = ggplot2::margin(-2 * scale, 0, 0, 0),
            panel.spacing.x = grid::unit(7 * scale, "pt"),
            panel.grid.major.y = ggplot2::element_line(
                color = "grey92", linewidth = 0.3 * scale
            ),
            plot.margin = ggplot2::margin(7, 4, 4, 4) * scale
        ) + combined_tag_theme(scale)

    ggplot2::ggplot(components, ggplot2::aes(x = statistic, y = taxon)) +
        ggplot2::geom_vline(
            xintercept = 0, color = "grey65", linewidth = 0.35 * scale
        ) +
        ggplot2::geom_vline(
            data = boundaries, ggplot2::aes(xintercept = boundary),
            color = "grey55", linetype = "dashed", linewidth = 0.45 * scale
        ) +
        ggplot2::geom_segment(
            ggplot2::aes(x = 0, xend = statistic, yend = taxon),
            color = "grey70", linewidth = 0.45 * scale
        ) +
        ggplot2::geom_point(
            ggplot2::aes(fill = selected), shape = 21, color = "grey15",
            size = 3 * scale, stroke = 0.7 * scale
        ) +
        ggplot2::facet_wrap(ggplot2::vars(component), nrow = 1) +
        ggplot2::scale_fill_manual(
            values = c("FALSE" = "white", "TRUE" = "grey15"),
            labels = c("BH-adjusted p > 0.05", "BH-adjusted p \u2264 0.05")
        ) +
        ggplot2::scale_y_discrete() +
        ggplot2::scale_x_continuous(
            limits = c(-7, 7), breaks = c(-6, -3, 0, 3, 6), expand = c(0, 0)
        ) +
        ggplot2::labs(
            x = "Signed component\nstatistic", y = NULL, fill = NULL, tag = "D"
        ) +
        combined_text_theme(panel_theme, scale) +
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
}

# CDI detection, zero sources, and abundance

make_cdi_panels <- function(components, allocations, observations, fitted, scale = 1) {
    taxa <- c("Coprococcus", "Enterococcus", "Faecalibacterium")
    groups <- c("Healthy", "CDI")
    group_colors <- c("Healthy" = "#4C78A8", "CDI" = "#E45756")
    position_rows <- function(data) {
        data$taxon <- factor(data$taxon, levels = rev(taxa))
        data$group <- factor(data$group, levels = groups)
        data$position <- as.numeric(data$taxon) +
            ifelse(data$group == "Healthy", 0.15, -0.15)
        data
    }
    allocations <- position_rows(allocations)
    positive <- observations[observations$relative_abundance > 0, , drop = FALSE]
    positive <- position_rows(positive)
    fitted <- position_rows(fitted)
    allocations$panel <- "Observed detection and\nestimated zero sources"
    positive$panel <- fitted$panel <- "Observed and fitted abundance"
    states <- c("Detected", "SA", "ND")
    segments <- do.call(rbind, lapply(seq_len(nrow(allocations)), function(i) {
        row <- allocations[i, ]
        fractions <- c(row$observed_detection, row$estimated_sa, row$estimated_nd)
        ends <- cumsum(100 * fractions)
        data.frame(
            taxon = row$taxon, group = row$group, position = row$position,
            state = factor(states, levels = states),
            start = c(0, head(ends, -1)), end = ends, panel = row$panel
        )
    }))
    square_key <- function(data, params, size) {
        grid::rectGrob(
            width = grid::unit(7 * scale, "pt"),
            height = grid::unit(7 * scale, "pt"),
            gp = grid::gpar(fill = data$fill, col = NA)
        )
    }
    panel_theme <- ggplot2::theme_classic(
        base_family = figure_font, base_size = 8 * scale
    ) + ggplot2::theme(
        axis.title = ggplot2::element_text(),
        axis.text = ggplot2::element_text(color = "#262626"),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(linewidth = 0.3 * scale),
        axis.ticks = ggplot2::element_line(linewidth = 0.3 * scale),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold"),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "bottom", legend.title = ggplot2::element_blank(),
        legend.key.width = grid::unit(8 * scale, "pt"),
        legend.key.height = grid::unit(8 * scale, "pt"),
        legend.spacing.x = grid::unit(3 * scale, "pt"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(-2 * scale, 0, 0, 0),
        plot.margin = ggplot2::margin(7, 4, 4, 4) * scale
    ) + combined_tag_theme(scale)
    panel_theme <- combined_text_theme(panel_theme, scale)

    allocation_panel <- ggplot2::ggplot(segments) +
        ggplot2::geom_hline(
            yintercept = c(1.5, 2.5), color = "#EEEEEE", linewidth = 0.3 * scale
        ) +
        ggplot2::geom_rect(
            ggplot2::aes(
                xmin = start, xmax = end, ymin = position - 0.105,
                ymax = position + 0.105, fill = state
            ),
            color = "white", linewidth = 0.2 * scale, key_glyph = square_key
        ) +
        ggplot2::geom_rect(
            data = allocations, ggplot2::aes(
                xmin = 0, xmax = 100, ymin = position - 0.105,
                ymax = position + 0.105
            ),
            fill = NA, color = "#B2B2B2", linewidth = 0.25 * scale
        ) +
        ggplot2::facet_wrap(ggplot2::vars(panel), nrow = 1) +
        ggplot2::scale_fill_manual(
            values = c(Detected = "#6A9DC0", SA = "#E5A36C", ND = "#E9E9E9"),
            labels = c(Detected = "Detected", SA = "Estimated\nSA", ND = "Estimated\nND"),
            drop = FALSE
        ) +
        ggplot2::scale_x_continuous(
            limits = c(-1, 101), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0.4, 3.6), breaks = allocations$position,
            labels = as.character(allocations$group), expand = c(0, 0)
        ) +
        ggplot2::labs(x = "Share of\nsamples (%)", y = NULL, tag = "E") +
        panel_theme + combined_text_theme(ggplot2::theme(
            axis.text.y = ggplot2::element_text(
                margin = ggplot2::margin(r = 3 * scale)
            ),
            legend.text = ggplot2::element_text(hjust = 0.5, vjust = 0.5)
        ), scale) +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

    fitted_label <- "Fitted abundance\nConditional on presence"
    abundance_panel <- ggplot2::ggplot(
        positive, ggplot2::aes(
            x = relative_abundance, y = position,
            group = interaction(taxon, group), color = group, fill = group
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 1:3, color = "#EEEEEE", linewidth = 0.3 * scale
        ) +
        ggplot2::geom_boxplot(
            orientation = "y", width = 0.19, outlier.shape = NA,
            alpha = 0.18, linewidth = 0.5 * scale, staplewidth = 0.45
        ) +
        ggplot2::geom_point(
            ggplot2::aes(shape = group),
            position = ggplot2::position_jitter(width = 0, height = 0.033, seed = 173),
            size = 0.6 * scale, alpha = 0.27, show.legend = FALSE
        ) +
        ggplot2::geom_point(
            data = fitted, ggplot2::aes(x = abundance, y = position, shape = fitted_label),
            inherit.aes = FALSE, fill = "white", color = "#252525",
            size = 2.15 * scale, stroke = 0.6 * scale
        ) +
        ggplot2::facet_wrap(ggplot2::vars(panel), nrow = 1) +
        ggplot2::scale_color_manual(values = group_colors) +
        ggplot2::scale_fill_manual(values = group_colors) +
        ggplot2::scale_shape_manual(
            values = stats::setNames(c(16, 17, 23), c(groups, fitted_label)),
            breaks = fitted_label
        ) +
        ggplot2::scale_x_log10(
            limits = c(0.005, 120), breaks = c(0.01, 0.1, 1, 10, 100),
            labels = c("0.01", "0.1", "1", "10", "100"), expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(0.4, 3.6), breaks = 1:3, expand = c(0, 0)
        ) +
        ggplot2::labs(x = "Relative abundance (%)\n(log scale)", y = NULL, tag = "F") +
        panel_theme + ggplot2::theme(
            axis.text.y = ggplot2::element_blank(), legend.box = "vertical"
        ) +
        ggplot2::guides(
            color = "none",
            fill = ggplot2::guide_legend(
                nrow = 1, order = 1,
                override.aes = list(alpha = 0.3, color = unname(group_colors))
            ),
            shape = ggplot2::guide_legend(
                nrow = 1, order = 2,
                theme = ggplot2::theme(
                    legend.text = ggplot2::element_text(hjust = 0.5, vjust = 0.5)
                ),
                override.aes = list(
                    size = 2.15 * scale, fill = "white", color = "#252525"
                )
            )
        )

    patchwork::wrap_plots(
        make_component_panel(components, scale), allocation_panel, abundance_panel,
        nrow = 1, widths = c(1.15, 1.25, 1.70)
    )
}

# Method-by-taxon discovery matrix.
columns <- data.frame(
    method = c("SACAT structural absence", "SACAT present-conditional abundance",
               "SACAT combined", "ZINQ prevalence", "ZINQ abundance", "ZINQ combined",
               "MaAsLin3 prevalence", "MaAsLin3 abundance", "MaAsLin3 combined",
               "ANCOM-BC2", "LinDA", "corncob", "edgeR", "DESeq2", "metagenomeSeq"),
    group = c(rep("SACAT", 3), rep("ZINQ", 3), rep("MaAsLin 3", 3),
              "ANCOM-BC2", "LinDA", "corncob", "edgeR", "DESeq2", "metagenomeSeq"),
    test = c("SA", "RA", "Omni.", "Prev.", "Abund.", "Joint",
             "Prev.", "Abund.", "Joint", rep("", 6)),
    width = c(rep(8.5, 9), 18.5, 11, 14, 11, 13, 24)
)
columns$right <- 35 + cumsum(columns$width)
columns$left <- columns$right - columns$width
columns$center <- (columns$left + columns$right) / 2
selection <- expand.grid(taxon = taxa, method = columns$method, stringsAsFactors = FALSE)
index <- match(paste(selection$taxon, selection$method), paste(results$taxon, results$method))
stopifnot(!anyNA(index))
selected <- results[index, ]
selected$x <- columns$center[match(selected$method, columns$method)]
selected$y <- c(30, 23, 16)[match(selected$taxon, taxa)]

matrix_grob <- function() {
    children <- list()
    add <- function(grob) children[[length(children) + 1L]] <<- grob
    mm <- function(value) unit(value * scale, "mm")
    text <- function(label, x, y, size = 7, face = "plain", just = "center") {
        textGrob(label, x = mm(x), y = mm(y), just = just,
                 gp = gpar(fontfamily = figure_font, fontsize = size,
                           fontface = face, col = "#262626"))
    }
    line <- function(x0, x1, y0, y1 = y0, color = "#D4D8DB", width = .55) {
        segmentsGrob(x0 = mm(x0), x1 = mm(x1), y0 = mm(y0), y1 = mm(y1),
                     gp = gpar(col = color, lwd = width * scale))
    }
    add(text("C", 3, 43, size = 8, face = "bold", just = "left"))
    add(rectGrob(x = mm(119), y = mm(23), width = mm(168), height = mm(7),
                 gp = gpar(fill = "#F3F5F6", col = NA)))
    for (group in unique(columns$group)) {
        block <- columns[columns$group == group, ]
        multiple <- nrow(block) > 1L
        add(text(group, mean(range(c(block$left, block$right))),
                 if (multiple) 42 else 39.6, face = "bold"))
        if (multiple) {
            add(line(min(block$left) + .6, max(block$right) - .6, 39.4))
            for (i in seq_len(nrow(block))) {
                add(text(block$test[i], block$center[i], 36.7))
            }
        }
    }
    for (i in seq_along(taxa)) {
        add(text(taxa[i], 32, c(30, 23, 16)[i], face = "italic", just = "right"))
    }
    for (y in c(33.5, 12.5)) add(line(35, 203, y, color = "#A5ADB3", width = .65))
    group_ends <- columns$right[!duplicated(columns$group, fromLast = TRUE)]
    for (x in head(group_ends, -1L)) add(line(x, x, 12.5, 34))
    add(pointsGrob(x = mm(selected$x), y = mm(selected$y), pch = 21, size = mm(2.7),
                   gp = gpar(fill = ifelse(selected$significant, "#282828", "white"),
                             col = "#282828", lwd = .8 * scale)))
    add(pointsGrob(x = mm(c(69, 121)), y = mm(c(5.5, 5.5)), pch = 21, size = mm(2.7),
                   gp = gpar(fill = c("#282828", "white"), col = "#282828", lwd = .8 * scale)))
    add(text("BH-adjusted p \u2264 0.05", 72, 5.5, just = "left"))
    add(text("BH-adjusted p > 0.05", 124, 5.5, just = "left"))
    gTree(children = do.call(gList, children))
}

draw_block <- function(grob, bottom, height) {
    pushViewport(viewport(x = unit(0, "mm"), y = unit(bottom * scale, "mm"),
                          width = unit(210 * scale, "mm"), height = unit(height * scale, "mm"),
                          just = c("left", "bottom")))
    grid.draw(grob)
    popViewport()
}

output_directory <- file.path(directory, "figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output <- file.path(output_directory, "main_real_data_combined.pdf")
if (capabilities("aqua")) {
    grDevices::quartz(type = "pdf", file = output, width = 183 / 25.4,
                     height = 181 * scale / 25.4, family = figure_font, bg = "white",
                     title = "Signal recovery and CDI associations")
} else {
    grDevices::cairo_pdf(output, width = 183 / 25.4, height = 181 * scale / 25.4,
                        family = figure_font, bg = "white")
}
grid.newpage()
draw_block(patchworkGrob(make_mock_panels(recovery, scale)), 128, 53)
draw_block(matrix_grob(), 83, 45)
draw_block(patchworkGrob(make_cdi_panels(
    components, allocations, observations, fitted, scale
)), 0, 83)
invisible(dev.off())
message("Saved ", output)
