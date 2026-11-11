#!/usr/bin/env Rscript

# Render the simulation figure from the two simulation summaries.
# Run with Rscript; input and output paths are relative to this script.
# Requires R >= 4.1, ggplot2 >= 3.5.0, and systemfonts.

required_packages <- c("ggplot2", "systemfonts")
missing_packages <- required_packages[!vapply(required_packages,
    requireNamespace, quietly = TRUE, FUN.VALUE = logical(1L))]
if (length(missing_packages)) {
    stop("Install the required R packages: ", paste(missing_packages, collapse = ", "),
         call. = FALSE)
}
if (packageVersion("ggplot2") < "3.5.0") {
    stop("ggplot2 >= 3.5.0 is required.", call. = FALSE)
}
available_fonts <- systemfonts::system_fonts()$family
font_family <- if ("Arial" %in% available_fonts) "Arial" else if (
    "Liberation Sans" %in% available_fonts) "Liberation Sans" else "sans"
suppressPackageStartupMessages(library(ggplot2))
library(grid)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
    stop("Run this generator with Rscript path/to/make_simulation_figure.R.", call. = FALSE)
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_argument), mustWork = TRUE))
split_simulations <- "--split" %in% commandArgs(trailingOnly = TRUE)
input_files <- file.path(dirname(script_dir),
    c("SimulationInternalValidation", "SimulationMetaSPARSim"), "tables", "main_plot_data.csv")
output_dir <- file.path(script_dir, "figures")
output_file <- file.path(output_dir, "main_simulation_calibrated_bh_power_combined.pdf")
if (!all(file.exists(input_files))) {
    stop("Missing simulation summaries: ",
         paste(input_files[!file.exists(input_files)], collapse = ", "), call. = FALSE)
}
plot_data <- do.call(rbind, lapply(input_files, read.csv,
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA")))
required_columns <- c("panel", "simulation", "setting_id", "scenario", "n_per_group",
    "covariate", "effect_parameter", "effect_index", "structural_difference",
    "abundance_log_effect", "method", "component", "plotted_metric", "plotted_value",
    "original_bh_value", "include_covariate")
stopifnot(all(required_columns %in% names(plot_data)), nrow(plot_data) > 0L,
          setequal(plot_data$panel, c("A", "B")),
          all(is.finite(plot_data$plotted_value)),
          all(plot_data$plotted_value >= 0 & plot_data$plotted_value <= 1),
          all(grepl("calibrated p + BH", plot_data$plotted_metric, fixed = TRUE)),
          !anyDuplicated(plot_data[c("panel", "setting_id", "method", "component")]))

core_methods <- c("SACAT", "ZINQ", "MaAsLin 3")
methods <- c(core_methods, "ANCOM-BC2", "DESeq2", "edgeR", "LinDA", "corncob", "metagenomeSeq")
internal <- plot_data[plot_data$panel == "A", ]
external <- plot_data[plot_data$panel == "B", ]
stopifnot(all(internal$n_per_group == 80L), setequal(plot_data$method, methods),
          all(plot_data$method[!plot_data$scenario %in% c("joint", "joint_effect")] %in% core_methods))
internal$effect_position <- ifelse(internal$scenario == "joint_effect",
                                  internal$effect_index, internal$effect_parameter)
internal_joint <- internal[internal$scenario == "joint_effect", ]

colors <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
                     "#882255", "#B79F00", "#222222", "#777777"), methods)
shapes <- setNames(c(16, 17, 15, 18, 3, 7, 8, 1, 4), methods)
linetypes <- setNames(c(rep("solid", 3), rep("dashed", 6)), methods)
covariate_labels <- c("Without covariate", "With covariate")
internal_scenarios <- c("observed_prevalence_only", "structural_matched_prevalence",
                        "present_conditional_abundance", "joint_effect")
internal_titles <- c("Different prevalence,\nsame structural absence",
                     "Same prevalence,\ndifferent structural absence",
                     "Abundance effects", "Joint effects")
external_scenarios <- c("structural_only", "abundance_only", "joint")
external_titles <- c("Structural effects", "Abundance effects", "Joint effects")
internal$method <- factor(internal$method, levels = methods)
internal$scenario_label <- factor(internal$scenario, levels = internal_scenarios,
                                  labels = internal_titles)
internal$covariate_label <- factor(internal$covariate,
    levels = c("without_covariate", "with_covariate"), labels = covariate_labels)
external$method <- factor(external$method, levels = methods)
internal <- internal[order(match(internal$method,
    c(core_methods, "edgeR", "DESeq2", "ANCOM-BC2", "LinDA", "corncob", "metagenomeSeq"))), ]
external <- external[order(external$method), ]
external$scenario_label <- factor(external$scenario, levels = external_scenarios,
                                  labels = external_titles)
external$covariate_label <- factor(external$include_covariate,
    levels = c(FALSE, TRUE), labels = covariate_labels)
stopifnot(!anyNA(internal[c("method", "scenario_label", "covariate_label")]),
          !anyNA(external[c("method", "scenario_label", "covariate_label")]))

internal_joint_labels <- unique(internal_joint[
    c("effect_index", "structural_difference", "abundance_log_effect")])
internal_joint_labels <- internal_joint_labels[order(internal_joint_labels$effect_index), ]
internal_joint_labels$label <- sprintf("%.2f\n%.2f",
    internal_joint_labels$structural_difference, internal_joint_labels$abundance_log_effect)
internal_breaks <- function(limits) {
    if (max(limits) > 3) return(internal_joint_labels$effect_index)
    if (max(limits) > .6) return(c(.20, .30, .40, .50, .60, .75))
    if (max(limits) > .4) return(c(0, .10, .20, .35, .50))
    c(0, .08, .16, .24, .32)
}
internal_axis_labels <- function(breaks) {
    if (max(breaks, na.rm = TRUE) > 3) {
        return(internal_joint_labels$label[match(breaks, internal_joint_labels$effect_index)])
    }
    ifelse(breaks == 0, "0", sub("^0", "", sprintf("%.2f", breaks)))
}
external$effect_key <- factor(paste(external$scenario, external$effect_index, sep = "_"),
    levels = unlist(lapply(external_scenarios, function(s) paste(s, 1:6, sep = "_"))))
external_labels <- unique(external[
    c("effect_key", "scenario", "structural_difference", "abundance_log_effect")])
external_labels$label <- ifelse(external_labels$scenario == "structural_only",
    sprintf("%.2f", external_labels$structural_difference),
    ifelse(external_labels$scenario == "abundance_only",
        sprintf("%.2f", external_labels$abundance_log_effect),
        sprintf("%.2f\n%.2f", external_labels$structural_difference, external_labels$abundance_log_effect)))

figure_theme <- theme_bw(base_size = 9, base_family = font_family) + theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E1E5E8", linewidth = .23),
    panel.border = element_rect(color = "#92989D", linewidth = .28),
    strip.background = element_rect(fill = "#F1F3F4", color = "#B0B5B8", linewidth = .25),
    strip.text.x = element_text(size = 9, face = "bold", lineheight = 1.1,
                                margin = margin(4, 2, 4, 2)),
    strip.text.y = element_text(size = 8, angle = -90, margin = margin(2, 3, 2, 3)),
    axis.text = element_text(size = 8, color = "#222222"),
    axis.text.x = element_text(lineheight = 1.05),
    axis.title = element_text(size = 9),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 4)),
    axis.ticks = element_line(linewidth = .26, color = "#555555"),
    axis.ticks.length = unit(1, "mm"),
    panel.spacing.x = unit(3.5, "mm"),
    panel.spacing.y = unit(3.0, "mm"),
    plot.margin = margin(1, 1, 1, 1),
    legend.position = "none", legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.width = unit(6, "mm"), legend.key.height = unit(3.5, "mm"),
    legend.key.spacing.x = unit(2.1, "mm"), legend.key.spacing.y = unit(.6, "mm"),
    legend.margin = margin(0, 0, 0, 0)
)
method_scales <- list(
    scale_color_manual(values = colors, breaks = methods, drop = FALSE),
    scale_shape_manual(values = shapes, breaks = methods, drop = FALSE),
    scale_linetype_manual(values = linetypes, breaks = methods, drop = FALSE),
    scale_y_continuous(limits = c(-.012, 1.025), breaks = c(0, .5, 1),
                       labels = c("0", "0.5", "1"), expand = expansion(mult = 0))
)
internal_plot <- ggplot(internal,
    aes(effect_position, plotted_value, color = method, shape = method,
        linetype = method, group = method)) +
    geom_line(linewidth = .50) + geom_point(size = 1.1, stroke = .3) +
    facet_grid(covariate_label ~ scenario_label, scales = "free_x") + method_scales +
    scale_x_continuous(breaks = internal_breaks, labels = internal_axis_labels,
                       expand = expansion(mult = c(.065, .065))) +
    labs(x = "Effect setting (SA / RA for joint effects)", y = "Rejection probability") + figure_theme
external_plot <- ggplot(external,
    aes(effect_key, plotted_value, color = method, shape = method,
        linetype = method, group = method)) +
    geom_line(linewidth = .50) + geom_point(size = 1.1, stroke = .3) +
    facet_grid(covariate_label ~ scenario_label, scales = "free_x") + method_scales +
    scale_x_discrete(labels = setNames(external_labels$label, as.character(external_labels$effect_key)),
                     expand = expansion(add = c(.32, .32))) +
    labs(x = "Effect setting (SA / RA for joint effects)", y = "Power") + figure_theme

plots <- list(internal_plot, external_plot)
for (plot in plots) {
    built <- ggplot_build(plot)
    stopifnot(nrow(built$data[[1L]]) == nrow(plot$data),
              nrow(built$data[[2L]]) == nrow(plot$data),
              identical(sort(built$data[[2L]]$y), sort(plot$data$plotted_value)))
}
make_legend <- function(selected_methods) {
    legend_data <- data.frame(method = rep(selected_methods, each = 2L),
                              x = rep(1:2, length(selected_methods)), y = 1)
    legend_plot <- ggplot(legend_data,
        aes(x, y, color = method, shape = method, linetype = method, group = method)) +
        geom_line(linewidth = .50) + geom_point(size = 1.35, stroke = .3) +
        scale_color_manual(values = colors, breaks = selected_methods) +
        scale_shape_manual(values = shapes, breaks = selected_methods) +
        scale_linetype_manual(values = linetypes, breaks = selected_methods) +
        figure_theme + theme(legend.position = "bottom") +
        guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1),
               linetype = guide_legend(nrow = 1))
    legend_gtable <- ggplotGrob(legend_plot)
    legend_gtable$grobs[[which(legend_gtable$layout$name == "guide-box-bottom")]]
}
draw_grob <- function(grob, x, y, width, height) {
    pushViewport(viewport(x = unit(x, "mm"), y = unit(y, "mm"),
                          width = unit(width, "mm"), height = unit(height, "mm"),
                          just = c("left", "bottom")))
    grid.draw(grob)
    popViewport()
}
draw_text <- function(label, x, y, size = 8, bold = FALSE, just = "left") {
    grid.text(label, x = unit(x, "mm"), y = unit(y, "mm"), just = just,
              gp = gpar(fontfamily = font_family, fontsize = size,
                         fontface = if (bold) "bold" else "plain", col = "#222222"))
}
draw_figure <- function(panel = NULL) {
    grid.newpage()
    legends <- lapply(list(core_methods, setdiff(methods, core_methods)), make_legend)
    grobs <- lapply(plots, ggplotGrob)
    if (is.null(panel)) {
        draw_text("A", 5, 235, 12, TRUE)
        draw_grob(grobs[[1L]], 2, 124, 216, 106)
        draw_text("B", 5, 117, 12, TRUE)
        draw_grob(grobs[[2L]], 2, 16, 216, 97)
    } else {
        draw_grob(grobs[[panel]], 2, 16, 216, c(106, 97)[panel])
    }
    draw_grob(legends[[1L]], 2, 8, 216, 5)
    draw_grob(legends[[2L]], 2, 2.5, 216, 5)
}

# Use Quartz on macOS and Cairo on other platforms.
write_figure <- function(path, panel = NULL) {
    height <- if (is.null(panel)) 242 else c(124, 115)[panel]
    if (isTRUE(capabilities("aqua"))) {
        grDevices::quartz(type = "pdf", file = path, width = 220 / 25.4,
            height = height / 25.4, family = font_family, bg = "white",
            title = "Simulation component tests and power")
    } else if (isTRUE(capabilities("cairo"))) {
        grDevices::cairo_pdf(filename = path, width = 220 / 25.4,
            height = height / 25.4, family = font_family, bg = "white", onefile = TRUE)
    } else {
        stop("An R installation with Quartz or Cairo PDF support is required.", call. = FALSE)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_figure(panel)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (split_simulations) {
    output_files <- file.path(output_dir,
        c("main_simulation_internal.pdf", "main_simulation_metasparsim.pdf"))
    for (panel in seq_along(output_files)) {
        write_figure(output_files[panel], panel)
        message("Wrote ", output_files[panel])
    }
} else {
    write_figure(output_file)
    message("Wrote ", output_file)
}
