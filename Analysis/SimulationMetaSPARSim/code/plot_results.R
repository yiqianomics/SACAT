#!/usr/bin/env Rscript

# Plot power, FDR, and global-null raw Type I error from the study tables.
# Requires data.table, ggplot2, and systemfonts. Only PDF files are written.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
})
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this script with Rscript.", call. = FALSE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)))
output <- file.path(root, "figures")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
summary <- fread(file.path(root, "tables", "setting_summary.csv"))
null <- fread(file.path(root, "tables", "type1_error.csv"))
fonts <- systemfonts::system_fonts()$family
font <- if ("Arial" %in% fonts) "Arial" else if ("Liberation Sans" %in% fonts) "Liberation Sans" else "sans"
methods <- c("SACAT", "ZINQ", "MaAsLin 3", "ANCOM-BC2", "DESeq2", "edgeR", "LinDA", "corncob", "metagenomeSeq")
colors <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
                     "#882255", "#B79F00", "#222222", "#777777"), methods)
shapes <- setNames(c(16, 17, 15, 18, 3, 7, 8, 1, 4), methods)
linetypes <- setNames(c(rep("solid", 3L), rep("dashed", 6L)), methods)
scenarios <- c("structural_only", "abundance_only", "joint")
scenario_titles <- c("Structural effects", "Abundance effects", "Joint effects")
covariate_titles <- c("Without covariate", "With covariate")
curves <- summary[
    (scenario == "structural_only" & component %in% c("structural_absence", "observed_prevalence")) |
    (scenario == "abundance_only" & component %in% c("present_conditional_abundance",
        "detected_quantile_abundance", "detected_log_abundance", "general_abundance")) |
    (scenario == "joint" & component %in% c("primary_omnibus", "native_joint", "general_abundance"))]
stopifnot(nrow(curves) == 252L, !anyDuplicated(curves[, .(setting_id, method)]),
          all(is.finite(curves$original_bh_power)), all(is.finite(curves$calibrated_bh_power)),
          all(is.finite(curves$original_bh_fdr)), all(is.finite(curves$calibrated_bh_fdr)))
curves[, `:=`(method = factor(method, levels = methods),
    scenario_label = factor(scenario, levels = scenarios, labels = scenario_titles),
    covariate_label = factor(include_covariate, levels = c(FALSE, TRUE), labels = covariate_titles),
    effect_key = factor(paste(scenario, effect_index, sep = "_"),
        levels = unlist(lapply(scenarios, function(s) paste(s, 1:6, sep = "_")))))]
axis_labels <- unique(curves[, .(effect_key, scenario, structural_difference, abundance_log_effect)])
axis_labels[, label := fifelse(scenario == "structural_only", sprintf("%.2f", structural_difference),
    fifelse(scenario == "abundance_only", sprintf("%.2f", abundance_log_effect),
            sprintf("%.2f\n%.2f", structural_difference, abundance_log_effect)))]

figure_theme <- theme_bw(base_size = 10, base_family = font) + theme(
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E1E5E8", linewidth = .25),
    panel.border = element_rect(color = "#92989D", linewidth = .3),
    strip.background = element_rect(fill = "#F1F3F4", color = "#B0B5B8", linewidth = .25),
    strip.text.x = element_text(size = 11, face = "bold", margin = margin(6, 3, 6, 3)),
    strip.text.y = element_text(size = 10, angle = -90, margin = margin(3, 5, 3, 5)),
    axis.text = element_text(size = 9, color = "#222222"),
    axis.text.x = element_text(lineheight = 1.1),
    axis.title = element_text(size = 10.5),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.ticks = element_line(linewidth = .3, color = "#555555"),
    axis.ticks.length = unit(1.2, "mm"),
    panel.spacing = unit(4, "mm"),
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10.5, color = "#444444", margin = margin(b = 10)),
    plot.title.position = "plot",
    plot.margin = margin(12, 14, 10, 12),
    legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 9.5), legend.key.width = unit(8, "mm"),
    legend.key.height = unit(4.5, "mm"), legend.key.spacing.x = unit(3, "mm"),
    legend.box.spacing = unit(3, "mm"), legend.margin = margin(1, 0, 1, 0)
)
method_scales <- list(
    scale_color_manual(values = colors, breaks = methods, drop = FALSE),
    scale_shape_manual(values = shapes, breaks = methods, drop = FALSE),
    scale_linetype_manual(values = linetypes, breaks = methods, drop = FALSE)
)
curve_plot <- function(metric, calibrated) {
    field <- paste0(if (calibrated) "calibrated" else "original", "_bh_", metric)
    value <- curves[[field]]
    metric_title <- if (metric == "power") "Power" else "False discovery rate"
    procedure <- if (calibrated) "Null-calibrated p-values + BH" else "Raw p-values + BH"
    figure <- ggplot(curves, aes(effect_key, .data[[field]], color = method,
        shape = method, linetype = method, group = method)) +
        geom_line(linewidth = .62) + geom_point(size = 1.7, stroke = .4) +
        facet_grid(covariate_label ~ scenario_label, scales = "free_x") +
        method_scales +
        scale_x_discrete(labels = setNames(axis_labels$label, as.character(axis_labels$effect_key)),
                         expand = expansion(add = c(.35, .35))) +
        labs(title = paste("metaSPARSim:", metric_title),
             subtitle = paste0(procedure, " (0.05)"),
             x = "Effect setting (SA / RA for joint effects)", y = metric_title) +
        guides(color = guide_legend(nrow = 2, byrow = TRUE),
               shape = guide_legend(nrow = 2, byrow = TRUE),
               linetype = guide_legend(nrow = 2, byrow = TRUE)) + figure_theme
    if (metric == "power") {
        figure <- figure + scale_y_continuous(limits = c(-.012, 1.025),
            breaks = c(0, .25, .5, .75, 1), expand = expansion(mult = 0))
    } else {
        upper <- if (calibrated) .10 else .75
        stopifnot(max(value) < upper)
        figure <- figure +
            geom_hline(yintercept = .05, color = "#B24335", linetype = "dotted", linewidth = .65) +
            scale_y_continuous(limits = c(-upper * .015, upper * 1.025),
                breaks = if (calibrated) seq(0, .1, .025) else seq(0, .75, .15),
                expand = expansion(mult = 0))
    }
    built <- ggplot_build(figure)
    stopifnot(nrow(built$data[[1L]]) == nrow(curves),
              nrow(built$data[[2L]]) == nrow(curves),
              identical(sort(built$data[[2L]]$y), sort(value)))
    figure
}

test_label <- function(method, component) {
    suffix <- fifelse(component == "structural_absence", "SA",
        fifelse(component == "present_conditional_abundance", "RA",
        fifelse(component == "primary_omnibus", "omnibus",
        fifelse(component == "observed_prevalence", "prevalence",
        fifelse(component %in% c("detected_quantile_abundance", "detected_log_abundance"), "abundance",
        fifelse(component == "native_joint", "joint", ""))))))
    trimws(paste(method, suffix))
}
test_order <- c("SACAT SA", "SACAT RA", "SACAT omnibus", "ZINQ prevalence", "ZINQ abundance", "ZINQ joint",
    "MaAsLin 3 prevalence", "MaAsLin 3 abundance", "MaAsLin 3 joint", "ANCOM-BC2", "DESeq2",
    "edgeR", "LinDA", "corncob", "metagenomeSeq")
stopifnot(nrow(null) == 30L, all(null$scenario == "global_null"))
null[, `:=`(label = factor(test_label(method, component), levels = rev(test_order)),
    covariate_label = factor(include_covariate, levels = c(FALSE, TRUE), labels = covariate_titles),
    lower = pmax(0, raw_type1 - 1.96 * raw_type1_mcse),
    upper = pmin(1, raw_type1 + 1.96 * raw_type1_mcse))]
stopifnot(!anyNA(null$label), all(null$upper < .25))
type1_plot <- ggplot(null, aes(raw_type1, label, color = method)) +
    geom_vline(xintercept = .05, color = "#B24335", linetype = "dotted", linewidth = .65) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .2, linewidth = .6) +
    geom_point(size = 2.1) +
    geom_text(aes(x = .258, label = sprintf("%.4f", raw_type1)), hjust = 0,
              family = font, size = 3.2, color = "#333333") +
    facet_grid(. ~ covariate_label) +
    scale_color_manual(values = colors, guide = "none") +
    scale_x_continuous(limits = c(0, .287), breaks = seq(0, .25, .05),
                       labels = c("0", "0.05", "0.10", "0.15", "0.20", "0.25"),
                       expand = expansion(mult = 0)) +
    scale_y_discrete(expand = expansion(add = c(.6, .6))) +
    labs(title = "metaSPARSim: global-null Type I error",
         subtitle = "Raw p < 0.05; 95% Monte Carlo intervals",
         x = "Mean taxonwise rejection probability", y = NULL) +
    figure_theme + theme(axis.ticks.y = element_blank(), panel.grid.major.x = element_blank(),
                         panel.grid.major.y = element_line(color = "#ECEFF1", linewidth = .25),
                         axis.text.y = element_text(size = 10), panel.spacing.x = unit(7, "mm"))

write_pdf <- function(filename, plots, width, height) {
    if (isTRUE(capabilities("aqua"))) {
        grDevices::quartz(type = "pdf", file = file.path(output, filename),
                         width = width / 25.4, height = height / 25.4,
                         family = font, bg = "white")
    } else {
        grDevices::cairo_pdf(file.path(output, filename), width = width / 25.4,
                            height = height / 25.4, family = font, bg = "white", onefile = TRUE)
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    for (plot in plots) print(plot)
}
write_pdf("power.pdf", list(curve_plot("power", FALSE), curve_plot("power", TRUE)), 270, 210)
write_pdf("fdr.pdf", list(curve_plot("fdr", FALSE), curve_plot("fdr", TRUE)), 270, 210)
write_pdf("type1_error.pdf", list(type1_plot), 270, 175)
message("Wrote power.pdf (2 pages), fdr.pdf (2 pages), and type1_error.pdf (1 page).")
