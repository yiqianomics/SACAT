#!/usr/bin/env Rscript

# Draw the manuscript comparisons from saved summaries without refitting tests.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]])
study_dir <- dirname(dirname(normalizePath(script)))
table_dir <- Sys.getenv("SACAT_SIMULATION_TABLES", file.path(study_dir, "tables"))
figure_dir <- Sys.getenv("SACAT_SIMULATION_FIGURES", file.path(study_dir, "figures"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
results <- fread(file.path(table_dir, "setting_summary.csv"))
methods <- c("SACAT", "ZINQ", "MaAsLin 3", "ANCOM-BC2", "DESeq2", "edgeR", "LinDA", "corncob", "metagenomeSeq")
colors <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#882255", "#B79F00", "#222222", "#777777"), methods)
shapes <- setNames(c(16, 17, 15, 18, 3, 7, 8, 1, 4), methods)
lines <- setNames(c(rep("solid", 3L), rep("dashed", 6L)), methods)
is_sa <- function(method, component) (method == "SACAT" & component == "structural_absence") |
    (method %in% c("ZINQ", "MaAsLin 3") & component == "observed_prevalence")
is_ra <- function(method, component) (method == "SACAT" & component == "present_conditional_abundance") |
    (method == "ZINQ" & component == "detected_quantile_abundance") |
    (method == "MaAsLin 3" & component == "detected_log_abundance") | component == "general_abundance"
decorate <- function(data) {
    data <- copy(data)
    data[, method := factor(method, levels = methods)]
    data[, covariate := factor(covariate, levels = c("without_covariate", "with_covariate"),
                               labels = c("Without covariate", "With covariate"))]
    data
}
page_theme <- theme_bw(base_size = 12, base_family = "sans") + theme(
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E1E5E8", linewidth = .3),
    strip.background = element_rect(fill = "#F1F3F4", color = "#B0B5B8"),
    strip.text = element_text(face = "bold", size = 11), plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11.5, margin = margin(b = 9)),
    legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 10.5), legend.key.width = grid::unit(1.1, "cm"),
    panel.spacing = grid::unit(.8, "lines"), plot.margin = margin(14, 17, 12, 12)
)
method_scales <- function(curves = TRUE, legend_rows = 2L) list(
    scale_color_manual(values = colors, drop = TRUE), scale_shape_manual(values = shapes, drop = TRUE),
    if (curves) scale_linetype_manual(values = lines, drop = TRUE),
    guides(color = guide_legend(nrow = legend_rows, byrow = TRUE), shape = guide_legend(nrow = legend_rows, byrow = TRUE),
           linetype = guide_legend(nrow = legend_rows, byrow = TRUE))
)
open_pdf <- function(filename, width = 13.5, height = 9) {
    grDevices::pdf(file.path(figure_dir, filename), width = width, height = height,
                   onefile = TRUE, useDingbats = FALSE, family = "Helvetica")
}
adjustments <- c(raw = "Raw p + BH (0.05)", calibrated = "Calibrated p + BH (0.05)")

mechanism <- decorate(results[run == "base" & n_per_group == 80L &
    scenario %in% c("observed_prevalence_only", "structural_matched_prevalence") & is_sa(method, component)])
mechanism[, story := factor(scenario,
    levels = c("observed_prevalence_only", "structural_matched_prevalence"),
    labels = c("Different prevalence, same structural absence", "Different structural absence, matched prevalence"))]
stopifnot(nrow(mechanism) == 60L)
open_pdf("mechanism_separation.pdf")
for (adjustment in names(adjustments)) {
    mechanism[, value := get(paste0(adjustment, "_bh_designated_rejection"))]
    plot <- ggplot(mechanism, aes(effect_parameter, value, color = method, shape = method,
                                  linetype = method, group = method)) +
        geom_line(linewidth = .75) + geom_point(size = 2.3) +
        facet_grid(covariate ~ story, scales = "free_x") + method_scales(legend_rows = 1L) + page_theme +
        scale_y_continuous(breaks = seq(0, 1, .2)) + coord_cartesian(ylim = c(0, 1)) +
        labs(title = "Mechanism separation", subtitle = paste("n = 80 per group |", adjustments[[adjustment]]),
             x = "Probability-difference effect parameter", y = "Designated-taxon rejection probability")
    print(plot)
}
invisible(dev.off())

sa <- results[run == "base" & n_per_group == 80L & scenario == "structural_only" & is_sa(method, component)]
ra <- results[run == "base" & n_per_group == 80L & scenario == "present_conditional_abundance" & is_ra(method, component)]
added <- results[run == "abundance"]
keys <- c("setting_id", "method", "component")
ra <- rbindlist(list(ra[!added, on = keys], added), use.names = TRUE)
joint <- results[run == "joint" & n_per_group == 80L & primary_comparison]
sa[, direction := "Structural-only"]
ra[, direction := "Abundance-only"]
joint[, direction := "Joint effects"]
comparison <- decorate(rbindlist(list(sa, ra, joint), use.names = TRUE))
comparison <- comparison[(direction == "Joint effects" & effect_index > 0L) |
    (direction != "Joint effects" & effect_parameter > 0)]
comparison[, direction := factor(direction, levels = c("Structural-only", "Abundance-only", "Joint effects"))]
comparison[, effect_label := ifelse(direction == "Joint effects",
    paste0(format(structural_difference, trim = TRUE), "/", format(abundance_log_effect, trim = TRUE)),
    format(effect_parameter, trim = TRUE))]
comparison[, effect_position := frank(if (direction[[1L]] == "Joint effects") effect_index else effect_parameter,
                                      ties.method = "dense"), by = direction]
comparison[, effect_key := paste(direction, effect_position, sep = "_")]
effect_order <- unique(comparison[order(direction, effect_position), .(effect_key, effect_label)])
comparison[, effect_key := factor(effect_key, levels = effect_order$effect_key)]
stopifnot(!anyDuplicated(comparison[, ..keys]), nrow(comparison) == 222L,
          all(is.finite(comparison$raw_bh_power)), all(is.finite(comparison$calibrated_bh_power)))
draw_comparison <- function(filename, metric, title, y_label) {
    open_pdf(filename)
    for (adjustment in names(adjustments)) {
        data <- copy(comparison)
        data[, value := get(paste0(adjustment, "_bh_", metric))]
        stopifnot(all(is.finite(data$value)), all(data$value >= 0 & data$value <= 1))
        # The 0.05 FDR reference is not an axis limit.
        ymax <- if (metric == "fdr") min(1, max(.10, max(data$value) * 1.08)) else 1
        plot <- ggplot(data, aes(effect_key, value, color = method, shape = method,
                                 linetype = method, group = method)) +
            geom_line(linewidth = .65) + geom_point(size = 2.1) +
            facet_grid(covariate ~ direction, scales = "free_x") + method_scales() + page_theme +
            scale_x_discrete(labels = setNames(effect_order$effect_label, effect_order$effect_key)) +
            scale_y_continuous(breaks = scales::breaks_pretty(5)) + coord_cartesian(ylim = c(0, ymax)) +
            theme(axis.text.x = element_text(size = 9, angle = 35, hjust = 1)) +
            labs(title = title, subtitle = paste("n = 80 per group |", adjustments[[adjustment]]),
                 x = "SA: probability difference | RA: absolute log effect | Joint: SA / RA effect", y = y_label)
        if (metric == "fdr") plot <- plot + geom_hline(yintercept = .05, linetype = "dotted", color = "#555555")
        print(plot)
        message(filename, " | ", adjustment, ": ", nrow(data), " values; range ", paste(range(data$value), collapse = " to "))
    }
    invisible(dev.off())
}
draw_comparison("power.pdf", "power", "Power across structural, abundance and joint effects", "Power")
draw_comparison("fdr.pdf", "fdr", "False discovery rate across the three effect directions", "False discovery rate")

nulls <- fread(file.path(table_dir, "type1_error.csv"))
sa_null <- nulls[run == "base" & scenario == "structural_only" & is_sa(method, component)]
ra_null <- nulls[run == "base" & scenario == "present_conditional_abundance" & is_ra(method, component)]
updated_null <- nulls[run == "abundance" & method == "ZINQ"]
ra_null <- rbindlist(list(ra_null[!updated_null, on = keys], updated_null), use.names = TRUE)
sa_null[, test := ifelse(method == "SACAT", "SACAT SA", paste(method, "prevalence"))]
ra_null[, test := ifelse(method == "SACAT", "SACAT RA", ifelse(
    method %in% c("ZINQ", "MaAsLin 3"), paste(method, "abundance"), method))]
component_null <- decorate(rbindlist(list(sa_null, ra_null), use.names = TRUE))
joint_keys <- results[run == "joint" & primary_comparison, ..keys]
joint_null <- decorate(merge(nulls[run == "joint"], unique(joint_keys), by = keys))
joint_null[, test := ifelse(method == "SACAT", "SACAT omnibus", ifelse(
    method %in% c("ZINQ", "MaAsLin 3"), paste(method, "joint"), as.character(method)))]
stopifnot(nrow(component_null) == 72L, nrow(joint_null) == 54L)
open_pdf("type1_error.pdf", width = 14, height = 12)
for (page in 1:2) {
    data <- copy(if (page == 1L) component_null else joint_null)
    data[, `:=`(lower = pmax(0, type1_error - 1.96 * monte_carlo_se),
                 upper = pmin(1, type1_error + 1.96 * monte_carlo_se))]
    test_order <- if (page == 1L) c("SACAT SA", "SACAT RA", "ZINQ prevalence", "ZINQ abundance",
        "MaAsLin 3 prevalence", "MaAsLin 3 abundance", methods[4:9]) else
        c("SACAT omnibus", "ZINQ joint", "MaAsLin 3 joint", methods[4:9])
    stopifnot(setequal(data$test, test_order))
    data[, test := factor(test, levels = rev(test_order))]
    data[, sample_size := factor(n_per_group, labels = paste("n =", sort(unique(n_per_group)), "per group"))]
    xmax <- min(1, max(.10, max(data$upper) * 1.05))
    plot <- ggplot(data, aes(type1_error, test, color = method, shape = method)) +
        geom_vline(xintercept = .05, linetype = "dotted", color = "#555555") +
        geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .2, linewidth = .45) +
        geom_point(size = 2.2) + facet_grid(sample_size ~ covariate) + method_scales(FALSE) + page_theme +
        theme(legend.position = "none", axis.text.y = element_text(size = 9), panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(color = "#E1E5E8", linewidth = .3)) +
        scale_x_continuous(breaks = scales::breaks_pretty(5)) + coord_cartesian(xlim = c(0, xmax)) +
        labs(title = "Global-null Type I error | raw p < 0.05",
             subtitle = if (page == 1L) "Component tests | 95% Monte Carlo intervals" else
                 "Joint-effect comparisons | 95% Monte Carlo intervals",
             x = "Per-taxon rejection probability", y = NULL)
    print(plot)
    message("type1_error.pdf | page ", page, ": ", nrow(data), " values; range ", paste(range(data$type1_error), collapse = " to "))
}
invisible(dev.off())
message("Created four two-page PDFs in ", figure_dir)
