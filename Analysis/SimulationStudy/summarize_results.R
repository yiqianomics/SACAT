#!/usr/bin/env Rscript

# Create simulation figures and tables from the aggregated results.

required_packages <- c("data.table", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
    stop(
        "Install the following R packages before running this script: ",
        paste(missing_packages, collapse = ", "),
        call. = FALSE
    )
}

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
    normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = TRUE)
} else {
    normalizePath("summarize_results.R", mustWork = TRUE)
}
study_dir <- dirname(script_path)
source_dir <- file.path(study_dir, "results_data", "source")
figure_dir <- file.path(study_dir, "figures")
table_dir <- file.path(study_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

figure_files <- c(
    "main_mechanism_separation.pdf",
    "main_abundance_power.pdf",
    "estimand_calibration.pdf",
    "component_specificity.pdf",
    "abundance_fdr_signal20.pdf",
    "abundance_fdr_signal40.pdf",
    "global_null_family_rejection_structural.pdf",
    "global_null_family_rejection_abundance.pdf",
    "correlated_community_abundance_power.pdf",
    "correlated_community_abundance_fdr.pdf",
    "correlated_community_structural_performance.pdf"
)
table_files <- c(
    "simulation_constants.csv",
    "simulation_design.csv",
    "method_targets.csv",
    "software_versions.csv",
    "confounding_design.csv",
    "test_availability.csv",
    "test_unavailability_reasons.csv",
    "estimand_calibration_summary.csv",
    "structural_estimand_summary.csv",
    "component_specificity_summary.csv",
    "abundance_performance_summary.csv",
    "global_null_family_rejection_summary.csv",
    "correlated_community_summary.csv"
)

input_files <- c(
    metrics = "combined_replication_metrics.csv",
    setting_summary = "combined_setting_summary.csv",
    truth_summary = "truth_summary.csv",
    failure_summary = "method_failure_summary.csv",
    simulation_design = "simulation_design.csv",
    confounding_summary = "confounding_design_summary.csv",
    package_versions = "package_versions.csv",
    completion_status = "setting_completion_status.csv"
)
input_paths <- stats::setNames(file.path(source_dir, input_files), names(input_files))
if (any(!file.exists(input_paths))) {
    stop(
        "Missing simulation result files: ",
        paste(names(input_paths)[!file.exists(input_paths)], collapse = ", "),
        call. = FALSE
    )
}

data.table::setDTthreads(0L)
metrics <- data.table::fread(input_paths[["metrics"]])
setting_summary <- data.table::fread(input_paths[["setting_summary"]])
truth_summary <- data.table::fread(input_paths[["truth_summary"]])
failure_summary <- data.table::fread(input_paths[["failure_summary"]])
simulation_design <- data.table::fread(input_paths[["simulation_design"]])
confounding_summary <- data.table::fread(input_paths[["confounding_summary"]])
package_versions <- data.table::fread(input_paths[["package_versions"]])
completion_status <- data.table::fread(input_paths[["completion_status"]])

expected_settings <- 540L
expected_replications <- 100L
if (
    data.table::uniqueN(completion_status$setting_id) != expected_settings ||
    data.table::uniqueN(completion_status$replication) != expected_replications ||
    nrow(completion_status) != expected_settings * expected_replications ||
    !all(completion_status$success)
) {
    stop("The simulation completion grid is incomplete.", call. = FALSE)
}

method_order <- c(
    "DASRA", "ZINQ", "MaAsLin 3", "ANCOM-BC2", "DESeq2",
    "edgeR", "LinDA", "corncob", "metagenomeSeq"
)
method_component_order <- c(
    "DASRA structural absence", "DASRA present-conditional abundance",
    "ZINQ Firth prevalence", "ZINQ quantile abundance",
    "MaAsLin 3 prevalence", "MaAsLin 3 abundance",
    "ANCOM-BC2", "DESeq2", "edgeR", "LinDA", "corncob",
    "metagenomeSeq"
)
method_colors <- c(
    "DASRA" = "#0072B2", "ZINQ" = "#D55E00",
    "MaAsLin 3" = "#009E73", "ANCOM-BC2" = "#CC79A7",
    "DESeq2" = "#56B4E9", "edgeR" = "#882255",
    "LinDA" = "#B79F00", "corncob" = "#000000",
    "metagenomeSeq" = "#777777"
)
method_shapes <- c(
    "DASRA" = 16, "ZINQ" = 17, "MaAsLin 3" = 15,
    "ANCOM-BC2" = 18, "DESeq2" = 3, "edgeR" = 7,
    "LinDA" = 8, "corncob" = 1, "metagenomeSeq" = 4
)
method_linetypes <- c(
    "DASRA" = 1, "ZINQ" = 2, "MaAsLin 3" = 4,
    "ANCOM-BC2" = 5, "DESeq2" = 6, "edgeR" = 3,
    "LinDA" = 7, "corncob" = 8, "metagenomeSeq" = 9
)
sample_colors <- c("60" = "#0072B2", "80" = "#D55E00", "120" = "#009E73")

short_method <- function(method) {
    data.table::fcase(
        method == "DASRA", "DASRA",
        method == "ZINQ", "ZINQ",
        method == "MaAsLin 3", "MaAsLin 3",
        default = method
    )
}

is_abundance_component <- function(component) {
    component %in% c(
        "present_conditional_abundance", "detected_quantile_abundance",
        "detected_log_abundance", "general_abundance"
    )
}

is_structural_component <- function(component) {
    component %in% c("structural_absence", "observed_prevalence")
}

short_method_label <- function(method_label) {
    data.table::fcase(
        grepl("^DASRA ", method_label), "DASRA",
        grepl("^ZINQ ", method_label), "ZINQ",
        grepl("^MaAsLin 3 ", method_label), "MaAsLin 3",
        default = method_label
    )
}

component_name <- function(component) {
    data.table::fcase(
        component == "structural_absence", "Structural absence",
        component == "observed_prevalence", "Observed prevalence",
        component == "present_conditional_abundance", "Present-conditional abundance",
        component == "detected_quantile_abundance", "Detected-sample quantile abundance",
        component == "detected_log_abundance", "Detected log abundance",
        component == "general_abundance", "General abundance",
        default = gsub("_", " ", component)
    )
}

prepare_design_fields <- function(data) {
    data <- data.table::copy(data)
    data[, method_short := factor(
        short_method(method), levels = method_order
    )]
    data[, n_label := factor(
        paste0("n = ", n_per_group),
        levels = paste0("n = ", c(60, 80, 120))
    )]
    data[, confounding_short := factor(
        confounding,
        levels = c("unconfounded", "confounded"),
        labels = c("Unconfounded", "Confounded")
    )]
    data[, design_column := factor(
        paste0(as.integer(round(100 * signal_fraction)), "%, ",
               ifelse(confounding == "unconfounded", "unconf.", "conf.")),
        levels = c("20%, unconf.", "20%, conf.",
                   "40%, unconf.", "40%, conf.")
    )]
    data[, design_tile := factor(
        paste0(n_per_group, " / ", as.integer(round(100 * signal_fraction)), "%"),
        levels = c("60 / 20%", "60 / 40%", "80 / 20%",
                   "80 / 40%", "120 / 20%", "120 / 40%")
    )]
    data
}

theme_simulation <- function(base_size = 8.5) {
    ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
        ggplot2::theme(
            plot.title = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(
                color = "#E5E5E5", linewidth = 0.28
            ),
            panel.border = ggplot2::element_rect(
                color = "#5A5A5A", linewidth = 0.35
            ),
            strip.background = ggplot2::element_rect(
                fill = "#F1F1F1", color = "#A0A0A0", linewidth = 0.3
            ),
            strip.text = ggplot2::element_text(size = 8, face = "plain"),
            axis.title = ggplot2::element_text(size = 9),
            axis.text = ggplot2::element_text(size = 7.5, color = "#222222"),
            legend.position = "bottom",
            legend.title = ggplot2::element_text(size = 8),
            legend.text = ggplot2::element_text(size = 7.5),
            legend.key.height = grid::unit(3.2, "mm"),
            legend.key.width = grid::unit(7.0, "mm"),
            plot.margin = ggplot2::margin(4, 5, 4, 5)
        )
}

method_scales <- function(nrow = 3L) {
    list(
        ggplot2::scale_color_manual(
            values = method_colors, breaks = method_order, drop = TRUE
        ),
        ggplot2::scale_shape_manual(
            values = method_shapes, breaks = method_order, drop = TRUE
        ),
        ggplot2::scale_linetype_manual(
            values = method_linetypes, breaks = method_order, drop = TRUE
        ),
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = nrow, byrow = TRUE),
            shape = ggplot2::guide_legend(nrow = nrow, byrow = TRUE),
            linetype = ggplot2::guide_legend(nrow = nrow, byrow = TRUE)
        )
    )
}

save_pdf <- function(plot, filename, width_mm, height_mm) {
    ggplot2::ggsave(
        filename = file.path(figure_dir, filename), plot = plot,
        device = grDevices::cairo_pdf, width = width_mm, height = height_mm,
        units = "mm", bg = "white"
    )
}

write_table <- function(data, filename) {
    data.table::fwrite(data, file.path(table_dir, filename), na = "")
}

metric_rows <- function(study_name, scenario_name, metric_name) {
    prepare_design_fields(setting_summary[
        study == study_name & scenario == scenario_name & metric == metric_name
    ])
}

calibration_long <- data.table::rbindlist(list(
    truth_summary[
        scenario == "observed_prevalence_only" & scenario_signal == TRUE,
        .(
            scenario = "Observed-prevalence perturbation",
            effect_parameter,
            quantity = "Observed prevalence",
            achieved = mean_standardized_expected_prevalence_difference
        )
    ],
    truth_summary[
        scenario == "observed_prevalence_only" & scenario_signal == TRUE,
        .(
            scenario = "Observed-prevalence perturbation",
            effect_parameter,
            quantity = "Structural absence",
            achieved = mean_structural_probability_difference
        )
    ],
    truth_summary[
        scenario == "structural_matched_prevalence" & scenario_signal == TRUE,
        .(
            scenario = "Structural perturbation with matched prevalence",
            effect_parameter,
            quantity = "Structural absence",
            achieved = mean_structural_probability_difference
        )
    ],
    truth_summary[
        scenario == "structural_matched_prevalence" & scenario_signal == TRUE,
        .(
            scenario = "Structural perturbation with matched prevalence",
            effect_parameter,
            quantity = "Observed prevalence",
            achieved = mean_standardized_expected_prevalence_difference
        )
    ]
))
calibration_summary <- calibration_long[, .(
    achieved_mean = mean(achieved),
    achieved_minimum = min(achieved),
    achieved_maximum = max(achieved),
    n_design_strata = .N
), by = .(scenario, effect_parameter, quantity)]
calibration_summary[, quantity := factor(
    quantity, levels = c("Observed prevalence", "Structural absence")
)]

calibration_plot <- ggplot2::ggplot(
    calibration_summary,
    ggplot2::aes(
        x = effect_parameter, y = achieved_mean,
        color = quantity, linetype = quantity, fill = quantity
    )
) +
    ggplot2::geom_abline(
        slope = 1, intercept = 0, color = "#777777",
        linetype = "dotted", linewidth = 0.45
    ) +
    ggplot2::geom_ribbon(
        ggplot2::aes(ymin = achieved_minimum, ymax = achieved_maximum),
        alpha = 0.10, color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.72) +
    ggplot2::geom_point(size = 1.25) +
    ggplot2::facet_wrap(~ scenario, nrow = 1, scales = "free_x") +
    ggplot2::scale_color_manual(values = c(
        "Observed prevalence" = "#0072B2", "Structural absence" = "#D55E00"
    )) +
    ggplot2::scale_fill_manual(values = c(
        "Observed prevalence" = "#0072B2", "Structural absence" = "#D55E00"
    )) +
    ggplot2::scale_linetype_manual(values = c(1, 2)) +
    ggplot2::labs(
        x = "Target difference",
        y = "Achieved absolute probability difference",
        color = NULL, fill = NULL, linetype = NULL
    ) +
    theme_simulation() +
    ggplot2::theme(legend.position = "bottom")
save_pdf(calibration_plot, "estimand_calibration.pdf", 178, 72)
calibration_output <- calibration_summary[, .(
    experiment = scenario,
    target_difference = effect_parameter,
    achieved_quantity = as.character(quantity),
    achieved_mean,
    achieved_minimum,
    achieved_maximum,
    design_strata = n_design_strata
)]
write_table(calibration_output, "estimand_calibration_summary.csv")

structural_rejection_data <- function(scenario_name) {
    data <- metric_rows(
        "structural_estimand", scenario_name,
        "designated_signal_bh_rejection"
    )
    data <- data[is_structural_component(component)]
    data[, method_short := factor(as.character(method_short), levels = method_order)]
    data
}

structural_observed <- structural_rejection_data("observed_prevalence_only")
structural_absence <- structural_rejection_data("structural_matched_prevalence")
structural_output <- data.table::rbindlist(
    list(structural_observed, structural_absence)
)[, .(
    experiment = ifelse(
        scenario == "observed_prevalence_only",
        "Observed-prevalence perturbation",
        "Structural perturbation with matched prevalence"
    ),
    samples_per_group = n_per_group,
    signal_taxa_percent = as.integer(round(100 * signal_fraction)),
    covariate_structure = as.character(confounding_short),
    target_difference = effect_parameter,
    method = as.character(method_short),
    tested_quantity = component_name(component),
    replications = n_replications,
    rejection_probability = mean,
    monte_carlo_sd = sd,
    monte_carlo_se = mcse,
    interval_lower = ci_lower,
    interval_upper = ci_upper
)]
write_table(structural_output, "structural_estimand_summary.csv")

main_mechanism <- data.table::rbindlist(list(
    data.table::copy(structural_observed)[, experiment :=
        "Observed-prevalence perturbation"],
    data.table::copy(structural_absence)[, experiment :=
        "Structural perturbation with matched prevalence"]
))
main_mechanism[, signal_label := factor(
    paste0(
        as.integer(round(50 * signal_fraction)), " taxa (",
        as.integer(round(100 * signal_fraction)), "%)"
    ),
    levels = c("10 taxa (20%)", "20 taxa (40%)")
)]
mechanism_method_labels <- c(
    "DASRA" = "DASRA: structural",
    "ZINQ" = "ZINQ: prevalence",
    "MaAsLin 3" = "MaAsLin 3: prevalence"
)

make_main_mechanism_panel <- function(experiment_name, x_label, panel_tag) {
    data <- main_mechanism[experiment == experiment_name]
    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = effect_parameter, y = mean, color = method_short,
            shape = method_short, linetype = signal_label,
            group = interaction(method_short, signal_label)
        )
    ) +
        ggplot2::geom_line(linewidth = 0.52) +
        ggplot2::geom_point(size = 0.90, stroke = 0.30) +
        ggplot2::facet_grid(confounding_short ~ n_label) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::scale_y_continuous(
            breaks = seq(0, 1, 0.25),
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(4)) +
        ggplot2::scale_color_manual(
            values = method_colors[c("DASRA", "ZINQ", "MaAsLin 3")],
            breaks = c("DASRA", "ZINQ", "MaAsLin 3"),
            labels = mechanism_method_labels
        ) +
        ggplot2::scale_shape_manual(
            values = method_shapes[c("DASRA", "ZINQ", "MaAsLin 3")],
            breaks = c("DASRA", "ZINQ", "MaAsLin 3"),
            labels = mechanism_method_labels
        ) +
        ggplot2::scale_linetype_manual(values = c(1, 2)) +
        ggplot2::labs(
            x = x_label, y = "Proportion of designated taxa rejected",
            color = NULL, shape = NULL, linetype = "Designated taxa",
            tag = panel_tag
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = 1, order = 1),
            shape = ggplot2::guide_legend(nrow = 1, order = 1),
            linetype = ggplot2::guide_legend(nrow = 1, order = 2)
        ) +
        theme_simulation() +
        ggplot2::theme(
            plot.tag = ggplot2::element_text(size = 10, face = "bold"),
            panel.spacing.x = grid::unit(2.6, "mm"),
            panel.spacing.y = grid::unit(2.0, "mm"),
            strip.text = ggplot2::element_text(size = 8.5),
            axis.text = ggplot2::element_text(size = 8),
            legend.text = ggplot2::element_text(size = 8),
            legend.title = ggplot2::element_text(size = 8)
        )
}

main_mechanism_plot <- make_main_mechanism_panel(
    "Observed-prevalence perturbation",
    "Target standardized observed-prevalence difference", "A"
) / make_main_mechanism_panel(
    "Structural perturbation with matched prevalence",
    "Target absolute standardized structural-absence probability difference", "B"
) + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(
        legend.position = "bottom", legend.box = "vertical",
        legend.box.just = "center"
    )
save_pdf(main_mechanism_plot, "main_mechanism_separation.pdf", 178, 184)

specificity_structural <- metric_rows(
        "abundance_comparison", "present_conditional_abundance",
        "designated_signal_bh_rejection"
    )[method == "DASRA" & component == "structural_absence"]
specificity_structural[, panel :=
    "Structural component under abundance perturbation"]
specificity_abundance <- metric_rows(
        "component_specificity", "structural_only",
        "designated_signal_bh_rejection"
    )[method == "DASRA" & component == "present_conditional_abundance"]
specificity_abundance[, panel :=
    "Abundance component under structural perturbation"]
specificity <- data.table::rbindlist(
    list(specificity_structural, specificity_abundance), fill = TRUE
)
specificity[, sample_size := factor(n_per_group, levels = c(60, 80, 120))]
specificity[, signal_label := factor(
    paste0(as.integer(round(50 * signal_fraction)), " designated taxa"),
    levels = c("10 designated taxa", "20 designated taxa")
)]
make_specificity_panel <- function(data, x_label, y_label, panel_tag) {
    ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = effect_parameter, y = mean, color = sample_size,
            linetype = signal_label,
            group = interaction(sample_size, signal_label)
        )
    ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = pmax(0, ci_lower), ymax = ci_upper),
            width = 0, linewidth = 0.24, alpha = 0.40
        ) +
        ggplot2::geom_line(linewidth = 0.58) +
        ggplot2::geom_point(size = 1.05) +
        ggplot2::facet_wrap(~ confounding_short, nrow = 1) +
        ggplot2::coord_cartesian(ylim = c(0, 0.015)) +
        ggplot2::scale_y_continuous(
            breaks = seq(0, 0.015, 0.005),
            labels = scales::label_number(accuracy = 0.001)
        ) +
        ggplot2::scale_color_manual(values = sample_colors) +
        ggplot2::scale_linetype_manual(values = c(1, 2)) +
        ggplot2::labs(
            x = x_label, y = y_label, tag = panel_tag,
            color = "Samples per group", linetype = NULL
        ) +
        theme_simulation() +
        ggplot2::theme(
            plot.tag = ggplot2::element_text(size = 10, face = "bold")
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = 1, order = 1),
            linetype = ggplot2::guide_legend(nrow = 1, order = 2)
        )
}

specificity_plot <- make_specificity_panel(
    specificity[panel == "Structural component under abundance perturbation"],
    "Target present-conditional mean log-relative-abundance difference",
    "Marginal Type I error", "A"
) / make_specificity_panel(
    specificity[panel == "Abundance component under structural perturbation"],
    "Target absolute standardized structural-absence probability difference",
    "Marginal Type I error", "B"
) + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
save_pdf(specificity_plot, "component_specificity.pdf", 178, 105)
specificity_output <- specificity[, .(
    perturbation = ifelse(
        panel == "Structural component under abundance perturbation",
        "Present-conditional abundance", "Structural absence"
    ),
    tested_component = ifelse(
        panel == "Structural component under abundance perturbation",
        "Structural absence", "Present-conditional abundance"
    ),
    samples_per_group = n_per_group,
    signal_taxa_percent = as.integer(round(100 * signal_fraction)),
    covariate_structure = as.character(confounding_short),
    target_effect = effect_parameter,
    replications = n_replications,
    marginal_type_i_error = mean,
    monte_carlo_sd = sd,
    monte_carlo_se = mcse,
    interval_lower = ci_lower,
    interval_upper = ci_upper
)]
write_table(specificity_output, "component_specificity_summary.csv")

abundance_metric_data <- function(signal_fraction_value, metric_name) {
    data <- metric_rows(
        "abundance_comparison", "present_conditional_abundance", metric_name
    )
    data <- data[
        is_abundance_component(component) & effect_parameter > 0 &
            abs(signal_fraction - signal_fraction_value) < 1e-12
    ]
    data
}

plot_abundance_fdr <- function(data, filename) {
    signal_text <- paste0(
        as.integer(round(100 * unique(data$signal_fraction))), "% signal taxa"
    )
    data[, confounding_signal := factor(
        paste(as.character(confounding_short), signal_text, sep = "\n"),
        levels = paste(c("Unconfounded", "Confounded"), signal_text, sep = "\n")
    )]
    plot <- ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = effect_parameter, y = mean, color = method_short,
            shape = method_short, linetype = method_short,
            group = method_short
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0.05, color = "#555555",
            linetype = "dotted", linewidth = 0.38
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = pmax(0, ci_lower), ymax = pmin(0.40, ci_upper)),
            width = 0, linewidth = 0.18, alpha = 0.32
        ) +
        ggplot2::geom_line(linewidth = 0.48) +
        ggplot2::geom_point(size = 0.82, stroke = 0.35) +
        ggplot2::facet_grid(n_label ~ confounding_signal) +
        ggplot2::coord_cartesian(ylim = c(0, 0.40)) +
        ggplot2::scale_y_continuous(
            breaks = seq(0, 0.4, 0.1),
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::scale_x_continuous(
            breaks = c(0.1, 0.4, 0.8, 1.2, 1.4)
        ) +
        ggplot2::labs(
            x = "Absolute present-conditional mean log-relative-abundance difference",
            y = "Empirical false discovery rate",
            color = NULL, shape = NULL, linetype = NULL
        ) +
        method_scales(nrow = 3L) +
        theme_simulation()
    save_pdf(plot, filename, 178, 112)
}

abundance_tables <- list(
    abundance_metric_data(0.20, "power"),
    abundance_metric_data(0.40, "power"),
    abundance_metric_data(0.20, "fdp"),
    abundance_metric_data(0.40, "fdp")
)
plot_abundance_fdr(
    data.table::copy(abundance_tables[[3L]]), "abundance_fdr_signal20.pdf"
)
plot_abundance_fdr(
    data.table::copy(abundance_tables[[4L]]), "abundance_fdr_signal40.pdf"
)
abundance_output <- data.table::rbindlist(abundance_tables)[, .(
    samples_per_group = n_per_group,
    signal_taxa_percent = as.integer(round(100 * signal_fraction)),
    covariate_structure = as.character(confounding_short),
    target_abundance_difference = effect_parameter,
    method = as.character(method_short),
    measure = ifelse(
        metric == "power",
        "Power after Benjamini-Hochberg adjustment",
        "Empirical false discovery rate"
    ),
    replications = n_replications,
    estimate = mean,
    monte_carlo_sd = sd,
    monte_carlo_se = mcse,
    interval_lower = ci_lower,
    interval_upper = ci_upper
)]
write_table(abundance_output, "abundance_performance_summary.csv")

main_abundance <- data.table::rbindlist(abundance_tables[1:2])
main_abundance[, design_row := factor(
    paste0(
        as.integer(round(50 * signal_fraction)), " taxa (",
        as.integer(round(100 * signal_fraction)), "%)\n",
        as.character(confounding_short)
    ),
    levels = c(
        "10 taxa (20%)\nUnconfounded",
        "10 taxa (20%)\nConfounded",
        "20 taxa (40%)\nUnconfounded",
        "20 taxa (40%)\nConfounded"
    )
)]

plot_main_abundance <- function(filename) {
    data <- main_abundance[metric == "power"]
    plot <- ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = effect_parameter, y = mean, color = method_short,
            shape = method_short, linetype = method_short,
            group = method_short
        )
    ) +
        ggplot2::geom_line(linewidth = 0.48) +
        ggplot2::geom_point(size = 0.80, stroke = 0.30) +
        ggplot2::facet_grid(design_row ~ n_label) +
        ggplot2::coord_cartesian(xlim = c(0.1, 1.4), ylim = c(0, 1)) +
        ggplot2::scale_x_continuous(
            breaks = c(0.1, 0.4, 0.8, 1.2, 1.4),
            expand = ggplot2::expansion(mult = c(0, 0.01))
        ) +
        ggplot2::scale_y_continuous(
            breaks = seq(0, 1, 0.25),
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
            x = "Absolute present-conditional mean log-relative-abundance difference",
            y = "Power", color = NULL, shape = NULL, linetype = NULL
        ) +
        method_scales(nrow = 3L) +
        theme_simulation() +
        ggplot2::theme(
            panel.spacing.x = grid::unit(5.5, "mm"),
            panel.spacing.y = grid::unit(1.8, "mm"),
            strip.text.x = ggplot2::element_text(size = 8.5),
            strip.text.y = ggplot2::element_text(angle = 0, size = 8),
            axis.text = ggplot2::element_text(size = 8),
            legend.text = ggplot2::element_text(size = 8)
        )
    save_pdf(plot, filename, 178, 172)
}

plot_main_abundance("main_abundance_power.pdf")

null_scenario_levels <- c(
    "Taxonwise global null", "Depth-imbalanced null",
    "Correlated-community null"
)
null_metrics <- prepare_design_fields(metrics[
    (scenario == "present_conditional_abundance" & effect_parameter == 0) |
        scenario %in% c("depth_imbalanced_null", "joint_global_null")
])
null_metrics[, null_scenario := factor(
    data.table::fcase(
        scenario == "present_conditional_abundance", "Taxonwise global null",
        scenario == "depth_imbalanced_null", "Depth-imbalanced null",
        default = "Correlated-community null"
    ),
    levels = null_scenario_levels
)]

summarize_null_family <- function(data) {
    per_replication <- data[, .(
        value = mean(fdp, na.rm = TRUE)
    ), by = .(
        replication, n_per_group, n_label, confounding_short,
        null_scenario, method_short, method_label, component
    )]
    per_replication[, .(
        n_replications = .N,
        mean = mean(value),
        sd = stats::sd(value),
        mcse = stats::sd(value) / sqrt(.N),
        ci_lower = pmax(0, mean(value) - 1.96 * stats::sd(value) / sqrt(.N)),
        ci_upper = pmin(
            1, mean(value) + 1.96 * stats::sd(value) / sqrt(.N)
        )
    ), by = .(
        n_per_group, n_label, confounding_short, null_scenario,
        method_short, method_label, component
    )]
}

plot_null_family_rejection <- function(family, filename, y_limit) {
    data <- if (family == "abundance") {
        null_metrics[is_abundance_component(component)]
    } else {
        null_metrics[is_structural_component(component)]
    }
    summary <- summarize_null_family(data)
    plot <- ggplot2::ggplot(
        summary,
        ggplot2::aes(
            x = factor(n_per_group), y = mean, color = method_short,
            shape = method_short, linetype = method_short,
            group = method_short
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0.05, color = "#555555",
            linetype = "dotted", linewidth = 0.38
        ) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = ci_lower, ymax = pmin(y_limit, ci_upper)),
            width = 0.08, linewidth = 0.20, alpha = 0.40
        ) +
        ggplot2::geom_line(linewidth = 0.50) +
        ggplot2::geom_point(size = 0.95) +
        ggplot2::facet_grid(confounding_short ~ null_scenario) +
        ggplot2::coord_cartesian(ylim = c(0, y_limit)) +
        ggplot2::scale_y_continuous(
            breaks = scales::breaks_width(if (family == "abundance") 0.25 else 0.02),
            labels = scales::label_number(accuracy = 0.01),
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
            x = "Samples per group",
            y = "Family-wise Type I error",
            color = NULL, shape = NULL, linetype = NULL
        ) +
        method_scales(nrow = if (family == "abundance") 3L else 1L) +
        theme_simulation()
    save_pdf(
        plot, filename, 178,
        if (family == "abundance") 92 else 78
    )
    summary[, component_family := family]
}

null_abundance <- plot_null_family_rejection(
    "abundance", "global_null_family_rejection_abundance.pdf", 1.00
)
null_structural <- plot_null_family_rejection(
    "structural", "global_null_family_rejection_structural.pdf", 0.12
)
null_output <- data.table::rbindlist(
    list(null_abundance, null_structural), fill = TRUE
)[, .(
    component_family = ifelse(
        component_family == "abundance", "Abundance", "Structural or prevalence"
    ),
    null_mechanism = as.character(null_scenario),
    samples_per_group = n_per_group,
    covariate_structure = as.character(confounding_short),
    method = as.character(method_short),
    method_component = method_label,
    tested_component = component_name(component),
    replications = n_replications,
    familywise_type_i_error = mean,
    monte_carlo_sd = sd,
    monte_carlo_se = mcse,
    interval_lower = ci_lower,
    interval_upper = ci_upper
)]
write_table(null_output, "global_null_family_rejection_summary.csv")

joint_abundance <- prepare_design_fields(setting_summary[
    study == "joint_robustness" & scenario == "joint_abundance" &
        metric %in% c("power", "fdp") & is_abundance_component(component)
])
joint_abundance[, effect_label := factor(
    paste0("Effect = ", format(effect_parameter, nsmall = 1)),
    levels = c("Effect = 0.4", "Effect = 0.8")
)]

plot_joint_abundance <- function(metric_name, filename, y_limit) {
    data <- joint_abundance[metric == metric_name]
    plot <- ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = n_per_group, y = mean, color = method_short,
            shape = method_short, linetype = method_short,
            group = method_short
        )
    )
    if (metric_name == "fdp") {
        plot <- plot + ggplot2::geom_hline(
            yintercept = 0.05, color = "#555555",
            linetype = "dotted", linewidth = 0.38
        )
    }
    plot <- plot +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = pmax(0, ci_lower),
                ymax = pmin(y_limit, ci_upper)
            ),
            width = 2.0, linewidth = 0.20, alpha = 0.38
        ) +
        ggplot2::geom_line(linewidth = 0.52) +
        ggplot2::geom_point(size = 0.90, stroke = 0.30) +
        ggplot2::facet_grid(effect_label ~ design_column) +
        ggplot2::coord_cartesian(ylim = c(0, y_limit)) +
        ggplot2::scale_x_continuous(breaks = c(60, 80, 120)) +
        ggplot2::scale_y_continuous(
            breaks = if (metric_name == "power") {
                seq(0, 1, 0.25)
            } else {
                seq(0, y_limit, 0.05)
            },
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::labs(
            x = "Samples per group",
            y = if (metric_name == "power") {
                "Power"
            } else {
                "Empirical false discovery rate"
            },
            color = NULL, shape = NULL, linetype = NULL
        ) +
        method_scales(nrow = 3L) +
        theme_simulation() +
        ggplot2::theme(
            panel.spacing.x = grid::unit(2.2, "mm"),
            panel.spacing.y = grid::unit(1.8, "mm")
        )
    save_pdf(plot, filename, 178, 100)
}

plot_joint_abundance(
    "power", "correlated_community_abundance_power.pdf", 1.00
)
plot_joint_abundance(
    "fdp", "correlated_community_abundance_fdr.pdf", 0.20
)

joint_structural <- prepare_design_fields(setting_summary[
    study == "joint_robustness" & scenario == "joint_structural" &
        metric %in% c("power", "fdp") &
        is_structural_component(component)
])
joint_structural[, effect_label := factor(
    paste0("Difference = ", format(effect_parameter, nsmall = 2)),
    levels = c("Difference = 0.20", "Difference = 0.35")
)]

joint_structural_key <- c(
    "n_per_group", "signal_fraction", "confounding", "effect_parameter",
    "method", "component", "metric"
)
valid_joint_structural_component <-
    (joint_structural$method == "DASRA" &
        joint_structural$component == "structural_absence") |
    (joint_structural$method %in% c("ZINQ", "MaAsLin 3") &
        joint_structural$component == "observed_prevalence")
if (
    nrow(joint_structural) != 144L ||
    data.table::uniqueN(joint_structural, by = joint_structural_key) != 144L ||
    !all(valid_joint_structural_component) ||
    data.table::uniqueN(
        joint_structural[, .(method, component)]
    ) != 3L
) {
    stop(
        "The correlated-community structural summary grid is incomplete.",
        call. = FALSE
    )
}

plot_joint_structural_metric <- function(
        metric_name, panel_tag, y_limit) {
    data <- joint_structural[metric == metric_name]
    method_breaks <- c("DASRA", "ZINQ", "MaAsLin 3")
    method_labels <- c(
        "DASRA: structural", "ZINQ: prevalence",
        "MaAsLin 3: prevalence"
    )
    dodge <- ggplot2::position_dodge(width = 3.2)
    plot <- ggplot2::ggplot(
        data,
        ggplot2::aes(
            x = n_per_group, y = mean, color = method_short,
            shape = method_short, linetype = method_short,
            group = method_short
        )
    )
    if (metric_name == "fdp") {
        plot <- plot + ggplot2::geom_hline(
            yintercept = 0.05, color = "#555555",
            linetype = "dotted", linewidth = 0.38
        )
    }
    plot +
        ggplot2::geom_errorbar(
            ggplot2::aes(
                ymin = pmax(0, ci_lower),
                ymax = pmin(y_limit, ci_upper)
            ),
            width = 2.0, linewidth = 0.20, alpha = 0.42,
            position = dodge
        ) +
        ggplot2::geom_line(linewidth = 0.54, position = dodge) +
        ggplot2::geom_point(
            size = 1.05, stroke = 0.32, position = dodge
        ) +
        ggplot2::facet_grid(effect_label ~ design_column) +
        ggplot2::coord_cartesian(ylim = c(0, y_limit)) +
        ggplot2::scale_x_continuous(breaks = c(60, 80, 120)) +
        ggplot2::scale_y_continuous(
            breaks = if (metric_name == "power") {
                seq(0, 1, 0.25)
            } else {
                seq(0, 0.10, 0.025)
            },
            expand = ggplot2::expansion(mult = c(0, 0.02))
        ) +
        ggplot2::scale_color_manual(
            values = method_colors[method_breaks], breaks = method_breaks,
            labels = method_labels
        ) +
        ggplot2::scale_shape_manual(
            values = method_shapes[method_breaks], breaks = method_breaks,
            labels = method_labels
        ) +
        ggplot2::scale_linetype_manual(
            values = method_linetypes[method_breaks], breaks = method_breaks,
            labels = method_labels
        ) +
        ggplot2::guides(
            color = ggplot2::guide_legend(nrow = 1),
            shape = ggplot2::guide_legend(nrow = 1),
            linetype = ggplot2::guide_legend(nrow = 1)
        ) +
        ggplot2::labs(
            x = if (metric_name == "power") NULL else "Samples per group",
            y = if (metric_name == "power") {
                "Power"
            } else {
                "Empirical false discovery rate"
            },
            color = NULL, shape = NULL, linetype = NULL, tag = panel_tag
        ) +
        theme_simulation() +
        ggplot2::theme(
            panel.spacing.x = grid::unit(2.2, "mm"),
            panel.spacing.y = grid::unit(1.8, "mm"),
            plot.tag = ggplot2::element_text(size = 10, face = "bold"),
            strip.text = ggplot2::element_text(size = 7.8),
            legend.text = ggplot2::element_text(size = 7.8)
        )
}

joint_structural_performance <-
    plot_joint_structural_metric("power", "A", 1.00) /
    plot_joint_structural_metric("fdp", "B", 0.10) +
    patchwork::plot_layout(guides = "collect", heights = c(1, 1)) &
    ggplot2::theme(
        legend.position = "bottom", legend.box = "horizontal",
        legend.justification = "center"
    )
save_pdf(
    joint_structural_performance,
    "correlated_community_structural_performance.pdf", 178, 170
)

joint_output <- data.table::rbindlist(
    list(joint_abundance, joint_structural), fill = TRUE
)[, .(
    experiment = ifelse(
        scenario == "joint_abundance", "Abundance perturbation",
        "Structural perturbation"
    ),
    samples_per_group = n_per_group,
    signal_taxa_percent = as.integer(round(100 * signal_fraction)),
    covariate_structure = as.character(confounding_short),
    target_effect = effect_parameter,
    method = as.character(method_short),
    measure = data.table::fcase(
        metric == "power", "Power after Benjamini-Hochberg adjustment",
        metric == "fdp", "Empirical false discovery rate",
        default = "Power after Benjamini-Hochberg adjustment"
    ),
    replications = n_replications,
    estimate = mean,
    monte_carlo_sd = sd,
    monte_carlo_se = mcse,
    interval_lower = ci_lower,
    interval_upper = ci_upper
)]
write_table(joint_output, "correlated_community_summary.csv")

scenario_descriptions <- data.table::data.table(
    study = c(
        "Structural estimand", "Structural estimand", "Abundance comparison",
        "Component specificity", "Depth calibration", "Joint robustness",
        "Joint robustness", "Joint robustness"
    ),
    scenario = c(
        "Observed-prevalence perturbation", "Structural perturbation with matched prevalence",
        "Present-conditional abundance", "Structural-only perturbation",
        "Depth-imbalanced global null", "Correlated-community global null",
        "Correlated-community abundance", "Correlated-community structural"
    ),
    data_generating_process = c(
        rep("Taxonwise two-part model", 5), rep("Correlated lognormal community", 3)
    ),
    target_quantity = c(
        "Standardized expected observed-prevalence difference",
        "Absolute standardized structural-absence probability difference",
        "Present-conditional mean log-relative-abundance difference",
        "Absolute standardized structural-absence probability difference",
        "Global null with fourfold median-depth difference",
        "Global null", "Direct log-absolute-abundance effect",
        "Absolute standardized structural-absence probability difference"
    ),
    source_scenario = c(
        "observed_prevalence_only", "structural_matched_prevalence",
        "present_conditional_abundance", "structural_only",
        "depth_imbalanced_null", "joint_global_null", "joint_abundance",
        "joint_structural"
    )
)
effect_grids <- simulation_design[, .(
    effect_grid = paste(
        format(sort(unique(effect_parameter)), trim = TRUE, scientific = FALSE),
        collapse = ", "
    )
), by = scenario]
data.table::setnames(effect_grids, "scenario", "source_scenario")
scenario_descriptions[, row_order := .I]
scenario_descriptions <- merge(
    scenario_descriptions, effect_grids,
    by = "source_scenario", all.x = TRUE, sort = FALSE
)
data.table::setorder(scenario_descriptions, row_order)
scenario_descriptions[, c("source_scenario", "row_order") := NULL]
data.table::setcolorder(
    scenario_descriptions,
    c("study", "scenario", "data_generating_process", "target_quantity", "effect_grid")
)
write_table(scenario_descriptions, "simulation_design.csv")

method_targets <- data.table::data.table(
    method = c(
        "DASRA", "DASRA", "ZINQ", "ZINQ", "MaAsLin 3", "MaAsLin 3",
        "ANCOM-BC2", "DESeq2", "edgeR", "LinDA", "corncob", "metagenomeSeq"
    ),
    component = c(
        "Structural absence", "Present-conditional abundance",
        "Observed prevalence", "Detected-sample quantile abundance",
        "Observed prevalence", "Detected log abundance",
        rep("General abundance", 6)
    ),
    simulation_target = c(
        "Structural-absence probability", "Present-conditional relative abundance",
        "Observed prevalence", "Detected-sample abundance quantiles",
        "Observed prevalence", "Detected log relative abundance",
        rep("Directly perturbed abundance taxa", 6)
    ),
    multiplicity_family = "Fifty focal taxa within each method-component setting",
    analysis_configuration = c(
        "group + z; original library size supplied",
        "group + z; original library size supplied",
        "group + z + standardized log depth; Firth model",
        "group + z + standardized log depth; three quantiles",
        "TSS and log transform; z and standardized log depth",
        "TSS and log transform; z and standardized log depth",
        "group + z; ANCOM-BC2 structural-zero and sensitivity procedures",
        "z + group; positive-count size factors; Wald test",
        "group + z; filterByExpr; TMM; quasi-likelihood test",
        "group + z; count input; adaptive zero handling",
        "group + z mean and dispersion models; robust Wald test",
        "group + z; CSS normalization; fitZig"
    )
)
write_table(method_targets, "method_targets.csv")

simulation_constants <- data.table::data.table(
    design_feature = c(
        "Focal taxa", "Samples per group", "Signal taxa",
        "Covariate structures", "Monte Carlo replications",
        "Base settings", "Complete design settings", "Nominal level",
        "Multiplicity adjustment"
    ),
    value = c(
        "50", "60, 80, 120", "10 (20%) or 20 (40%)",
        "Unconfounded and confounded", "100", "45", "540", "0.05",
        "Benjamini-Hochberg within each 50-taxon method-component family"
    )
)
write_table(simulation_constants, "simulation_constants.csv")

status_labels <- c(
    "conditional_present_information_nonpositive" = "Nonpositive conditional information",
    "conditional_present_root_polish_not_closed" = "Root refinement did not close",
    "conditional_present_root_polish_no_descent" = "Root refinement found no descent",
    "conditional_present_gradient" = "Gradient criterion not met",
    "conditional_present_persistent_boundary" = "Persistent conditional boundary",
    "background_mode_not_converged" = "Background mode did not converge",
    "structural_absence_nonoptimal_nuisance_fit" = "Structural nuisance fit not optimal",
    "structural_absence_persistent_boundary" = "Persistent structural boundary",
    "inference_derivative_unstable" = "Unstable inference derivative",
    "filterByExpr_excluded" = "Excluded by edgeR filterByExpr"
)
availability <- failure_summary[, .(
    eligible_tests = sum(N),
    available_tests = sum(N[available]),
    unavailable_tests = sum(N[!available])
), by = .(method_label, component)]
availability[, availability_percent := 100 * available_tests / eligible_tests]
availability[, method_short := short_method_label(method_label)]
availability[, method_order_index := match(method_short, method_order)]
data.table::setorder(availability, method_order_index, component)
availability_output <- availability[, .(
    method = method_label,
    component = component_name(component),
    eligible_tests,
    available_tests,
    unavailable_tests,
    availability_percent = round(availability_percent, 3)
)]
availability_output[, method_order_index := match(method, method_component_order)]
data.table::setorder(availability_output, method_order_index)
availability_output[, method_order_index := NULL]
write_table(availability_output, "test_availability.csv")

unavailability <- failure_summary[available == FALSE, .(
    unavailable_tests = sum(N)
), by = .(method_label, component, status)]
unavailability[, reason := unname(status_labels[status])]
unavailability[is.na(reason), reason := gsub("_", " ", status)]
unavailability[, percent_within_component :=
    100 * unavailable_tests / sum(unavailable_tests),
    by = .(method_label, component)
]
data.table::setorder(unavailability, method_label, -unavailable_tests)
unavailability_output <- unavailability[, .(
    method = method_label,
    component = component_name(component),
    reason,
    unavailable_tests,
    percent_of_component_unavailability = round(percent_within_component, 3)
)]
unavailability_output[, method_order_index := match(method, method_component_order)]
data.table::setorder(
    unavailability_output, method_order_index, -unavailable_tests
)
unavailability_output[, method_order_index := NULL]
write_table(unavailability_output, "test_unavailability_reasons.csv")

confounding_table <- confounding_summary[, .(
    samples_per_group = n_per_group,
    covariate_structure = confounding_label,
    replications = n_replications,
    mean_group_z_correlation = round(mean_group_z_correlation, 6),
    sd_group_z_correlation = round(sd_group_z_correlation, 6),
    mean_z_mean_difference = round(mean_z_mean_difference, 6),
    sd_z_mean_difference = round(sd_z_mean_difference, 6)
)]
confounding_table[, confounding_order := match(
    covariate_structure, c("Unconfounded", "Confounded")
)]
data.table::setorder(confounding_table, confounding_order, samples_per_group)
confounding_table[, confounding_order := NULL]
write_table(confounding_table, "confounding_design.csv")

software_versions <- package_versions[, .(package, version)]
write_table(software_versions, "software_versions.csv")

output_paths <- c(
    file.path(figure_dir, figure_files),
    file.path(table_dir, table_files)
)
if (any(!file.exists(output_paths))) {
    stop(
        "Expected outputs were not created: ",
        paste(basename(output_paths[!file.exists(output_paths)]), collapse = ", "),
        call. = FALSE
    )
}
message(
    "Created ", length(figure_files), " figures and ",
    length(table_files), " CSV tables in ", study_dir
)
