# Plot-data preparation ----------------------------------------------------

.dasra_plot_profile_failure <- function(status) {
    list(
        reference = NA_real_,
        comparison = NA_real_,
        status = as.character(status)[1L]
    )
}

.dasra_structural_plot_profile <- function(y, N, fit) {
    if (is.null(fit) || !isTRUE(fit$tested) || !isTRUE(fit$regular)) {
        reason <- if (is.null(fit$reason)) "unavailable" else fit$reason
        return(.dasra_plot_profile_failure(reason))
    }
    stored <- fit$diagnostics$fit
    if (is.null(stored) || is.null(stored$conditional_present$par) ||
        is.null(stored$X_eta) || is.null(stored$quadrature)) {
        return(.dasra_plot_profile_failure("fit_details_unavailable"))
    }

    X_eta <- stored$X_eta
    if (!("Group" %in% colnames(X_eta))) {
        return(.dasra_plot_profile_failure("group_design_unavailable"))
    }

    companion <- tryCatch(
        suppressWarnings(
            zt_fit_alpha(
                beta = stored$conditional_present$par,
                y = y,
                N = N,
                X_rho = X_eta,
                X_eta = X_eta,
                gh = stored$quadrature
            )
        ),
        error = function(e) NULL
    )
    if (is.null(companion) || !isTRUE(companion$ok)) {
        reason <- if (is.null(companion$reason)) {
            "companion_fit_unavailable"
        } else {
            companion$reason
        }
        return(.dasra_plot_profile_failure(reason))
    }

    components <- zt_beta_detection_components(
        stored$conditional_present$par,
        N,
        X_eta,
        stored$quadrature
    )
    if (is.null(components)) {
        return(.dasra_plot_profile_failure(
            "companion_detection_components_unavailable"
        ))
    }
    finite_nll <- tryCatch(
        companion$nll(companion$par),
        error = function(e) NA_real_
    )
    if (!is.finite(finite_nll)) {
        return(.dasra_plot_profile_failure(
            "companion_objective_unavailable"
        ))
    }
    audit <- zt_structural_zero_limit_audit(
        finite_nll = finite_nll,
        log_r = components$log_r,
        y = y
    )
    if (isTRUE(audit$dominated)) {
        return(.dasra_plot_profile_failure(
            "companion_zero_limit_dominated"
        ))
    }

    X_reference <- X_comparison <- X_eta
    X_reference[, "Group"] <- 0
    X_comparison[, "Group"] <- 1
    reference <- mean(stats::plogis(
        as.numeric(X_reference %*% companion$par)
    ))
    comparison <- mean(stats::plogis(
        as.numeric(X_comparison %*% companion$par)
    ))
    if (!is.finite(reference) || !is.finite(comparison)) {
        return(.dasra_plot_profile_failure(
            "companion_standardization_unavailable"
        ))
    }

    list(reference = reference, comparison = comparison, status = "ok")
}

.dasra_abundance_plot_profile <- function(taxon, abundance_fits) {
    if (is.null(abundance_fits$taxon) ||
        is.null(abundance_fits$raw_fits)) {
        return(.dasra_plot_profile_failure("fit_details_unavailable"))
    }
    detail_index <- match(taxon, abundance_fits$taxon$taxon)
    raw_fit <- abundance_fits$raw_fits[[taxon]]
    if (is.na(detail_index) || is.null(raw_fit) ||
        !isTRUE(raw_fit$available) ||
        !isTRUE(abundance_fits$taxon$formed[[detail_index]])) {
        reason <- if (is.na(detail_index)) {
            "fit_details_unavailable"
        } else {
            abundance_fits$taxon$reason[[detail_index]]
        }
        return(.dasra_plot_profile_failure(reason))
    }

    theta <- raw_fit$theta
    X_eta <- raw_fit$solver_diagnostics$X_eta
    if (is.null(theta) || is.null(X_eta) ||
        !("Group" %in% colnames(X_eta))) {
        return(.dasra_plot_profile_failure("group_design_unavailable"))
    }
    p_eta <- ncol(X_eta)
    if (length(theta) != p_eta + 1L) {
        return(.dasra_plot_profile_failure("parameter_layout_unavailable"))
    }

    sigma <- exp(theta[[p_eta + 1L]])
    gh <- .dasra_make_abundance_gh_rule(raw_fit$quadrature_Q)
    X_reference <- X_comparison <- X_eta
    X_reference[, "Group"] <- 0
    X_comparison[, "Group"] <- 1
    location_reference <- as.numeric(
        X_reference %*% theta[seq_len(p_eta)]
    )
    location_comparison <- as.numeric(
        X_comparison %*% theta[seq_len(p_eta)]
    )
    reference <- 100 * exp(mean(.dasra_abundance_mean_log_relative(
        location_reference, sigma, gh
    )))
    comparison <- 100 * exp(mean(.dasra_abundance_mean_log_relative(
        location_comparison, sigma, gh
    )))
    if (!is.finite(reference) || !is.finite(comparison) ||
        reference <= 0 || comparison <= 0) {
        return(.dasra_plot_profile_failure(
            "abundance_standardization_unavailable"
        ))
    }

    list(reference = reference, comparison = comparison, status = "ok")
}

.dasra_prepare_plot_data <- function(Y, N, results, structural_fits,
                                     abundance_fits, contrast) {
    taxa <- results$taxon
    structural <- lapply(seq_along(taxa), function(j) {
        tryCatch(
            .dasra_structural_plot_profile(
                y = Y[, j],
                N = N,
                fit = structural_fits[[taxa[[j]]]]
            ),
            error = function(e) .dasra_plot_profile_failure(
                paste0("plot_profile_error: ", conditionMessage(e))
            )
        )
    })
    abundance <- lapply(taxa, function(taxon) {
        tryCatch(
            .dasra_abundance_plot_profile(taxon, abundance_fits),
            error = function(e) .dasra_plot_profile_failure(
                paste0("plot_profile_error: ", conditionMessage(e))
            )
        )
    })

    profiles <- data.frame(
        feature = taxa,
        structural_reference = vapply(
            structural, `[[`, numeric(1), "reference"
        ),
        structural_comparison = vapply(
            structural, `[[`, numeric(1), "comparison"
        ),
        structural_status = vapply(
            structural, `[[`, character(1), "status"
        ),
        abundance_reference = vapply(
            abundance, `[[`, numeric(1), "reference"
        ),
        abundance_comparison = vapply(
            abundance, `[[`, numeric(1), "comparison"
        ),
        abundance_status = vapply(
            abundance, `[[`, character(1), "status"
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    list(
        version = 1L,
        contrast = contrast,
        profiles = profiles
    )
}

# Plot specification -------------------------------------------------------

.dasra_plot_stars <- function(p) {
    if (!is.finite(p) || p > 0.05) return("")
    if (p <= 0.001) return("***")
    if (p <= 0.01) return("**")
    "*"
}

.dasra_plot_adjustment_label <- function(method) {
    if (identical(method, "none")) "p" else paste(method, "adjusted p")
}

.dasra_plot_validate_integer <- function(value, name) {
    if (length(value) != 1L || is.na(value) || !is.numeric(value) ||
        !is.finite(value) || value < 1 || value != as.integer(value)) {
        stop(sprintf("`%s` must be one positive integer.", name),
             call. = FALSE)
    }
    as.integer(value)
}

.dasra_plot_select <- function(x, features, selection, max_features,
                               alpha) {
    results <- x$results
    diagnostics <- x$diagnostics
    eligible <- diagnostics$retained & diagnostics$formed_omnibus &
        (is.finite(results$z_structural_absence) |
             is.finite(results$z_relative_abundance)) &
        is.finite(results$p_adj_omnibus)

    if (!is.null(features)) {
        if (!is.character(features) || anyNA(features) ||
            any(!nzchar(features)) || anyDuplicated(features)) {
            stop(
                "`features` must contain unique, non-empty feature names.",
                call. = FALSE
            )
        }
        missing_features <- setdiff(features, results$taxon)
        if (length(missing_features)) {
            stop(
                sprintf(
                    "Unknown feature name(s): %s.",
                    paste(missing_features, collapse = ", ")
                ),
                call. = FALSE
            )
        }
        indices <- match(features, results$taxon)
        invalid <- !eligible[indices]
        if (any(invalid)) {
            stop(
                sprintf(
                    "Feature(s) without a formed omnibus result: %s.",
                    paste(features[invalid], collapse = ", ")
                ),
                call. = FALSE
            )
        }
        return(indices)
    }

    indices <- which(eligible)
    if (identical(selection, "significant")) {
        indices <- indices[results$p_adj_omnibus[indices] <= alpha]
        if (!length(indices)) {
            stop(
                paste(
                    "No features meet the adjusted-p threshold; use",
                    "`selection = \"top\"` or supply `features`."
                ),
                call. = FALSE
            )
        }
    }
    if (!length(indices)) {
        stop("No formed omnibus results are available to plot.",
             call. = FALSE)
    }
    indices <- indices[order(
        results$p_adj_omnibus[indices],
        results$p_omnibus[indices],
        indices,
        method = "radix"
    )]
    if (!identical(selection, "all")) {
        indices <- utils::head(indices, max_features)
    }
    indices
}

.dasra_plot_color_limits <- function(data, p_color_limits) {
    if (is.character(p_color_limits)) {
        if (length(p_color_limits) != 1L ||
            !identical(p_color_limits, "adaptive")) {
            stop(
                "`p_color_limits` must be \"adaptive\" or numeric c(lower, upper).",
                call. = FALSE
            )
        }
        component_p <- c(
            data$p_structural[
                is.finite(data$z_structural)
            ],
            data$p_abundance[
                is.finite(data$z_abundance)
            ]
        )
        component_p <- component_p[
            is.finite(component_p) & component_p >= 0 & component_p <= 1
        ]
        positive <- component_p[component_p > 0]
        if (length(positive)) {
            lower <- max(
                .Machine$double.xmin,
                10 ^ floor(log10(min(positive)))
            )
            if (lower >= 1) lower <- 0.1
        } else if (length(component_p)) {
            lower <- 1e-16
        } else {
            stop("No finite component adjusted p-values are available.",
                 call. = FALSE)
        }
        return(c(lower = lower, upper = 1))
    }

    if (!is.numeric(p_color_limits) || length(p_color_limits) != 2L ||
        anyNA(p_color_limits) || any(!is.finite(p_color_limits)) ||
        p_color_limits[[1L]] <= 0 || p_color_limits[[1L]] >=
            p_color_limits[[2L]] || p_color_limits[[2L]] > 1) {
        stop(
            paste(
                "Numeric `p_color_limits` must satisfy",
                "0 < lower < upper <= 1."
            ),
            call. = FALSE
        )
    }
    setNames(as.numeric(p_color_limits), c("lower", "upper"))
}

.dasra_plot_palette_function <- function(palette, limits) {
    if (!is.character(palette) || length(palette) < 2L || anyNA(palette)) {
        stop("`evidence_palette` must contain at least two colors.",
             call. = FALSE)
    }
    valid <- tryCatch(
        grDevices::col2rgb(palette),
        error = function(e) NULL
    )
    if (is.null(valid)) {
        stop("`evidence_palette` contains an invalid color.",
             call. = FALSE)
    }
    ramp <- grDevices::colorRamp(palette, space = "Lab")
    lower_evidence <- -log10(limits[["upper"]])
    upper_evidence <- -log10(limits[["lower"]])

    function(p) {
        p <- pmin(limits[["upper"]], pmax(limits[["lower"]], p))
        fraction <- (-log10(p) - lower_evidence) /
            (upper_evidence - lower_evidence)
        rgb <- round(ramp(pmin(1, pmax(0, fraction))))
        grDevices::rgb(rgb[, 1L], rgb[, 2L], rgb[, 3L], maxColorValue = 255)
    }
}

.dasra_plot_build_spec <- function(x, features, selection, max_features,
                                   alpha, p_color_limits,
                                   evidence_palette, group_colors) {
    if (!inherits(x, "dasra")) {
        stop("`x` must be a `dasra` result.", call. = FALSE)
    }
    required_results <- c(
        "taxon", "p_omnibus", "p_adj_omnibus",
        "p_adj_structural_absence", "z_structural_absence",
        "p_adj_relative_abundance", "z_relative_abundance"
    )
    if (!all(required_results %in% names(x$results)) ||
        !identical(x$settings$component, "all")) {
        stop(
            "The dual-component profile requires `component = \"all\"`.",
            call. = FALSE
        )
    }
    if (is.null(x$plot_data) || is.null(x$plot_data$profiles)) {
        stop(
            paste(
                "This fit does not contain plotting summaries; rerun",
                "`dasra(..., component = \"all\", store_plot_data = TRUE)`."
            ),
            call. = FALSE
        )
    }
    if (!identical(x$plot_data$version, 1L)) {
        stop("The stored plotting-summary version is not supported.",
             call. = FALSE)
    }
    fitted_contrast <- unname(
        x$settings$contrast[c("reference", "comparison")]
    )
    if (!identical(unname(x$plot_data$contrast), fitted_contrast)) {
        stop("Stored plotting summaries do not match the fitted contrast.",
             call. = FALSE)
    }
    required_profiles <- c(
        "feature", "structural_reference", "structural_comparison",
        "structural_status", "abundance_reference",
        "abundance_comparison", "abundance_status"
    )
    profiles <- x$plot_data$profiles
    if (!is.data.frame(profiles) ||
        !all(required_profiles %in% names(profiles)) ||
        anyNA(profiles$feature) || anyDuplicated(profiles$feature)) {
        stop("The stored plotting summaries are malformed.",
             call. = FALSE)
    }
    selection <- match.arg(selection, c("significant", "top", "all"))
    max_features <- .dasra_plot_validate_integer(
        max_features, "max_features"
    )
    if (length(alpha) != 1L || is.na(alpha) || !is.numeric(alpha) ||
        !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
        stop("`alpha` must be one number strictly between zero and one.",
             call. = FALSE)
    }
    if (!is.character(group_colors) || length(group_colors) != 2L ||
        anyNA(group_colors) ||
        is.null(tryCatch(grDevices::col2rgb(group_colors),
                         error = function(e) NULL))) {
        stop("`group_colors` must contain two valid colors.",
             call. = FALSE)
    }

    indices <- .dasra_plot_select(
        x, features, selection, max_features, alpha
    )
    profile_index <- match(x$results$taxon[indices], profiles$feature)
    if (anyNA(profile_index)) {
        stop("Stored plotting summaries do not match `x$results`.",
             call. = FALSE)
    }

    results <- x$results[indices, , drop = FALSE]
    data <- cbind(
        data.frame(
            feature = results$taxon,
            feature_label = gsub("_", " ", results$taxon, fixed = TRUE),
            stars = vapply(
                results$p_adj_omnibus,
                .dasra_plot_stars,
                character(1)
            ),
            z_structural = results$z_structural_absence,
            p_structural = results$p_adj_structural_absence,
            z_abundance = results$z_relative_abundance,
            p_abundance = results$p_adj_relative_abundance,
            p_omnibus = results$p_omnibus,
            p_adjusted_omnibus = results$p_adj_omnibus,
            stringsAsFactors = FALSE,
            check.names = FALSE
        ),
        profiles[profile_index, setdiff(
            names(profiles), "feature"
        ), drop = FALSE]
    )
    limits <- .dasra_plot_color_limits(data, p_color_limits)
    color_function <- .dasra_plot_palette_function(
        evidence_palette, limits
    )
    data$color_structural <- rep(NA_character_, nrow(data))
    structural_color <- is.finite(data$z_structural) &
        is.finite(data$p_structural)
    data$color_structural[structural_color] <- color_function(
        data$p_structural[structural_color]
    )
    data$color_abundance <- rep(NA_character_, nrow(data))
    abundance_color <- is.finite(data$z_abundance) &
        is.finite(data$p_abundance)
    data$color_abundance[abundance_color] <- color_function(
        data$p_abundance[abundance_color]
    )

    list(
        data = data,
        contrast = fitted_contrast,
        adjustment_label = .dasra_plot_adjustment_label(
            x$settings$p_adjust_method
        ),
        color_limits = limits,
        evidence_palette = evidence_palette,
        group_colors = unname(group_colors),
        selection = selection,
        alpha = alpha
    )
}

# Grid drawing -------------------------------------------------------------

.dasra_plot_format_abundance_tick <- function(exponent) {
    if (exponent >= -2L && exponent <= 2L) {
        return(paste0(
            format(10 ^ exponent, scientific = FALSE, trim = TRUE),
            "%"
        ))
    }
    paste0("10^", exponent, "*'%'")
}

.dasra_plot_abundance_scale <- function(values, max_ticks = 4L) {
    values <- values[is.finite(values) & values > 0]
    limits <- if (length(values)) {
        c(floor(log10(min(values))), ceiling(log10(max(values))))
    } else {
        c(-6, 0)
    }
    if (limits[[1L]] == limits[[2L]]) {
        limits[[1L]] <- limits[[1L]] - 1
    }
    tick_step <- max(
        1L,
        ceiling(diff(limits) / (max_ticks - 1L))
    )
    ticks <- seq(limits[[1L]], limits[[2L]], by = tick_step)
    domain <- limits + c(-1, 1) * 0.035 * diff(limits)

    list(
        limits = limits,
        domain = domain,
        ticks = ticks,
        labels = vapply(
            ticks, .dasra_plot_format_abundance_tick, character(1)
        )
    )
}

.dasra_plot_axis <- function(at, labels, reversed = FALSE,
                             domain = range(at), gp_text, gp_line) {
    if (length(domain) != 2L || any(!is.finite(domain)) ||
        domain[[1L]] >= domain[[2L]]) {
        stop("Internal plotting axis has an invalid domain.", call. = FALSE)
    }
    plotted_at <- if (reversed) sum(domain) - at else at
    grid::pushViewport(grid::viewport(xscale = domain))
    grid::grid.segments(
        x0 = grid::unit(domain[[1L]], "native"),
        x1 = grid::unit(domain[[2L]], "native"),
        y0 = grid::unit(0.08, "npc"),
        y1 = grid::unit(0.08, "npc"),
        gp = gp_line
    )
    for (i in seq_along(at)) {
        grid::grid.segments(
            x0 = grid::unit(plotted_at[[i]], "native"),
            x1 = grid::unit(plotted_at[[i]], "native"),
            y0 = grid::unit(0.08, "npc"),
            y1 = grid::unit(0.23, "npc"),
            gp = gp_line
        )
        label <- labels[[i]]
        if (is.character(label) && grepl("^10\\^", label)) {
            label <- parse(text = label)
        }
        grid::grid.text(
            label,
            x = grid::unit(plotted_at[[i]], "native"),
            y = grid::unit(0.48, "npc"),
            just = if (abs(plotted_at[[i]] - domain[[1L]]) < 1e-12) {
                "left"
            } else if (abs(plotted_at[[i]] - domain[[2L]]) < 1e-12) {
                "right"
            } else {
                "centre"
            },
            gp = gp_text
        )
    }
    grid::popViewport()
}

.dasra_plot_header_groups <- function(title, contrast, colors, gp_title,
                                      gp_small) {
    label_widths <- vapply(contrast, function(label) {
        grid::convertWidth(
            grid::grobWidth(grid::textGrob(label, gp = gp_small)),
            "npc",
            valueOnly = TRUE
        )
    }, numeric(1))
    key_width <- 0.07
    key_gap <- 0.025
    entry_gap <- 0.06
    fixed_width <- length(contrast) * (key_width + key_gap) +
        (length(contrast) - 1L) * entry_gap
    group_gp <- gp_small
    available_label_width <- 0.94 - fixed_width
    if (is.finite(sum(label_widths)) &&
        sum(label_widths) > available_label_width) {
        scale <- available_label_width / sum(label_widths)
        group_gp <- grid::gpar(
            col = gp_small$col,
            fontsize = gp_small$fontsize * scale
        )
        label_widths <- label_widths * scale
    }
    grid::grid.text(
        title, x = grid::unit(0, "npc"), y = grid::unit(0.78, "npc"),
        just = "left", gp = gp_title
    )
    cursor <- 0
    for (i in seq_along(contrast)) {
        grid::grid.segments(
            x0 = grid::unit(cursor, "npc"),
            x1 = grid::unit(cursor + key_width, "npc"),
            y0 = grid::unit(0.20, "npc"),
            y1 = grid::unit(0.20, "npc"),
            gp = grid::gpar(
                col = colors[[i]], lwd = 2.6, lineend = "butt"
            )
        )
        label_x <- cursor + key_width + key_gap
        grid::grid.text(
            contrast[[i]],
            x = grid::unit(label_x, "npc"),
            y = grid::unit(0.20, "npc"),
            just = "left",
            gp = group_gp
        )
        cursor <- label_x + label_widths[[i]] + entry_gap
    }
}

.dasra_plot_draw_header <- function(spec, columns, gp_title, gp_small,
                                    gp_tiny) {
    data <- spec$data
    grid::pushViewport(grid::viewport(
        layout.pos.row = 1L, layout.pos.col = columns[[1L]]
    ))
    grid::grid.text(
        "Feature", x = 0, y = 0.78, just = "left", gp = gp_title
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(
        layout.pos.row = 1L, layout.pos.col = columns[[2L]]
    ))
    .dasra_plot_header_groups(
        "Structural-absence\nprobability", spec$contrast, spec$group_colors,
        gp_title, gp_small
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(
        layout.pos.row = 1L, layout.pos.col = columns[[3L]]
    ))
    grid::grid.text("Signed evidence", y = 0.82, gp = gp_title)
    grid::grid.points(
        x = grid::unit(c(0.24, 0.55), "npc"),
        y = grid::unit(c(0.54, 0.54), "npc"),
        pch = c(24, 21),
        size = grid::unit(c(1.55, 1.45), "mm"),
        gp = grid::gpar(col = "#1A1C1F", fill = "white", lwd = 0.7)
    )
    grid::grid.text(
        parse(text = c("Z[j]^SA", "Z[j]^RA")),
        x = grid::unit(c(0.28, 0.59), "npc"),
        y = grid::unit(c(0.54, 0.54), "npc"),
        just = "left",
        gp = gp_small
    )
    ramp_left <- 0.43
    ramp_right <- 0.78
    ramp_colors <- grDevices::colorRampPalette(
        spec$evidence_palette, space = "Lab"
    )(48L)
    ramp_x <- seq(
        ramp_left, ramp_right, length.out = length(ramp_colors) + 1L
    )
    for (i in seq_along(ramp_colors)) {
        grid::grid.rect(
            x = mean(ramp_x[c(i, i + 1L)]),
            y = 0.18,
            width = diff(ramp_x[c(i, i + 1L)]),
            height = 0.075,
            gp = grid::gpar(
                col = NA, fill = ramp_colors[[i]]
            )
        )
    }
    grid::grid.text(
        paste("Component", spec$adjustment_label),
        x = 0.39, y = 0.18, just = "right", gp = gp_tiny
    )
    legend_p <- c(
        spec$color_limits[["upper"]],
        if (spec$color_limits[["lower"]] < 0.05 &&
            spec$color_limits[["upper"]] > 0.05) 0.05 else numeric(),
        spec$color_limits[["lower"]]
    )
    legend_fraction <- (
        -log10(legend_p) + log10(spec$color_limits[["upper"]])
    ) / (
        -log10(spec$color_limits[["lower"]]) +
            log10(spec$color_limits[["upper"]])
    )
    legend_x <- ramp_left + legend_fraction * (ramp_right - ramp_left)
    legend_labels <- vapply(legend_p, function(p) {
        if (abs(p - 0.05) < .Machine$double.eps ^ 0.5) {
            ".05"
        } else {
            format.pval(p, digits = 1, eps = 0)
        }
    }, character(1))
    grid::grid.segments(
        x0 = grid::unit(legend_x, "npc"),
        x1 = grid::unit(legend_x, "npc"),
        y0 = grid::unit(0.125, "npc"),
        y1 = grid::unit(0.235, "npc"),
        gp = grid::gpar(col = "#6F7479", lwd = 0.35)
    )
    for (i in seq_along(legend_x)) {
        horizontal_justification <- if (i == 1L) {
            "left"
        } else if (i == length(legend_x)) {
            "right"
        } else {
            "centre"
        }
        grid::grid.text(
            legend_labels[[i]],
            x = grid::unit(legend_x[[i]], "npc"),
            y = grid::unit(0.035, "npc"),
            just = horizontal_justification,
            gp = gp_tiny
        )
    }
    grid::popViewport()

    grid::pushViewport(grid::viewport(
        layout.pos.row = 1L, layout.pos.col = columns[[4L]]
    ))
    .dasra_plot_header_groups(
        "Present-conditional\nrelative abundance", spec$contrast,
        spec$group_colors, gp_title, gp_small
    )
    grid::popViewport()
    invisible(data)
}

.dasra_plot_draw <- function(spec) {
    data <- spec$data
    n <- nrow(data)
    ink <- "#1A1C1F"
    muted <- "#73777C"
    rule <- "#D9DDE0"
    faint <- "#ECEFF1"
    gp_title <- grid::gpar(col = ink, fontsize = 8.3, fontface = "bold")
    gp_label <- grid::gpar(col = ink, fontsize = 7.2)
    gp_small <- grid::gpar(col = ink, fontsize = 6.6)
    gp_tiny <- grid::gpar(col = muted, fontsize = 6.6)
    gp_axis <- grid::gpar(col = muted, fontsize = 6.6)
    gp_rule <- grid::gpar(col = rule, lwd = 0.45)

    star_gp <- grid::gpar(col = ink, fontsize = 6.4, fontface = "bold")
    label_widths <- lapply(seq_len(n), function(i) {
        width <- grid::grobWidth(grid::textGrob(
            data$feature_label[[i]], gp = gp_label
        ))
        if (nzchar(data$stars[[i]])) {
            width <- width + grid::unit(1, "pt") +
                grid::grobWidth(grid::textGrob(
                    data$stars[[i]], gp = star_gp
                ))
        }
        width
    })
    feature_width <- grid::unit.pmax(
        do.call(grid::unit.pmax, label_widths),
        grid::grobWidth(grid::textGrob("Feature", gp = gp_title))
    ) + grid::unit(7, "pt")
    layout <- grid::grid.layout(
        nrow = n + 3L,
        ncol = 7L,
        widths = grid::unit.c(
            feature_width,
            grid::unit(2.7, "mm"),
            grid::unit(1.05, "null"),
            grid::unit(2.7, "mm"),
            grid::unit(1.72, "null"),
            grid::unit(2.7, "mm"),
            grid::unit(1.18, "null")
        ),
        heights = grid::unit.c(
            grid::unit(14, "mm"),
            grid::unit(7.5, "mm"),
            rep(grid::unit(1, "null"), n),
            grid::unit(3.8, "mm")
        )
    )
    columns <- c(1L, 3L, 5L, 7L)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
        width = grid::unit(1, "npc") - grid::unit(4, "mm"),
        height = grid::unit(1, "npc") - grid::unit(4, "mm"),
        layout = layout,
        name = "dasra-association-profile"
    ))
    .dasra_plot_draw_header(
        spec, columns, gp_title, gp_small, gp_tiny
    )

    structural_ticks <- seq(0, 1, by = 0.25)
    structural_domain <- c(-0.035, 1.035)
    z_values <- c(data$z_structural, data$z_abundance)
    z_values <- z_values[is.finite(z_values)]
    z_step <- max(1, ceiling(max(abs(z_values)) / 2))
    z_limit <- 2 * z_step
    z_ticks <- seq(-z_limit, z_limit, by = z_step)
    z_domain <- c(-1.035, 1.035) * z_limit
    abundance_scale <- .dasra_plot_abundance_scale(c(
        data$abundance_reference,
        data$abundance_comparison
    ))
    abundance_domain <- abundance_scale$domain

    grid::pushViewport(grid::viewport(
        layout.pos.row = 2L, layout.pos.col = columns[[2L]]
    ))
    .dasra_plot_axis(
        structural_ticks,
        paste0(round(100 * structural_ticks), "%"),
        reversed = TRUE,
        domain = structural_domain,
        gp_text = gp_axis,
        gp_line = gp_rule
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(
        layout.pos.row = 2L, layout.pos.col = columns[[3L]]
    ))
    .dasra_plot_axis(
        z_ticks,
        format(z_ticks, trim = TRUE, scientific = FALSE),
        domain = z_domain,
        gp_text = gp_axis,
        gp_line = gp_rule
    )
    grid::grid.text(
        paste(spec$contrast[[2L]], "lower"),
        x = 0, y = 0.88, just = "left", gp = gp_tiny
    )
    grid::grid.text(
        paste(spec$contrast[[2L]], "higher"),
        x = 1, y = 0.88, just = "right", gp = gp_tiny
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(
        layout.pos.row = 2L, layout.pos.col = columns[[4L]]
    ))
    .dasra_plot_axis(
        abundance_scale$ticks,
        abundance_scale$labels,
        domain = abundance_domain,
        gp_text = gp_axis,
        gp_line = gp_rule
    )
    grid::grid.text(
        "log scale", x = 1, y = 0.88, just = "right", gp = gp_tiny
    )
    grid::popViewport()

    for (i in seq_len(n)) {
        row <- i + 2L
        for (column in columns) {
            grid::pushViewport(grid::viewport(
                layout.pos.row = row, layout.pos.col = column
            ))
            grid::grid.segments(
                x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                gp = grid::gpar(col = faint, lwd = 0.35)
            )
            grid::popViewport()
        }

        grid::pushViewport(grid::viewport(
            layout.pos.row = row, layout.pos.col = columns[[1L]]
        ))
        label <- grid::textGrob(
            data$feature_label[[i]],
            x = 0, y = 0.50, just = "left", gp = gp_label
        )
        grid::grid.draw(label)
        if (nzchar(data$stars[[i]])) {
            grid::grid.text(
                data$stars[[i]],
                x = grid::grobWidth(label) + grid::unit(1, "pt"),
                y = 0.68,
                just = "left",
                gp = star_gp
            )
        }
        grid::popViewport()

        grid::pushViewport(grid::viewport(
            layout.pos.row = row,
            layout.pos.col = columns[[2L]],
            xscale = structural_domain
        ))
        structural_profile <- c(
            data$structural_reference[[i]],
            data$structural_comparison[[i]]
        )
        if (all(is.finite(structural_profile))) {
            y <- c(0.64, 0.36)
            plotted_profile <- 1 - structural_profile
            grid::grid.segments(
                x0 = grid::unit(1, "native"),
                x1 = grid::unit(plotted_profile, "native"),
                y0 = y, y1 = y,
                gp = grid::gpar(
                    col = spec$group_colors,
                    lwd = 2.8,
                    lineend = "round"
                )
            )
            grid::grid.points(
                x = grid::unit(plotted_profile, "native"),
                y = y,
                pch = 21,
                size = grid::unit(1.0, "mm"),
                gp = grid::gpar(
                    col = spec$group_colors,
                    fill = spec$group_colors,
                    lwd = 0.5
                )
            )
        } else {
            grid::grid.text("--", x = 0.98, just = "right", gp = gp_axis)
        }
        grid::popViewport()

        grid::pushViewport(grid::viewport(
            layout.pos.row = row,
            layout.pos.col = columns[[3L]],
            xscale = z_domain
        ))
        for (tick in z_ticks) {
            grid::grid.segments(
                x0 = grid::unit(tick, "native"),
                x1 = grid::unit(tick, "native"),
                y0 = 0, y1 = 1,
                gp = grid::gpar(
                    col = if (tick == 0) "#8A8D91" else faint,
                    lwd = if (tick == 0) 0.65 else 0.30,
                    lty = if (tick == 0) 2 else 1
                )
            )
        }
        if (is.finite(data$z_structural[[i]])) {
            grid::grid.segments(
                x0 = grid::unit(0, "native"),
                x1 = grid::unit(data$z_structural[[i]], "native"),
                y0 = 0.64, y1 = 0.64,
                gp = grid::gpar(col = "#B8BDC1", lwd = 0.48)
            )
            grid::grid.points(
                x = grid::unit(data$z_structural[[i]], "native"),
                y = 0.64,
                pch = 24,
                size = grid::unit(1.65, "mm"),
                gp = grid::gpar(
                    col = ink,
                    fill = data$color_structural[[i]],
                    lwd = 0.65
                )
            )
        }
        if (is.finite(data$z_abundance[[i]])) {
            grid::grid.segments(
                x0 = grid::unit(0, "native"),
                x1 = grid::unit(data$z_abundance[[i]], "native"),
                y0 = 0.36, y1 = 0.36,
                gp = grid::gpar(col = "#B8BDC1", lwd = 0.48)
            )
            grid::grid.points(
                x = grid::unit(data$z_abundance[[i]], "native"),
                y = 0.36,
                pch = 21,
                size = grid::unit(1.52, "mm"),
                gp = grid::gpar(
                    col = ink,
                    fill = data$color_abundance[[i]],
                    lwd = 0.65
                )
            )
        }
        grid::popViewport()

        grid::pushViewport(grid::viewport(
            layout.pos.row = row,
            layout.pos.col = columns[[4L]],
            xscale = abundance_domain
        ))
        abundance_profile <- c(
            data$abundance_reference[[i]],
            data$abundance_comparison[[i]]
        )
        if (all(is.finite(abundance_profile)) &&
            all(abundance_profile > 0)) {
            transformed <- log10(abundance_profile)
            y <- c(0.56, 0.44)
            grid::grid.segments(
                x0 = grid::unit(transformed[[1L]], "native"),
                x1 = grid::unit(transformed[[2L]], "native"),
                y0 = y[[1L]], y1 = y[[2L]],
                gp = grid::gpar(
                    col = "#B8BDC1",
                    lwd = 0.62,
                    lineend = "round"
                )
            )
            grid::grid.points(
                x = grid::unit(transformed, "native"),
                y = y,
                pch = 21,
                size = grid::unit(if (n > 12L) 1.08 else 1.25, "mm"),
                gp = grid::gpar(
                    col = "#50555A",
                    fill = spec$group_colors,
                    lwd = 0.45
                )
            )
        } else {
            grid::grid.text("--", x = 0.02, just = "left", gp = gp_axis)
        }
        grid::popViewport()
    }

    grid::pushViewport(grid::viewport(
        layout.pos.row = n + 3L,
        layout.pos.col = 1L:7L
    ))
    grid::grid.text(
        paste0(
            "Omnibus ", spec$adjustment_label,
            ":  * <= .05    ** <= .01    *** <= .001"
        ),
        x = 0, y = 0.55, just = "left", gp = gp_tiny
    )
    grid::popViewport()
    grid::popViewport()
    invisible(spec)
}

.dasra_plot_open_device <- function(file, width, height, dpi) {
    extension <- tolower(sub("^.*\\.", "", file))
    if (identical(extension, "pdf")) {
        grDevices::pdf(
            file, width = width, height = height,
            onefile = FALSE, useDingbats = FALSE, bg = "white"
        )
    } else if (identical(extension, "png")) {
        grDevices::png(
            file, width = width, height = height,
            units = "in", res = dpi, bg = "white"
        )
    } else if (identical(extension, "svg")) {
        if (!isTRUE(capabilities("cairo"))) {
            stop("SVG output requires an R build with Cairo support.",
                 call. = FALSE)
        }
        grDevices::svg(file, width = width, height = height, bg = "white")
    } else {
        stop("`file` must end in .pdf, .png, or .svg.", call. = FALSE)
    }
    invisible(NULL)
}

#' Plot a DASRA dual-component association profile
#'
#' Draws a publication-oriented, taxon-aligned summary of a DASRA analysis.
#' The center displays the structural-absence and relative-abundance signed
#' statistics on one standardized scale. Their colors encode the corresponding
#' component adjusted p-values. Feature superscripts encode the primary omnibus
#' adjusted p-value. The side panels show group-standardized model summaries
#' prepared when [dasra()] was called with `store_plot_data = TRUE`.
#'
#' The structural side panel is a descriptive unrestricted companion fit; the
#' prespecified structural inference remains the central score statistic. The
#' abundance side panel displays standardized present-conditional geometric
#' mean relative abundance, expressed as a percentage on a log scale, from the
#' fitted mark model. Its paired points avoid implying a zero baseline on the
#' log scale. The central abundance statistic retains the target-excluded
#' reference correction used by DASRA.
#' Side-panel segments are descriptive fitted group summaries, not effect
#' estimates or confidence intervals; the abundance connector only links its
#' two fitted group means. Component inference is carried by the central signed
#' statistics. An omitted central marker or `--` denotes an unavailable
#' component or side summary, not a zero effect.
#'
#' @param x A `dasra` object fitted with `component = "all"` and
#'   `store_plot_data = TRUE`.
#' @param features Optional character vector of feature names to display. The
#'   supplied order is preserved and overrides `selection`, `max_features`, and
#'   `alpha`.
#' @param selection Feature-selection rule when `features` is `NULL`.
#'   `"significant"` displays primary omnibus adjusted p-values no greater than
#'   `alpha`; `"top"` displays the smallest available adjusted p-values; and
#'   `"all"` displays every formed omnibus result.
#' @param max_features Maximum number of rows for `"significant"` or `"top"`.
#'   The default is `24L`. Significant results are never padded with
#'   nonsignificant features.
#' @param alpha Adjusted-p threshold used by `selection = "significant"`.
#' @param p_color_limits Either `"adaptive"` or numeric `c(lower, upper)` with
#'   `0 < lower < upper <= 1`. The adaptive strong-evidence endpoint is derived
#'   from the smallest positive finite component adjusted p-value among the
#'   displayed markers and rounded down to a power of ten; the weak-evidence
#'   endpoint is one. Degenerate all-zero or all-one displays use a finite
#'   fallback. Colors are interpolated on the `-log10(p)` scale and values
#'   outside the range are clamped.
#' @param evidence_palette Colors from weak to strong component evidence. The
#'   default is a sequential gray-plum scale chosen to remain distinct from
#'   the reference/comparison group colors.
#' @param group_colors Two colors for the reference and comparison groups.
#' @param group_labels Optional two-element character vector used only as the
#'   displayed reference and comparison labels. The fitted contrast is used by
#'   default. This does not change the analysis or group ordering.
#' @param file Optional output filename ending in `.pdf`, `.png`, or `.svg`.
#'   When `NULL`, the current graphics device is used.
#' @param width,height Output dimensions in inches when `file` is supplied.
#'   Defaults use a 180-mm manuscript width and a height chosen from the number
#'   of displayed features.
#' @param dpi Resolution for PNG output.
#' @param ... Reserved for future graphical options.
#'
#' @return Invisibly, the data and resolved graphical settings used to draw the
#'   figure.
#'
#' @examples
#' \dontrun{
#' fit <- dasra(
#'   counts, metadata, ~ group + age, "group", "library_size",
#'   component = "all", store_plot_data = TRUE
#' )
#' plot(fit)
#' plot(fit, p_color_limits = c(1e-8, 1))
#' plot(fit, file = "dasra-profile.pdf", width = 7.2, height = 6.5)
#' }
#'
#' @method plot dasra
#' @export
#' @importFrom graphics plot
plot.dasra <- function(
        x,
        features = NULL,
        selection = c("significant", "top", "all"),
        max_features = 24L,
        alpha = 0.05,
        p_color_limits = "adaptive",
        evidence_palette = c(
            "#e8ecef", "#d6cdd5", "#c0acba",
            "#a7879c", "#89617c", "#673b5c"
        ),
        group_colors = c("#6090c1", "#f28e4b"),
        file = NULL,
        width = NULL,
        height = NULL,
        dpi = 300,
        group_labels = NULL,
        ...) {
    selection <- match.arg(selection)
    spec <- .dasra_plot_build_spec(
        x = x,
        features = features,
        selection = selection,
        max_features = max_features,
        alpha = alpha,
        p_color_limits = p_color_limits,
        evidence_palette = evidence_palette,
        group_colors = group_colors
    )

    if (!is.null(group_labels)) {
        if (!is.character(group_labels) || length(group_labels) != 2L ||
            anyNA(group_labels) || any(!nzchar(group_labels))) {
            stop("`group_labels` must contain two non-empty labels.",
                 call. = FALSE)
        }
        spec$contrast <- unname(group_labels)
    }

    if (!is.null(file)) {
        if (!is.character(file) || length(file) != 1L || is.na(file) ||
            !nzchar(file)) {
            stop("`file` must be one non-empty filename.", call. = FALSE)
        }
        if (is.null(width)) width <- 180 / 25.4
        if (is.null(height)) {
            height <- max(2.8, 1.15 + 0.215 * nrow(spec$data))
        }
        if (!is.numeric(width) || length(width) != 1L ||
            !is.finite(width) || width <= 0 ||
            !is.numeric(height) || length(height) != 1L ||
            !is.finite(height) || height <= 0) {
            stop("`width` and `height` must be positive numbers.",
                 call. = FALSE)
        }
        if (!is.numeric(dpi) || length(dpi) != 1L || !is.finite(dpi) ||
            dpi <= 0) {
            stop("`dpi` must be one positive number.", call. = FALSE)
        }
        .dasra_plot_open_device(file, width, height, dpi)
        on.exit(grDevices::dev.off(), add = TRUE)
    }
    .dasra_plot_draw(spec)
    invisible(spec)
}
