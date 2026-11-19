# Plot-data preparation ----------------------------------------------------

.sacat_plot_profile_failure <- function(status) {
    list(
        reference = NA_real_,
        comparison = NA_real_,
        status = as.character(status)[1L]
    )
}

.sacat_structural_plot_profile <- function(y, N, fit) {
    if (is.null(fit) || !isTRUE(fit$tested) || !isTRUE(fit$regular)) {
        reason <- if (is.null(fit$reason)) "unavailable" else fit$reason
        return(.sacat_plot_profile_failure(reason))
    }
    stored <- fit$diagnostics$fit
    if (is.null(stored) || is.null(stored$conditional_present$par) ||
        is.null(stored$X_eta) || is.null(stored$quadrature)) {
        return(.sacat_plot_profile_failure("fit_details_unavailable"))
    }

    X_eta <- stored$X_eta
    if (!("Group" %in% colnames(X_eta))) {
        return(.sacat_plot_profile_failure("group_design_unavailable"))
    }

    companion_state <- tryCatch(
        suppressWarnings({
            conditional_present <- zt_beta_detection_components(
                stored$conditional_present$par,
                N,
                X_eta,
                stored$quadrature
            )
            companion <- if (is.null(conditional_present)) {
                list(
                    ok = FALSE,
                    reason = "structural_absence_no_finite_fit"
                )
            } else {
                zt_fit_alpha(
                    beta = stored$conditional_present$par,
                    y = y,
                    N = N,
                    X_rho = X_eta,
                    X_eta = X_eta,
                    gh = stored$quadrature,
                    conditional_present = conditional_present
                )
            }
            list(
                companion = companion,
                conditional_present = conditional_present
            )
        }),
        error = function(e) NULL
    )
    companion <- companion_state$companion
    if (is.null(companion) || !isTRUE(companion$ok)) {
        reason <- if (is.null(companion$reason)) {
            "companion_fit_unavailable"
        } else {
            companion$reason
        }
        return(.sacat_plot_profile_failure(reason))
    }

    components <- companion_state$conditional_present
    if (is.null(components)) {
        return(.sacat_plot_profile_failure(
            "companion_detection_components_unavailable"
        ))
    }
    finite_nll <- tryCatch(
        companion$nll(companion$par),
        error = function(e) NA_real_
    )
    if (!is.finite(finite_nll)) {
        return(.sacat_plot_profile_failure(
            "companion_objective_unavailable"
        ))
    }
    zero_limit_comparison <- zt_structural_zero_limit_comparison(
        finite_nll = finite_nll,
        log_r = components$log_r,
        y = y
    )
    if (isTRUE(zero_limit_comparison$dominated)) {
        return(.sacat_plot_profile_failure(
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
        return(.sacat_plot_profile_failure(
            "companion_standardization_unavailable"
        ))
    }

    state <- zt_detection_state_from_components(
        companion$par, y, X_eta, components
    )
    if (is.null(state) || any(!is.finite(state$gamma))) {
        return(.sacat_plot_profile_failure(
            "companion_posterior_unavailable"
        ))
    }

    list(
        reference = reference, comparison = comparison,
        posterior_sa = state$gamma, status = "ok"
    )
}

.sacat_abundance_plot_profile <- function(taxon, abundance_fits) {
    if (is.null(abundance_fits$taxon) ||
        is.null(abundance_fits$raw_fits)) {
        return(.sacat_plot_profile_failure("fit_details_unavailable"))
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
        return(.sacat_plot_profile_failure(reason))
    }

    theta <- raw_fit$theta
    X_eta <- raw_fit$solver_diagnostics$X_eta
    if (is.null(theta) || is.null(X_eta) ||
        !("Group" %in% colnames(X_eta))) {
        return(.sacat_plot_profile_failure("group_design_unavailable"))
    }
    p_eta <- ncol(X_eta)
    if (length(theta) != p_eta + 1L) {
        return(.sacat_plot_profile_failure("parameter_layout_unavailable"))
    }

    sigma <- exp(theta[[p_eta + 1L]])
    gh <- .sacat_make_abundance_gh_rule(raw_fit$quadrature_Q)
    X_reference <- X_comparison <- X_eta
    X_reference[, "Group"] <- 0
    X_comparison[, "Group"] <- 1
    location_reference <- as.numeric(
        X_reference %*% theta[seq_len(p_eta)]
    )
    location_comparison <- as.numeric(
        X_comparison %*% theta[seq_len(p_eta)]
    )
    reference <- 100 * exp(mean(.sacat_abundance_mean_log_relative(
        location_reference, sigma, gh
    )))
    comparison <- 100 * exp(mean(.sacat_abundance_mean_log_relative(
        location_comparison, sigma, gh
    )))
    if (!is.finite(reference) || !is.finite(comparison) ||
        reference <= 0 || comparison <= 0) {
        return(.sacat_plot_profile_failure(
            "abundance_standardization_unavailable"
        ))
    }

    list(reference = reference, comparison = comparison, status = "ok")
}

.sacat_prepare_plot_data <- function(Y, N, group, results, structural_fits,
                                     abundance_fits, contrast,
                                     cluster = NULL) {
    taxa <- results$taxon
    structural_one <- function(j) {
        tryCatch(
            .sacat_structural_plot_profile(
                y = Y[, j],
                N = N,
                fit = structural_fits[[taxa[[j]]]]
            ),
            error = function(e) .sacat_plot_profile_failure(
                paste0("plot_profile_error: ", conditionMessage(e))
            )
        )
    }
    structural <- .sacat_taxon_lapply(
        seq_along(taxa), structural_one, cluster = cluster
    )
    abundance <- lapply(taxa, function(taxon) {
        tryCatch(
            .sacat_abundance_plot_profile(taxon, abundance_fits),
            error = function(e) .sacat_plot_profile_failure(
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

    zero_sources <- do.call(rbind, lapply(seq_along(taxa), function(j) {
        posterior <- structural[[j]]$posterior_sa
        if (is.null(posterior)) posterior <- rep(NA_real_, nrow(Y))
        posterior[Y[, j] > 0] <- 0
        do.call(rbind, lapply(0:1, function(g) {
            selected <- group == g
            available <- all(is.finite(posterior[selected]))
            data.frame(
                feature = taxa[[j]], group = g,
                detected = mean(Y[selected, j] > 0),
                estimated_sa = if (available) mean(posterior[selected]) else NA_real_,
                estimated_nd = if (available) {
                    mean((Y[selected, j] == 0) - posterior[selected])
                } else NA_real_,
                n = sum(selected),
                status = if (available) "ok" else structural[[j]]$status,
                stringsAsFactors = FALSE
            )
        }))
    }))
    positive_abundance <- do.call(rbind, lapply(seq_along(taxa), function(j) {
        selected <- Y[, j] > 0
        data.frame(
            feature = rep(taxa[[j]], sum(selected)),
            group = group[selected],
            abundance = 100 * Y[selected, j] / N[selected],
            stringsAsFactors = FALSE
        )
    }))

    list(
        version = 2L,
        contrast = contrast,
        profiles = profiles,
        zero_sources = zero_sources,
        positive_abundance = positive_abundance
    )
}

# Plot specification -------------------------------------------------------

.sacat_plot_stars <- function(p) {
    if (!is.finite(p) || p > 0.05) return("")
    if (p <= 0.001) return("***")
    if (p <= 0.01) return("**")
    "*"
}

.sacat_plot_adjustment_label <- function(method) {
    if (identical(method, "none")) "p" else paste(method, "adjusted p")
}

.sacat_plot_component_bh_guide <- function(
        retained, raw_p, adjusted_p, z, alpha) {
    n_taxa <- length(retained)
    if (!is.numeric(raw_p) || !is.numeric(adjusted_p) || !is.numeric(z) ||
        length(raw_p) != n_taxa || length(adjusted_p) != n_taxa ||
        length(z) != n_taxa) {
        stop("the component p-values and statistics are malformed",
             call. = FALSE)
    }

    family <- retained & is.finite(raw_p)
    family_size <- sum(family)
    guide <- list(
        z = NA_real_, p = NA_real_, discoveries = 0L,
        family_size = family_size
    )
    if (!family_size) {
        return(guide)
    }

    raw_family <- raw_p[family]
    adjusted_family <- adjusted_p[family]
    z_family <- z[family]
    recomputed <- stats::p.adjust(raw_family, method = "BH")
    tolerance <- 1e-10
    finite_z <- is.finite(z_family)
    expected_p <- 2 * stats::pnorm(-abs(z_family[finite_z]))
    if (any(raw_family < 0 | raw_family > 1)) {
        stop("raw p-values fall outside [0, 1]", call. = FALSE)
    }
    if (any(!is.finite(adjusted_family)) ||
        max(abs(adjusted_family - recomputed)) > tolerance) {
        stop("stored adjusted p-values do not match BH adjustment",
             call. = FALSE)
    }
    if (any(abs(raw_family[finite_z] - expected_p) > tolerance)) {
        stop("raw p-values do not match the plotted statistics",
             call. = FALSE)
    }

    rejected <- recomputed <= alpha
    guide$discoveries <- sum(rejected)
    if (!guide$discoveries) {
        return(guide)
    }
    critical_p <- alpha * guide$discoveries / family_size
    critical_z <- stats::qnorm(
        critical_p / 2, lower.tail = FALSE
    )
    if (any(rejected & !finite_z)) {
        stop("a BH discovery has no finite statistic", call. = FALSE)
    }
    if (!is.finite(critical_z) || critical_z <= 0) {
        stop("the BH boundary has no finite statistic", call. = FALSE)
    }

    boundary_membership <- raw_family <= critical_p
    boundary_mismatch <- boundary_membership != rejected
    boundary_tolerance <- 16 * .Machine$double.eps * pmax(
        abs(raw_family), abs(critical_p), .Machine$double.xmin
    )
    if (any(boundary_mismatch &
            abs(raw_family - critical_p) > boundary_tolerance)) {
        stop("the BH decisions do not match the step-up boundary",
             call. = FALSE)
    }

    guide$z <- critical_z
    guide$p <- critical_p
    guide
}

.sacat_plot_component_guides <- function(x, alpha, show) {
    if (length(show) != 1L || is.na(show) || !is.logical(show)) {
        stop("`show_component_guides` must be TRUE or FALSE.",
             call. = FALSE)
    }
    method <- x$settings$p_adjust_method
    supported <- length(method) == 1L && !is.na(method) &&
        method %in% c("BH", "fdr")
    empty <- list(
        z = NA_real_, p = NA_real_, discoveries = NA_integer_,
        family_size = NA_integer_
    )
    if (!show || !supported) {
        return(list(
            enabled = FALSE, method = method, alpha = alpha,
            structural = empty, abundance = empty
        ))
    }

    retained <- x$diagnostics$retained
    if (!is.logical(retained) || length(retained) != nrow(x$results) ||
        anyNA(retained)) {
        stop("The retained-taxon indicator is malformed.", call. = FALSE)
    }
    check_component <- function(raw_p, adjusted_p, z) {
        tryCatch(
            list(
                guide = .sacat_plot_component_bh_guide(
                    retained, raw_p, adjusted_p, z, alpha
                ),
                problem = NULL
            ),
            error = function(e) list(guide = empty,
                                     problem = conditionMessage(e))
        )
    }
    structural <- check_component(
        x$results$p_structural_absence,
        x$results$p_adj_structural_absence,
        x$results$z_structural_absence
    )
    abundance <- check_component(
        x$results$p_relative_abundance,
        x$results$p_adj_relative_abundance,
        x$results$z_relative_abundance
    )
    problems <- c(
        if (!is.null(structural$problem)) {
            paste0("structural absence (", structural$problem, ")")
        },
        if (!is.null(abundance$problem)) {
            paste0("present-conditional abundance (",
                   abundance$problem, ")")
        }
    )
    if (length(problems)) {
        warning(
            paste0(
                "Omitted component BH gate",
                if (length(problems) > 1L) "s" else "",
                ": ", paste(problems, collapse = "; "), "."
            ),
            call. = FALSE
        )
    }
    list(
        enabled = TRUE, method = method, alpha = alpha,
        structural = structural$guide, abundance = abundance$guide
    )
}

.sacat_plot_validate_integer <- function(value, name) {
    if (length(value) != 1L || is.na(value) || !is.numeric(value) ||
        !is.finite(value) || value < 1 || value != as.integer(value)) {
        stop(sprintf("`%s` must be one positive integer.", name),
             call. = FALSE)
    }
    as.integer(value)
}

.sacat_plot_select <- function(x, features, selection, max_features,
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

.sacat_plot_build_spec <- function(x, features, selection, max_features,
                                   alpha, group_colors, zero_colors,
                                   show_component_guides = TRUE) {
    if (!inherits(x, "sacat")) {
        stop("`x` must be a `sacat` result.", call. = FALSE)
    }
    required_results <- c(
        "taxon", "p_omnibus", "p_adj_omnibus",
        "p_structural_absence", "p_adj_structural_absence",
        "z_structural_absence", "p_relative_abundance",
        "p_adj_relative_abundance", "z_relative_abundance"
    )
    if (!all(required_results %in% names(x$results)) ||
        !identical(x$settings$component, "all")) {
        stop('The profile requires `component = "all"`.', call. = FALSE)
    }
    payload <- x$plot_data
    if (is.null(payload) || !identical(payload$version, 2L)) {
        stop(paste(
            "This fit does not contain the required plotting data; refit with",
            "`store_plot_data = TRUE`."
        ), call. = FALSE)
    }
    fitted_contrast <- unname(
        x$settings$contrast[c("reference", "comparison")]
    )
    if (!identical(unname(payload$contrast), fitted_contrast)) {
        stop("Stored plotting summaries do not match the fitted contrast.",
             call. = FALSE)
    }
    required_profiles <- c(
        "feature", "abundance_reference", "abundance_comparison",
        "abundance_status"
    )
    profiles <- payload$profiles
    zeros <- payload$zero_sources
    positive <- payload$positive_abundance
    if (!is.data.frame(profiles) ||
        !all(required_profiles %in% names(profiles)) ||
        anyNA(profiles$feature) || anyDuplicated(profiles$feature) ||
        !is.data.frame(zeros) ||
        !all(c("feature", "group", "detected", "estimated_sa",
               "estimated_nd", "status") %in% names(zeros)) ||
        !is.data.frame(positive) ||
        !all(c("feature", "group", "abundance") %in% names(positive))) {
        stop("The stored plotting summaries are malformed.", call. = FALSE)
    }
    selection <- match.arg(selection, c("significant", "top", "all"))
    max_features <- .sacat_plot_validate_integer(max_features, "max_features")
    if (length(alpha) != 1L || is.na(alpha) || !is.numeric(alpha) ||
        !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
        stop("`alpha` must be one number strictly between zero and one.",
             call. = FALSE)
    }
    check_colors <- function(value, size, name) {
        if (!is.character(value) || length(value) != size || anyNA(value) ||
            is.null(tryCatch(grDevices::col2rgb(value),
                            error = function(e) NULL))) {
            stop(sprintf("`%s` must contain %d valid colors.", name, size),
                 call. = FALSE)
        }
        unname(value)
    }
    group_colors <- check_colors(group_colors, 2L, "group_colors")
    zero_colors <- check_colors(zero_colors, 3L, "zero_colors")

    indices <- .sacat_plot_select(x, features, selection, max_features, alpha)
    results <- x$results[indices, , drop = FALSE]
    profile_index <- match(results$taxon, profiles$feature)
    if (anyNA(profile_index)) {
        stop("Stored plotting summaries do not match `x$results`.",
             call. = FALSE)
    }
    zero_index <- match(
        paste(rep(results$taxon, each = 2L), rep(0:1, nrow(results))),
        paste(zeros$feature, zeros$group)
    )
    if (anyNA(zero_index)) {
        stop("Stored zero-source summaries do not match `x$results`.",
             call. = FALSE)
    }
    data <- cbind(
        data.frame(
            feature = results$taxon,
            feature_label = gsub("_", " ", results$taxon, fixed = TRUE),
            stars = vapply(results$p_adj_omnibus, .sacat_plot_stars,
                           character(1)),
            z_structural = results$z_structural_absence,
            p_structural = results$p_adj_structural_absence,
            z_abundance = results$z_relative_abundance,
            p_abundance = results$p_adj_relative_abundance,
            p_omnibus = results$p_omnibus,
            p_adjusted_omnibus = results$p_adj_omnibus,
            stringsAsFactors = FALSE
        ),
        profiles[profile_index, setdiff(names(profiles), "feature"),
                 drop = FALSE]
    )
    list(
        data = data,
        zero_sources = zeros[zero_index, , drop = FALSE],
        positive_abundance = positive[
            positive$feature %in% results$taxon, , drop = FALSE
        ],
        contrast = fitted_contrast,
        adjustment_label = .sacat_plot_adjustment_label(
            x$settings$p_adjust_method
        ),
        group_colors = group_colors,
        zero_colors = zero_colors,
        selection = selection,
        alpha = alpha,
        component_guides = .sacat_plot_component_guides(
            x, alpha, show_component_guides
        )
    )
}

# Grid drawing -------------------------------------------------------------

.sacat_plot_abundance_scale <- function(values, max_ticks = 5L) {
    values <- values[is.finite(values) & values > 0]
    limits <- if (length(values)) {
        c(floor(log10(min(values))), ceiling(log10(max(values))))
    } else {
        c(-3, 2)
    }
    if (limits[[1L]] == limits[[2L]]) limits[[1L]] <- limits[[1L]] - 1
    tick_step <- max(1L, ceiling(diff(limits) / (max_ticks - 1L)))
    ticks <- seq(limits[[1L]], limits[[2L]], by = tick_step)
    list(
        limits = limits,
        domain = limits + c(-1, 1) * 0.04 * diff(limits),
        ticks = ticks,
        labels = vapply(ticks, function(z) {
            format(10 ^ z, scientific = z < -3, trim = TRUE)
        }, character(1))
    )
}

.sacat_plot_box_summary <- function(values) {
    values <- log10(values[is.finite(values) & values > 0])
    if (!length(values)) return(rep(NA_real_, 5L))
    quartiles <- stats::quantile(values, c(0.25, 0.5, 0.75), names = FALSE)
    iqr <- quartiles[[3L]] - quartiles[[1L]]
    inside <- values[values >= quartiles[[1L]] - 1.5 * iqr &
                         values <= quartiles[[3L]] + 1.5 * iqr]
    c(min(inside), quartiles, max(inside))
}

.sacat_plot_axis <- function(at, labels, domain, title) {
    grid::pushViewport(grid::viewport(xscale = domain))
    grid::grid.segments(
        x0 = grid::unit(min(at), "native"),
        x1 = grid::unit(max(at), "native"), y0 = 0.91, y1 = 0.91,
        gp = grid::gpar(col = "#777777", lwd = 0.6)
    )
    grid::grid.segments(
        x0 = grid::unit(at, "native"), x1 = grid::unit(at, "native"),
        y0 = 0.91, y1 = 0.82, gp = grid::gpar(col = "#777777", lwd = 0.6)
    )
    grid::grid.text(
        labels, x = grid::unit(at, "native"), y = 0.64,
        gp = grid::gpar(col = "#333333", fontsize = 8)
    )
    if (nzchar(title)) {
        grid::grid.text(title, y = 0.17,
                        gp = grid::gpar(fontsize = 8.5, lineheight = 1))
    }
    grid::popViewport()
}

.sacat_plot_draw_box <- function(values, center, color, shape) {
    if (!length(values)) return(invisible(NULL))
    x <- log10(values)
    jitter <- 0.045 * sin(seq_along(x) * pi * (3 - sqrt(5)))
    grid::grid.points(
        x = grid::unit(x, "native"), y = center + jitter, pch = shape,
        size = grid::unit(0.7, "mm"),
        gp = grid::gpar(col = grDevices::adjustcolor(color, alpha.f = 0.28))
    )
    box <- .sacat_plot_box_summary(values)
    grid::grid.segments(
        x0 = grid::unit(box[[1L]], "native"),
        x1 = grid::unit(box[[5L]], "native"),
        y0 = center, y1 = center, gp = grid::gpar(col = color, lwd = 0.8)
    )
    grid::grid.segments(
        x0 = grid::unit(box[c(1L, 5L)], "native"),
        x1 = grid::unit(box[c(1L, 5L)], "native"),
        y0 = center - 0.05, y1 = center + 0.05,
        gp = grid::gpar(col = color, lwd = 0.8)
    )
    grid::grid.rect(
        x = grid::unit(mean(box[c(2L, 4L)]), "native"), y = center,
        width = grid::unit(box[[4L]] - box[[2L]], "native"), height = 0.19,
        gp = grid::gpar(col = color, lwd = 1,
                       fill = grDevices::adjustcolor(color, alpha.f = 0.18))
    )
    grid::grid.segments(
        x0 = grid::unit(box[[3L]], "native"),
        x1 = grid::unit(box[[3L]], "native"),
        y0 = center - 0.095, y1 = center + 0.095,
        gp = grid::gpar(col = color, lwd = 1.3)
    )
    invisible(box)
}

.sacat_plot_draw <- function(spec) {
    data <- spec$data
    n <- nrow(data)
    ink <- "#252525"
    muted <- "#6F6F6F"
    gp_title <- grid::gpar(col = ink, fontsize = 9, fontface = "bold",
                          lineheight = 1.05)
    gp_label <- grid::gpar(col = ink, fontsize = 8.5, fontface = "italic",
                          lineheight = 1)
    gp_small <- grid::gpar(col = ink, fontsize = 8, lineheight = 1)
    labels <- lapply(seq_len(n), function(i) {
        grid::textGrob(
            paste(strwrap(data$feature_label[[i]], width = 20), collapse = "\n"),
            rot = 45, gp = gp_label,
            name = paste0("sacat-feature-label-", i)
        )
    })
    group_labels <- vapply(spec$contrast, function(label) {
        paste(strwrap(label, width = 15), collapse = "\n")
    }, character(1))
    feature_width <- do.call(grid::unit.pmax, lapply(labels, grid::grobWidth)) +
        grid::unit(10, "pt")
    group_width <- do.call(grid::unit.pmax, lapply(group_labels, function(label) {
        grid::grobWidth(grid::textGrob(label, gp = gp_small))
    })) + grid::unit(4, "pt")
    layout <- grid::grid.layout(
        nrow = n + 5L, ncol = 10L,
        widths = grid::unit.c(
            feature_width, grid::unit(2, "mm"),
            grid::unit(0.74, "null"), grid::unit(2.5, "mm"),
            grid::unit(0.74, "null"), grid::unit(5, "mm"),
            group_width, grid::unit(1.3, "null"), grid::unit(5, "mm"),
            grid::unit(1.8, "null")
        ),
        heights = grid::unit.c(
            grid::unit(11, "mm"), grid::unit(5, "mm"),
            rep(grid::unit(1, "null"), n),
            grid::unit(16, "mm"), grid::unit(19, "mm"),
            grid::unit(5, "mm")
        )
    )
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
        width = grid::unit(1, "npc") - grid::unit(6, "mm"),
        height = grid::unit(1, "npc") - grid::unit(5, "mm"),
        layout = layout, gp = grid::gpar(fontfamily = "sans"),
        name = "sacat-association-profile"
    ))
    in_cell <- function(row, column, draw, ...) {
        grid::pushViewport(grid::viewport(
            layout.pos.row = row, layout.pos.col = column, ...
        ))
        draw()
        grid::popViewport()
    }
    headings <- list(
        list(cols = 3:5, tag = "A", text = "Component evidence"),
        list(cols = 7:8, tag = "B",
             text = "Observed detection and\nestimated zero sources"),
        list(cols = 10L, tag = "C", text = "Observed and fitted\nabundance")
    )
    for (heading in headings) {
        in_cell(1L, heading$cols, function() {
            grid::grid.text(heading$tag, x = 0, y = 0.9, just = c("left", "top"),
                            gp = grid::gpar(fontsize = 10, fontface = "bold"))
            grid::grid.text(heading$text, y = 0.45, gp = gp_title)
        })
    }
    in_cell(2L, 3L, function() grid::grid.text("SA", gp = gp_title))
    in_cell(2L, 5L, function() grid::grid.text("RA", gp = gp_title))

    z <- c(data$z_structural, data$z_abundance)
    if (isTRUE(spec$component_guides$enabled)) {
        z <- c(z, spec$component_guides$structural$z,
               spec$component_guides$abundance$z)
    }
    z <- z[is.finite(z)]
    step <- max(1, ceiling(max(abs(z), 1) / 2))
    z_ticks <- seq(-2 * step, 2 * step, by = step)
    z_domain <- range(z_ticks) * 1.07
    abundance <- .sacat_plot_abundance_scale(c(
        spec$positive_abundance$abundance,
        data$abundance_reference, data$abundance_comparison
    ))
    group_y <- c(0.68, 0.32)

    for (i in seq_len(n)) {
        row <- i + 2L
        for (column in list(1L, 3:5, 7:8, 10L)) {
            in_cell(row, column, function() {
                grid::grid.segments(
                    x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                    gp = grid::gpar(col = "#EAEAEA", lwd = 0.4)
                )
            })
        }
        in_cell(row, 1L, function() {
            grid::grid.draw(grid::editGrob(
                labels[[i]],
                x = grid::unit(1, "npc") - grid::unit(5, "pt") -
                    grid::grobWidth(labels[[i]]) / 2,
                y = grid::unit(0.5, "npc")
            ))
            if (nzchar(data$stars[[i]])) {
                grid::grid.text(
                    data$stars[[i]], x = 1, y = 0.78, just = "left",
                    gp = grid::gpar(fontsize = 7, fontface = "bold")
                )
            }
        })
        for (k in 1:2) {
            statistic <- if (k == 1L) data$z_structural[[i]] else
                data$z_abundance[[i]]
            adjusted_p <- if (k == 1L) data$p_structural[[i]] else
                data$p_abundance[[i]]
            guide <- if (k == 1L) spec$component_guides$structural$z else
                spec$component_guides$abundance$z
            in_cell(row, c(3L, 5L)[[k]], function() {
                grid::grid.segments(
                    x0 = grid::unit(0, "native"), x1 = grid::unit(0, "native"),
                    y0 = 0, y1 = 1, gp = grid::gpar(col = "#B6B6B6", lwd = 0.6)
                )
                if (isTRUE(spec$component_guides$enabled) && is.finite(guide)) {
                    grid::grid.segments(
                        x0 = grid::unit(c(-guide, guide), "native"),
                        x1 = grid::unit(c(-guide, guide), "native"),
                        y0 = 0, y1 = 1,
                        gp = grid::gpar(col = "#888888", lty = 2, lwd = 0.7)
                    )
                }
                if (is.finite(statistic)) {
                    grid::grid.segments(
                        x0 = grid::unit(0, "native"),
                        x1 = grid::unit(statistic, "native"),
                        y0 = 0.5, y1 = 0.5,
                        gp = grid::gpar(col = "#B5B5B5", lwd = 0.7)
                    )
                    grid::grid.points(
                        x = grid::unit(statistic, "native"), y = 0.5,
                        pch = 21, size = grid::unit(2.2, "mm"),
                        gp = grid::gpar(
                            col = ink, lwd = 0.9,
                            fill = if (is.finite(adjusted_p) &&
                                       adjusted_p <= spec$alpha) ink else "white"
                        )
                    )
                } else {
                    grid::grid.text("--", gp = gp_small)
                }
            }, xscale = z_domain)
        }

        in_cell(row, 7L, function() {
            grid::grid.text(group_labels, x = 0.95, y = group_y,
                            just = "right", gp = gp_small)
        })
        zero_rows <- spec$zero_sources[(2L * i - 1L):(2L * i), , drop = FALSE]
        in_cell(row, 8L, function() {
            for (g in 1:2) {
                fractions <- unlist(zero_rows[g, c(
                    "detected", "estimated_sa", "estimated_nd"
                )], use.names = FALSE)
                if (all(is.finite(fractions))) {
                    ends <- c(0, cumsum(fractions))
                    for (j in 1:3) {
                        grid::grid.rect(
                            x = mean(ends[j:(j + 1L)]), y = group_y[[g]],
                            width = fractions[[j]], height = 0.22,
                            gp = grid::gpar(fill = spec$zero_colors[[j]],
                                           col = "white", lwd = 0.3)
                        )
                    }
                } else {
                    detected <- fractions[[1L]]
                    grid::grid.rect(
                        x = detected / 2, y = group_y[[g]],
                        width = detected, height = 0.22,
                        gp = grid::gpar(fill = spec$zero_colors[[1L]],
                                       col = NA)
                    )
                    if (detected < 1) {
                        grid::grid.rect(
                            x = (1 + detected) / 2, y = group_y[[g]],
                            width = 1 - detected, height = 0.22,
                            gp = grid::gpar(fill = "white", col = muted, lty = 2)
                        )
                        grid::grid.text("--", x = (1 + detected) / 2,
                                        y = group_y[[g]], gp = gp_small)
                    }
                }
                grid::grid.rect(
                    x = 0.5, y = group_y[[g]], width = 1, height = 0.22,
                    gp = grid::gpar(fill = NA, col = "#A5A5A5", lwd = 0.4)
                )
            }
        })
        in_cell(row, 10L, function() {
            for (g in 0:1) {
                values <- spec$positive_abundance$abundance[
                    spec$positive_abundance$feature == data$feature[[i]] &
                    spec$positive_abundance$group == g
                ]
                .sacat_plot_draw_box(
                    values, group_y[[g + 1L]], spec$group_colors[[g + 1L]],
                    c(16, 17)[[g + 1L]]
                )
                fitted <- c(data$abundance_reference[[i]],
                            data$abundance_comparison[[i]])[[g + 1L]]
                if (is.finite(fitted) && fitted > 0) {
                    grid::grid.points(
                        x = grid::unit(log10(fitted), "native"),
                        y = group_y[[g + 1L]], pch = 23,
                        size = grid::unit(2, "mm"),
                        gp = grid::gpar(fill = "white", col = ink, lwd = 0.8)
                    )
                }
            }
        }, xscale = abundance$domain)
    }
    axis_row <- n + 3L
    for (column in c(3L, 5L)) {
        in_cell(axis_row, column, function() {
            .sacat_plot_axis(z_ticks, format(z_ticks, trim = TRUE), z_domain, "")
        })
    }
    in_cell(axis_row, 3:5, function() {
        grid::grid.text("Signed component statistic", y = 0.17,
                        gp = grid::gpar(fontsize = 8.5))
    })
    in_cell(axis_row, 8L, function() {
        .sacat_plot_axis(seq(0, 1, by = 0.25), c("0", "25", "50", "75", "100"),
                         c(0, 1), "Share of samples (%)")
    })
    in_cell(axis_row, 10L, function() {
        .sacat_plot_axis(abundance$ticks, abundance$labels, abundance$domain,
                         "Relative abundance (%)\n(log scale)")
    })

    legend_row <- n + 4L
    in_cell(legend_row, 3:5, function() {
        for (j in 1:2) {
            y <- c(0.72, 0.37)[[j]]
            grid::grid.points(
                x = 0.05, y = y, pch = 21, size = grid::unit(1.9, "mm"),
                gp = grid::gpar(fill = c(ink, "white")[[j]], col = ink)
            )
            grid::grid.text(
                paste(spec$adjustment_label, c("<=", ">")[[j]],
                      format(spec$alpha, trim = TRUE)),
                x = 0.14, y = y, just = "left", gp = gp_small
            )
        }
    })
    in_cell(legend_row, 7:8, function() {
        x <- c(0.16, 0.5, 0.84)
        grid::grid.rect(x = x, y = 0.73, width = grid::unit(2.1, "mm"),
                         height = grid::unit(2.1, "mm"),
                         gp = grid::gpar(fill = spec$zero_colors, col = NA))
        grid::grid.text(c("Detected", "Estimated\nSA", "Estimated\nND"),
                        x = x, y = 0.35, gp = gp_small)
    })
    in_cell(legend_row, 10L, function() {
        group_text <- lapply(group_labels, function(label) {
            grid::textGrob(label, gp = gp_small)
        })
        item_widths <- lapply(group_text, function(label) {
            grid::grobWidth(label) + grid::unit(4, "mm")
        })
        gap <- grid::unit(6, "mm")
        left <- grid::unit(0.5, "npc") -
            (item_widths[[1L]] + item_widths[[2L]] + gap) / 2
        for (g in 1:2) {
            grid::grid.points(
                x = left + grid::unit(1, "mm"), y = 0.82,
                pch = c(16, 17)[[g]],
                size = grid::unit(1.4, "mm"),
                gp = grid::gpar(col = spec$group_colors[[g]]),
                name = paste0("sacat-abundance-group-key-", g)
            )
            grid::grid.text(
                group_labels[[g]], x = left + grid::unit(4, "mm"), y = 0.82,
                just = "left", gp = gp_small,
                name = paste0("sacat-abundance-group-label-", g)
            )
            left <- left + item_widths[[g]] + gap
        }
        fitted_label <- grid::textGrob(
            "Fitted abundance\nConditional on presence", gp = gp_small
        )
        label_width <- grid::grobWidth(fitted_label)
        left <- grid::unit(0.5, "npc") -
            (label_width + grid::unit(4, "mm")) / 2
        grid::grid.points(
            x = left + grid::unit(1, "mm"), y = 0.30, pch = 23,
            size = grid::unit(2, "mm"),
            gp = grid::gpar(fill = "white", col = ink)
        )
        grid::grid.text(
            fitted_label$label,
            x = left + grid::unit(4, "mm") + label_width / 2,
            y = 0.30, gp = gp_small
        )
    })
    in_cell(n + 5L, 1:10, function() {
        grid::grid.text(
            paste0("Omnibus ", spec$adjustment_label,
                   ":  * <= 0.05   ** <= 0.01   *** <= 0.001"),
            x = 0, just = "left", gp = grid::gpar(col = muted, fontsize = 7)
        )
    })
    grid::popViewport()
    invisible(spec)
}

.sacat_plot_open_device <- function(file, width, height, dpi) {
    extension <- tolower(sub("^.*\\.", "", file))
    if (identical(extension, "pdf")) {
        grDevices::pdf(file, width = width, height = height,
                       onefile = FALSE, useDingbats = FALSE, bg = "white")
    } else if (identical(extension, "png")) {
        grDevices::png(file, width = width, height = height,
                       units = "in", res = dpi, bg = "white")
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

#' Plot a SACAT association profile
#'
#' Aligns component evidence, observed detection with estimated zero sources,
#' and observed and fitted abundance for selected taxa.
#'
#' Panel A shows the structural-absence (SA) score statistic and the signed
#' square root of the adjusted relative-abundance (RA) likelihood-ratio
#' statistic. Positive values indicate higher SA or reference-corrected RA
#' in the comparison group. Filled circles mark component adjusted p-values
#' at or below `alpha`; open circles mark values above it. Dashed lines show
#' component-specific BH discovery boundaries computed from the complete
#' retained families. Taxon superscripts summarize omnibus adjusted p-values.
#'
#' Panel B shows each group's observed detection fraction and estimated
#' contributions of structural absence and nondetection (ND) to observed zeros.
#' For each zero, the SA contribution is its posterior SA probability from the
#' group-inclusive structural plotting fit; the ND contribution is its
#' complement. Positive counts contribute to detected samples. These quantities
#' are averaged within the observed groups and sum to one. An undivided,
#' unfilled zero segment marked `--` denotes an unavailable zero-source fit.
#'
#' Panel C shows positive observed relative abundances as points and boxplots,
#' with library sizes from the fitted analysis. Boxplots summarize log10
#' relative abundance: the box spans the quartiles, the line marks the median,
#' and whiskers extend to observations within 1.5 interquartile ranges.
#' Diamonds show the exponentiated mean fitted log relative abundance
#' conditional on presence, standardized over all observed covariate profiles.
#' Both are displayed as percentages on a log scale. The natural log of the
#' ratio of fitted abundances (comparison group over reference group) gives
#' the raw abundance contrast. Panel A shows the RA test after cross-taxon
#' reference correction.
#'
#' @param x A `sacat` object fitted with `component = "all"` and
#'   `store_plot_data = TRUE`.
#' @param features Optional unique taxon names, displayed in the supplied order.
#' @param selection Select omnibus-significant taxa, the top-ranked taxa,
#'   or all eligible taxa. Ignored when `features` is supplied.
#' @param max_features Maximum number of automatically selected taxa.
#'   Ignored for `selection = "all"` or explicit `features`.
#' @param alpha Adjusted-p threshold for selection and component markers.
#' @param group_colors Two colors, in reference/comparison order, for observed
#'   abundance distributions.
#' @param zero_colors Three colors for detected samples, estimated SA, and
#'   estimated ND, respectively.
#' @param file Optional PDF, PNG, or SVG filename. Without it, draw on the
#'   current graphics device.
#' @param width,height Output dimensions in inches when `file` is supplied.
#'   Defaults to 10.5 inches wide and a height based on the number of taxa.
#' @param dpi Resolution of PNG output.
#' @param group_labels Optional reference/comparison labels for display.
#' @param show_component_guides Show BH component discovery boundaries when
#'   `p_adjust_method` is `"BH"` or `"fdr"`. A component with no discoveries
#'   has no boundary.
#' @param ... Reserved for graphical extensions.
#'
#' @return Invisibly, a list containing the displayed data and graphical settings.
#'
#' @examples
#' \dontrun{
#' example("sacat", run.donttest = TRUE)
#' plot(
#'     fit, features = c("Coprococcus", "Enterococcus", "Faecalibacterium"),
#'     group_labels = c("Healthy", "CDI"), file = "sacat-profile.pdf"
#' )
#' }
#'
#' @method plot sacat
#' @export
#' @importFrom graphics plot
plot.sacat <- function(
        x,
        features = NULL,
        selection = c("significant", "top", "all"),
        max_features = 24L,
        alpha = 0.05,
        group_colors = c("#4C78A8", "#E45756"),
        zero_colors = c("#6A9DC0", "#E5A36C", "#E9E9E9"),
        file = NULL,
        width = NULL,
        height = NULL,
        dpi = 300,
        group_labels = NULL,
        show_component_guides = TRUE,
        ...) {
    selection <- match.arg(selection)
    spec <- .sacat_plot_build_spec(
        x, features, selection, max_features, alpha, group_colors,
        zero_colors, show_component_guides
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
        if (is.null(width)) width <- 10.5
        if (is.null(height)) height <- max(3.5, 2.35 + 0.75 * nrow(spec$data))
        if (!is.numeric(width) || length(width) != 1L ||
            !is.finite(width) || width <= 0 ||
            !is.numeric(height) || length(height) != 1L ||
            !is.finite(height) || height <= 0) {
            stop("`width` and `height` must be positive numbers.",
                 call. = FALSE)
        }
        if (!is.numeric(dpi) || length(dpi) != 1L || !is.finite(dpi) || dpi <= 0) {
            stop("`dpi` must be one positive number.", call. = FALSE)
        }
        .sacat_plot_open_device(file, width, height, dpi)
        on.exit(grDevices::dev.off(), add = TRUE)
    }
    .sacat_plot_draw(spec)
    invisible(spec)
}
