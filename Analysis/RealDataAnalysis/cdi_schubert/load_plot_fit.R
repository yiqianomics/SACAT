# Reuse the CDI fit only when its plot data and tests match the completed analysis.
load_cdi_plot_fit <- function(input, results, fit_path, refit = FALSE) {
    same_numbers <- function(x, y) {
        is.numeric(x) && is.numeric(y) && length(x) == length(y) &&
            identical(is.na(x), is.na(y)) &&
            all(abs(x[!is.na(x)] - y[!is.na(y)]) <=
                1e-10 * pmax(abs(x[!is.na(x)]), abs(y[!is.na(y)]),
                            .Machine$double.xmin))
    }
    matches_analysis <- function(fit) {
        if (!inherits(fit, "sacat") ||
            !identical(fit$plot_data$version, 2L) ||
            !all(c("profiles", "zero_sources", "positive_abundance") %in%
                 names(fit$plot_data)) ||
            !identical(unname(fit$settings$contrast), c("H", "CDI")) ||
            fit$settings$n_samples != ncol(input$counts) ||
            fit$settings$n_taxa != nrow(input$counts) ||
            !identical(fit$settings$p_adjust_method, "BH") ||
            !setequal(fit$results$taxon, rownames(input$counts)) ||
            anyDuplicated(fit$results$taxon)) {
            return(FALSE)
        }
        definitions <- list(
            "SACAT structural absence" = c(
                p_value = "p_structural_absence",
                statistic = "z_structural_absence"
            ),
            "SACAT present-conditional abundance" = c(
                p_value = "p_relative_abundance",
                statistic = "z_relative_abundance",
                estimate = "estimate_relative_abundance"
            ),
            "SACAT combined" = c(p_value = "p_omnibus")
        )
        availability_columns <- c(
            "SACAT structural absence" = "formed_structural_absence",
            "SACAT present-conditional abundance" = "formed_relative_abundance",
            "SACAT combined" = "formed_omnibus"
        )
        all(vapply(names(definitions), function(method) {
            rows <- results[results$method == method, , drop = FALSE]
            if (nrow(rows) != nrow(input$counts) ||
                !setequal(rows$taxon, fit$results$taxon) ||
                anyDuplicated(rows$taxon)) return(FALSE)
            fitted_rows <- fit$results[
                match(rows$taxon, fit$results$taxon), , drop = FALSE
            ]
            diagnostics <- fit$diagnostics[
                match(rows$taxon, fit$diagnostics$taxon), , drop = FALSE
            ]
            available <- diagnostics$retained &
                diagnostics[[availability_columns[[method]]]]
            if (!identical(as.logical(rows$available), available)) return(FALSE)
            columns <- definitions[[method]]
            all(vapply(names(columns), function(column) {
                same_numbers(rows[[column]], fitted_rows[[columns[[column]]]])
            }, logical(1)))
        }, logical(1)))
    }
    fit <- if (file.exists(fit_path)) readRDS(fit_path) else NULL
    valid <- isTRUE(tryCatch(matches_analysis(fit), error = function(e) FALSE))
    if (!valid && !refit) {
        stop(
            "A saved CDI fit with current plot data matching the completed ",
            "method-results table is required: ", fit_path,
            call. = FALSE
        )
    }
    if (!valid) {
        message("Fitting CDI SACAT to prepare current plot data.")
        fit <- SACAT::sacat(
            counts = input$counts,
            metadata = input$metadata,
            formula = ~ age_z + sex + antibiotics_3mo + group,
            group = "group", reference = "H",
            library_size = "library_size", store_plot_data = TRUE
        )
        if (!isTRUE(matches_analysis(fit))) {
            stop(
                "The current SACAT fit differs from the completed CDI results. ",
                "Rerun the CDI analysis before making figures.",
                call. = FALSE
            )
        }
        dir.create(dirname(fit_path), recursive = TRUE, showWarnings = FALSE)
        saveRDS(fit, fit_path)
    }
    fit
}
