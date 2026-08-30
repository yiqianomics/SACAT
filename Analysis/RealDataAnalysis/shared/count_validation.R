validate_identifiers <- function(values, label) {
    values <- as.character(values)
    if (!length(values) || anyNA(values) || any(!nzchar(values)) ||
        anyDuplicated(values)) {
        stop(sprintf("%s must be nonempty and unique.", label), call. = FALSE)
    }
    invisible(values)
}

as_count_matrix <- function(values, label = "The count table") {
    values <- as.matrix(values)
    numeric_values <- suppressWarnings(as.numeric(values))
    if (!length(numeric_values) || any(!is.finite(numeric_values)) ||
        any(numeric_values < 0) ||
        any(abs(numeric_values - round(numeric_values)) > 1e-8) ||
        any(numeric_values > .Machine$integer.max)) {
        stop(
            sprintf("%s must contain nonnegative integer counts.", label),
            call. = FALSE
        )
    }
    counts <- matrix(
        round(numeric_values),
        nrow = nrow(values),
        ncol = ncol(values),
        dimnames = dimnames(values)
    )
    storage.mode(counts) <- "integer"
    invisible(counts)
}
