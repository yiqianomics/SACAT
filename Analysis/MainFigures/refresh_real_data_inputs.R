#!/usr/bin/env Rscript

# Export manuscript plot data from completed analyses without fitting models.
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
directory <- dirname(normalizePath(script))
analysis_directory <- file.path(directory, "..", "RealDataAnalysis")
cdi_directory <- file.path(analysis_directory, "cdi_schubert")
data_directory <- file.path(directory, "data", "real_data")
taxa <- c("Coprococcus", "Enterococcus", "Faecalibacterium")
groups <- c("Healthy", "CDI")
depths <- c("native", "1000", "750", "500", "250", "100", "50")

source(file.path(cdi_directory, "load_plot_fit.R"))
input <- readRDS(file.path(
    cdi_directory, "processed", "cdi_schubert_sacat_input.rds"
))
results <- read.csv(file.path(
    cdi_directory, "table", "schubert_cdi_method_results_all_taxa.csv"
), stringsAsFactors = FALSE)
fit <- load_cdi_plot_fit(
    input, results, file.path(cdi_directory, "work", "sacat_fit.rds")
)
stopifnot(length(unique(results$method)) == 15L,
          nrow(input$counts) == 88L, ncol(input$counts) == 243L,
          identical(colnames(input$counts), rownames(input$metadata)),
          !anyDuplicated(results[, c("method", "taxon")]))
for (method in unique(results$method)) {
    rows <- results[results$method == method, ]
    stopifnot(nrow(rows) == 88L,
              setequal(rows$taxon, rownames(input$counts)),
              !anyNA(rows$available),
              all(is.finite(rows$p_for_bh)),
              all(rows$p_for_bh >= 0 & rows$p_for_bh <= 1),
              max(abs(p.adjust(rows$p_for_bh, "BH") - rows$q_value)) < 1e-12,
              all(rows$significant == (rows$available & rows$q_value <= .05)))
}

zero_sources <- fit$plot_data$zero_sources
zero_sources <- zero_sources[zero_sources$feature %in% taxa, , drop = FALSE]
stopifnot(nrow(zero_sources) == 6L, all(zero_sources$status == "ok"),
          !anyDuplicated(zero_sources[, c("feature", "group")]),
          setequal(zero_sources$group, 0:1))
zero_sources$taxon <- zero_sources$feature
zero_sources$group <- groups[zero_sources$group + 1L]
names(zero_sources)[names(zero_sources) == "detected"] <- "observed_detection"
zero_sources <- zero_sources[, c(
    "taxon", "group", "observed_detection", "estimated_sa", "estimated_nd", "n"
)]
shares <- as.matrix(zero_sources[, c(
    "observed_detection", "estimated_sa", "estimated_nd"
)])
stopifnot(all(is.finite(shares)), all(shares >= -1e-12 & shares <= 1 + 1e-12),
          max(abs(rowSums(shares) - 1)) < 1e-12)
for (i in seq_len(nrow(zero_sources))) {
    row <- zero_sources[i, ]
    selected <- as.character(input$metadata$group) ==
        c(Healthy = "H", CDI = "CDI")[[row$group]]
    stopifnot(row$n == sum(selected),
              abs(row$observed_detection -
                  mean(input$counts[row$taxon, selected] > 0)) < 1e-12)
}

positive <- fit$plot_data$positive_abundance
positive <- positive[positive$feature %in% taxa, , drop = FALSE]
positive$taxon <- positive$feature
positive$group <- groups[positive$group + 1L]
names(positive)[names(positive) == "abundance"] <- "relative_abundance"
positive <- positive[, c("taxon", "group", "relative_abundance")]
stopifnot(!anyNA(positive), all(is.finite(positive$relative_abundance)),
          all(positive$relative_abundance > 0 & positive$relative_abundance <= 100))
for (taxon in taxa) for (group in groups) {
    selected <- as.character(input$metadata$group) ==
        c(Healthy = "H", CDI = "CDI")[[group]]
    observed <- 100 * input$counts[taxon, selected] /
        input$metadata$library_size[selected]
    stored <- positive$relative_abundance[
        positive$taxon == taxon & positive$group == group
    ]
    stopifnot(isTRUE(all.equal(sort(stored), sort(observed[observed > 0]),
                              check.attributes = FALSE, tolerance = 1e-12)))
}

profiles <- fit$plot_data$profiles[
    match(taxa, fit$plot_data$profiles$feature), , drop = FALSE
]
stopifnot(identical(profiles$feature, taxa),
          all(profiles$abundance_status == "ok"))
fitted <- rbind(
    data.frame(taxon = taxa, group = "Healthy",
               abundance_percent = profiles$abundance_reference),
    data.frame(taxon = taxa, group = "CDI",
               abundance_percent = profiles$abundance_comparison)
)
stopifnot(all(is.finite(fitted$abundance_percent)),
          all(fitted$abundance_percent > 0 & fitted$abundance_percent <= 100))

nonnull <- read.csv(file.path(
    analysis_directory, "rauer_mock", "table", "rauer_mock_nonnull_results.csv"
), stringsAsFactors = FALSE)
definitions <- data.frame(
    target = rep(c("Structural", "Abundance"), each = 3L),
    method = rep(c("SACAT", "MaAsLin3", "ZINQ"), 2L),
    source_method = c(
        "SACAT structural absence", "MaAsLin3 prevalence", "ZINQ prevalence",
        "SACAT present-conditional abundance", "MaAsLin3 abundance", "ZINQ abundance"
    ),
    n_targets = rep(c(8L, 3L), each = 3L)
)
stopifnot(nrow(nonnull) == 52L * 7L * 6L * 11L,
          identical(unique(nonnull$design), "nonnull"),
          setequal(nonnull$method, definitions$source_method),
          setequal(nonnull$depth, depths),
          setequal(nonnull$allocation_id, seq_len(52L)),
          !anyDuplicated(nonnull[, c("method", "depth", "allocation_id", "taxon")]),
          !anyNA(nonnull$available), !anyNA(nonnull$significant),
          all(nonnull$significant == (nonnull$available & nonnull$q_value <= .05)))
blocks <- split(seq_len(nrow(nonnull)), interaction(
    nonnull$allocation_id, nonnull$depth, nonnull$method, drop = TRUE
))
stopifnot(all(vapply(blocks, function(index) {
    rows <- nonnull[index, ]
    expected <- p.adjust(ifelse(rows$available, rows$p_value, 1), "BH")
    nrow(rows) == 11L && all(is.finite(expected)) &&
        max(abs(expected - rows$q_value)) < 1e-12
}, logical(1))))
recovery <- do.call(rbind, lapply(seq_len(nrow(definitions)), function(i) {
    definition <- definitions[i, ]
    do.call(rbind, lapply(depths, function(depth) {
        rows <- nonnull[
            nonnull$method == definition$source_method &
                nonnull$target_set == paste(definition$target, "targets") &
                nonnull$depth == depth, , drop = FALSE
        ]
        stopifnot(nrow(rows) == 52L * definition$n_targets,
                  length(unique(rows$taxon)) == definition$n_targets,
                  all(table(rows$allocation_id) == definition$n_targets))
        discoveries <- sum(rows$significant)
        data.frame(target = definition$target, method = definition$method,
                   depth = depth, n_allocations = 52L,
                   n_targets = definition$n_targets, n_discoveries = discoveries,
                   recovery_percent = 100 * discoveries / (52L * definition$n_targets))
    }))
}))

outputs <- list(
    schubert_cdi_figure_method_results = results[, c(
        "method", "taxon", "available", "p_value", "p_for_bh", "q_value",
        "significant", "estimate", "statistic"
    )],
    schubert_cdi_figure_zero_sources = zero_sources,
    schubert_cdi_figure_positive_abundance = positive,
    schubert_cdi_figure_fitted_abundance = fitted,
    rauer_mock_recovery = recovery
)
dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)
staging_directory <- tempfile("real-data-inputs-", tmpdir = data_directory)
dir.create(staging_directory)
for (name in names(outputs)) {
    write.csv(outputs[[name]], file.path(staging_directory, paste0(name, ".csv")),
              row.names = FALSE, na = "")
}
staged <- list.files(staging_directory, full.names = TRUE)
stopifnot(all(file.copy(staged, data_directory, overwrite = TRUE)))
unlink(staging_directory, recursive = TRUE)
message("Refreshed manuscript inputs from the completed CDI and spike-in analyses.")
