# Shared preparation for public Qiita count tables
#
# Each dataset script supplies its cohort definition and adjustment variables.
# This function reads an already downloaded BIOM table and metadata file,
# retains the unfiltered sample totals as library sizes, aggregates features to
# genus or the deepest resolved higher rank, and writes the same four-part RDS
# object used by the existing real-data analyses.

read_qiita_biom <- function(path) {
    if (!requireNamespace("biomformat", quietly = TRUE)) {
        stop("Package 'biomformat' is required to read the BIOM table.",
             call. = FALSE)
    }

    biom <- suppressWarnings(biomformat::read_biom(path))
    counts <- as.matrix(biomformat::biom_data(biom))
    taxonomy <- biomformat::observation_metadata(biom)

    if (!is.numeric(counts) || anyNA(counts) || any(counts < 0) ||
        any(counts != round(counts))) {
        stop("The BIOM table must contain non-negative integer counts.",
             call. = FALSE)
    }
    if (is.null(rownames(counts)) || is.null(colnames(counts))) {
        stop("The BIOM table must contain feature and sample identifiers.",
             call. = FALSE)
    }
    if (is.null(taxonomy) || nrow(taxonomy) != nrow(counts)) {
        stop(
            paste(
                "Taxonomy is absent from the BIOM table.",
                "Fetch the Qiita study with taxonomy included."
            ),
            call. = FALSE
        )
    }

    if (!is.null(rownames(taxonomy)) &&
        setequal(rownames(taxonomy), rownames(counts))) {
        taxonomy <- taxonomy[rownames(counts), , drop = FALSE]
    }

    taxonomy_columns <- grep(
        "taxon|taxonomy", names(taxonomy), ignore.case = TRUE, value = TRUE
    )
    if (!length(taxonomy_columns)) {
        taxonomy_columns <- names(taxonomy)
    }

    taxonomy_entries <- lapply(seq_len(nrow(taxonomy)), function(index) {
        unlist(taxonomy[index, taxonomy_columns, drop = FALSE],
               use.names = FALSE)
    })
    storage.mode(counts) <- "integer"
    list(counts = counts, taxonomy = taxonomy_entries)
}

taxonomy_bin <- function(entry) {
    pieces <- unlist(lapply(as.character(entry), function(value) {
        strsplit(value, ";", fixed = TRUE)[[1L]]
    }), use.names = FALSE)
    pieces <- trimws(pieces)
    pieces <- gsub("^['\"]+|['\"]+$", "", pieces)
    pieces <- sub("^\\[", "", pieces)
    pieces <- sub("\\]$", "", pieces)

    rank_names <- c(
        "kingdom", "phylum", "class", "order", "family", "genus",
        "species"
    )
    resolved <- stats::setNames(rep(NA_character_, length(rank_names)),
                                rank_names)

    for (index in seq_along(pieces)) {
        piece <- pieces[[index]]
        if (!nzchar(piece)) next

        rank <- NA_character_
        value <- piece
        if (grepl("^D_[0-6]__", piece)) {
            rank_index <- as.integer(sub("^D_([0-6])__.*", "\\1", piece)) + 1L
            rank <- rank_names[[rank_index]]
            value <- sub("^D_[0-6]__", "", piece)
        } else if (grepl("^[dkpcofgs]__", piece, ignore.case = TRUE)) {
            prefix <- tolower(substr(piece, 1L, 1L))
            rank <- c(
                d = "kingdom", k = "kingdom", p = "phylum", c = "class",
                o = "order", f = "family", g = "genus", s = "species"
            )[[prefix]]
            value <- sub("^[dkpcofgs]__", "", piece, ignore.case = TRUE)
        } else if (index <= length(rank_names)) {
            rank <- rank_names[[index]]
        }

        value <- trimws(value)
        value <- sub("^\\[", "", value)
        value <- sub("\\]$", "", value)
        unresolved <- !nzchar(value) || grepl(
            "^(unclassified|uncultured|unknown|unidentified|na|none)$",
            value,
            ignore.case = TRUE
        )
        if (!is.na(rank) && !unresolved) resolved[[rank]] <- value
    }

    if (!is.na(resolved[["genus"]])) {
        return(c(taxon = resolved[["genus"]], taxonomic_level = "genus"))
    }

    fallback <- c("family", "order", "class", "phylum", "kingdom")
    fallback <- fallback[!is.na(resolved[fallback])]
    if (!length(fallback)) {
        return(c(
            taxon = "Unclassified_Bacteria",
            taxonomic_level = "unclassified"
        ))
    }
    level <- fallback[[1L]]
    c(
        taxon = paste0("Unclassified_", resolved[[level]]),
        taxonomic_level = level
    )
}

read_qiita_metadata <- function(path, sample_id_column) {
    metadata <- utils::read.delim(
        path,
        check.names = FALSE,
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE
    )
    if (!(sample_id_column %in% names(metadata))) {
        stop(sprintf("Metadata column '%s' is missing.", sample_id_column),
             call. = FALSE)
    }
    names(metadata)[names(metadata) == sample_id_column] <- "sample_id"
    metadata$sample_id <- trimws(as.character(metadata$sample_id))
    if (anyNA(metadata$sample_id) || any(!nzchar(metadata$sample_id)) ||
        anyDuplicated(metadata$sample_id)) {
        stop("Qiita metadata sample identifiers must be unique and non-empty.",
             call. = FALSE)
    }

    character_columns <- vapply(metadata, is.character, logical(1))
    metadata[character_columns] <- lapply(
        metadata[character_columns],
        function(value) {
            value <- trimws(value)
            missing <- is.na(value) | !nzchar(value) | tolower(value) %in%
                c("na", "n/a", "not provided", "not applicable")
            value[missing] <- NA_character_
            value
        }
    )
    metadata
}

write_qiita_count_csv <- function(counts, path) {
    output <- data.frame(
        taxon = rownames(counts),
        counts,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    utils::write.csv(output, path, row.names = FALSE, quote = TRUE)
}

prepare_qiita_dataset <- function(dataset_directory, configuration) {
    required_configuration <- c(
        "dataset_id", "dataset_title", "biom_file", "metadata_file",
        "select_samples", "group_column", "reference", "comparison",
        "covariates", "subject_column", "source_publication",
        "source_publication_doi", "primary_comparison"
    )
    missing_configuration <- setdiff(
        required_configuration, names(configuration)
    )
    if (length(missing_configuration)) {
        stop(sprintf(
            "Dataset configuration is missing: %s",
            paste(missing_configuration, collapse = ", ")
        ), call. = FALSE)
    }

    biom_path <- file.path(dataset_directory, configuration$biom_file)
    metadata_path <- file.path(dataset_directory, configuration$metadata_file)
    if (!file.exists(biom_path) || !file.exists(metadata_path)) {
        stop(
            paste(
                "The prepared source files are missing.",
                sprintf("Expected %s and %s.", biom_path, metadata_path)
            ),
            call. = FALSE
        )
    }

    message(sprintf("Reading Qiita study for %s", configuration$dataset_id))
    biom_input <- read_qiita_biom(biom_path)
    all_counts <- biom_input$counts
    all_library_sizes <- colSums(all_counts)

    sample_id_column <- if (!is.null(configuration$sample_id_column)) {
        configuration$sample_id_column
    } else {
        "#SampleID"
    }
    metadata_raw <- read_qiita_metadata(metadata_path, sample_id_column)

    if (!is.null(configuration$prepare_metadata)) {
        original_sample_ids <- metadata_raw$sample_id
        metadata_raw <- configuration$prepare_metadata(metadata_raw)
        if (!is.data.frame(metadata_raw) ||
            !identical(metadata_raw$sample_id, original_sample_ids)) {
            stop(
                paste(
                    "prepare_metadata must return a data frame with sample",
                    "identifiers in their original order."
                ),
                call. = FALSE
            )
        }
    }

    numeric_covariates <- if (is.null(configuration$numeric_covariates)) {
        character()
    } else {
        configuration$numeric_covariates
    }
    if (!all(numeric_covariates %in% names(configuration$covariates))) {
        stop(
            "numeric_covariates must name output variables in covariates.",
            call. = FALSE
        )
    }

    required_metadata <- unique(c(
        configuration$group_column,
        unname(configuration$covariates),
        configuration$subject_column,
        configuration$selection_columns
    ))
    missing_metadata <- setdiff(required_metadata, names(metadata_raw))
    if (length(missing_metadata)) {
        stop(sprintf(
            "Qiita metadata is missing: %s",
            paste(missing_metadata, collapse = ", ")
        ), call. = FALSE)
    }
    numeric_source_columns <- unname(
        configuration$covariates[numeric_covariates]
    )
    nonnumeric_sources <- numeric_source_columns[!vapply(
        metadata_raw[, numeric_source_columns, drop = FALSE],
        is.numeric,
        logical(1)
    )]
    if (length(nonnumeric_sources)) {
        stop(sprintf(
            paste(
                "Numeric covariate source columns must be numeric after",
                "prepare_metadata: %s"
            ),
            paste(nonnumeric_sources, collapse = ", ")
        ), call. = FALSE)
    }

    common_samples <- intersect(colnames(all_counts), metadata_raw$sample_id)
    if (!length(common_samples)) {
        stop("No sample identifiers are shared by the BIOM table and metadata.",
             call. = FALSE)
    }
    metadata <- metadata_raw[match(common_samples, metadata_raw$sample_id), ,
                             drop = FALSE]
    metadata$library_size <- as.numeric(all_library_sizes[common_samples])
    positive_depth <- is.finite(metadata$library_size) &
        metadata$library_size > 0

    selection <- configuration$select_samples(metadata)
    if (!is.list(selection) || is.null(selection$keep) ||
        length(selection$keep) != nrow(metadata)) {
        stop("select_samples must return a list containing one keep value per sample.",
             call. = FALSE)
    }
    selected_by_definition <- !is.na(selection$keep) & selection$keep

    group_value <- as.character(metadata[[configuration$group_column]])
    group_complete <- group_value %in%
        c(configuration$reference, configuration$comparison)
    adjustment_source_columns <- unname(configuration$covariates)
    adjustment_complete <- stats::complete.cases(
        metadata[, adjustment_source_columns, drop = FALSE]
    )
    if (length(numeric_source_columns)) {
        numeric_finite <- vapply(
            seq_len(nrow(metadata)),
            function(index) all(vapply(
                metadata[index, numeric_source_columns, drop = FALSE],
                is.finite,
                logical(1)
            )),
            logical(1)
        )
        adjustment_complete <- adjustment_complete & numeric_finite
    }
    subject_complete <- !is.na(metadata[[configuration$subject_column]]) &
        nzchar(as.character(metadata[[configuration$subject_column]]))
    eligible <- positive_depth & selected_by_definition & group_complete &
        adjustment_complete & subject_complete

    candidate_metadata <- metadata[eligible, , drop = FALSE]
    candidate_metadata$group_value <- group_value[eligible]
    candidate_metadata <- candidate_metadata[order(
        as.character(candidate_metadata[[configuration$subject_column]]),
        -candidate_metadata$library_size,
        candidate_metadata$sample_id
    ), , drop = FALSE]

    one_sample_per_subject <- !identical(
        configuration$one_sample_per_subject, FALSE
    )
    if (one_sample_per_subject) {
        duplicate_subject <- duplicated(
            as.character(candidate_metadata[[configuration$subject_column]])
        )
        duplicate_sample_ids <- candidate_metadata$sample_id[duplicate_subject]
        candidate_metadata <- candidate_metadata[!duplicate_subject, ,
                                                 drop = FALSE]
    } else {
        duplicate_sample_ids <- character()
    }

    candidate_metadata <- candidate_metadata[order(
        factor(candidate_metadata$group_value,
               levels = c(configuration$reference,
                          configuration$comparison)),
        candidate_metadata$sample_id
    ), , drop = FALSE]

    analysis_metadata <- data.frame(
        sample_id = candidate_metadata$sample_id,
        group = factor(
            candidate_metadata$group_value,
            levels = c(configuration$reference, configuration$comparison)
        ),
        library_size = as.numeric(candidate_metadata$library_size),
        subject_id = as.character(
            candidate_metadata[[configuration$subject_column]]
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    numeric_covariate_standardization <- list()
    for (output_name in names(configuration$covariates)) {
        source_name <- configuration$covariates[[output_name]]
        if (output_name %in% numeric_covariates) {
            value <- candidate_metadata[[source_name]]
            center <- mean(value)
            spread <- stats::sd(value)
            if (!is.finite(spread) || spread <= 0) {
                stop(sprintf(
                    "Numeric covariate '%s' has no variation in the final cohort.",
                    output_name
                ), call. = FALSE)
            }
            analysis_metadata[[output_name]] <- as.numeric(
                (value - center) / spread
            )
            numeric_covariate_standardization[[output_name]] <- list(
                source_column = source_name,
                center = center,
                scale = spread
            )
        } else {
            value <- as.character(candidate_metadata[[source_name]])
            analysis_metadata[[output_name]] <- factor(
                value, levels = sort(unique(value))
            )
        }
    }
    source_metadata_columns <- intersect(
        configuration$source_metadata_columns, names(candidate_metadata)
    )
    for (column in setdiff(
            source_metadata_columns, names(analysis_metadata))) {
        analysis_metadata[[column]] <- candidate_metadata[[column]]
    }
    rownames(analysis_metadata) <- analysis_metadata$sample_id

    primary_counts <- all_counts[, analysis_metadata$sample_id, drop = FALSE]
    taxon_map_matrix <- t(vapply(
        biom_input$taxonomy,
        taxonomy_bin,
        character(2)
    ))
    taxon_map <- data.frame(
        feature_id = rownames(all_counts),
        taxon = taxon_map_matrix[, "taxon"],
        taxonomic_level = taxon_map_matrix[, "taxonomic_level"],
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    taxon_counts <- rowsum(
        primary_counts,
        group = taxon_map$taxon,
        reorder = TRUE
    )
    storage.mode(taxon_counts) <- "integer"
    minimum_positive_samples <- ceiling(0.05 * ncol(taxon_counts))
    prevalence_keep <- rowSums(taxon_counts > 0) >= minimum_positive_samples
    taxon_counts <- taxon_counts[prevalence_keep, , drop = FALSE]

    reference_samples <- analysis_metadata$group == configuration$reference
    comparison_samples <- analysis_metadata$group == configuration$comparison
    reference_positive <- rowSums(
        taxon_counts[, reference_samples, drop = FALSE] > 0
    )
    comparison_positive <- rowSums(
        taxon_counts[, comparison_samples, drop = FALSE] > 0
    )
    group_support <- reference_positive > 0L & comparison_positive > 0L
    group_support_exclusions <- data.frame(
        taxon = rownames(taxon_counts)[!group_support],
        reference_positive_samples = as.integer(
            reference_positive[!group_support]
        ),
        comparison_positive_samples = as.integer(
            comparison_positive[!group_support]
        ),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    taxon_bins_before_group_support <- nrow(taxon_counts)
    taxon_counts <- taxon_counts[group_support, , drop = FALSE]

    taxon_levels <- tapply(
        taxon_map$taxonomic_level,
        taxon_map$taxon,
        function(value) paste(sort(unique(value)), collapse = ";")
    )
    source_features_by_taxon <- table(taxon_map$taxon)
    relative_abundance <- sweep(
        taxon_counts, 2L, analysis_metadata$library_size, "/"
    )
    taxon_metadata <- data.frame(
        taxon = rownames(taxon_counts),
        taxonomic_level = unname(taxon_levels[rownames(taxon_counts)]),
        source_feature_count = as.integer(
            source_features_by_taxon[rownames(taxon_counts)]
        ),
        total_count = as.numeric(rowSums(taxon_counts)),
        positive_samples = as.integer(rowSums(taxon_counts > 0)),
        reference_group = configuration$reference,
        comparison_group = configuration$comparison,
        positive_samples_reference = as.integer(rowSums(
            taxon_counts[, reference_samples, drop = FALSE] > 0
        )),
        positive_samples_comparison = as.integer(rowSums(
            taxon_counts[, comparison_samples, drop = FALSE] > 0
        )),
        prevalence = rowMeans(taxon_counts > 0),
        prevalence_reference = rowMeans(
            taxon_counts[, reference_samples, drop = FALSE] > 0
        ),
        prevalence_comparison = rowMeans(
            taxon_counts[, comparison_samples, drop = FALSE] > 0
        ),
        mean_relative_abundance = rowMeans(relative_abundance),
        mean_relative_abundance_reference = rowMeans(
            relative_abundance[, reference_samples, drop = FALSE]
        ),
        mean_relative_abundance_comparison = rowMeans(
            relative_abundance[, comparison_samples, drop = FALSE]
        ),
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    metadata_ids <- metadata_raw$sample_id
    disposition <- data.frame(
        sample_id = metadata_ids,
        in_biom = metadata_ids %in% colnames(all_counts),
        library_size = as.numeric(all_library_sizes[metadata_ids]),
        analysis_status = "not retained",
        exclusion_reason = "not represented in the BIOM table",
        stringsAsFactors = FALSE
    )
    common_index <- match(metadata$sample_id, disposition$sample_id)
    disposition$exclusion_reason[common_index] <- "outside the prespecified cohort"
    disposition$exclusion_reason[common_index[!positive_depth]] <-
        "non-positive library size"
    disposition$exclusion_reason[common_index[
        positive_depth & selected_by_definition & !group_complete
    ]] <- "not in the prespecified comparison groups"
    disposition$exclusion_reason[common_index[
        positive_depth & selected_by_definition & group_complete &
            !adjustment_complete
    ]] <- "missing a required adjustment variable"
    disposition$exclusion_reason[common_index[
        positive_depth & selected_by_definition & group_complete &
            adjustment_complete & !subject_complete
    ]] <- "missing subject identifier"
    duplicate_index <- match(duplicate_sample_ids, disposition$sample_id)
    disposition$exclusion_reason[duplicate_index] <-
        "additional sample from a retained subject"
    retained_index <- match(analysis_metadata$sample_id, disposition$sample_id)
    disposition$analysis_status[retained_index] <- "retained"
    disposition$exclusion_reason[retained_index] <- ""

    observed_group_counts <- table(analysis_metadata$group)
    preprocessing <- list(
        dataset_id = configuration$dataset_id,
        dataset_alias = configuration$dataset_title,
        source_publication = configuration$source_publication,
        source_publication_doi = configuration$source_publication_doi,
        data_source = configuration$data_source,
        primary_comparison = configuration$primary_comparison,
        count_source = configuration$count_source,
        cohort_selection = selection$note,
        cohort_selection_details = selection$details,
        repeated_sample_rule = if (!is.null(
            configuration$repeated_sample_rule
        )) {
            configuration$repeated_sample_rule
        } else if (one_sample_per_subject) {
            paste(
                "one sample per subject; the sample with the largest original",
                "library size was retained, with sample identifier breaking ties"
            )
        } else {
            "all eligible samples retained"
        },
        current_count_processing = paste(
            "unrarefied BIOM features aggregated to genus or deepest resolved",
            "higher rank; taxa retained at 5% prevalence and required to have",
            "a positive count in both comparison groups"
        ),
        prevalence_filter_scope = "final primary cohort",
        minimum_positive_samples = minimum_positive_samples,
        group_support_filter =
            "at least one positive-count sample in each comparison group",
        group_support_exclusions = group_support_exclusions,
        rarefaction = "not applied",
        reference_group = configuration$reference,
        comparison_group = configuration$comparison,
        planned_adjustment_variables = names(configuration$covariates),
        numeric_covariate_standardization =
            numeric_covariate_standardization,
        dasra_input_orientation = "taxa by samples",
        dasra_group_column = "group",
        dasra_library_size_column = "library_size",
        library_size_definition =
            "sum of the full downloaded BIOM feature table before filtering",
        sample_counts = as.list(observed_group_counts),
        feature_counts = list(
            source_features = nrow(all_counts),
            taxon_bins_before_prevalence_filter = length(unique(taxon_map$taxon)),
            taxon_bins_after_prevalence_filter =
                taxon_bins_before_group_support,
            taxa_excluded_without_two_group_positive_support =
                nrow(group_support_exclusions),
            retained_taxon_bins = nrow(taxon_counts)
        )
    )

    summary_table <- data.frame(
        item = c(
            "dataset", "source DOI", "data source", "primary comparison",
            "cohort selection", "repeated-sample rule", "reference samples",
            "comparison samples", "source features", "prevalence rule",
            "taxonomic aggregation", "two-group positive-count support rule",
            "retained taxon bins", "DASRA count orientation",
            "DASRA library-size definition"
        ),
        value = c(
            configuration$dataset_title,
            configuration$source_publication_doi,
            configuration$data_source,
            configuration$primary_comparison,
            selection$note,
            preprocessing$repeated_sample_rule,
            as.character(observed_group_counts[[configuration$reference]]),
            as.character(observed_group_counts[[configuration$comparison]]),
            as.character(nrow(all_counts)),
            sprintf(
                "present in at least %d of %d primary samples (5%%)",
                minimum_positive_samples, ncol(taxon_counts)
            ),
            "genus when resolved; otherwise deepest resolved higher-rank bin",
            preprocessing$group_support_filter,
            as.character(nrow(taxon_counts)),
            preprocessing$dasra_input_orientation,
            preprocessing$library_size_definition
        ),
        stringsAsFactors = FALSE
    )

    output_directory <- file.path(dataset_directory, "processed")
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    prefix <- configuration$dataset_id
    analysis_input <- list(
        counts = taxon_counts,
        metadata = analysis_metadata,
        taxon_metadata = taxon_metadata,
        preprocessing = preprocessing
    )

    saveRDS(
        analysis_input,
        file.path(output_directory, paste0(prefix, "_dasra_input.rds")),
        compress = "xz"
    )
    write_qiita_count_csv(
        taxon_counts,
        file.path(output_directory, paste0(prefix, "_taxon_counts.csv"))
    )
    utils::write.csv(
        analysis_metadata,
        file.path(output_directory, paste0(prefix, "_metadata.csv")),
        row.names = FALSE,
        quote = TRUE
    )
    utils::write.csv(
        taxon_metadata,
        file.path(output_directory, paste0(prefix, "_taxon_metadata.csv")),
        row.names = FALSE,
        quote = TRUE
    )
    utils::write.csv(
        disposition,
        file.path(output_directory, paste0(prefix, "_sample_disposition.csv")),
        row.names = FALSE,
        quote = TRUE,
        na = ""
    )
    utils::write.csv(
        summary_table,
        file.path(
            output_directory,
            paste0(prefix, "_preprocessing_summary.csv")
        ),
        row.names = FALSE,
        quote = TRUE
    )

    message(sprintf(
        "%s preparation complete: %d %s, %d %s, %d retained taxon bins.",
        configuration$dataset_id,
        observed_group_counts[[configuration$reference]],
        configuration$reference,
        observed_group_counts[[configuration$comparison]],
        configuration$comparison,
        nrow(taxon_counts)
    ))
    invisible(analysis_input)
}
