# Prepare the Zupancic Old Order Amish obesity cohort from MicrobiomeHD

options(stringsAsFactors = FALSE)

dataset_directory <- dirname(normalizePath(sub(
    "^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]
)))
source_directory <- file.path(dataset_directory, "raw", "ob_zupancic_results")
count_file <- file.path(
    source_directory, "RDP", "ob_zupancic.otu_table.100.denovo.rdp_assigned"
)
metadata_file <- file.path(source_directory, "ob_zupancic.metadata.txt")
if (!file.exists(count_file) || !file.exists(metadata_file)) {
    stop("Place the extracted ob_zupancic_results directory under raw/.",
         call. = FALSE)
}
if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
}
source(file.path(
    dirname(dataset_directory), "shared", "count_validation.R"
))

read_count_table <- function(path) {
    table <- data.table::fread(
        path, sep = "\t", data.table = FALSE, check.names = FALSE,
        showProgress = interactive()
    )
    feature_ids <- validate_identifiers(
        table[[1L]], "Feature identifiers"
    )
    table[[1L]] <- NULL
    counts <- as_count_matrix(table, "The MicrobiomeHD count table")
    validate_identifiers(colnames(counts), "Sample identifiers")
    rownames(counts) <- feature_ids
    counts
}

parse_rdp_taxon <- function(lineage) {
    pieces <- strsplit(lineage, ";", fixed = TRUE)[[1L]]
    prefixes <- c(
        kingdom = "k__", phylum = "p__", class = "c__",
        order = "o__", family = "f__", genus = "g__"
    )
    values <- vapply(prefixes, function(prefix) {
        match <- pieces[startsWith(pieces, prefix)]
        if (!length(match)) return(NA_character_)
        value <- sub(prefix, "", match[[1L]], fixed = TRUE)
        if (nzchar(value)) value else NA_character_
    }, character(1))
    if (!is.na(values[["genus"]])) {
        return(c(taxon = values[["genus"]], taxonomic_level = "genus"))
    }
    fallback <- c("family", "order", "class", "phylum", "kingdom")
    fallback <- fallback[!is.na(values[fallback])]
    if (!length(fallback)) {
        return(c(taxon = "Unclassified_Bacteria",
                 taxonomic_level = "unclassified"))
    }
    level <- fallback[[1L]]
    c(taxon = paste0("Unclassified_", values[[level]]),
      taxonomic_level = level)
}

write_count_csv <- function(counts, path) {
    utils::write.csv(
        data.frame(taxon = rownames(counts), counts, check.names = FALSE),
        path, row.names = FALSE, quote = TRUE
    )
}

message("Reading the Zupancic MicrobiomeHD count table")
all_counts <- read_count_table(count_file)
all_library_sizes <- colSums(all_counts)
metadata_raw <- data.table::fread(
    metadata_file, sep = "\t", data.table = FALSE, check.names = FALSE,
    showProgress = FALSE
)
metadata_raw$sample_id <- as.character(metadata_raw[[1L]])
rownames(metadata_raw) <- metadata_raw$sample_id

common_samples <- intersect(colnames(all_counts), metadata_raw$sample_id)
counts <- all_counts[, common_samples, drop = FALSE]
metadata <- metadata_raw[common_samples, , drop = FALSE]
baseline_group_counts_before_filter <- table(factor(
    metadata$DiseaseState[
        metadata$visit_number == 1 &
            metadata$DiseaseState %in% c("H", "OB")
    ],
    levels = c("H", "OB")
))

sample_keep <- colSums(counts) > 100
counts <- counts[, sample_keep, drop = FALSE]
metadata <- metadata[colnames(counts), , drop = FALSE]
baseline_group_counts_after_initial_depth <- table(factor(
    metadata$DiseaseState[
        metadata$visit_number == 1 &
            metadata$DiseaseState %in% c("H", "OB")
    ],
    levels = c("H", "OB")
))
source_feature_count <- nrow(counts)
feature_keep <- rowSums(counts) >= 10 & rowMeans(counts > 0) > 0.01
counts <- counts[feature_keep, , drop = FALSE]
sample_keep <- colSums(counts) > 100
counts <- counts[, sample_keep, drop = FALSE]
metadata <- metadata[colnames(counts), , drop = FALSE]

subject_id <- trimws(as.character(metadata$submitted_subject_id_s))
eligible <- metadata$DiseaseState %in% c("H", "OB") &
    metadata$visit_number == 1 & metadata$sex_s %in% c("female", "male") &
    !is.na(subject_id) & nzchar(subject_id)
candidate_metadata <- metadata[eligible, , drop = FALSE]
candidate_metadata$subject_id <- subject_id[eligible]
candidate_metadata$library_size <- as.numeric(
    all_library_sizes[candidate_metadata$sample_id]
)
candidate_metadata <- candidate_metadata[order(
    candidate_metadata$subject_id,
    -candidate_metadata$library_size,
    candidate_metadata$sample_id
), , drop = FALSE]
duplicate_subject_sample <- duplicated(candidate_metadata$subject_id)
duplicate_sample_ids <- candidate_metadata$sample_id[duplicate_subject_sample]
candidate_metadata <- candidate_metadata[!duplicate_subject_sample, ,
                                         drop = FALSE]
candidate_metadata <- candidate_metadata[order(
    factor(candidate_metadata$DiseaseState, levels = c("H", "OB")),
    candidate_metadata$sample_id
), , drop = FALSE]

analysis_metadata <- data.frame(
    sample_id = candidate_metadata$sample_id,
    group = factor(candidate_metadata$DiseaseState, levels = c("H", "OB")),
    library_size = candidate_metadata$library_size,
    subject_id = candidate_metadata$subject_id,
    sex = factor(candidate_metadata$sex_s, levels = c("female", "male")),
    visit = as.integer(candidate_metadata$visit_number),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
rownames(analysis_metadata) <- analysis_metadata$sample_id

primary_counts <- counts[, analysis_metadata$sample_id, drop = FALSE]
taxon_map_matrix <- t(vapply(
    rownames(primary_counts), parse_rdp_taxon, character(2)
))
taxon_map <- data.frame(
    feature_id = rownames(primary_counts),
    taxon = taxon_map_matrix[, "taxon"],
    taxonomic_level = taxon_map_matrix[, "taxonomic_level"],
    stringsAsFactors = FALSE
)
taxon_counts <- rowsum(primary_counts, taxon_map$taxon, reorder = TRUE)
storage.mode(taxon_counts) <- "integer"
minimum_positive_samples <- ceiling(0.05 * ncol(taxon_counts))
taxon_counts <- taxon_counts[
    rowSums(taxon_counts > 0) >= minimum_positive_samples, , drop = FALSE
]
reference_samples <- analysis_metadata$group == "H"
comparison_samples <- analysis_metadata$group == "OB"
group_support <-
    rowSums(taxon_counts[, reference_samples, drop = FALSE] > 0) > 0 &
    rowSums(taxon_counts[, comparison_samples, drop = FALSE] > 0) > 0
taxon_counts <- taxon_counts[group_support, , drop = FALSE]

taxon_levels <- tapply(
    taxon_map$taxonomic_level, taxon_map$taxon,
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
    positive_samples_reference = as.integer(rowSums(
        taxon_counts[, reference_samples, drop = FALSE] > 0
    )),
    positive_samples_comparison = as.integer(rowSums(
        taxon_counts[, comparison_samples, drop = FALSE] > 0
    )),
    prevalence = rowMeans(taxon_counts > 0),
    mean_relative_abundance = rowMeans(relative_abundance),
    stringsAsFactors = FALSE,
    row.names = NULL
)

group_counts <- table(analysis_metadata$group)
preprocessing <- list(
    dataset_id = "microbiomehd_zupancic_obesity",
    dataset_alias = "Zupancic Old Order Amish obesity",
    source_publication = paste(
        "Zupancic et al. (2012), Analysis of the gut microbiota in the Old",
        "Order Amish and its relation to the metabolic syndrome"
    ),
    source_publication_doi = "10.1371/journal.pone.0043052",
    data_source = "MicrobiomeHD v3, Zenodo record 1146764",
    primary_comparison = "obesity versus healthy category at visit 1",
    count_source = "RDP-assigned 100% de novo MicrobiomeHD feature table",
    microbiomehd_count_filter = paste(
        "in the complete cohort: sample reads >100; feature total reads >=10;",
        "feature prevalence >1%; sample reads rechecked after feature filtering"
    ),
    baseline_rule = "visit_number equals 1",
    subject_identifier = "submitted_subject_id_s",
    repeated_sample_rule = paste(
        "one visit-1 sample per subject; retain the sample with the largest",
        "original library size, with sample identifier breaking ties; the",
        "rule does not use disease group or microbiome association results"
    ),
    duplicate_subject_samples_excluded = duplicate_sample_ids,
    current_count_processing = paste(
        "aggregate to genus or deepest resolved higher rank; retain taxa",
        "present in at least 5% of the final cohort and in both groups"
    ),
    rarefaction = "not applied",
    reference_group = "H",
    comparison_group = "OB",
    planned_adjustment_variables = "sex",
    sample_counts = as.list(group_counts),
    baseline_group_counts_before_microbiomehd_filter =
        as.list(baseline_group_counts_before_filter),
    baseline_group_counts_after_initial_depth_filter =
        as.list(baseline_group_counts_after_initial_depth),
    source_features_before_microbiomehd_filter = source_feature_count,
    source_features_after_microbiomehd_filter = nrow(counts),
    retained_taxa = nrow(taxon_counts),
    minimum_positive_samples = minimum_positive_samples,
    library_size_definition = "full count-table total before taxon filtering"
)

output_directory <- file.path(dataset_directory, "processed")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_input <- list(
    counts = taxon_counts,
    metadata = analysis_metadata,
    taxon_metadata = taxon_metadata,
    preprocessing = preprocessing
)
saveRDS(
    analysis_input,
    file.path(output_directory,
              "microbiomehd_zupancic_obesity_sacat_input.rds"),
    compress = "xz"
)
write_count_csv(
    taxon_counts,
    file.path(output_directory,
              "microbiomehd_zupancic_obesity_taxon_counts.csv")
)
utils::write.csv(
    analysis_metadata,
    file.path(output_directory,
              "microbiomehd_zupancic_obesity_metadata.csv"),
    row.names = FALSE, quote = TRUE
)
utils::write.csv(
    taxon_metadata,
    file.path(output_directory,
              "microbiomehd_zupancic_obesity_taxon_metadata.csv"),
    row.names = FALSE, quote = TRUE
)
utils::write.csv(
    data.frame(
        item = c(
            "source DOI", "visit-1 healthy before MicrobiomeHD filter",
            "visit-1 obesity before MicrobiomeHD filter",
            "visit-1 healthy after initial depth filter",
            "visit-1 obesity after initial depth filter",
            "stable subject identifier", "repeated-sample rule",
            "duplicate subject samples excluded",
            "final healthy samples", "final obesity samples",
            "retained taxa", "minimum positive samples", "rarefaction"
        ),
        value = c(preprocessing$source_publication_doi,
                  baseline_group_counts_before_filter[["H"]],
                  baseline_group_counts_before_filter[["OB"]],
                  baseline_group_counts_after_initial_depth[["H"]],
                  baseline_group_counts_after_initial_depth[["OB"]],
                  preprocessing$subject_identifier,
                  preprocessing$repeated_sample_rule,
                  paste(preprocessing$duplicate_subject_samples_excluded,
                        collapse = ";"),
                  group_counts[["H"]], group_counts[["OB"]],
                  nrow(taxon_counts), minimum_positive_samples,
                  preprocessing$rarefaction)
    ),
    file.path(output_directory,
              "microbiomehd_zupancic_obesity_preprocessing_summary.csv"),
    row.names = FALSE, quote = TRUE
)

message(sprintf(
    "Zupancic preparation complete: %d H, %d OB, %d retained taxa.",
    group_counts[["H"]], group_counts[["OB"]], nrow(taxon_counts)
))
