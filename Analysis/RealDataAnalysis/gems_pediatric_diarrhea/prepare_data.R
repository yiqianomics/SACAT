# Prepare the GEMS pediatric diarrhea cohort

options(stringsAsFactors = FALSE)

dataset_directory <- dirname(normalizePath(sub(
    "^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]
)))
source_file <- file.path(dataset_directory, "raw", "forserveroptim.rdata")
if (!file.exists(source_file)) {
    stop("Place forserveroptim.rdata in the dataset raw directory.",
         call. = FALSE)
}
if (!requireNamespace("Biobase", quietly = TRUE) ||
    !requireNamespace("metagenomeSeq", quietly = TRUE)) {
    stop("Packages 'Biobase' and 'metagenomeSeq' are required.", call. = FALSE)
}
source(file.path(
    dirname(dataset_directory), "shared", "count_validation.R"
))

standardize <- function(value) {
    as.numeric((value - mean(value)) / stats::sd(value))
}

write_count_csv <- function(counts, path) {
    utils::write.csv(
        data.frame(taxon = rownames(counts), counts, check.names = FALSE),
        path, row.names = FALSE, quote = TRUE
    )
}

source_environment <- new.env(parent = emptyenv())
loaded_objects <- load(source_file, envir = source_environment)
required_objects <- c("gates", "graw", "totalCounts")
if (!all(required_objects %in% loaded_objects)) {
    stop(
        "The GEMS source file lacks gates, graw, or totalCounts.",
        call. = FALSE
    )
}
phenotype <- Biobase::pData(source_environment$gates)
counts <- as_count_matrix(
    source_environment$graw, "The GEMS count table"
)
validate_identifiers(rownames(counts), "Taxon identifiers")
validate_identifiers(colnames(counts), "Sample identifiers")
library_sizes <- as.numeric(source_environment$totalCounts)
names(library_sizes) <- names(source_environment$totalCounts)
validate_identifiers(names(library_sizes), "Library-size sample identifiers")
if (any(!is.finite(library_sizes)) || any(library_sizes <= 0)) {
    stop("GEMS library sizes must be positive and finite.", call. = FALSE)
}

common_samples <- Reduce(intersect, list(
    colnames(counts), rownames(phenotype), names(library_sizes)
))
phenotype <- phenotype[common_samples, , drop = FALSE]
counts <- counts[, common_samples, drop = FALSE]

group_value <- ifelse(
    phenotype$Type == "Control", "Control",
    ifelse(phenotype$Type == "Case", "MSD", NA_character_)
)
age <- suppressWarnings(as.numeric(phenotype$Age))
country <- trimws(as.character(phenotype$Country))
eligible <- group_value %in% c("Control", "MSD") &
    is.finite(age) & nzchar(country) &
    is.finite(library_sizes[common_samples]) &
    library_sizes[common_samples] > 0

sample_ids <- common_samples[eligible]
sample_order <- order(
    factor(group_value[eligible], levels = c("Control", "MSD")),
    sample_ids
)
sample_ids <- sample_ids[sample_order]
analysis_metadata <- data.frame(
    sample_id = sample_ids,
    group = factor(group_value[eligible][sample_order],
                   levels = c("Control", "MSD")),
    library_size = library_sizes[sample_ids],
    age_months = age[eligible][sample_order],
    age_z = standardize(age[eligible][sample_order]),
    country = factor(country[eligible][sample_order],
                     levels = sort(unique(country[eligible]))),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
rownames(analysis_metadata) <- analysis_metadata$sample_id

taxon_counts <- counts[, sample_ids, drop = FALSE]
minimum_positive_samples <- ceiling(0.05 * ncol(taxon_counts))
prevalence_keep <- rowSums(taxon_counts > 0) >= minimum_positive_samples
taxon_counts <- taxon_counts[prevalence_keep, , drop = FALSE]
reference_samples <- analysis_metadata$group == "Control"
comparison_samples <- analysis_metadata$group == "MSD"
group_support <-
    rowSums(taxon_counts[, reference_samples, drop = FALSE] > 0) > 0 &
    rowSums(taxon_counts[, comparison_samples, drop = FALSE] > 0) > 0
taxon_counts <- taxon_counts[group_support, , drop = FALSE]

relative_abundance <- sweep(
    taxon_counts, 2L, analysis_metadata$library_size, "/"
)
taxon_metadata <- data.frame(
    taxon = rownames(taxon_counts),
    taxonomic_level = "genus",
    source_feature_count = 1L,
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
    dataset_id = "gems_pediatric_diarrhea",
    dataset_alias = "GEMS pediatric diarrhea",
    source_publication = paste(
        "Pop et al. (2014), Diarrhea in young children from low-income",
        "countries: a systematic study of microbial etiology using",
        "multiplex sequencing"
    ),
    source_publication_doi = "10.1186/gb-2014-15-6-r76",
    data_source = "Public MSD1000 analysis data accompanying the GEMS study",
    primary_comparison = "moderate-to-severe diarrhea versus matched control",
    count_source = "published unrarefied genus count matrix graw",
    current_count_processing = paste(
        "retain genera present in at least 5% of the final cohort and with",
        "a positive count in both groups"
    ),
    rarefaction = "not applied",
    reference_group = "Control",
    comparison_group = "MSD",
    planned_adjustment_variables = c("age_z", "country"),
    sample_counts = as.list(group_counts),
    retained_taxa = nrow(taxon_counts),
    minimum_positive_samples = minimum_positive_samples,
    library_size_definition = "published totalCounts before genus filtering"
)

output_directory <- file.path(dataset_directory, "processed")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_input <- list(
    counts = taxon_counts,
    metadata = analysis_metadata,
    taxon_metadata = taxon_metadata,
    preprocessing = preprocessing
)
saveRDS(analysis_input,
        file.path(output_directory, "gems_pediatric_diarrhea_sacat_input.rds"),
        compress = "xz")
write_count_csv(
    taxon_counts,
    file.path(output_directory, "gems_pediatric_diarrhea_taxon_counts.csv")
)
utils::write.csv(
    analysis_metadata,
    file.path(output_directory, "gems_pediatric_diarrhea_metadata.csv"),
    row.names = FALSE, quote = TRUE
)
utils::write.csv(
    taxon_metadata,
    file.path(output_directory, "gems_pediatric_diarrhea_taxon_metadata.csv"),
    row.names = FALSE, quote = TRUE
)
utils::write.csv(
    data.frame(
        item = c("source DOI", "control samples", "MSD samples",
                 "retained taxa", "minimum positive samples", "rarefaction"),
        value = c(preprocessing$source_publication_doi,
                  group_counts[["Control"]], group_counts[["MSD"]],
                  nrow(taxon_counts), minimum_positive_samples,
                  preprocessing$rarefaction)
    ),
    file.path(output_directory,
              "gems_pediatric_diarrhea_preprocessing_summary.csv"),
    row.names = FALSE, quote = TRUE
)

message(sprintf(
    "GEMS preparation complete: %d controls, %d MSD cases, %d genera.",
    group_counts[["Control"]], group_counts[["MSD"]], nrow(taxon_counts)
))
