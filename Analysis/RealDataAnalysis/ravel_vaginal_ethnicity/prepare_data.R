# Prepare the Ravel cross-sectional vaginal microbiome cohort
#
# The public table contains one vaginal sample from each participant.  This
# analysis compares Black and White participants, the two largest groups for
# the prespecified binary comparison.  The published, unrarefied taxon counts
# and their exact column totals are retained.

options(stringsAsFactors = FALSE)

locate_dataset_directory <- function() {
    file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                          value = TRUE)
    if (length(file_argument) == 1L) {
        return(dirname(normalizePath(sub("^--file=", "", file_argument))))
    }
    file.path(getwd(), "Analysis", "RealDataAnalysis",
              "ravel_vaginal_ethnicity")
}

write_count_csv <- function(counts, path) {
    utils::write.csv(
        data.frame(taxon = rownames(counts), counts, check.names = FALSE),
        path, row.names = FALSE, quote = TRUE
    )
}

dataset_directory <- locate_dataset_directory()
source(file.path(
    dirname(dataset_directory), "shared", "count_validation.R"
))
raw_directory <- file.path(dataset_directory, "raw")
output_directory <- file.path(dataset_directory, "processed")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

message("Reading the Ravel count and participant tables")
counts_all <- as_count_matrix(read.table(
    file.path(raw_directory, "counts.tsv"), header = TRUE,
    check.names = FALSE
), "The published count table")
validate_identifiers(rownames(counts_all), "Taxon identifiers")
validate_identifiers(colnames(counts_all), "Sample identifiers")
metadata_all <- read.table(
    file.path(raw_directory, "metadata.tsv"), header = TRUE,
    check.names = FALSE, stringsAsFactors = FALSE
)
metadata_all$sample_id <- rownames(metadata_all)

sample_ids <- intersect(colnames(counts_all), metadata_all$sample_id)
metadata_all <- metadata_all[match(sample_ids, metadata_all$sample_id), ,
                             drop = FALSE]
counts_all <- counts_all[, sample_ids, drop = FALSE]

eligible <- metadata_all$Ethnic_Group %in% c("White", "Black")
selected_ids <- metadata_all$sample_id[eligible]
selected_order <- order(
    factor(metadata_all$Ethnic_Group[eligible], levels = c("White", "Black")),
    selected_ids
)
selected_ids <- selected_ids[selected_order]
selected_metadata <- metadata_all[match(selected_ids, metadata_all$sample_id), ,
                                  drop = FALSE]

analysis_metadata <- data.frame(
    sample_id = selected_ids,
    group = factor(selected_metadata$Ethnic_Group,
                   levels = c("White", "Black")),
    library_size = as.numeric(selected_metadata$Depth),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
rownames(analysis_metadata) <- analysis_metadata$sample_id

taxon_counts <- counts_all[, selected_ids, drop = FALSE]
taxon_counts <- rowsum(taxon_counts, rownames(taxon_counts), reorder = TRUE)
storage.mode(taxon_counts) <- "integer"

if (!identical(as.numeric(colSums(taxon_counts)),
               analysis_metadata$library_size)) {
    stop("Published library sizes do not match the taxon count totals.",
         call. = FALSE)
}

minimum_positive_samples <- ceiling(0.05 * ncol(taxon_counts))
prevalence_keep <- rowSums(taxon_counts > 0) >= minimum_positive_samples
taxon_counts <- taxon_counts[prevalence_keep, , drop = FALSE]
reference_samples <- analysis_metadata$group == "White"
comparison_samples <- analysis_metadata$group == "Black"
group_support <-
    rowSums(taxon_counts[, reference_samples, drop = FALSE] > 0) > 0 &
    rowSums(taxon_counts[, comparison_samples, drop = FALSE] > 0) > 0
taxon_counts <- taxon_counts[group_support, , drop = FALSE]

relative_abundance <- sweep(
    taxon_counts, 2L, analysis_metadata$library_size, "/"
)
taxon_metadata <- data.frame(
    taxon = rownames(taxon_counts),
    taxonomic_level = "published species-level label",
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
    dataset_id = "ravel_vaginal_ethnicity",
    dataset_alias = "Ravel cross-sectional vaginal microbiome",
    source_publication = paste(
        "Ravel et al. (2011), Vaginal microbiome of reproductive-age",
        "women"
    ),
    source_publication_doi = "10.1073/pnas.1002611107",
    data_source = paste(
        "Published supplementary count table, mirrored by the ETH Zurich",
        "Microbiome Data Analysis workshop"
    ),
    primary_comparison = "Black versus White participants",
    count_source = "published unrarefied taxon count table",
    rarefaction = "not applied",
    reference_group = "White",
    comparison_group = "Black",
    planned_adjustment_variables = character(),
    sample_counts = as.list(group_counts),
    retained_taxa = nrow(taxon_counts),
    minimum_positive_samples = minimum_positive_samples,
    library_size_definition = "published count-table column total"
)

analysis_input <- list(
    counts = taxon_counts,
    metadata = analysis_metadata,
    taxon_metadata = taxon_metadata,
    preprocessing = preprocessing
)
saveRDS(
    analysis_input,
    file.path(output_directory,
              "ravel_vaginal_ethnicity_sacat_input.rds"),
    compress = "xz"
)
write_count_csv(
    taxon_counts,
    file.path(output_directory, "ravel_vaginal_ethnicity_taxon_counts.csv")
)
utils::write.csv(
    analysis_metadata,
    file.path(output_directory, "ravel_vaginal_ethnicity_metadata.csv"),
    row.names = FALSE, quote = TRUE
)
utils::write.csv(
    taxon_metadata,
    file.path(output_directory,
              "ravel_vaginal_ethnicity_taxon_metadata.csv"),
    row.names = FALSE, quote = TRUE
)
utils::write.csv(
    data.frame(
        item = c(
            "source DOI", "White participants", "Black participants",
            "retained taxa", "prevalence rule", "rarefaction",
            "library size"
        ),
        value = c(
            preprocessing$source_publication_doi,
            group_counts[["White"]], group_counts[["Black"]],
            nrow(taxon_counts),
            sprintf("positive in at least %d participants and in both groups",
                    minimum_positive_samples),
            preprocessing$rarefaction,
            preprocessing$library_size_definition
        ),
        stringsAsFactors = FALSE
    ),
    file.path(output_directory,
              "ravel_vaginal_ethnicity_preprocessing_summary.csv"),
    row.names = FALSE, quote = TRUE
)

message(sprintf(
    "Ravel preparation complete: %d White, %d Black, %d retained taxa.",
    group_counts[["White"]], group_counts[["Black"]], nrow(taxon_counts)
))
