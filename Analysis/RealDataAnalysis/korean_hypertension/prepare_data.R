# Prepare the Korean hypertension gut microbiome cohort
#
# The public workbook contains sample-level ASV counts, taxonomic assignments,
# and clinical metadata for 120 adults with hypertension and 503 adults with
# normotension.  Systolic and diastolic blood pressure define the study group
# and are therefore not used as adjustment variables.

options(stringsAsFactors = FALSE)

locate_dataset_directory <- function() {
    file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                          value = TRUE)
    if (length(file_argument) == 1L) {
        return(dirname(normalizePath(sub("^--file=", "", file_argument))))
    }
    file.path(getwd(), "Analysis", "RealDataAnalysis",
              "korean_hypertension")
}

standardize <- function(value) {
    as.numeric((value - mean(value)) / stats::sd(value))
}

write_count_csv <- function(counts, path) {
    utils::write.csv(
        data.frame(taxon = rownames(counts), counts, check.names = FALSE),
        path, row.names = FALSE, quote = TRUE
    )
}

dataset_directory <- locate_dataset_directory()
source(file.path(
    dirname(dataset_directory), "shared", "prepare_qiita_dataset.R"
))
source(file.path(
    dirname(dataset_directory), "shared", "count_validation.R"
))
source_file <- file.path(dataset_directory, "raw", "Table_1.xlsx")
output_directory <- file.path(dataset_directory, "processed")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file)) {
    stop("Place Table_1.xlsx in the dataset raw directory.", call. = FALSE)
}
if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required to prepare this dataset.",
         call. = FALSE)
}

message("Reading the Korean hypertension public workbook")
metadata_raw <- as.data.frame(readxl::read_excel(
    source_file, sheet = "Metadata"
), check.names = FALSE)
count_table <- as.data.frame(readxl::read_excel(
    source_file, sheet = "ASV table"
), check.names = FALSE)
taxonomy_raw <- as.data.frame(readxl::read_excel(
    source_file, sheet = "Taxonomy"
), check.names = FALSE)

sample_ids <- validate_identifiers(count_table[[1L]], "Sample identifiers")
feature_ids <- validate_identifiers(
    names(count_table)[-1L], "Feature identifiers"
)
counts_by_sample <- as_count_matrix(
    count_table[, -1L, drop = FALSE], "The ASV count table"
)
rownames(counts_by_sample) <- sample_ids
colnames(counts_by_sample) <- feature_ids
all_counts <- t(counts_by_sample)
all_library_sizes <- colSums(all_counts)

metadata_raw <- metadata_raw[match(sample_ids, metadata_raw$Sample), ,
                             drop = FALSE]
taxonomy_raw <- taxonomy_raw[
    match(feature_ids, taxonomy_raw$`Feature ID`), , drop = FALSE
]
if (anyNA(metadata_raw$Sample) || anyNA(taxonomy_raw$`Feature ID`)) {
    stop("Sample or feature identifiers do not align across the workbook.",
         call. = FALSE)
}

group_value <- as.character(metadata_raw$group)
age <- as.numeric(metadata_raw$Age)
bmi <- as.numeric(metadata_raw$BMI)
sex_value <- as.character(metadata_raw$Sex)
eligible <- group_value %in% c("Normotension", "Hypertension") &
    is.finite(age) & is.finite(bmi) &
    sex_value %in% c("Female", "Male") & all_library_sizes > 0

selected_ids <- sample_ids[eligible]
selected_order <- order(
    factor(group_value[eligible],
           levels = c("Normotension", "Hypertension")),
    selected_ids
)
selected_ids <- selected_ids[selected_order]
analysis_metadata <- data.frame(
    sample_id = selected_ids,
    group = factor(group_value[eligible][selected_order],
                   levels = c("Normotension", "Hypertension")),
    library_size = as.numeric(all_library_sizes[selected_ids]),
    age = age[eligible][selected_order],
    age_z = standardize(age[eligible][selected_order]),
    sex = factor(sex_value[eligible][selected_order],
                 levels = c("Female", "Male")),
    bmi = bmi[eligible][selected_order],
    bmi_z = standardize(bmi[eligible][selected_order]),
    stringsAsFactors = FALSE,
    check.names = FALSE
)
rownames(analysis_metadata) <- analysis_metadata$sample_id

taxon_map_matrix <- t(vapply(
    taxonomy_raw$Taxon, taxonomy_bin, character(2)
))
taxon_map <- data.frame(
    feature_id = feature_ids,
    taxon = taxon_map_matrix[, "taxon"],
    taxonomic_level = taxon_map_matrix[, "taxonomic_level"],
    stringsAsFactors = FALSE
)
taxon_counts <- rowsum(
    all_counts[, selected_ids, drop = FALSE], taxon_map$taxon,
    reorder = TRUE
)
storage.mode(taxon_counts) <- "integer"

minimum_positive_samples <- ceiling(0.05 * ncol(taxon_counts))
taxon_counts <- taxon_counts[
    rowSums(taxon_counts > 0) >= minimum_positive_samples, , drop = FALSE
]
reference_samples <- analysis_metadata$group == "Normotension"
comparison_samples <- analysis_metadata$group == "Hypertension"
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
    dataset_id = "korean_hypertension",
    dataset_alias = "Korean adult hypertension gut microbiome cohort",
    source_publication = paste(
        "Song et al. (2023), The association between gut microbiome and",
        "hypertension varies according to enterotypes: a Korean study"
    ),
    source_publication_doi = "10.3389/frmbi.2023.1072059",
    data_source = "Public article supplementary workbook Table 1",
    primary_comparison = "hypertension versus normotension",
    count_source = "published unrarefied ASV count table",
    rarefaction = "not applied",
    reference_group = "Normotension",
    comparison_group = "Hypertension",
    planned_adjustment_variables = c("age_z", "sex", "bmi_z"),
    variables_not_adjusted = c(
        "systolic blood pressure and diastolic blood pressure because they",
        "define the hypertension group"
    ),
    sample_counts = as.list(group_counts),
    retained_taxa = nrow(taxon_counts),
    minimum_positive_samples = minimum_positive_samples,
    library_size_definition =
        "sum of all ASV counts in the published table before taxon filtering"
)

analysis_input <- list(
    counts = taxon_counts,
    metadata = analysis_metadata,
    taxon_metadata = taxon_metadata,
    preprocessing = preprocessing
)
saveRDS(
    analysis_input,
    file.path(output_directory, "korean_hypertension_dasra_input.rds"),
    compress = "xz"
)
write_count_csv(
    taxon_counts,
    file.path(output_directory, "korean_hypertension_taxon_counts.csv")
)
utils::write.csv(
    analysis_metadata,
    file.path(output_directory, "korean_hypertension_metadata.csv"),
    row.names = FALSE, quote = TRUE
)
utils::write.csv(
    taxon_metadata,
    file.path(output_directory, "korean_hypertension_taxon_metadata.csv"),
    row.names = FALSE, quote = TRUE
)
utils::write.csv(
    data.frame(
        item = c(
            "source DOI", "normotension", "hypertension", "retained taxa",
            "prevalence rule", "rarefaction", "library size"
        ),
        value = c(
            preprocessing$source_publication_doi,
            group_counts[["Normotension"]],
            group_counts[["Hypertension"]], nrow(taxon_counts),
            sprintf("positive in at least %d samples and in both groups",
                    minimum_positive_samples),
            preprocessing$rarefaction,
            preprocessing$library_size_definition
        ),
        stringsAsFactors = FALSE
    ),
    file.path(
        output_directory, "korean_hypertension_preprocessing_summary.csv"
    ),
    row.names = FALSE, quote = TRUE
)

message(sprintf(
    paste(
        "Korean hypertension preparation complete: %d normotensive,",
        "%d hypertensive, %d retained taxa."
    ),
    group_counts[["Normotension"]], group_counts[["Hypertension"]],
    nrow(taxon_counts)
))
