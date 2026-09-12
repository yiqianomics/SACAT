# Prepare the Baxter colorectal cancer cohort for SACAT
#
# Source publication:
# Baxter NT, Ruffin MT IV, Rogers MAM, Schloss PD. Genome Medicine (2016).
# https://doi.org/10.1186/s13073-016-0290-3
#
# Source data:
# MicrobiomeHD v3, Zenodo record 1146764, crc_baxter_results.tar.gz.
# The archive contains the RDP-assigned 100% de novo feature table produced by
# the MicrobiomeHD pipeline. The source publication analyzed a 97% OTU table.
#
# The publication rarefied samples to 10,000 reads and retained OTUs observed
# in at least 5% of samples. This analysis retains the complete MicrobiomeHD
# counts and original pre-filtering library sizes. The standard MicrobiomeHD
# count filters are applied to the full 490-sample cohort before the primary
# comparison is selected. The publication's 5% feature-prevalence rule is then
# evaluated in the CRC-versus-no-lesion cohort.

options(stringsAsFactors = FALSE)

locate_dataset_directory <- function() {
    file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                          value = TRUE)
    candidates <- character()
    if (length(file_argument) == 1L) {
        script_path <- sub("^--file=", "", file_argument)
        candidates <- c(candidates, dirname(normalizePath(script_path)))
    }
    candidates <- c(
        candidates,
        getwd(),
        file.path(getwd(), "Analysis", "RealDataAnalysis", "crc_baxter")
    )
    candidates <- unique(normalizePath(candidates, mustWork = FALSE))
    expected_file <- file.path(
        "raw", "crc_baxter_results", "RDP",
        "crc_baxter.otu_table.100.denovo.rdp_assigned"
    )
    matches <- candidates[file.exists(file.path(candidates, expected_file))]
    if (length(matches) != 1L) {
        stop(
            "Could not identify the crc_baxter data directory. Run this script ",
            "with Rscript or from the SACAT repository root.",
            call. = FALSE
        )
    }
    matches[[1L]]
}

read_count_table <- function(path) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required to read the count table.",
             call. = FALSE)
    }
    table <- data.table::fread(
        path,
        sep = "\t",
        header = TRUE,
        data.table = FALSE,
        check.names = FALSE,
        showProgress = interactive()
    )
    if (ncol(table) < 3L || nrow(table) < 2L) {
        stop("The RDP-assigned count table has an unexpected shape.",
             call. = FALSE)
    }
    feature_ids <- as.character(table[[1L]])
    table[[1L]] <- NULL
    if (anyNA(feature_ids) || any(!nzchar(feature_ids)) ||
        anyDuplicated(feature_ids)) {
        stop("Feature identifiers must be unique and non-empty.",
             call. = FALSE)
    }
    if (is.null(names(table)) || anyNA(names(table)) ||
        any(!nzchar(names(table))) || anyDuplicated(names(table))) {
        stop("Sample identifiers in the count table must be unique and non-empty.",
             call. = FALSE)
    }
    if (!all(vapply(table, is.numeric, logical(1)))) {
        stop("All count columns must be numeric.", call. = FALSE)
    }
    counts <- as.matrix(table)
    if (any(!is.finite(counts)) || any(counts < 0) ||
        any(counts != round(counts))) {
        stop("The count table must contain non-negative integer counts.",
             call. = FALSE)
    }
    storage.mode(counts) <- "integer"
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
        if (!nzchar(value)) NA_character_ else value
    }, character(1))

    if (!is.na(values[["genus"]])) {
        return(c(taxon = values[["genus"]], taxonomic_level = "genus"))
    }

    fallback_levels <- c("family", "order", "class", "phylum", "kingdom")
    resolved_level <- fallback_levels[!is.na(values[fallback_levels])][1L]
    if (is.na(resolved_level)) {
        return(c(
            taxon = "Unclassified_Bacteria",
            taxonomic_level = "unclassified"
        ))
    }
    c(
        taxon = paste0("Unclassified_", values[[resolved_level]]),
        taxonomic_level = resolved_level
    )
}

standardize <- function(x, variable_name) {
    if (any(!is.finite(x))) {
        stop(sprintf("%s contains missing or non-finite values.", variable_name),
             call. = FALSE)
    }
    scale_value <- stats::sd(x)
    if (!is.finite(scale_value) || scale_value == 0) {
        stop(sprintf("%s has no usable variation.", variable_name),
             call. = FALSE)
    }
    (x - mean(x)) / scale_value
}

write_count_csv <- function(counts, path) {
    output <- data.frame(
        taxon = rownames(counts),
        counts,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    utils::write.csv(output, path, row.names = FALSE, quote = TRUE)
}

dataset_directory <- locate_dataset_directory()
raw_directory <- file.path(dataset_directory, "raw", "crc_baxter_results")
output_directory <- file.path(dataset_directory, "processed")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(
    raw_directory, "RDP", "crc_baxter.otu_table.100.denovo.rdp_assigned"
)
metadata_file <- file.path(raw_directory, "crc_baxter.metadata.txt")
if (!file.exists(metadata_file)) {
    stop("The Baxter metadata file is missing.", call. = FALSE)
}

message("Reading the Baxter RDP-assigned feature table")
all_counts <- read_count_table(count_file)
all_library_sizes <- colSums(all_counts)

metadata_raw <- data.table::fread(
    metadata_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    check.names = FALSE,
    encoding = "Latin-1",
    showProgress = FALSE
)
required_metadata <- c(
    "Sample_Name_s", "DiseaseState", "diagnosis_s", "Age_s", "Gender_s",
    "BMI_s"
)
missing_metadata <- setdiff(required_metadata, names(metadata_raw))
if (length(missing_metadata)) {
    stop(
        sprintf(
            "The Baxter metadata file is missing required columns: %s",
            paste(missing_metadata, collapse = ", ")
        ),
        call. = FALSE
    )
}

metadata_raw$sample_id <- as.character(metadata_raw$Sample_Name_s)
if (anyNA(metadata_raw$sample_id) || any(!nzchar(metadata_raw$sample_id)) ||
    anyDuplicated(metadata_raw$sample_id)) {
    stop("Baxter metadata sample identifiers must be unique and non-empty.",
         call. = FALSE)
}
rownames(metadata_raw) <- metadata_raw$sample_id

all_sample_ids <- union(colnames(all_counts), metadata_raw$sample_id)
metadata_index <- match(all_sample_ids, metadata_raw$sample_id)
count_index <- match(all_sample_ids, colnames(all_counts))
sample_disposition <- data.frame(
    sample_id = all_sample_ids,
    in_count_table = !is.na(count_index),
    in_metadata = !is.na(metadata_index),
    disease_state = metadata_raw$DiseaseState[metadata_index],
    diagnosis = metadata_raw$diagnosis_s[metadata_index],
    library_size = as.numeric(all_library_sizes[count_index]),
    stringsAsFactors = FALSE
)

common_samples <- intersect(colnames(all_counts), metadata_raw$sample_id)
if (!length(common_samples)) {
    stop("No sample identifiers are shared by the count table and metadata.",
         call. = FALSE)
}
counts <- all_counts[, common_samples, drop = FALSE]
metadata <- metadata_raw[common_samples, , drop = FALSE]
metadata$library_size <- as.numeric(all_library_sizes[common_samples])

microbiomehd_minimum_sample_reads <- 100L
microbiomehd_minimum_feature_reads <- 10L
microbiomehd_minimum_feature_prevalence <- 0.01

sample_keep_microbiomehd <- colSums(counts) >
    microbiomehd_minimum_sample_reads
counts <- counts[, sample_keep_microbiomehd, drop = FALSE]
metadata <- metadata[colnames(counts), , drop = FALSE]

source_feature_count_before_microbiomehd <- nrow(counts)
feature_keep_microbiomehd <-
    rowSums(counts) >= microbiomehd_minimum_feature_reads &
    rowMeans(counts > 0) > microbiomehd_minimum_feature_prevalence
counts <- counts[feature_keep_microbiomehd, , drop = FALSE]

sample_keep_after_feature_filter <- colSums(counts) >
    microbiomehd_minimum_sample_reads
counts <- counts[, sample_keep_after_feature_filter, drop = FALSE]
metadata <- metadata[colnames(counts), , drop = FALSE]

if (ncol(counts) != 490L || nrow(counts) != 18448L) {
    stop(
        sprintf(
            paste(
                "Unexpected MicrobiomeHD cleaning result:",
                "%d samples and %d features."
            ),
            ncol(counts), nrow(counts)
        ),
        call. = FALSE
    )
}

metadata$group_clean <- ifelse(
    metadata$DiseaseState %in% c("H", "CRC"),
    metadata$DiseaseState,
    NA_character_
)
metadata$age_clean <- suppressWarnings(as.numeric(metadata$Age_s))
metadata$bmi_clean <- suppressWarnings(as.numeric(metadata$BMI_s))
sex_value <- tolower(trimws(as.character(metadata$Gender_s)))
unexpected_sex <- setdiff(unique(sex_value), c("female", "male"))
if (length(unexpected_sex)) {
    stop(
        sprintf(
            "Unexpected values in Gender_s: %s",
            paste(unexpected_sex, collapse = ", ")
        ),
        call. = FALSE
    )
}
metadata$sex_clean <- ifelse(
    sex_value == "female", "F",
    ifelse(sex_value == "male", "M", NA_character_)
)

positive_depth <- is.finite(metadata$library_size) & metadata$library_size > 0
primary_group <- metadata$group_clean %in% c("H", "CRC")
complete_adjustment <- is.finite(metadata$age_clean) &
    metadata$sex_clean %in% c("F", "M")
analysis_keep <- positive_depth & primary_group & complete_adjustment

analysis_metadata_source <- metadata[analysis_keep, , drop = FALSE]
analysis_metadata_source <- analysis_metadata_source[
    order(
        factor(analysis_metadata_source$group_clean,
               levels = c("H", "CRC")),
        analysis_metadata_source$sample_id
    ),
    ,
    drop = FALSE
]

analysis_metadata <- data.frame(
    sample_id = analysis_metadata_source$sample_id,
    group = factor(
        analysis_metadata_source$group_clean,
        levels = c("H", "CRC")
    ),
    library_size = as.numeric(analysis_metadata_source$library_size),
    age = as.numeric(analysis_metadata_source$age_clean),
    age_z = standardize(
        as.numeric(analysis_metadata_source$age_clean), "age"
    ),
    sex = factor(
        analysis_metadata_source$sex_clean,
        levels = c("F", "M")
    ),
    bmi = as.numeric(analysis_metadata_source$bmi_clean),
    diagnosis = analysis_metadata_source$diagnosis_s,
    disease_state = analysis_metadata_source$DiseaseState,
    stringsAsFactors = FALSE,
    check.names = FALSE
)
original_columns <- setdiff(
    names(analysis_metadata_source),
    c(
        names(analysis_metadata), "group_clean", "age_clean", "bmi_clean",
        "sex_clean"
    )
)
analysis_metadata <- cbind(
    analysis_metadata,
    analysis_metadata_source[, original_columns, drop = FALSE]
)
rownames(analysis_metadata) <- analysis_metadata$sample_id

expected_group_counts <- c(H = 172L, CRC = 120L)
observed_group_counts <- table(analysis_metadata$group)
if (!identical(as.integer(observed_group_counts),
               as.integer(expected_group_counts))) {
    stop(
        sprintf(
            "Unexpected final group counts: %s",
            paste(
                paste(names(observed_group_counts), observed_group_counts,
                      sep = "="),
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

primary_counts <- counts[, analysis_metadata$sample_id, drop = FALSE]
minimum_positive_samples <- ceiling(0.05 * ncol(primary_counts))
feature_keep <- rowSums(primary_counts > 0) >= minimum_positive_samples
filtered_feature_counts <- primary_counts[feature_keep, , drop = FALSE]

taxon_map <- t(vapply(
    rownames(filtered_feature_counts),
    parse_rdp_taxon,
    character(2)
))
taxon_map <- data.frame(
    feature_id = rownames(filtered_feature_counts),
    taxon = taxon_map[, "taxon"],
    taxonomic_level = taxon_map[, "taxonomic_level"],
    stringsAsFactors = FALSE,
    row.names = NULL
)

taxon_levels <- tapply(
    taxon_map$taxonomic_level,
    taxon_map$taxon,
    function(x) paste(sort(unique(x)), collapse = ";")
)
if (any(grepl(";", taxon_levels, fixed = TRUE))) {
    stop("At least one taxon label maps to more than one taxonomic level.",
         call. = FALSE)
}

taxon_counts <- rowsum(
    filtered_feature_counts,
    group = taxon_map$taxon,
    reorder = TRUE
)
storage.mode(taxon_counts) <- "integer"
taxon_keep <- rowSums(taxon_counts > 0) >= minimum_positive_samples
taxon_counts <- taxon_counts[taxon_keep, , drop = FALSE]

reference_samples <- analysis_metadata$group == "H"
comparison_samples <- analysis_metadata$group == "CRC"
positive_samples_H <- rowSums(
    taxon_counts[, reference_samples, drop = FALSE] > 0
)
positive_samples_CRC <- rowSums(
    taxon_counts[, comparison_samples, drop = FALSE] > 0
)
two_group_support <- positive_samples_H > 0L & positive_samples_CRC > 0L
group_support_exclusions <- data.frame(
    taxon = rownames(taxon_counts)[!two_group_support],
    positive_samples_H = as.integer(positive_samples_H[!two_group_support]),
    positive_samples_CRC = as.integer(
        positive_samples_CRC[!two_group_support]
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
)
taxon_bins_before_group_support_filter <- nrow(taxon_counts)
taxon_counts <- taxon_counts[two_group_support, , drop = FALSE]

if (!identical(colnames(taxon_counts), rownames(analysis_metadata))) {
    stop("Counts and metadata are not in the same sample order.",
         call. = FALSE)
}
if (any(analysis_metadata$library_size < colSums(taxon_counts))) {
    stop("A filtered taxon total exceeds its original sample library size.",
         call. = FALSE)
}

relative_abundance <- sweep(
    taxon_counts,
    2L,
    analysis_metadata$library_size,
    "/"
)
source_features_by_taxon <- table(taxon_map$taxon)
taxon_metadata <- data.frame(
    taxon = rownames(taxon_counts),
    taxonomic_level = unname(taxon_levels[rownames(taxon_counts)]),
    source_feature_count = as.integer(
        source_features_by_taxon[rownames(taxon_counts)]
    ),
    total_count = as.numeric(rowSums(taxon_counts)),
    positive_samples = as.integer(rowSums(taxon_counts > 0)),
    positive_samples_H = as.integer(
        rowSums(taxon_counts[, reference_samples, drop = FALSE] > 0)
    ),
    positive_samples_CRC = as.integer(
        rowSums(taxon_counts[, comparison_samples, drop = FALSE] > 0)
    ),
    prevalence = rowMeans(taxon_counts > 0),
    prevalence_H = rowMeans(taxon_counts[, reference_samples, drop = FALSE] > 0),
    prevalence_CRC = rowMeans(taxon_counts[, comparison_samples, drop = FALSE] > 0),
    mean_relative_abundance = rowMeans(relative_abundance),
    mean_relative_abundance_H = rowMeans(
        relative_abundance[, reference_samples, drop = FALSE]
    ),
    mean_relative_abundance_CRC = rowMeans(
        relative_abundance[, comparison_samples, drop = FALSE]
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
)

if (nrow(taxon_metadata) != 126L) {
    stop(
        sprintf(
            "Expected 126 retained taxon bins, but obtained %d.",
            nrow(taxon_metadata)
        ),
        call. = FALSE
    )
}

sample_disposition$analysis_status <- "not retained"
sample_disposition$exclusion_reason <- ifelse(
    !sample_disposition$in_count_table,
    "not represented in count table",
    ifelse(
        !sample_disposition$in_metadata,
        "no matching metadata",
        ifelse(
            !is.finite(sample_disposition$library_size) |
                sample_disposition$library_size <= 0,
            "non-positive library size",
            ifelse(
                !(sample_disposition$disease_state %in% c("H", "CRC")),
                "not in primary CRC-versus-H comparison",
                "missing required age or sex covariate"
            )
        )
    )
)
retained_index <- match(analysis_metadata$sample_id,
                        sample_disposition$sample_id)
sample_disposition$analysis_status[retained_index] <- "retained"
sample_disposition$exclusion_reason[retained_index] <- ""

preprocessing <- list(
    dataset_id = "crc_baxter",
    dataset_alias = "BaxterCRC16S",
    source_publication = paste(
        "Baxter et al. (2016), Microbiota-based model improves the",
        "sensitivity of fecal immunochemical test for detecting colonic lesions"
    ),
    source_publication_doi = "10.1186/s13073-016-0290-3",
    data_source = "MicrobiomeHD v3, Zenodo record 1146764",
    primary_comparison = "CRC versus no-lesion healthy control",
    excluded_group = "adenoma (MicrobiomeHD DiseaseState nonCRC)",
    count_source = "RDP-assigned 100% de novo MicrobiomeHD feature table",
    publication_sequence_processing = paste(
        "mothur 1.36; 97% average-neighbor OTUs; RDP taxonomy;",
        "10,000-read rarefaction; OTUs present in at least 5% of samples"
    ),
    microbiomehd_count_filter = paste(
        "in the complete cohort: sample reads >100; feature total reads >=10;",
        "feature prevalence >1%; sample reads rechecked after feature filtering"
    ),
    current_count_processing = paste(
        "MicrobiomeHD count filter followed by unrarefied upstream features",
        "present in at least 5% of the primary cohort; aggregation to genus",
        "or deepest resolved higher rank; taxon bins required to have a",
        "positive count in both comparison groups"
    ),
    prevalence_filter_scope = "CRC-versus-H primary cohort",
    minimum_positive_samples = minimum_positive_samples,
    group_support_filter = paste(
        "at least one positive-count sample in each comparison group;",
        "required for the two-arm conditional-abundance estimand"
    ),
    group_support_exclusions = group_support_exclusions,
    sample_depth_rule = "positive count-table library size",
    rarefaction = "not applied",
    reference_group = "H",
    comparison_group = "CRC",
    planned_adjustment_variables = c("age_z", "sex"),
    sacat_input_orientation = "taxa by samples",
    sacat_group_column = "group",
    sacat_library_size_column = "library_size",
    sample_counts = as.list(observed_group_counts),
    feature_counts = list(
        source_features_before_microbiomehd_filter =
            source_feature_count_before_microbiomehd,
        source_features_after_microbiomehd_filter = nrow(counts),
        features_after_primary_prevalence_filter =
            nrow(filtered_feature_counts),
        taxon_bins_before_group_support_filter =
            taxon_bins_before_group_support_filter,
        taxa_excluded_without_two_group_positive_support =
            nrow(group_support_exclusions),
        retained_taxon_bins = nrow(taxon_counts)
    )
)

summary_table <- data.frame(
    item = c(
        "dataset", "source DOI", "primary comparison", "excluded group",
        "source count table", "original publication processing",
        "MicrobiomeHD count filter",
        "source features before MicrobiomeHD filter",
        "features after MicrobiomeHD filter",
        "current sample-depth rule", "current rarefaction",
        "current feature-prevalence rule", "taxonomic aggregation",
        "two-group positive-count support rule",
        "taxa excluded without two-group positive-count support",
        "excluded taxa and positive samples by group",
        "healthy samples", "CRC samples", "retained upstream features",
        "retained taxon bins", "SACAT count orientation",
        "SACAT library-size definition"
    ),
    value = c(
        "Baxter colorectal cancer 16S stool cohort",
        preprocessing$source_publication_doi,
        preprocessing$primary_comparison,
        preprocessing$excluded_group,
        preprocessing$count_source,
        preprocessing$publication_sequence_processing,
        preprocessing$microbiomehd_count_filter,
        as.character(source_feature_count_before_microbiomehd),
        as.character(nrow(counts)),
        preprocessing$sample_depth_rule,
        preprocessing$rarefaction,
        sprintf(
            "present in at least %d of %d primary samples (5%%)",
            minimum_positive_samples, ncol(primary_counts)
        ),
        "genus when resolved; otherwise deepest resolved higher-rank bin",
        preprocessing$group_support_filter,
        as.character(nrow(group_support_exclusions)),
        if (nrow(group_support_exclusions)) {
            paste(
                sprintf(
                    "%s (H=%d, CRC=%d)",
                    group_support_exclusions$taxon,
                    group_support_exclusions$positive_samples_H,
                    group_support_exclusions$positive_samples_CRC
                ),
                collapse = "; "
            )
        } else {
            "none"
        },
        as.character(observed_group_counts[["H"]]),
        as.character(observed_group_counts[["CRC"]]),
        as.character(nrow(filtered_feature_counts)),
        as.character(nrow(taxon_counts)),
        preprocessing$sacat_input_orientation,
        "sum of the full MicrobiomeHD feature table before taxon filtering"
    ),
    stringsAsFactors = FALSE
)

analysis_input <- list(
    counts = taxon_counts,
    metadata = analysis_metadata,
    taxon_metadata = taxon_metadata,
    preprocessing = preprocessing
)

saveRDS(
    analysis_input,
    file.path(output_directory, "crc_baxter_sacat_input.rds"),
    compress = "xz"
)
write_count_csv(
    taxon_counts,
    file.path(output_directory, "crc_baxter_taxon_counts.csv")
)
utils::write.csv(
    analysis_metadata,
    file.path(output_directory, "crc_baxter_metadata.csv"),
    row.names = FALSE,
    quote = TRUE
)
utils::write.csv(
    taxon_metadata,
    file.path(output_directory, "crc_baxter_taxon_metadata.csv"),
    row.names = FALSE,
    quote = TRUE
)
utils::write.csv(
    sample_disposition,
    file.path(output_directory, "crc_baxter_sample_disposition.csv"),
    row.names = FALSE,
    quote = TRUE,
    na = ""
)
utils::write.csv(
    summary_table,
    file.path(output_directory, "crc_baxter_preprocessing_summary.csv"),
    row.names = FALSE,
    quote = TRUE
)

message(sprintf(
    "Baxter preparation complete: %d H, %d CRC, %d retained taxon bins.",
    observed_group_counts[["H"]], observed_group_counts[["CRC"]],
    nrow(taxon_counts)
))
