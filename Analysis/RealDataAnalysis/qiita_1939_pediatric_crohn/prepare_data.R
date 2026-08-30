# Prepare the Qiita 1939 pediatric Crohn's disease cohort
#
# The comparison uses rectal mucosal samples from the RISK collection. Children
# with newly diagnosed Crohn's disease are compared with non-IBD controls after
# excluding samples collected during antibiotic exposure. One sample per child
# is retained using the original library size.

options(stringsAsFactors = FALSE)

locate_dataset_directory <- function() {
    file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                          value = TRUE)
    if (length(file_argument) == 1L) {
        return(dirname(normalizePath(sub("^--file=", "", file_argument))))
    }
    file.path(getwd(), "Analysis", "RealDataAnalysis",
              "qiita_1939_pediatric_crohn")
}

dataset_directory <- locate_dataset_directory()
source(file.path(
    dirname(dataset_directory), "shared", "prepare_qiita_dataset.R"
))

prepare_crohn_metadata <- function(metadata) {
    metadata$age_numeric <- as.numeric(metadata$host_age)
    metadata$crohn_status <- ifelse(
        metadata$diagnosis == "no", "Control",
        ifelse(metadata$diagnosis == "CD", "Crohn", NA_character_)
    )
    metadata
}

select_crohn_samples <- function(metadata) {
    keep <- metadata$collection == "RISK" &
        metadata$sample_type_qiita == "rectum mucosa" &
        metadata$crohn_status %in% c("Control", "Crohn") &
        metadata$antibiotics == "FALSE"

    list(
        keep = keep,
        note = paste(
            "RISK collection rectal mucosal samples from children with",
            "Crohn's disease or non-IBD controls, excluding samples recorded",
            "as collected during antibiotic exposure"
        ),
        details = list(
            collection = "RISK",
            specimen = "rectum mucosa",
            diagnoses = c("no", "CD"),
            antibiotics = "FALSE"
        )
    )
}

configuration <- list(
    dataset_id = "qiita_1939_pediatric_crohn",
    dataset_title = "Qiita 1939 pediatric Crohn's disease rectal cohort",
    biom_file = file.path("raw", "qiita_1939.biom"),
    metadata_file = file.path("raw", "qiita_1939.tsv"),
    sample_id_column = "#SampleID",
    prepare_metadata = prepare_crohn_metadata,
    select_samples = select_crohn_samples,
    selection_columns = c(
        "collection", "sample_type_qiita", "diagnosis", "antibiotics",
        "crohn_status"
    ),
    group_column = "crohn_status",
    reference = "Control",
    comparison = "Crohn",
    covariates = c(
        age_z = "age_numeric",
        sex = "sex"
    ),
    numeric_covariates = "age_z",
    subject_column = "host_subject_id",
    one_sample_per_subject = TRUE,
    source_metadata_columns = c(
        "collection", "sample_type_qiita", "biopsy_location",
        "diagnosis", "antibiotics", "host_age", "sex"
    ),
    source_publication = paste(
        "Gevers et al. (2014), The treatment-naive microbiome in",
        "new-onset Crohn's disease"
    ),
    source_publication_doi = "10.1016/j.chom.2014.02.005",
    data_source = paste(
        "Qiita study 1939; European Nucleotide Archive PRJEB13679;",
        "standardized public redbiom count table"
    ),
    count_source = paste(
        "public Qiita BIOM count table with taxonomy fetched through redbiom"
    ),
    primary_comparison = paste(
        "Crohn's disease versus non-IBD control in RISK rectal mucosa,",
        "adjusted for age and sex"
    )
)

prepare_qiita_dataset(dataset_directory, configuration)
