# Prepare the Rauer defined mock-community benchmark.

options(stringsAsFactors = FALSE, warn = 1)

locate_dataset_directory <- function() {
    file_argument <- grep(
        "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
    )
    if (length(file_argument) == 1L) {
        script_path <- sub("^--file=", "", file_argument)
        return(dirname(normalizePath(script_path)))
    }
    candidates <- unique(c(
        getwd(),
        file.path(getwd(), "Analysis", "RealDataAnalysis", "rauer_mock")
    ))
    candidates <- normalizePath(candidates, mustWork = FALSE)
    matches <- candidates[file.exists(file.path(candidates, "prepare_data.R"))]
    if (length(matches) != 1L) {
        stop("Could not identify the rauer_mock directory.", call. = FALSE)
    }
    matches[[1L]]
}

required_packages <- c("dada2", "Biostrings", "stringdist", "ape")
missing_packages <- required_packages[!vapply(
    required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
    stop(
        sprintf(
            "Required R packages are unavailable: %s",
            paste(missing_packages, collapse = ", ")
        ),
        call. = FALSE
    )
}

workers <- suppressWarnings(as.integer(Sys.getenv("SACAT_WORKERS", "1")))
if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("SACAT_WORKERS must be a positive integer.", call. = FALSE)
}

dataset_directory <- locate_dataset_directory()
raw_directory <- file.path(dataset_directory, "raw")
fastq_directory <- file.path(raw_directory, "fastq")
reference_directory <- file.path(raw_directory, "reference")
metadata_file <- file.path(raw_directory, "Metafile.csv")
filtered_directory <- file.path(dataset_directory, "work", "filtered_fastq")
processed_directory <- file.path(dataset_directory, "processed")
dir.create(filtered_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_directory, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(metadata_file)) {
    stop("The source metadata file raw/Metafile.csv is missing.", call. = FALSE)
}

metadata <- utils::read.csv2(metadata_file, fileEncoding = "latin1")
required_metadata <- c(
    "Match_ID", "buffer", "kit", "lysis.prot", "Sample_ID",
    "Broad_type", "Dil"
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata)) {
    stop(
        sprintf(
            "The source metadata is missing required columns: %s",
            paste(missing_metadata, collapse = ", ")
        ),
        call. = FALSE
    )
}
metadata$Sample_ID <- as.character(metadata$Sample_ID)
if (anyNA(metadata$Sample_ID) || any(!nzchar(metadata$Sample_ID)) ||
    anyDuplicated(metadata$Sample_ID)) {
    stop("Sample identifiers must be unique and non-empty.", call. = FALSE)
}

forward_files <- sort(list.files(
    fastq_directory,
    pattern = "_R1_001[.]fastq[.]gz$",
    full.names = TRUE
))
sample_ids <- sub("_.*$", "", basename(forward_files))
if (length(forward_files) != 94L || anyDuplicated(sample_ids) ||
    !setequal(sample_ids, metadata$Sample_ID)) {
    stop(
        "Expected 94 forward-read FASTQ files matching the source metadata.",
        call. = FALSE
    )
}
names(forward_files) <- sample_ids
filtered_files <- file.path(
    filtered_directory, paste0(sample_ids, "_filtered.fastq.gz")
)
names(filtered_files) <- sample_ids

message("Filtering forward reads")
filter_tracking <- dada2::filterAndTrim(
    fwd = forward_files,
    filt = filtered_files,
    truncLen = 299L,
    trimLeft = 20L,
    truncQ = 2L,
    maxEE = Inf,
    maxLen = Inf,
    minLen = 20L,
    rm.phix = TRUE,
    multithread = workers,
    verbose = TRUE
)

message("Learning the forward-read error model")
error_model <- dada2::learnErrors(
    filtered_files,
    multithread = workers,
    verbose = 1L,
    nbases = 1e9,
    randomize = FALSE
)

message("Inferring amplicon sequence variants")
dada_fits <- dada2::dada(
    filtered_files[file.exists(filtered_files)],
    err = error_model,
    multithread = workers,
    verbose = TRUE
)
sequence_table <- dada2::makeSequenceTable(dada_fits)
if (nrow(sequence_table) != 94L ||
    any(nchar(colnames(sequence_table)) != 279L)) {
    stop("The inferred sequence table has an unexpected shape.", call. = FALSE)
}

denoised_reads <- vapply(
    dada_fits, function(fit) sum(dada2::getUniques(fit)), numeric(1)
)
read_tracking <- data.frame(
    sample_id = rownames(filter_tracking),
    input_reads = filter_tracking[, "reads.in"],
    filtered_reads = filter_tracking[, "reads.out"],
    denoised_reads = denoised_reads[rownames(filter_tracking)],
    row.names = NULL,
    check.names = FALSE
)

d6300_files <- sort(list.files(
    file.path(
        reference_directory,
        "ZymoBIOMICS.STD.refseq.v2",
        "ssrRNAs"
    ),
    pattern = "_16S_.*[.]fasta$",
    full.names = TRUE
))
d6300_files <- d6300_files[!grepl(
    "Cryptococcus|Saccharomyces", basename(d6300_files)
)]
d6321_files <- sort(list.files(
    file.path(reference_directory, "D6321.refseq", "16S"),
    pattern = "[.]16S[.]fasta$",
    full.names = TRUE
))
if (length(d6300_files) != 8L || length(d6321_files) != 3L) {
    stop(
        "Expected eight bacterial D6300 and three D6321 16S references.",
        call. = FALSE
    )
}

read_reference <- function(path, source_standard) {
    records <- Biostrings::readDNAStringSet(path)
    species <- if (source_standard == "D6300") {
        sub("_16S_.*$", "", basename(path))
    } else {
        gsub("[.]", "_", sub("[.]16S[.]fasta$", "", basename(path)))
    }
    sequences <- as.character(records)
    primer <- "AGTTTGAT(C|T)(A|C)TGGCTCAG"
    if (!all(grepl(primer, sequences))) {
        stop(sprintf("The forward primer was not found in %s.", basename(path)),
             call. = FALSE)
    }
    sequences <- sub(paste0(".*", primer), "", sequences)
    sequences <- substring(sequences, 1L, 279L)
    if (any(nchar(sequences) != 279L)) {
        stop(sprintf("Reference sequences are too short in %s.", basename(path)),
             call. = FALSE)
    }
    data.frame(
        taxon = species,
        source_standard = source_standard,
        sequence = sequences,
        stringsAsFactors = FALSE
    )
}

reference_records <- do.call(rbind, c(
    lapply(d6300_files, read_reference, source_standard = "D6300"),
    lapply(d6321_files, read_reference, source_standard = "D6321")
))
reference_records <- unique(reference_records)

asv_sequences <- colnames(sequence_table)
distance_matrix <- stringdist::stringdistmatrix(
    asv_sequences,
    reference_records$sequence,
    method = "lv",
    useNames = "none"
)
minimum_distance <- apply(distance_matrix, 1L, min)

resolve_tie <- function(asv_index, candidate_columns) {
    sequences <- c(
        reference_records$sequence[candidate_columns],
        asv_sequences[[asv_index]]
    )
    labels <- c(
        paste0(
            reference_records$taxon[candidate_columns],
            "_reference_",
            seq_along(candidate_columns)
        ),
        "observed_ASV"
    )
    alignment <- ape::as.DNAbin(strsplit(stats::setNames(sequences, labels), ""))
    distances <- as.matrix(ape::dist.dna(alignment))
    observed_distance <- distances["observed_ASV", -length(labels)]
    reference_records$taxon[candidate_columns[which.min(observed_distance)]]
}

assigned_taxa <- vapply(seq_along(asv_sequences), function(index) {
    candidates <- which(distance_matrix[index, ] == minimum_distance[[index]])
    candidate_taxa <- unique(reference_records$taxon[candidates])
    if (minimum_distance[[index]] > 4L) return(NA_character_)
    if (length(candidate_taxa) == 1L) return(candidate_taxa)
    resolve_tie(index, candidates)
}, character(1))

d6300_taxa <- sort(unique(
    reference_records$taxon[reference_records$source_standard == "D6300"]
))
d6321_taxa <- sort(unique(
    reference_records$taxon[reference_records$source_standard == "D6321"]
))
taxa <- c(d6300_taxa, d6321_taxa)
taxon_counts_all <- vapply(taxa, function(taxon) {
    selected_asvs <- which(!is.na(assigned_taxa) & assigned_taxa == taxon)
    if (!length(selected_asvs)) return(integer(nrow(sequence_table)))
    as.integer(rowSums(sequence_table[, selected_asvs, drop = FALSE]))
}, integer(nrow(sequence_table)))
rownames(taxon_counts_all) <- rownames(sequence_table)

metadata <- metadata[match(rownames(sequence_table), metadata$Sample_ID), ]
if (!identical(metadata$Sample_ID, rownames(sequence_table))) {
    stop("The metadata and sequence table could not be aligned.", call. = FALSE)
}

selected_match_ids <- c(917:932, 965:980)
selected <- metadata$Match_ID %in% selected_match_ids
analysis_metadata <- metadata[selected, , drop = FALSE]
analysis_counts <- t(taxon_counts_all[selected, , drop = FALSE])
storage.mode(analysis_counts) <- "integer"
colnames(analysis_counts) <- analysis_metadata$Sample_ID
rownames(analysis_counts) <- taxa

analysis_metadata$sample_id <- analysis_metadata$Sample_ID
analysis_metadata$library_size <- as.integer(rowSums(sequence_table)[selected])
analysis_metadata$source_standard <- ifelse(
    analysis_metadata$Broad_type == "Even mock", "D6300", "D6321"
)
analysis_metadata$input_stratum <- ifelse(
    analysis_metadata$Match_ID %in% c(917:924, 965:972),
    "higher",
    "lower"
)
analysis_metadata$input_cells <- ifelse(
    analysis_metadata$Match_ID %in% c(917:924, 965:972),
    "100000",
    ifelse(analysis_metadata$Match_ID %in% 925:932, "10000", "6000")
)
analysis_metadata$extraction_kit <- analysis_metadata$kit
analysis_metadata$lysis <- analysis_metadata$lysis.prot
analysis_metadata$extraction_buffer <- analysis_metadata$buffer
pair_index <- match(
    analysis_metadata$Match_ID,
    c(917:932, 965:980)
)
analysis_metadata$pair_cell <- ifelse(
    pair_index <= 16L, pair_index, pair_index - 16L
)
analysis_metadata$match_id <- analysis_metadata$Match_ID
analysis_metadata <- analysis_metadata[, c(
    "sample_id", "match_id", "library_size", "source_standard",
    "input_stratum", "input_cells", "extraction_kit", "lysis",
    "extraction_buffer", "pair_cell"
)]
rownames(analysis_metadata) <- analysis_metadata$sample_id

if (ncol(analysis_counts) != 32L || nrow(analysis_counts) != 11L ||
    any(colSums(analysis_counts) > analysis_metadata$library_size) ||
    !all(table(analysis_metadata$source_standard) == 16L)) {
    stop("The prepared benchmark input failed validation.", call. = FALSE)
}

relative_abundance <- sweep(
    analysis_counts,
    2L,
    analysis_metadata$library_size,
    "/"
)
taxon_metadata <- data.frame(
    taxon = taxa,
    source_standard = c(rep("D6300", 8L), rep("D6321", 3L)),
    target_set = c(
        rep("structural absence", 8L),
        rep("present-conditional abundance", 3L)
    ),
    total_count = rowSums(analysis_counts),
    positive_samples = rowSums(analysis_counts > 0),
    positive_samples_D6300 = rowSums(
        analysis_counts[, analysis_metadata$source_standard == "D6300",
                        drop = FALSE] > 0
    ),
    positive_samples_D6321 = rowSums(
        analysis_counts[, analysis_metadata$source_standard == "D6321",
                        drop = FALSE] > 0
    ),
    mean_relative_abundance_D6300 = rowMeans(
        relative_abundance[, analysis_metadata$source_standard == "D6300",
                           drop = FALSE]
    ),
    mean_relative_abundance_D6321 = rowMeans(
        relative_abundance[, analysis_metadata$source_standard == "D6321",
                           drop = FALSE]
    ),
    stringsAsFactors = FALSE
)

preprocessing <- list(
    dataset_id = "rauer_mock",
    source_publication_doi = "10.1186/s40168-024-01998-4",
    raw_sequence_accession = "PRJEB67827",
    forward_read_files = length(forward_files),
    prepared_libraries = ncol(analysis_counts),
    tested_taxa = nrow(analysis_counts),
    selected_match_ids = selected_match_ids,
    dada2_forward_read_settings = list(
        trim_left = 20L,
        trunc_length = 299L,
        trunc_quality = 2L,
        max_expected_errors = Inf,
        minimum_length = 20L,
        remove_phi_x = TRUE,
        remove_chimeras = FALSE
    ),
    reference_assignment = paste(
        "minimum Levenshtein distance no greater than four;",
        "evolutionary distance resolves tied taxa"
    ),
    library_size_definition =
        "total DADA2-denoised forward-read ASV counts before reference-taxon selection",
    count_orientation = "taxa by samples"
)

prepared <- list(
    counts = analysis_counts,
    metadata = analysis_metadata,
    taxon_metadata = taxon_metadata,
    preprocessing = preprocessing
)
saveRDS(
    prepared,
    file.path(processed_directory, "rauer_mock_sacat_input.rds"),
    compress = "xz"
)

count_output <- data.frame(
    taxon = rownames(analysis_counts),
    analysis_counts,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
utils::write.csv(
    count_output,
    file.path(processed_directory, "rauer_mock_taxon_counts.csv"),
    row.names = FALSE,
    quote = TRUE
)
utils::write.csv(
    analysis_metadata,
    file.path(processed_directory, "rauer_mock_metadata.csv"),
    row.names = FALSE,
    quote = TRUE
)
utils::write.csv(
    taxon_metadata,
    file.path(processed_directory, "rauer_mock_taxon_metadata.csv"),
    row.names = FALSE,
    quote = TRUE
)
utils::write.csv(
    read_tracking,
    file.path(processed_directory, "rauer_mock_read_tracking.csv"),
    row.names = FALSE,
    quote = TRUE
)
preprocessing_summary <- data.frame(
    item = c(
        "forward_read_files",
        "prepared_libraries",
        "tested_taxa",
        "native_library_size_median",
        "native_library_size_minimum",
        "native_library_size_maximum"
    ),
    value = c(
        length(forward_files),
        ncol(analysis_counts),
        nrow(analysis_counts),
        stats::median(analysis_metadata$library_size),
        min(analysis_metadata$library_size),
        max(analysis_metadata$library_size)
    ),
    stringsAsFactors = FALSE
)
utils::write.csv(
    preprocessing_summary,
    file.path(processed_directory, "rauer_mock_preprocessing_summary.csv"),
    row.names = FALSE,
    quote = TRUE
)

message(sprintf(
    "Prepared %d libraries and %d reference taxa.",
    ncol(analysis_counts), nrow(analysis_counts)
))
