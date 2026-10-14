# Prepare the cohort with Analysis/RealDataAnalysis/cdi_schubert/prepare_data.R.
input <- readRDS(
    "Analysis/RealDataAnalysis/cdi_schubert/processed/cdi_schubert_sacat_input.rds"
)

schubert_cdi <- list(
    counts = input$counts,
    metadata = input$metadata[, c(
        "group", "age_z", "sex", "antibiotics_3mo", "library_size"
    )]
)

dir.create("data", showWarnings = FALSE)
save(schubert_cdi, file = "data/schubert_cdi.rda", compress = "xz")
