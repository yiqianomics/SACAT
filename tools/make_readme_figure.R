library(SACAT)
data("schubert_cdi")

fit <- sacat(
    counts = schubert_cdi$counts,
    metadata = schubert_cdi$metadata,
    formula = ~ age_z + sex + antibiotics_3mo + group,
    group = "group",
    reference = "H",
    library_size = "library_size",
    store_plot_data = TRUE
)

dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
for (extension in c("svg", "pdf")) {
    plot(
        fit, features = c("Coprococcus", "Enterococcus", "Faecalibacterium"),
        group_labels = c("Healthy", "CDI"),
        file = paste0("man/figures/schubert-cdi-example.", extension),
        width = 9.5, height = 3.8
    )
}
