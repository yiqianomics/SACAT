test_that("the CDI example has aligned counts and analysis covariates", {
    example <- new.env()
    data("schubert_cdi", package = "SACAT", envir = example)
    x <- example$schubert_cdi
    expect_identical(dim(x$counts), c(88L, 243L))
    expect_identical(colnames(x$counts), rownames(x$metadata))
    expect_type(x$counts, "integer")
    expect_true(all(x$counts >= 0))
    expect_false(anyDuplicated(rownames(x$counts)) > 0L)
    expect_identical(names(x$metadata), c(
        "group", "age_z", "sex", "antibiotics_3mo", "library_size"
    ))
    expect_identical(levels(x$metadata$group), c("H", "CDI"))
    expect_identical(as.integer(table(x$metadata$group)), c(152L, 91L))
    expect_false(anyNA(x$metadata))
    expect_true(all(x$metadata$library_size >= colSums(x$counts)))
    expect_equal(mean(x$metadata$age_z), 0, tolerance = 1e-12)
    expect_equal(stats::sd(x$metadata$age_z), 1, tolerance = 1e-12)
    expect_true(all(c("Coprococcus", "Enterococcus", "Faecalibacterium") %in%
                        rownames(x$counts)))
    design <- stats::model.matrix(
        ~ age_z + sex + antibiotics_3mo + group, data = x$metadata
    )
    expect_identical(qr(design)$rank, ncol(design))
})
