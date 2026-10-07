#!/usr/bin/env Rscript

# Install metaSPARSim 1.1.2 with its read sampler seeded by R.
install_metasparsim <- function() {
    directory <- tempfile("metasparsim-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    archive <- file.path(directory, "source.tar.gz")
    url <- paste0("https://gitlab.com/sysbiobig/metasparsim/-/archive/",
                  "master/metasparsim-master.tar.gz")
    download.file(url, archive, mode = "wb", quiet = TRUE)
    untar(archive, exdir = directory)
    source <- file.path(directory, "metasparsim-master")
    stopifnot(read.dcf(file.path(source, "DESCRIPTION"))[1L, "Version"] == "1.1.2")
    sampler <- file.path(source, "src", "r_unif.cpp")
    lines <- readLines(sampler)
    original <- "std::default_random_engine rng(time(0));"
    target <- grep(original, lines, fixed = TRUE)
    stopifnot(length(target) == 1L)
    lines[target] <- paste0(
        "  // Modified October 7, 2026: obtain the sampler seed from R's RNG.\n",
        "  std::default_random_engine rng(",
        "static_cast<unsigned int>(R::runif(1.0, 2147483647.0)));")
    writeLines(lines, sampler)
    install.packages(source, repos = NULL, type = "source")
}

install_metasparsim()
