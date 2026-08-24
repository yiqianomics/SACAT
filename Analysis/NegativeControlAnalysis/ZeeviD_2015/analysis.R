#!/usr/bin/env Rscript

script_arguments <- commandArgs(trailingOnly = FALSE)
script_file <- grep("^--file=", script_arguments, value = TRUE)
dataset_directory <- if (length(script_file)) {
    normalizePath(dirname(sub("^--file=", "", script_file[[1L]])))
} else {
    normalizePath(getwd())
}
source(file.path(dataset_directory, "..", "run_negative_control.R"))
run_negative_control(dataset_directory, workers = 7L, replicates = 100L)
