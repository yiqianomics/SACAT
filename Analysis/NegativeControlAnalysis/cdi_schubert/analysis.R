#!/usr/bin/env Rscript

script_file <- grep("^--file=", commandArgs(), value = TRUE)
dataset_directory <- if (length(script_file)) {
    dirname(normalizePath(sub("^--file=", "", script_file[[1L]])))
} else {
    normalizePath(getwd())
}
source(file.path(dataset_directory, "..", "run_negative_control.R"))
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) <= 1L)
workers <- if (length(arguments)) as.integer(arguments[[1L]]) else 7L
run_negative_control(dataset_directory, workers = workers, replicates = 100L)
