#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(optparse)
    library(tidyverse)
    library(data.table)
})

# Parse command line arguments
option_list <- list(
    make_option(c("--count_files"), type="character", help="Path to count files"),
    make_option(c("--sample_info"), type="character", help="Path to sample info file"),
    make_option(c("--count_matrix"), type="character", help="Output count matrix file"),
    make_option(c("--design_matrix"), type="character", help="Output design matrix file")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read and process count files
count_files <- strsplit(opt$count_files, " ")[[1]]
count_data <- lapply(count_files, function(f) {
    sample_id <- gsub("_counts.txt", "", basename(f))
    counts <- fread(f, header=TRUE)
    setnames(counts, "count", sample_id)
    return(counts)
})

# Merge count files
count_matrix <- Reduce(function(x, y) {
    merge(x, y, by="transcript_id", all=TRUE)
}, count_data)

# Fill NA with 0
count_matrix[is.na(count_matrix)] <- 0

# Write count matrix
write.table(count_matrix, opt$count_matrix, sep="\t", quote=FALSE, row.names=FALSE)

# Create design matrix
sample_info <- read.table(opt$sample_info, header=TRUE, sep="\t")
design_matrix <- data.frame(
    sample_id = sample_info$SampleId,
    condition = sample_info$treatment
)

# Write design matrix
write.table(design_matrix, opt$design_matrix, sep="\t", quote=FALSE, row.names=FALSE) 