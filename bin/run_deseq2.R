#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript run_deseq2.R <count_matrix.txt> <design_matrix.txt>")
}

count_matrix_file <- args[1]
design_matrix_file <- args[2]

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# Load required packages
if (!require("DESeq2")) install.packages("DESeq2", repos="https://cloud.r-project.org")
if (!require("ggplot2")) install.packages("ggplot2", repos="https://cloud.r-project.org")

# Read input files
count_matrix <- read.table(count_matrix_file, header=TRUE, sep="\t", row.names=1)
design_matrix <- read.table(design_matrix_file, header=TRUE, sep="\t")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = design_matrix,
    design = ~ condition
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Save results
write.table(res, "results/differential_expression.tsv", sep="\t", quote=FALSE)

# Create plots
pdf("plots/MA_plot.pdf")
plotMA(res)
dev.off()

pdf("plots/PCA_plot.pdf")
plotPCA(vst(dds))
dev.off() 