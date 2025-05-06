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

# Convert counts to integers
count_matrix <- round(count_matrix)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = design_matrix,
    design = ~ condition
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds, independentFiltering=FALSE)

# Save results
write.table(res, "results/differential_expression.tsv", sep="\t", quote=FALSE)

# Create MA plot
pdf("plots/MA_plot.pdf")
plotMA(res)
dev.off()

# Create PCA plot
pdf("plots/PCA_plot.pdf")
# Use variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds)
# Calculate PCA
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Create PCA plot
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA Plot") +
    theme_bw()
dev.off() 