#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript deseq2_analysis.R <count_matrix.txt> <design_matrix.txt> <output_prefix>")
}

count_file <- args[1]
design_file <- args[2]
output_prefix <- args[3]

# Load required libraries
suppressPackageStartupMessages({
    library(DESeq2)
    library(tidyverse)
    library(pheatmap)
})

# Read data
counts <- read.table(count_file, header=TRUE, row.names=1)
design <- read.table(design_file, header=TRUE, row.names=1)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = design,
    design = ~ treatment
)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Save results
write.csv(res, paste0(output_prefix, "_results.csv"))

# Generate MA plot
pdf(paste0(output_prefix, "_MA_plot.pdf"))
plotMA(res)
dev.off()

# Generate PCA plot
vsd <- vst(dds)
pdf(paste0(output_prefix, "_PCA_plot.pdf"))
plotPCA(vsd, intgroup="treatment")
dev.off()

# Generate heatmap of top genes
top_genes <- head(order(res$padj), 50)
pdf(paste0(output_prefix, "_heatmap.pdf"))
pheatmap(assay(vsd)[top_genes,], 
         scale="row",
         show_rownames=TRUE,
         show_colnames=TRUE)
dev.off() 