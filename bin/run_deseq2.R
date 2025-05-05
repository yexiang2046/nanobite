#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(tidyverse)
    library(pheatmap)
    library(RColorBrewer)
})

# Parse command line arguments
option_list <- list(
    make_option(c("--count_matrix"), type="character", help="Path to count matrix file"),
    make_option(c("--design_matrix"), type="character", help="Path to design matrix file"),
    make_option(c("--output_dir"), type="character", help="Output directory for results"),
    make_option(c("--plot_dir"), type="character", help="Output directory for plots")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Read input files
count_matrix <- read.table(opt$count_matrix, header=TRUE, sep="\t", row.names=1)
design_matrix <- read.table(opt$design_matrix, header=TRUE, sep="\t")

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
res_df <- as.data.frame(res)
res_df$transcript_id <- rownames(res_df)

# Write results
write.table(
    res_df,
    file.path(opt$output_dir, "deseq2_results.txt"),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# Generate MA plot
pdf(file.path(opt$plot_dir, "MA_plot.pdf"))
plotMA(res, main="MA Plot")
dev.off()

# Generate PCA plot
vsd <- vst(dds)
pdf(file.path(opt$plot_dir, "PCA_plot.pdf"))
plotPCA(vsd, intgroup="condition") +
    theme_bw() +
    ggtitle("PCA Plot")
dev.off()

# Generate heatmap of top differentially expressed genes
top_genes <- head(order(res$padj), 50)
top_counts <- assay(vsd)[top_genes,]
pdf(file.path(opt$plot_dir, "heatmap.pdf"))
pheatmap(
    top_counts,
    scale="row",
    clustering_distance_rows="correlation",
    clustering_distance_cols="correlation",
    show_rownames=TRUE,
    show_colnames=TRUE,
    main="Top 50 Differentially Expressed Genes"
)
dev.off() 