// FeatureCounts for transcript quantification
process featurecounts {
    container 'biocontainers/subread:v1.6.3dfsg-1-deb_cv1'
    publishDir "${params.output_dir}/featurecounts", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file)
    path gtf_file
    path reference_fasta
    
    output:
    tuple val(sample_id), path("${sample_id}_featurecounts.txt"), emit: counts
    path("${sample_id}_featurecounts.txt.summary"), emit: summary
    
    script:
    """
    # Install featureCounts if not available
    if ! command -v featureCounts &> /dev/null; then
        R -e "if (!require('Rsubread')) BiocManager::install('Rsubread', update = FALSE)"
    fi
    
    # Run featureCounts
    featureCounts \
        -a ${gtf_file} \
        -o ${sample_id}_featurecounts.txt \
        -g gene_id \
        -t exon \
        -p \
        -B \
        -C \
        -T ${task.cpus} \
        ${bam_file}
    
    # Create a clean count matrix
    tail -n +3 ${sample_id}_featurecounts.txt | cut -f1,7 > ${sample_id}_counts_clean.txt
    """
}

// Aggregate featureCounts results
process aggregate_featurecounts {
    publishDir "${params.output_dir}/featurecounts", mode: 'copy'
    
    input:
    path count_files
    path sample_info
    
    output:
    path "count_matrix.txt", emit: count_matrix
    path "count_matrix_clean.txt", emit: count_matrix_clean
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    if (!require("dplyr")) install.packages("dplyr", repos="https://cloud.r-project.org")
    if (!require("readr")) install.packages("readr", repos="https://cloud.r-project.org")
    
    library(dplyr)
    library(readr)
    
    # Read sample info
    sample_info <- read_tsv("${sample_info}")
    
    # Function to read featureCounts output
    read_featurecounts <- function(file_path) {
        # Skip the first two lines (header and command info)
        counts <- read_tsv(file_path, skip = 1, col_names = TRUE)
        # Extract gene_id and count columns
        counts <- counts[, c(1, 7)]
        colnames(counts) <- c("Geneid", "Count")
        return(counts)
    }
    
    # Read all count files
    count_files <- list.files(pattern = "*_featurecounts.txt")
    count_data <- list()
    
    for (file in count_files) {
        sample_name <- gsub("_featurecounts.txt", "", file)
        counts <- read_featurecounts(file)
        count_data[[sample_name]] <- counts
    }
    
    # Merge all count data
    merged_counts <- count_data[[1]]
    for (i in 2:length(count_data)) {
        merged_counts <- merge(merged_counts, count_data[[i]], by = "Geneid", all = TRUE)
    }
    
    # Set column names
    colnames(merged_counts) <- c("Geneid", names(count_data))
    
    # Replace NA with 0
    merged_counts[is.na(merged_counts)] <- 0
    
    # Write output
    write_tsv(merged_counts, "count_matrix.txt")
    
    # Create clean version without gene info
    clean_counts <- merged_counts[, -1]
    rownames(clean_counts) <- merged_counts$Geneid
    write_tsv(clean_counts, "count_matrix_clean.txt")
    """
}

// Generate featureCounts statistics
process featurecounts_stats {
    publishDir "${params.output_dir}/featurecounts/stats", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file)
    path gtf_file
    
    output:
    tuple val(sample_id), path("${sample_id}_stats.txt"), emit: stats
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    if (!require("Rsubread")) BiocManager::install("Rsubread", update = FALSE)
    library(Rsubread)
    
    # Get BAM statistics
    bam_stats <- propmapped("${bam_file}")
    
    # Write statistics
    cat("Sample: ${sample_id}\\n", file = "${sample_id}_stats.txt")
    cat("Total reads: ", sum(bam_stats), "\\n", file = "${sample_id}_stats.txt", append = TRUE)
    cat("Mapped reads: ", sum(bam_stats[bam_stats > 0]), "\\n", file = "${sample_id}_stats.txt", append = TRUE)
    cat("Mapping rate: ", round(sum(bam_stats[bam_stats > 0]) / sum(bam_stats) * 100, 2), "%\\n", file = "${sample_id}_stats.txt", append = TRUE)
    """
} 