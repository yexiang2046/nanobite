#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript merge_counts.R <sample_info.txt> <count_files_dir>")
}

sample_info <- args[1]
count_dir <- args[2]

# Read sample info
samples <- read.table(sample_info, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Initialize empty list to store count data
count_data <- list()

# Process each sample
for (i in 1:nrow(samples)) {
    sample_id <- trimws(samples$SampleId[i])
    count_file <- file.path(count_dir, paste0(sample_id, "_aligned_transcript_counts.tsv"))
    
    if (file.exists(count_file)) {
        cat("Processing", sample_id, "...\n")
        # Read count file
        counts <- read.table(count_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        # Keep only transcript_name and est_count columns
        counts <- counts[, c("transcript_name", "est_count")]
        # Rename est_count column to sample_id
        names(counts)[2] <- sample_id
        # Store in list
        count_data[[sample_id]] <- counts
    } else {
        warning("No count file found for ", sample_id, "\nLooked for: ", count_file)
    }
}

# Check if we have any data to merge
if (length(count_data) == 0) {
    stop("Error: No count files were processed. Please check your sample info file and count files.")
}

# Merge all count files
cat("Merging count files...\n")
# Start with first sample
merged_counts <- count_data[[1]]

# Merge remaining samples
for (i in 2:length(count_data)) {
    merged_counts <- merge(merged_counts, count_data[[i]], by = "transcript_name", all = TRUE)
}

# Replace NA with 0
merged_counts[is.na(merged_counts)] <- 0

# Write merged count matrix
write.table(merged_counts, 
            file = "count_matrix.txt",
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# Create design matrix
cat("Creating design matrix...\n")
design_matrix <- data.frame(
    sample_id = samples$SampleId,
    condition = samples$treatment
)
write.table(design_matrix,
            file = "design_matrix.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Debug output
cat("\nFinal count matrix:\n")
print(head(merged_counts, 5))
cat("\nNumber of columns in count matrix:", ncol(merged_counts), "\n")

cat("\nDone! Output files are in current directory\n")
cat("- count_matrix.txt: Merged count matrix\n")
cat("- design_matrix.txt: Sample design matrix\n") 