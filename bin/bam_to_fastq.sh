#!/bin/bash

# Function to print usage
usage() {
    echo "Usage: $0 <input_directory> <output_directory>"
    echo "Convert BAM files to FASTQ.gz format"
    echo ""
    echo "Arguments:"
    echo "  input_directory   Directory containing BAM files"
    echo "  output_directory  Directory where FASTQ.gz files will be saved"
    exit 1
}

# Check if samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is not installed. Please install samtools first."
    exit 1
fi

# Check arguments
if [ "$#" -ne 2 ]; then
    usage
fi

input_dir="$1"
output_dir="$2"

# Check if directories exist
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory does not exist: $input_dir"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Process each BAM file
for bam_file in "$input_dir"/*.bam; do
    # Check if there are any BAM files
    if [ ! -e "$bam_file" ]; then
        echo "No BAM files found in $input_dir"
        exit 1
    fi

    # Get base name for output files
    base_name=$(basename "$bam_file" .bam)
    echo "Processing: $base_name"

    # Convert BAM to FASTQ and compress with gzip
    echo "Converting $bam_file to FASTQ.gz..."
    samtools bam2fq "$bam_file" | gzip > "$output_dir/${base_name}.fastq.gz" 2> "$output_dir/${base_name}_conversion.log"

    # Check if conversion was successful
    if [ $? -eq 0 ]; then
        echo "Successfully converted $base_name"
    else
        echo "Error converting $base_name. Check log file: $output_dir/${base_name}_conversion.log"
    fi
done

echo "All conversions completed!" 