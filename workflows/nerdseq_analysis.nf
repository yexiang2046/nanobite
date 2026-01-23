#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { bam_to_fastq_mod; minimap2_align_mod } from '../modules/align.nf'
include { samtools_sort } from '../modules/align.nf'
include { CHOPPER_QC } from '../modules/chopper_filter.nf'
include { modkit_pileup; modkit_extract; modkit_summary } from '../modules/modkit.nf'

// Default parameters
params.help = false
params.reference = null
params.bam_dir = null
params.output_dir = "results"
params.min_coverage = 5
params.prob_threshold = 0.8
params.mm2opts_nerd = "-ax sr"

// Print help message
def helpMessage() {
    log.info"""
    ========================================
    NERD-seq Style RNA Modification Analysis Pipeline
    ========================================

    This pipeline processes basecalled BAM files (with MM/ML tags) using minimap2
    alignment with short read RNA parameters while preserving modification tags,
    followed by comprehensive modification analysis using modkit tools.

    Usage:
    nextflow run workflows/nerdseq_analysis.nf --reference reference.fa --bam_dir /path/to/bams

    Required Arguments:
        --reference         Reference genome FASTA file
        --bam_dir           Directory containing basecalled BAM files with MM/ML tags

    Optional Arguments:
        --output_dir        Output directory (default: "results")
        --mm2opts_nerd      Minimap2 alignment options for NERD-seq (default: "-ax sr")
        --min_coverage      Minimum coverage for modification calling (default: 5)
        --prob_threshold    Probability threshold for modification calls (default: 0.8)

    Output Structure:
        results/
        ├── alignment/            # Aligned BAM files (MM/ML tags preserved)
        ├── bam/                  # Sorted BAM files and indices
        ├── modkit_pileup/        # Per-site modification frequencies
        ├── modkit_extract/       # Read-level modification calls
        └── modkit_summary/       # Summary statistics

    Example:
    # Process basecalled BAMs with RNA modifications
    nextflow run workflows/nerdseq_analysis.nf \\
        --reference genome.fa \\
        --bam_dir basecalling_output/ \\
        --mm2opts_nerd "-ax sr" \\
        --prob_threshold 0.8
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (!params.reference || !params.bam_dir) {
    error "Both --reference and --bam_dir are required. Use --help for usage information."
}

// Validate reference file exists
reference_file = file(params.reference)
if (!reference_file.exists()) {
    error "Reference file not found: ${params.reference}"
}

// Validate BAM directory exists
def bam_directory = file(params.bam_dir)
if (!bam_directory.exists() || !bam_directory.isDirectory()) {
    error "BAM directory not found or is not a directory: ${params.bam_dir}"
}

// Create BAM channel from directory
def create_bam_channel() {
    return Channel
        .fromPath("${params.bam_dir}/*.bam")
        .map { bam ->
            tuple(bam.getSimpleName(), bam)
        }
}

// Main workflow
workflow {
    // Create BAM input channel
    bam_input_ch = create_bam_channel()

    // Step 1: Convert BAM to FASTQ while preserving MM/ML tags
    bam_to_fastq_mod(bam_input_ch)

    // Step 2: Filter FASTQ with Chopper
    CHOPPER_QC(bam_to_fastq_mod.out)

    // Step 3: Align FASTQ with minimap2 preserving MM/ML tags
    fastq_ref_ch = CHOPPER_QC.out.filtered_fastq
        .map { fastq ->
            def sample_id = fastq.getSimpleName().replaceAll(/_filtered$/, '')
            tuple(sample_id, fastq, reference_file)
        }
    minimap2_align_mod(fastq_ref_ch)

    // Step 4: Sort BAM and create index
    samtools_sort(minimap2_align_mod.out)

    // Create channel with sample ID, BAM, and BAI files for modkit pileup
    bam_bai_ch = samtools_sort.out
        .map { bam, bai ->
            tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam, bai, reference_file)
        }

    // Step 5: Run modkit pileup to get per-site modification frequencies
    modkit_pileup(bam_bai_ch)

    // Step 6: Extract read-level modification calls
    modkit_extract(
        samtools_sort.out
            .map { bam, bai ->
                tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam)
            }
    )

    // Step 7: Generate modification summary statistics
    modkit_summary(
        samtools_sort.out
            .map { bam, bai ->
                tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam)
            }
    )
}

workflow.onComplete {
    log.info"""
    ========================================
    Pipeline execution completed!
    ========================================
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Output directory: ${params.output_dir}

    Alignment: Minimap2 with ${params.mm2opts_nerd} (MM/ML tags preserved)
    Modification probability threshold: ${params.prob_threshold}
    Minimum coverage: ${params.min_coverage}

    Results:
    - Aligned BAMs: ${params.output_dir}/alignment/
    - Sorted BAMs: ${params.output_dir}/bam/
    - Modification pileup: ${params.output_dir}/modkit_pileup/
    - Modification calls: ${params.output_dir}/modkit_extract/
    - Summary statistics: ${params.output_dir}/modkit_summary/
    """.stripIndent()
}
