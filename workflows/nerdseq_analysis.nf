#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { bam_to_fastq_mod; minimap2_align_mod } from '../modules/align.nf'
include { samtools_sort } from '../modules/align.nf'
include { NANOFILT_QC } from '../modules/nanofilt_filter.nf'
include { modkit_pileup; modkit_summary } from '../modules/modkit.nf'

// Default parameters
params.help = false
params.reference = null
params.bam_dir = null
params.output_dir = "results"
params.min_coverage = 5
params.prob_threshold = 0.8
params.mm2opts_nerd = "-ax sr"
params.run_nanofilt = false
params.pileup_only = false

params.nanofilt_min_qual   = 10
params.nanofilt_min_len    = 500
params.nanofilt_headcrop   = 0
params.nanofilt_tailcrop   = 0

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
        --pileup_only       Skip alignment; run only modkit pileup on aligned BAMs in --bam_dir (default: false).
                            BAMs must be coordinate-sorted and indexed (each .bam must have a .bam.bai).
        --mm2opts_nerd      Minimap2 alignment options for NERD-seq (default: "-ax sr")
        --min_coverage      Minimum coverage for modification calling (default: 5)
        --prob_threshold    Probability threshold for modification calls (default: 0.8)
        --run_nanofilt      Enable NanoFilt filtering step (default: false)
        --nanofilt_min_qual Minimum quality score for NanoFilt filtering (default: 10)
        --nanofilt_min_len  Minimum read length for NanoFilt filtering (default: 500)
        --nanofilt_headcrop Trim N bases from the start of reads (default: 0)
        --nanofilt_tailcrop Trim N bases from the end of reads (default: 0)

    Output Structure:
        results/
        ├── nanofilt_qc/          # NanoFilt filtered FASTQ files (if --run_nanofilt)
        ├── alignment/            # Aligned BAM files (MM/ML tags preserved)
        ├── bam/                  # Sorted BAM files and indices
        ├── modkit_pileup/        # Per-site modification frequencies
        └── modkit_summary/       # Summary statistics

    Example:
    # Process basecalled BAMs with RNA modifications
    nextflow run workflows/nerdseq_analysis.nf \\
        --reference genome.fa \\
        --bam_dir basecalling_output/ \\
        --mm2opts_nerd "-ax sr" \\
        --prob_threshold 0.8

    # Only run modkit pileup on already-aligned BAMs (BAMs must be sorted and indexed)
    nextflow run workflows/nerdseq_analysis.nf \\
        --reference genome.fa \\
        --bam_dir /path/to/aligned_bams/ \\
        --pileup_only
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
            tuple(bam.getSimpleName().replaceAll(/\.bam$/, ''), bam)
        }
}

// Create channel for pileup-only mode: (sample_id, bam, bai, reference). BAMs must be indexed.
def create_pileup_only_channel() {
    return Channel
        .fromPath("${params.bam_dir}/*.bam")
        .map { bam ->
            def sample_id = bam.getSimpleName().replaceAll(/\.bam$/, '')
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                error "With --pileup_only, each BAM must have an index (e.g. ${bai}). Create one with: samtools index ${bam}"
            }
            tuple(sample_id, bam, bai, reference_file)
        }
}

// Main workflow
workflow {
    // Create BAM input channel (always needed)
    bam_input_ch = create_bam_channel()

    if (!params.pileup_only) {
        // ── Full mode: basecalled BAM → FASTQ → (optional NanoFilt) → minimap2 → sort → modkit ──
        // Step 1: Convert BAM to FASTQ while preserving MM/ML tags
        bam_to_fastq_mod(bam_input_ch)

        // Step 2: Optionally filter FASTQ with NanoFilt
        if (params.run_nanofilt) {
            NANOFILT_QC(bam_to_fastq_mod.out)
            fastq_for_align = NANOFILT_QC.out.filtered_fastq
                .map { fastq ->
                    def sample_id = fastq.getSimpleName().replaceAll(/_nanofilt$/, '')
                    tuple(sample_id, fastq)
                }
        } else {
            fastq_for_align = bam_to_fastq_mod.out
        }

        // Step 3: Align with minimap2 (MM/ML tags preserved)
        fastq_ref_ch = fastq_for_align
            .map { sample_id, fastq ->
                tuple(sample_id, fastq, reference_file)
            }
        minimap2_align_mod(fastq_ref_ch)

        // Step 4: Sort BAM + index
        samtools_sort(minimap2_align_mod.out)

        // Prepare channel for modkit (sample_id, bam, bai, reference)
        bam_bai_ch = samtools_sort.out
            .map { bam, bai ->
                tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam, bai, reference_file)
            }

        // Step 5: modkit pileup
        modkit_pileup(bam_bai_ch)

        // Step 6: modkit summary
        modkit_summary(
            samtools_sort.out
                .map { bam, bai ->
                    tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam)
                }
        )
    }
    else {
        // ── Pileup-only mode (already aligned + indexed BAMs) ──
        pileup_only_ch = create_pileup_only_channel()
        modkit_pileup(pileup_only_ch)
        modkit_summary(bam_input_ch)
    }
}

workflow.onComplete {
    def nanofilt_result = params.run_nanofilt ? "- Filtered FASTQ: ${params.output_dir}/nanofilt_qc/" : ""
    def pileup_only_result = params.pileup_only ? "- Pileup-only mode (aligned BAMs used as input)\n    " : ""
    log.info"""
    ========================================
    Pipeline execution completed!
    ========================================
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Output directory: ${params.output_dir}

    Mode: ${params.pileup_only ? 'Pileup only (aligned BAMs)' : 'Full (basecalled BAM → align → pileup)'}
    ${params.pileup_only ? '' : "Alignment: Minimap2 with ${params.mm2opts_nerd} (MM/ML tags preserved)"}
    ${params.pileup_only ? '' : "NanoFilt filtering: ${params.run_nanofilt ? 'ENABLED' : 'DISABLED'}"}
    Modification probability threshold: ${params.prob_threshold}
    Minimum coverage: ${params.min_coverage}

    Results:
    ${params.pileup_only ? pileup_only_result : nanofilt_result}
    ${params.pileup_only ? '' : "- Aligned BAMs: ${params.output_dir}/alignment/\n    - Sorted BAMs: ${params.output_dir}/bam/"}
    - Modification pileup: ${params.output_dir}/modkit_pileup/
    - Summary statistics: ${params.output_dir}/modkit_summary/
    """.stripIndent()
}
