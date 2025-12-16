#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { modkit_pileup; modkit_extract; modkit_summary; filter_pseU } from '../modules/modkit.nf'
include { samtools_sort } from '../modules/align.nf'

// Default parameters
params.help = false
params.bam_dir = null
params.reference = null
params.output_dir = "results"
params.min_coverage = 5
params.prob_threshold = 0.8
params.filter_pseu = false
params.cpus = 16

// Print help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run workflows/modification_analysis.nf --bam_dir /path/to/bams --reference reference.fa

    Description:
        Analyzes RNA modifications from aligned BAM files containing modification calls.
        Input BAM files should have base modifications already called (e.g., from dorado
        basecaller with --modified-bases flag).

    Required Arguments:
        --bam_dir          Directory containing aligned BAM files with modification calls
        --reference        Reference genome FASTA file

    Optional Arguments:
        --output_dir       Output directory (default: ${params.output_dir})
        --min_coverage     Minimum coverage for modkit pileup (default: ${params.min_coverage})
        --prob_threshold   Probability threshold for modifications (default: ${params.prob_threshold})
        --filter_pseu      Filter for pseudouridine sites (default: ${params.filter_pseu})
        --cpus             Number of CPUs (default: ${params.cpus})

    Output:
        results/modkit_pileup/    - BED files with genome-wide modification pileup
        results/modkit_extract/   - TSV files with per-read modification calls
        results/modkit_summary/   - Summary statistics for all modifications in BAM
        results/bam/              - Sorted BAM files and indices (if input was unsorted)
        results/pseU_sites/       - Pseudouridine sites (if --filter_pseu is enabled)
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (!params.bam_dir || !params.reference) {
    error "Both --bam_dir and --reference are required"
}

if (!file(params.bam_dir).exists()) {
    error "BAM directory not found: ${params.bam_dir}"
}

if (!file(params.reference).exists()) {
    error "Reference file not found: ${params.reference}"
}

// Create channel from BAM files in directory
Channel
    .fromPath("${params.bam_dir}/*.bam")
    .map { bam_file ->
        def sample_id = bam_file.getSimpleName()
        // Check if sorted BAI exists
        def bai_file = file("${bam_file}.bai")
        def sorted_bai = file("${bam_file.getParent()}/${sample_id}.srt.bam.bai")

        // Return tuple with indication if already sorted/indexed
        if (bai_file.exists()) {
            return tuple(sample_id, bam_file, bai_file, true)
        } else if (sorted_bai.exists()) {
            return tuple(sample_id, bam_file, sorted_bai, true)
        } else {
            return tuple(sample_id, bam_file, null, false)
        }
    }
    .set { bam_input_ch }

// Main workflow
workflow {
    // Split channel based on whether BAMs are already sorted/indexed
    bam_input_ch
        .branch {
            sorted: it[3] == true      // Already has index
            unsorted: it[3] == false   // Needs sorting/indexing
        }
        .set { branched_bams }

    // Process unsorted BAMs
    samtools_sort(
        branched_bams.unsorted.map { sample_id, bam_file, bai_file, is_sorted -> bam_file }
    )
    
    // Extract sample_id from sorted BAM filename (sorted filename is based on original filename)
    sorted_bams_ch = samtools_sort.out
        .map { sorted_bam, bai ->
            // Remove .srt extension to get original filename base, which matches sample_id
            def sample_id = sorted_bam.getSimpleName().replaceAll(/\.srt$/, '')
            tuple(sample_id, sorted_bam, bai)
        }

    // Combine already sorted BAMs with newly sorted ones
    all_sorted_bams = branched_bams.sorted
        .map { sample_id, bam, bai, is_sorted -> tuple(sample_id, bam, bai) }
        .mix(sorted_bams_ch)

    // Add reference to each tuple for modkit_pileup
    pileup_input = all_sorted_bams
        .map { sample_id, bam, bai ->
            tuple(sample_id, bam, bai, file(params.reference))
        }

    // Run modkit pileup
    modkit_pileup(pileup_input)

    // Run modkit extract (doesn't need reference)
    modkit_extract(
        all_sorted_bams.map { sample_id, bam, bai -> tuple(sample_id, bam) }
    )

    // Run modkit summary
    modkit_summary(
        all_sorted_bams.map { sample_id, bam, bai -> tuple(sample_id, bam) }
    )

    // Optionally filter for pseudouridine sites
    if (params.filter_pseu) {
        filter_pseU(modkit_pileup.out)
    }
}
