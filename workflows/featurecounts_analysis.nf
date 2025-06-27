#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import featureCounts modules
include { featurecounts_stats } from '../modules/featurecounts.nf'
include { featurecounts } from '../modules/featurecounts.nf'
include { aggregate_featurecounts } from '../modules/featurecounts.nf'

// Default parameters
params.help = false
params.bam_dir = null
params.gtf_file = null
params.reference = null
params.sample_info = null
params.output_dir = "results/featurecounts"

// Print help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run featurecounts_analysis.nf --bam_dir /path/to/bams --gtf_file annotation.gtf --sample_info sample_info.txt

    Required Arguments:
        --bam_dir         Directory containing sorted BAM files
        --gtf_file        GTF annotation file
        --sample_info     Sample info file with SampleId and treatment columns

    Optional Arguments:
        --reference       Reference genome FASTA file (optional)
        --output_dir      Output directory (default: results/featurecounts)
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (!params.bam_dir || !params.gtf_file || !params.sample_info) {
    error "Required parameters: --bam_dir, --gtf_file, and --sample_info"
}

// Create BAM file channel
Channel
    .fromPath("${params.bam_dir}/*.bam")
    .ifEmpty { error "No BAM files found in ${params.bam_dir}" }
    .map { bam ->
        def sample_id = bam.getSimpleName().replace('_aligned.srt', '').replace('.bam', '')
        tuple(sample_id, bam)
    }
    .set { bam_ch }

// Main workflow
workflow FEATURECOUNTS_ANALYSIS {
    take:
        bam_files
        gtf_file
        reference_file
        sample_info

    main:
        // Run featureCounts for each sample
        featurecounts(bam_files, gtf_file, reference_file)

        // Generate statistics
        featurecounts_stats(bam_files, gtf_file)

        // Aggregate results
        aggregate_featurecounts(
            featurecounts.out.counts.collect(),
            sample_info
        )

    emit:
        counts = featurecounts.out.counts
        summary = featurecounts.out.summary
        stats = featurecounts_stats.out.stats
        count_matrix = aggregate_featurecounts.out.count_matrix
        count_matrix_clean = aggregate_featurecounts.out.count_matrix_clean
}

// Entry point
workflow {
    // Reference file (optional)
    def reference_ch = params.reference ? Channel.fromPath(params.reference) : Channel.value(file('NO_FILE'))

    FEATURECOUNTS_ANALYSIS(
        bam_ch,
        file(params.gtf_file),
        reference_ch,
        file(params.sample_info)
    )
} 