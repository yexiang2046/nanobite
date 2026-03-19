#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { featurecounts_genome_aligned } from '../modules/featurecounts.nf'
include { aggregate_featurecounts } from '../modules/featurecounts.nf'
include { prepare_deseq2; deseq2_analysis } from '../modules/deseq2.nf'

// Default parameters
params.help = false
params.bam_dir = null
params.gtf_file = null
params.sample_info = null
params.output_dir = "results/transcript_counts"

// Help message
def helpMessage() {
    log.info"""
    ========================================
    Genome Transcript Quantification Pipeline
    ========================================

    Quantify all transcripts from genome-aligned DRS BAM files using featureCounts
    with optimized parameters for Nanopore long-read data.

    Usage:
        nextflow run genome_transcript_quantification.nf --bam_dir <dir> --gtf_file <gtf> --sample_info <file>

    Required Arguments:
        --bam_dir           Directory containing genome-aligned BAM files (.bam + .bai indices)
        --gtf_file          GTF annotation file (GENCODE, Ensembl, etc.)
        --sample_info       Sample info file (tab-separated with SampleId, treatment columns)

    Optional Arguments:
        --output_dir        Output directory (default: results/transcript_counts)
        --stranded          Strandedness: 0=unstranded, 1=stranded, 2=reverse (default: 1 for DRS)
        --min_mapq          Minimum mapping quality (default: 10)
        --feature_type      Feature type to count from GTF (default: exon)
        --gene_attribute    GTF attribute for grouping (default: gene_id)
        --run_deseq2        Run differential expression analysis (default: false)

    Multi-mapper Handling Options:
        --count_multimappers        Count multi-mapping reads with -M flag (default: false)
        --use_fractional_counting   Use fractional counting for multi-mappers (default: false, requires -M)
        --count_overlapping         Count reads overlapping multiple features (default: false)

    Examples:

        1. Exclude multi-mappers (conservative, default):
            nextflow run genome_transcript_quantification.nf \\
                --bam_dir bams/ \\
                --gtf_file gencode.v44.annotation.gtf \\
                --sample_info sample_info.txt

        2. Count multi-mappers with fractional counting (recommended for small RNAs):
            nextflow run genome_transcript_quantification.nf \\
                --bam_dir bams/ \\
                --gtf_file gencode.v44.annotation.gtf \\
                --sample_info sample_info.txt \\
                --count_multimappers true \\
                --use_fractional_counting true

        3. Count all multi-mapper alignments (permissive):
            nextflow run genome_transcript_quantification.nf \\
                --bam_dir bams/ \\
                --gtf_file gencode.v44.annotation.gtf \\
                --sample_info sample_info.txt \\
                --count_multimappers true

        4. With differential expression analysis:
            nextflow run genome_transcript_quantification.nf \\
                --bam_dir bams/ \\
                --gtf_file gencode.v44.annotation.gtf \\
                --sample_info sample_info.txt \\
                --run_deseq2 true

        5. Custom feature counting (use gene-level instead of exon):
            nextflow run genome_transcript_quantification.nf \\
                --bam_dir bams/ \\
                --gtf_file gencode.v44.annotation.gtf \\
                --sample_info sample_info.txt \\
                --feature_type gene \\
                --gene_attribute gene_name

    Multi-mapper Handling Modes:

        Mode 1 - Exclude multi-mappers (Default, Conservative):
            --count_multimappers false
            • Multi-mapping reads are NOT counted
            • Use for: Unique gene quantification
            • Downside: Underestimates gene families, small RNAs with multiple copies

        Mode 2 - Count all multi-mapper alignments (Permissive):
            --count_multimappers true --use_fractional_counting false
            • Each alignment of multi-mapper is counted separately
            • Read mapping to 5 locations → adds 5 to counts
            • Warning: Inflates counts for repetitive genes

        Mode 3 - Fractional counting (Recommended for small RNAs):
            --count_multimappers true --use_fractional_counting true
            • Multi-mapper counts distributed proportionally
            • Read mapping to 5 locations → each gets 0.2 counts
            • Best for: Small RNAs (tRNAs, snRNAs with multiple genomic copies)

    Sample Info File Format (tab-separated):
        SampleId    treatment
        sample1     control
        sample2     treated
        sample3     control
        sample4     treated

    Output Structure:
        results/transcript_counts/
        ├── sample1_featurecounts.txt              # Per-sample counts
        ├── sample1_featurecounts.txt.summary      # Assignment statistics
        ├── count_matrix.txt                       # Aggregated counts (all columns)
        └── count_matrix_clean.txt                 # Counts only (for DESeq2)

    Notes:
        • BAM files must be sorted and indexed (.bai)
        • This workflow uses corrected featureCounts parameters for single-end DRS data
        • Removes incorrect paired-end flags (-p, -B, -C) from original implementation
        • Adds long-read mode (-L) for Nanopore data
        • Adds strand specification (-s) for DRS data
        • Uses --primary flag to avoid double-counting secondary alignments

    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required inputs
if (!params.bam_dir || !params.gtf_file || !params.sample_info) {
    error """
    ERROR: Missing required parameters!

    Required:
        --bam_dir       Directory with BAM files
        --gtf_file      GTF annotation file
        --sample_info   Sample information file

    Run with --help for usage information.
    """
}

// Validate files exist
def gtf = file(params.gtf_file)
if (!gtf.exists()) {
    error "GTF file not found: ${params.gtf_file}"
}

def sample_info = file(params.sample_info)
if (!sample_info.exists()) {
    error "Sample info file not found: ${params.sample_info}"
}

def bam_dir = file(params.bam_dir)
if (!bam_dir.exists()) {
    error "BAM directory not found: ${params.bam_dir}"
}

// Create BAM file channel
Channel
    .fromPath("${params.bam_dir}/*.bam")
    .ifEmpty { error "No BAM files found in ${params.bam_dir}" }
    .map { bam ->
        // Extract sample ID from filename
        def sample_id = bam.getSimpleName()
            .replace('_aligned.srt', '')
            .replace('_aligned', '')
            .replace('.srt', '')
            .replace('.sorted', '')

        // Check if BAI index exists
        def bai = file("${bam}.bai")
        if (!bai.exists()) {
            log.warn "BAI index not found for ${bam.name}, may cause issues"
        }

        tuple(sample_id, bam)
    }
    .set { bam_ch }

// Main workflow
workflow {
    log.info """
    ========================================
    Genome Transcript Quantification Pipeline
    ========================================
    BAM directory    : ${params.bam_dir}
    GTF file         : ${params.gtf_file}
    Sample info      : ${params.sample_info}
    Output directory : ${params.output_dir}

    featureCounts Parameters:
    -------------------------
    Stranded         : ${params.stranded}
    Min MAPQ         : ${params.min_mapq}
    Feature type     : ${params.feature_type}
    Gene attribute   : ${params.gene_attribute}

    Multi-mapper Handling:
    ----------------------
    Count multi-mappers    : ${params.count_multimappers}
    Fractional counting    : ${params.use_fractional_counting}
    Count overlapping      : ${params.count_overlapping}

    Analysis:
    ---------
    Run DESeq2       : ${params.run_deseq2}
    ========================================
    """.stripIndent()

    // Run featureCounts on each BAM file
    featurecounts_genome_aligned(bam_ch, gtf)

    // Aggregate counts across all samples
    aggregate_featurecounts(
        featurecounts_genome_aligned.out.counts.collect(),
        sample_info
    )

    // Optional: Run DESeq2 differential expression analysis
    if (params.run_deseq2) {
        log.info "Running DESeq2 differential expression analysis..."

        prepare_deseq2(
            aggregate_featurecounts.out.count_matrix_clean,
            sample_info
        )

        deseq2_analysis(
            prepare_deseq2.out.count_matrix,
            prepare_deseq2.out.design_matrix
        )
    }
}

workflow.onComplete {
    log.info """
    ========================================
    Pipeline completed!
    ========================================
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration  : ${workflow.duration}
    Output    : ${params.output_dir}
    ========================================
    """.stripIndent()
}
