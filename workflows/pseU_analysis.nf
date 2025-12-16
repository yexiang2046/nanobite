#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { basecalling_pseU } from '../modules/basecalling.nf'
include { dorado_align } from '../modules/align.nf'
include { samtools_sort } from '../modules/align.nf'
include { modkit_pileup; modkit_extract; modkit_summary; filter_pseU } from '../modules/modkit.nf'

// Default parameters
params.help = false
params.reference = null
params.sample_info = null
params.output_dir = "results"
params.min_coverage = 5
params.prob_threshold = 0.8
params.min_qscore = 10
params.mm2opts = '--mm2-opts \'ax sr\''
params.skip_basecalling = false
params.bam_dir = null

// Print help message
def helpMessage() {
    log.info"""
    ========================================
    Pseudouridine (ψ) Modification Analysis Pipeline
    ========================================

    This pipeline performs basecalling with pseudouridine modification detection,
    alignment with Q score filtering using dorado, and comprehensive modification
    analysis using modkit.

    Usage:
    nextflow run workflows/pseU_analysis.nf --reference reference.fa --sample_info sample_info.txt

    Required Arguments:
        --reference         Reference genome FASTA file
        --sample_info      Tab-delimited sample info file with columns: SampleId, pod5_path
                          (or SampleId, bam_path when using --skip_basecalling)

    Optional Arguments:
        --output_dir        Output directory (default: "results")
        --skip_basecalling  Skip basecalling step and start with existing BAM files (default: false)
        --bam_dir           Directory containing basecalled BAM files (required when --skip_basecalling)
        --min_coverage      Minimum coverage for modification calling (default: 5)
        --prob_threshold    Probability threshold for modification calls (default: 0.8)
        --min_qscore        Minimum Q score for read filtering during alignment (default: 10)
        --mm2opts_sr        Minimap2 short read alignment options (default: "--mm2-opts '-ax sr'")

    Output Structure:
        results/
        ├── basecalling_pseU/     # Basecalled BAM files with ψ modifications
        ├── alignment/            # Dorado aligned BAM files (Q score filtered)
        ├── bam/                  # Sorted BAM files and indices
        ├── modkit_pileup/        # Per-site modification frequencies
        ├── modkit_extract/       # Read-level modification calls
        ├── modkit_summary/       # Summary statistics
        └── pseU_sites/           # Filtered pseudouridine sites
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (!params.reference || !params.sample_info) {
    error "Both --reference and --sample_info are required. Use --help for usage information."
}

// Validate skip_basecalling options
if (params.skip_basecalling && !params.bam_dir) {
    error "When using --skip_basecalling, --bam_dir must be specified. Use --help for usage information."
}

// Validate reference file exists
def reference_file = file(params.reference)
if (!reference_file.exists()) {
    error "Reference file not found: ${params.reference}"
}

// Validate BAM directory exists when skip_basecalling is enabled
if (params.skip_basecalling) {
    def bam_directory = file(params.bam_dir)
    if (!bam_directory.exists() || !bam_directory.isDirectory()) {
        error "BAM directory not found or is not a directory: ${params.bam_dir}"
    }
}

// Create sample channel from sample info file
def create_sample_channel() {
    def required_column = params.skip_basecalling ? 'bam_path' : 'pod5_path'

    return Channel
        .fromPath(params.sample_info)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            // Validate required columns
            if (!row.containsKey('SampleId') || !row.containsKey(required_column)) {
                error "Sample info file must contain 'SampleId' and '${required_column}' columns. Found columns: ${row.keySet()}"
            }

            def sample_id = row.SampleId?.trim()
            def input_path = file(row[required_column]?.trim())

            // Validate sample data
            if (!sample_id) {
                error "SampleId cannot be empty"
            }
            if (!input_path.exists()) {
                error "${params.skip_basecalling ? 'BAM' : 'POD5'} file not found: ${input_path}"
            }

            return tuple(sample_id, input_path)
        }
}

// Create BAM channel from directory when skipping basecalling
def create_bam_channel() {
    return Channel
        .fromPath("${params.bam_dir}/*.bam")
        .map { bam ->
            tuple(bam.getSimpleName(), bam)
        }
}

// Main workflow
workflow {
    // Conditional workflow based on skip_basecalling parameter
    if (params.skip_basecalling) {
        // Start from existing BAM files
        bam_input_ch = create_bam_channel()

        // Step 1: Align BAM to reference with Q score filtering using dorado
        dorado_align(
            bam_input_ch
                .map { sampleId, bam ->
                    tuple(params.mm2opts, sampleId, bam, reference_file)
                }
        )

        // Step 2: Sort BAM and create index
        samtools_sort(dorado_align.out)

        // Create channel with sample ID, BAM, and BAI files
        bam_bai_ch = samtools_sort.out
            .map { bam, bai ->
                tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam, bai, reference_file)
            }

        // Step 3: Run modkit pileup to get per-site modification frequencies
        modkit_pileup(bam_bai_ch)

        // Step 4: Filter for pseudouridine sites
        filter_pseU(modkit_pileup.out)

        // Step 5: Extract read-level modification calls
        modkit_extract(
            samtools_sort.out
                .map { bam, bai ->
                    tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam)
                }
        )

        // Step 6: Generate modification summary statistics
        modkit_summary(
            samtools_sort.out
                .map { bam, bai ->
                    tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam)
                }
        )
    } else {
        // Standard workflow with basecalling
        // Create sample channel
        samples_ch = create_sample_channel()

        // Step 1: Basecalling with pseudouridine modification detection
        basecalling_pseU(samples_ch)

        // Step 2: Align BAM to reference with Q score filtering using dorado
        dorado_align(
            basecalling_pseU.out
                .map { bam ->
                    tuple(params.mm2opts, bam.getSimpleName(), bam, reference_file)
                }
        )

        // Step 3: Sort BAM and create index
        samtools_sort(dorado_align.out)

        // Create channel with sample ID, BAM, and BAI files
        bam_bai_ch = samtools_sort.out
            .map { bam, bai ->
                tuple(bam.getSimpleName().replaceAll(/\.srt$/, ''), bam, bai, reference_file)
            }

        // Step 4: Run modkit pileup to get per-site modification frequencies
        modkit_pileup(bam_bai_ch)

        // Step 5: Filter for pseudouridine sites
        filter_pseU(modkit_pileup.out)

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
}

workflow.onComplete {
    def resultsMessage = params.skip_basecalling ?
        """
    - Aligned BAMs (Q>${params.min_qscore}): ${params.output_dir}/alignment/
    - Sorted BAMs: ${params.output_dir}/bam/
    - Modification pileup: ${params.output_dir}/modkit_pileup/
    - Pseudouridine sites: ${params.output_dir}/pseU_sites/
    - Modification calls: ${params.output_dir}/modkit_extract/
    - Summary statistics: ${params.output_dir}/modkit_summary/
        """ :
        """
    - Basecalled BAMs: ${params.output_dir}/basecalling_pseU/
    - Aligned BAMs (Q>${params.min_qscore}): ${params.output_dir}/alignment/
    - Sorted BAMs: ${params.output_dir}/bam/
    - Modification pileup: ${params.output_dir}/modkit_pileup/
    - Pseudouridine sites: ${params.output_dir}/pseU_sites/
    - Modification calls: ${params.output_dir}/modkit_extract/
    - Summary statistics: ${params.output_dir}/modkit_summary/
        """

    log.info"""
    ========================================
    Pipeline execution completed!
    ========================================
    Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration: ${workflow.duration}
    Output directory: ${params.output_dir}
    Mode: ${params.skip_basecalling ? 'Skipped basecalling (started from BAM files)' : 'Full pipeline (basecalling included)'}
    Q score filtering: ${params.min_qscore}
    Alignment: Dorado aligner with ${params.mm2opts_sr} options

    Results:${resultsMessage}
    """.stripIndent()
}
