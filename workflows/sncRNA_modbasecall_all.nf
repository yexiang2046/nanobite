#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { RNA_MOD_ALL } from '../modules/rna_modbasecall.nf'
include { samtools_fastq } from '../modules/align.nf'
include { nanofilt_filter } from '../modules/align.nf'
include { minimap2_align_sr } from '../modules/align.nf'
include { sam_to_bam } from '../modules/align.nf'
include { samtools_sort } from '../modules/align.nf'
include { nanocount } from '../modules/nanocount.nf'

// Default parameters
params.help = false
params.reference = null
params.sample_info = null
params.output_dir = "results"
params.model_sup = "rna004_130bps_sup@v5.2.0"
params.mm2opts_sr = "-ax sr"
params.min_qscore = 10

// Print help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run sncRNA_modbasecall_all.nf --reference reference.fa --sample_info sample_info.txt

    Required Arguments:
        --reference         Reference genome FASTA file
        --sample_info      Tab-delimited sample info file with columns: SampleId, pod5_path

    Optional Arguments:
        --output_dir       Output directory (default: "results")
        --model_sup        Dorado basecalling model (default: "rna004_130bps_sup@v5.2.0")
        --mm2opts_sr       Minimap2 short read alignment options (default: "-ax sr")
        --min_qscore       Minimum Q score for filtering (default: 10)
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate inputs
if (!params.reference || !params.sample_info) {
    error "Both --reference and --sample_info are required"
}

// Validate reference file exists
def ref_file = file(params.reference)
if (!ref_file.exists()) {
    error "Reference file does not exist: ${params.reference}"
}

// Create sample channel from sample info file
def create_sample_channel() {
    Channel.fromPath(params.sample_info)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            // Validate required columns
            if (!row.containsKey('SampleId') || !row.containsKey('pod5_path')) {
                error "Sample info file must contain 'SampleId' and 'pod5_path' columns. Found columns: ${row.keySet()}"
            }

            def sample_id = row.SampleId?.trim()
            def pod5_path = file(row.pod5_path?.trim())

            // Validate pod5 path exists
            if (!pod5_path.exists()) {
                error "POD5 path does not exist for sample ${sample_id}: ${pod5_path}"
            }

            return tuple(sample_id, pod5_path)
        }
}

// Main workflow
workflow {
    // Create sample channel
    samples_ch = create_sample_channel()

    // Step 1: Basecalling with all modifications
    RNA_MOD_ALL(samples_ch)

    // Step 2: Convert BAM to FASTQ for alignment with minimap2
    // RNA_MOD_ALL outputs: tuple val(sample_id), path("${sample_id}_all_mods.bam")
    samtools_fastq(RNA_MOD_ALL.out)

    // Step 3: Filter reads by quality score
    nanofilt_filter(samtools_fastq.out)

    // Step 4: Align with minimap2
    minimap2_align_sr(
        nanofilt_filter.out
            .map { sample_id, fastq ->
                tuple(sample_id, fastq, ref_file)
            }
    )

    // Step 5: Convert SAM to BAM
    sam_to_bam(minimap2_align_sr.out)

    // Step 6: Sort and index BAM files
    samtools_sort(sam_to_bam.out)

    // Step 7: Quantify with nanocount
    // samtools_sort outputs: tuple path("${bam_file.getSimpleName()}.srt.bam"), path("${bam_file.getSimpleName()}.srt.bam.bai")
    // nanocount expects: tuple path(bam_file), val(sample_id)
    sorted_bam_for_nanocount = samtools_sort.out
        .map { bam, bai ->
            def sample_id = bam.getSimpleName().replaceAll('_aligned\\.srt', '').replaceAll('_all_mods', '')
            tuple(bam, sample_id)
        }

    nanocount(sorted_bam_for_nanocount)
}