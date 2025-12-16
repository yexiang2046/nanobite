#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { RNA_MOD_M6A_DRACH } from '../modules/rna_modbasecall.nf'
include { RNA_MOD_M5C } from '../modules/rna_modbasecall.nf'
include { RNA_MOD_PSEU } from '../modules/rna_modbasecall.nf'
include { RNA_MOD_INOSINE } from '../modules/rna_modbasecall.nf'
include { RNA_MOD_C_2OME } from '../modules/rna_modbasecall.nf'
include { RNA_MOD_A_2OME } from '../modules/rna_modbasecall.nf'
include { RNA_MOD_U_2OME } from '../modules/rna_modbasecall.nf'
include { RNA_MOD_G_2OME } from '../modules/rna_modbasecall.nf'
include { modkit_combine } from '../modules/modkit.nf'
include { align } from '../modules/align.nf'
include { nanocount } from '../modules/nanocount.nf'

// Default parameters
params.help = false
params.reference = null
params.sample_info = null
params.output_dir = "results"
params.model_sup = "rna004_130bps_sup@v5.2.0"
params.model_hac = "rna004_130bps_hac@v5.2.0"
params.gpu_device = "cuda:0"
params.use_gpu = true
params.modkit_threads = 16
params.modkit_threshold = 0.0 // Minimum probability threshold for modifications

// Options for minimap2 alignment on short reads 
params.mm2opts = "-ax sr" // Used for snRNA profiling

// Print help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run sncRNA_basecalling_align.nf --reference reference.fa --sample_info sample_info.txt

    Required Arguments:
        --reference         Reference genome FASTA file
        --sample_info      Tab-delimited sample info file with columns: SampleId, pod5_path

    Optional Arguments:
        --modkit_threads    Number of threads for modkit (default: ${params.modkit_threads})
        --modkit_threshold  Minimum probability threshold for modkit modifications (default: ${params.modkit_threshold})
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

// Create sample channel from sample info file
Channel
    .fromPath(params.sample_info)
    .splitCsv(header: true, sep: '\t')
    .map { row -> 
        // Print row contents for debugging
        println "Processing row: ${row}"
        
        if (!row.containsKey('SampleId') || !row.containsKey('pod5_path')) {
            error "Sample info file must contain 'SampleId' and 'pod5_path' columns. Found columns: ${row.keySet()}"
        }
        
        def sample_id = row.SampleId?.trim()
        def pod5_path = file(row.pod5_path?.trim())
        
        if (!sample_id) {
            error "SampleId cannot be empty"
        }
        if (!pod5_path.exists()) {
            error "POD5 file not found: ${pod5_path}"
        }
        
        println "Created tuple: [${sample_id}, ${pod5_path}]"
        return tuple(sample_id, pod5_path)
    }
    .set { samples_ch }



// Main workflow
workflow {
    // Run basecalling
    RNA_MOD_M6A_DRACH(samples_ch)
    RNA_MOD_M5C(samples_ch)
    RNA_MOD_PSEU(samples_ch)
    RNA_MOD_INOSINE(samples_ch)
    RNA_MOD_C_2OME(samples_ch)
    RNA_MOD_A_2OME(samples_ch)
    RNA_MOD_U_2OME(samples_ch)
    RNA_MOD_G_2OME(samples_ch)
    
    // Combine basecalling outputs by sample_id
    all_basecalled_bams = RNA_MOD_M6A_DRACH.out
        .mix(RNA_MOD_M5C.out)
        .mix(RNA_MOD_PSEU.out)
        .mix(RNA_MOD_INOSINE.out)
        .mix(RNA_MOD_C_2OME.out)
        .mix(RNA_MOD_A_2OME.out)
        .mix(RNA_MOD_U_2OME.out)
        .mix(RNA_MOD_G_2OME.out)
        .groupTuple(by: 0) // Group BAMs by sample_id
        .map { sample_id, bams -> tuple(sample_id, bams, file(params.reference)) } // Add reference
    
    // Run modkit to combine modification information
    modkit_combine(all_basecalled_bams)
    
    // Run alignment on original BAM files
    align(
        modkit_combine.out.bam_files
            .flatMap { sample_id, bams ->
                bams.collect { bam -> tuple(params.mm2opts, sample_id, bam, file(params.reference)) }
            }
    )
    
    // Run nanocount on alignment output
    nanocount(
        align.out
            .map { bam ->
                tuple(bam, bam.getSimpleName())
            }
    )
}