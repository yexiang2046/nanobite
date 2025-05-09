#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { mod_basecalling_rna } from '../modules/basecalling.nf'
include { align } from '../modules/align.nf'

// Default parameters
params.help = false
params.reference = null
params.sample_info = null
params.output_dir = "results"
params.model = "rna004_130bps_hac@v5.1.0"
params.gpu_device = "cuda:0"
params.use_gpu = true

// Print help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run rna_mod_basecalling.nf --reference reference.fa --sample_info sample_info.txt

    Required Arguments:
        --reference         Reference genome FASTA file
        --sample_info      Tab-delimited sample info file with columns: SampleId, pod5_path
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
    mod_basecalling_rna(samples_ch)
    // Run alignment
    // Run alignment for each sample
    align(
        mod_basecalling_rna.out
            .map { bam ->
                tuple(params.mm2opts, bam.getSimpleName(), bam, file(params.reference))
            }
    )
} 