#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { basecalling_rna } from '../modules/basecalling.nf'
include { align } from '../modules/align.nf'

// Default parameters
params.help = false
params.reference = null
params.sample_info = null
params.output_dir = "results"
params.model = "rna004_130bps_hac@v5.1.0"
params.modified_bases = "m,a,u"
params.batch_size = 100
params.chunks_per_runner = 208
params.gpu_device = "cuda:0"
params.use_gpu = true

// Print help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run rna_mod_basecalling.nf --reference reference.fa --sample_info sample_info.txt

    Required Arguments:
        --reference         Reference genome FASTA file
        --sample_info      Tab-delimited sample info file with columns: sample_id, pod5_path

    Optional Arguments:
        --output_dir       Output directory (default: results)
        --model           Dorado model (default: rna004_130bps_hac@v5.1.0)
        --modified_bases  Modifications to detect (default: m,a,u)
        --batch_size     Batch size for processing (default: 100)
        --gpu_device     GPU device to use (default: cuda:0)
        --use_gpu        Use GPU for processing (default: true)
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
        def sample_id = row.sample_id
        def pod5_path = file(row.pod5_path)
        if (!pod5_path.exists()) {
            error "POD5 file not found: ${pod5_path}"
        }
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
        basecalling_rna.out
            .map { bam ->
                tuple(params.mm2opts, bam.getSimpleName(), bam, file(params.reference))
            }
    )
} 