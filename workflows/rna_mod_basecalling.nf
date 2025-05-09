#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Default parameters
params.help = false
params.reference = null
params.sample_info = null
params.output_dir = "results"
params.model = "dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
params.modified_bases = "m,a,u"
params.batch_size = 100
params.chunks_per_runner = 208
params.gpu_device = "cuda:all"

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
        --model           Dorado model (default: dna_r10.4.1_e8.2_400bps_sup@v4.2.0)
        --modified_bases  Modifications to detect (default: m,a,u)
        --batch_size     Batch size for processing (default: 100)
        --gpu_device     GPU device to use (default: cuda:all)
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

// Process for RNA modification basecalling
process MOD_BASECALLING {
    tag "basecalling_${sample_id}"
    
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    
    publishDir "${params.output_dir}/basecalls", mode: 'copy'
    
    input:
    tuple val(sample_id), path(pod5_file)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.bam")
    
    script:
    """
    dorado basecaller \
        ${params.model} \
        ${pod5_file} \
        --reference ${reference} \
        --modified-bases ${params.modified_bases} \
        --batch-size ${params.batch_size} \
        --chunks-per-runner ${params.chunks_per_runner} \
        --device ${params.gpu_device} \
        > ${sample_id}.bam
    """
}

// Main workflow
workflow {
    // Reference genome channel
    reference_ch = Channel.fromPath(params.reference)
    
    // Run basecalling
    MOD_BASECALLING(samples_ch, reference_ch)
} 