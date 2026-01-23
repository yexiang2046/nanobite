process CHOPPER_QC {
    publishDir "${params.output_dir}/chopper_qc", mode: 'copy'
    container "quay.io/biocontainers/chopper:1.1.0--h9f5ecd1_0"
    time '30m'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    path("${sample_id}_filtered.fastq"), emit: filtered_fastq

    script:
    """
    chopper filter \
        --input ${fastq_file} \
        --output ${sample_id}_filtered.fastq \
        --minlength 1000 \
        --minqual 7 \
        --maxqual 40
    """
}