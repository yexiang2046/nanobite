process CHOPPER_QC {
    publishDir "${params.output_dir}/chopper_qc", mode: 'copy'
    container "xiang2019/chopper:v0.12.0b"
    time '30m'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    path("${sample_id}_chopper.fastq"), emit: filtered_fastq

    script:
    """
    cat ${fastq_file} | chopper --trim-approach best-read-segment --cutoff 15 > ${sample_id}_chopper.fastq
    """
}