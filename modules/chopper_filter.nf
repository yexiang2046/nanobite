process CHOPPER_QC {
    publishDir "${params.output_dir}/chopper_qc", mode: 'copy'
    container "xiang2019/chopper:v0.12.0b"
    time '30m'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    path("${sample_id}_chopper.fastq"), emit: filtered_fastq

    script:
    def cutoff = params.chopper_cutoff ?: 10
    """
    chopper --input ${fastq_file} --output ${sample_id}_chopper.fastq --trim-approach best-read-segment --cutoff ${cutoff}
    """
}