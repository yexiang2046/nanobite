process NANOFILT_QC {
    publishDir "${params.output_dir}/nanofilt_qc", mode: 'copy'
    container "nanozoo/nanofilt:2.8.0--f3f0c60"  // Available container; or build one with Python/Biopython
    time '30m'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    path("${sample_id}_nanofilt.fastq"), emit: filtered_fastq

    script:
    def min_qual = params.nanofilt_min_qual ?: 10
    def min_len = params.nanofilt_min_len ?: 500
    def headcrop = params.nanofilt_headcrop ?: 0
    def tailcrop = params.nanofilt_tailcrop ?: 0
    """
    gunzip -c ${fastq_file} | NanoFilt -q ${min_qual} -l ${min_len} --headcrop ${headcrop} --tailcrop ${tailcrop} > ${sample_id}_nanofilt.fastq
    """
}