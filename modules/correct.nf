process correct_dna {
    input:
    val sampleId
    path fastq_file

    output:
    path '*.fastq'

    script:
    """
    dorado correct ${fastq_file} > ${sampleId}.corrected.fastq
    """
}