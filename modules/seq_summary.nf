process seq_summary {
    input:
    val sampleId
    path bam_file

    output:
    path '*.tsv'

    script:
    """
    dorado summary ${bam_file} > ${sampleId}.summary.tsv
    """
}