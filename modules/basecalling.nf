process basecalling_dna {
    input:
    val sampleId
    path pod5_folder

    output:
    path '*.bam'
    path '*.fastq'

    script:
    """
    dorado basecaller hac  ${pod5_folder} > ${sampleId}.bam
    samtools fastq ${sampleId}.bam > ${sampleId}.fastq
    """
}