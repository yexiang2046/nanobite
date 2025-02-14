process align_dna {
    input:
    val sampleId
    path fasta_file
    path bam_file

    output:
    path '*.bam'
    path '*.bam.bai'

    script:
    """
    dorado aligner ${fasta_file} ${bam_file} > ${sampleId}.sam | samtools sort --threads <num_threads> > ${sampleId}.srt.bam 
    samtools index ${sampleId}.srt.bam
    """
}