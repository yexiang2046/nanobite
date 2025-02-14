process polish_assembly {
    input:
    val sampleId
    path fasta_file 
    path srt_bam_file
    
    output:
    path '*.fasta'

    script:
    """
    dorado polish ${srt_bam_file} ${fasta_file} > ${sampleId}.polished.fasta
    """
}   