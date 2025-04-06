process align {
    input:
    val sampleId
    val mm2opts
    path fasta_file
    path fastq_file

    output:
    path '*.bam'
    path '*.bam.bai'

    script:
    """
    dorado aligner -t ${cpus} ${mm2opts} ${fasta_file} ${fastq_file} > ${sampleId}.sam | samtools sort --threads ${cpus} > ${sampleId}.srt.bam 
    samtools index ${sampleId}.srt.bam
    """
}