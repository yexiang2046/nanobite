process align {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    input:
    tuple val(mm2opts), val(sampleId), path(bam_file), path(fasta_file)

    output:
    path '*.bam'
    path '*.bam.bai'

    script:
    """
    dorado aligner -t ${cpus} ${mm2opts} ${fasta_file} ${bam_file} > ${sampleId}.sam | samtools sort --threads ${cpus} > ${sampleId}.srt.bam 
    samtools index ${sampleId}.srt.bam
    """
}