process align {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    cpus 4  // Default to 4 CPUs, can be overridden in config
    publishDir "${params.output_dir}/alignment", mode: 'copy'

    input:
    tuple val(mm2opts), val(sampleId), path(bam_file), path(fasta_file)

    output:
    path "${sampleId}.sam"

    script:
    """
    dorado aligner -t ${task.cpus} ${mm2opts} ${fasta_file} ${bam_file} > ${sampleId}.sam
    """
}

process process_sam {
    container 'staphb/samtools:1.21'
    cpus 4  // Default to 4 CPUs, can be overridden in config
    publishDir "${params.output_dir}/bam", mode: 'copy'

    input:
    path sam_file

    output:
    path "${sam_file.getSimpleName()}.srt.bam"
    path "${sam_file.getSimpleName()}.srt.bam.bai"

    script:
    """
    samtools sort --threads ${task.cpus} ${sam_file} > ${sam_file.getSimpleName()}.srt.bam
    samtools index ${sam_file.getSimpleName()}.srt.bam
    """
}