process align {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    cpus 4  // Default to 4 CPUs, can be overridden in config
    publishDir "${params.output_dir}/alignment", mode: 'copy'

    input:
    tuple val(mm2opts), val(sampleId), path(bam_file), path(fasta_file)

    output:
    path "${sampleId}_aligned.bam"

    script:
    """
    dorado aligner -t ${task.cpus} ${mm2opts} ${fasta_file} ${bam_file} > ${sampleId}_aligned.bam
    """
}

process process_sam {
    container 'staphb/samtools:1.21'
    cpus 4  // Default to 4 CPUs, can be overridden in config
    publishDir "${params.output_dir}/bam", mode: 'copy'

    input:
    path bam_file

    output:
    path "${bam_file.getSimpleName()}.srt.bam"
    path "${bam_file.getSimpleName()}.srt.bam.bai"

    script:
    """
    samtools sort --threads ${task.cpus} ${bam_file} > ${bam_file.getSimpleName()}.srt.bam
    samtools index ${bam_file.getSimpleName()}.srt.bam
    """
}