// Alignment processes for Nanopore sequencing data

process align {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    cpus 4
    publishDir "${params.output_dir}/alignment", mode: 'copy'

    input:
    tuple val(mm2opts), val(sample_id), path(bam_file), path(fasta_file)

    output:
    path "${sample_id}_aligned.bam"

    script:
    """
    dorado aligner -t ${task.cpus} ${mm2opts} ${fasta_file} ${bam_file} > ${sample_id}_aligned.bam
    """
}

process process_sam {
    container 'staphb/samtools:1.21'
    cpus 4
    publishDir "${params.output_dir}/bam", mode: 'copy'

    input:
    path bam_file

    output:
    tuple path("${bam_file.getSimpleName()}.srt.bam"), path("${bam_file.getSimpleName()}.srt.bam.bai")

    script:
    """
    samtools sort --threads ${task.cpus} ${bam_file} > ${bam_file.getSimpleName()}.srt.bam
    samtools index ${bam_file.getSimpleName()}.srt.bam
    """
}

process bam_to_fastq {
    container 'staphb/samtools:1.21'
    cpus 2

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.fastq")

    script:
    """
    samtools fastq ${bam_file} > ${sample_id}.fastq
    """
}

process nanofilt {
    container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'
    cpus 2
    publishDir "${params.output_dir}/filtered", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_file)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz")

    script:
    def min_qscore = params.min_qscore ?: 10
    """
    cat ${fastq_file} | \
        NanoFilt -q ${min_qscore} | \
        gzip > ${sample_id}_filtered.fastq.gz
    """
}

process bwa_index {
    container 'staphb/bwa:0.7.17'
    cpus 2

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*")

    script:
    """
    bwa index ${reference}
    """
}

process bwa_align {
    container 'staphb/bwa:0.7.17'
    cpus 4

    input:
    tuple val(sample_id), path(fastq_file), path(reference), path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}_aligned.sam")

    script:
    """
    bwa mem -C -W 13 -k 6 -x ont2d -t ${task.cpus} ${reference} ${fastq_file} > ${sample_id}_aligned.sam
    """
}

process sam_to_bam {
    container 'staphb/samtools:1.21'
    cpus 2
    publishDir "${params.output_dir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_file)

    output:
    path "${sample_id}_aligned.bam"

    script:
    """
    samtools view -bS ${sam_file} > ${sample_id}_aligned.bam
    """
}

process minimap2_align {
    container 'staphb/minimap2:2.28'
    cpus 4
    publishDir "${params.output_dir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_file), path(reference)

    output:
    path "${sample_id}_aligned.bam"

    script:
    def mm2opts = params.mm2opts ?: "-ax splice -uf -k14"
    """
    minimap2 ${mm2opts} -t ${task.cpus} ${reference} ${fastq_file} | \
        samtools view -bS - > ${sample_id}_aligned.bam
    """
}

process minimap2_align_sr {
    container 'staphb/minimap2:2.28'
    cpus 4
    publishDir "${params.output_dir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_file), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_aligned.sam")

    script:
    def mm2opts = params.mm2opts_sr ?: "-ax sr"
    """
    minimap2 ${mm2opts} -t ${task.cpus} ${reference} ${fastq_file} > ${sample_id}_aligned.sam
    """
}
