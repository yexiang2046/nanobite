// Basecalling processes for Nanopore sequencing data using Dorado

process mod_basecalling_rna {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_rna', mode: 'copy'

    input:
    tuple val(sample_id), path(pod5_folder)

    output:
    path '*.bam'

    script:
    """
    dorado basecaller -r \
        --device "cuda:all" \
        --emit-moves \
        --estimate-poly-a \
        sup,m5C,m6A_DRACH,pseU \
        ${pod5_folder} > ${sample_id}.bam
    """
}

process basecalling_rna {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_rna', mode: 'copy'

    input:
    tuple val(sample_id), path(pod5_folder)

    output:
    path '*.bam'

    script:
    """
    dorado basecaller -r \
        --device "cuda:all" \
        --emit-moves \
        hac \
        ${pod5_folder} > ${sample_id}.bam
    """
}

process basecalling_dna {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_dna', mode: 'copy'

    input:
    tuple val(sample_id), path(pod5_folder)

    output:
    path '*.bam'

    script:
    """
    dorado basecaller -r \
        --device "cuda:all" \
        --emit-moves \
        hac \
        ${pod5_folder} > ${sample_id}.bam
    """
}

process basecalling_small_rna {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_small_rna', mode: 'copy'

    input:
    tuple val(sample_id), path(pod5_folder)

    output:
    path '*.bam'

    script:
    """
    dorado basecaller -r \
        --device "cuda:all" \
        --emit-moves \
        hac \
        ${pod5_folder} > ${sample_id}.bam
    """
}

process basecalling_small_rna_fastq {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_small_rna_fastq', mode: 'copy'

    input:
    tuple val(sample_id), path(pod5_folder)

    output:
    tuple val(sample_id), path('*.fastq.gz')

    script:
    """
    dorado basecaller -r \
        --device "cuda:all" \
        --emit-moves \
        hac \
        ${pod5_folder} | samtools fastq -T "*" - | gzip > ${sample_id}.fastq.gz
    """
}

process basecalling_pseU {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_pseU', mode: 'copy'

    input:
    tuple val(sample_id), path(pod5_folder)

    output:
    path '*.bam'

    script:
    """
    dorado basecaller -r \
        --device "cuda:all" \
        --emit-moves \
        --estimate-poly-a \
        sup,pseU \
        ${pod5_folder} > ${sample_id}.bam
    """
}