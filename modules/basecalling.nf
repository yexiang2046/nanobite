// Basecalling processes for Nanopore sequencing data using Dorado

process mod_basecalling_rna {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
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
        --min-qscore ${params.min_qscore} \
        sup,m5C,m6A_DRACH,pseU \
        ${pod5_folder} > ${sample_id}_q${params.min_qscore}.bam
    """
}

process basecalling_rna {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
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
        --min-qscore ${params.min_qscore} \
        hac \
        ${pod5_folder} > ${sample_id}_q${params.min_qscore}.bam
    """
}

process basecalling_dna {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
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
        --min-qscore ${params.min_qscore} \
        hac \
        ${pod5_folder} > ${sample_id}_q${params.min_qscore}.bam
    """
}

process basecalling_small_rna {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
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
        --min-qscore ${params.min_qscore} \
        hac \
        ${pod5_folder} > ${sample_id}_q${params.min_qscore}.bam
    """
}

process basecalling_small_rna_fastq {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
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
        --min-qscore ${params.min_qscore} \
        hac \
        ${pod5_folder} | samtools fastq -T "*" - | gzip > ${sample_id}_q${params.min_qscore}.fastq.gz
    """
}

process basecalling_pseU {
    container 'ontresearch/dorado:shae423e761540b9d08b526a1eb32faf498f32e8f22'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_pseU', mode: 'copy'

    input:
    tuple val(sample_id), path(pod5_folder)

    output:
    path "${sample_id}_q${params.min_qscore}.bam"

    script:
    def model = params.model ?: "rna004_130bps_sup@v5.1.0"
    """
    dorado download --model ${model}
    dorado basecaller -r \
        --device "cuda:all" \
        --emit-moves \
        --estimate-poly-a \
        --modified-bases pseU \
        --min-qscore ${params.min_qscore} \
        ${model} \
        ${pod5_folder} > ${sample_id}_q${params.min_qscore}.bam
    """
}