process mod_basecalling_rna {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_rna', mode: 'copy'
    input:
    val sampleId
    path pod5_folder

    output:
    path '*.bam'

    script:
    """
    dorado basecaller -r \
    --device "cuda:all" \
    --emit-moves \
    sup,m5C,m6A_DRACH,pseU \
    ${pod5_folder} > ${sampleId}.bam
    """
}

process basecalling_rna {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_rna', mode: 'copy'
    input:
    tuple val(sampleId), path(pod5_folder)

    output:
    path '*.bam'

    script:
    """
    dorado basecaller -r \
    --device "cuda:all" \
    --emit-moves \
    hac \
    ${pod5_folder} > ${sampleId}.bam
    """
}

process basecalling_dna {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
    publishDir 'results/basecalling_dna', mode: 'copy'
    input:
    val sampleId
    path pod5_folder

    output:
    path '*.bam'

    script:
    """
    dorado basecaller -r \
    --device "cuda:all" \
    --emit-moves \
    hac \
    ${pod5_folder} > ${sampleId}.bam
    """
}