process mod_basecalling_rna {
    container 'staphb/dorado:0.9.0-cuda12.2.0'
    containerOptions '--gpus all'
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