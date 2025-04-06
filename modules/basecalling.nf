process mod_basecalling_rna {
    input:
    val sampleId
    path pod5_folder

    output:
    path '*.fastq'

    script:
    """
    dorado basecaller -r \
    --device "cuda:all" \
    --emit-moves \
    --emit-fastq \
    sup,m5C,m6A_DRACH,pseU \
    ${pod5_folder} > ${sampleId}.fastq
    """
}

process basecalling_rna {
    input:
    val sampleId
    path pod5_folder

    output:
    path '*.fastq'

    script:
    """
    dorado basecaller -r \
    --device "cuda:all" \
    --emit-fastq \
    --emit-moves \
    hac \
    ${pod5_folder} > ${sampleId}.fastq
    """
}

process basecalling_dna {
    input:
    val sampleId
    path pod5_folder

    output:
    path '*.fastq'

    script:
    """
    dorado basecaller -r \
    --device "cuda:all" \
    --emit-fastq \
    --emit-moves \
    hac \
    ${pod5_folder} > ${sampleId}.fastq
    """
}