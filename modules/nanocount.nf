process nanocount {
    container 'biocorecrg/nanocount:1.0.0.post6'
    publishDir "${params.output_dir}/transcript_counts", mode: 'copy'

    input:
    tuple path(bam_file), val(sample_id)

    output:
    path "${sample_id}_transcript_counts.tsv", emit: counts

    script:
    """
    NanoCount -i ${bam_file} -o ${sample_id}_transcript_counts.tsv
    """
} 