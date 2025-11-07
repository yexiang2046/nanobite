process modkit_pileup {
    container 'xiang2019/modkit:v1.0.0'
    cpus 8
    publishDir "${params.output_dir}/modkit_pileup", mode: 'copy'

    input:
    tuple val(sampleId), path(bam_file), path(bai_file), path(reference)

    output:
    tuple val(sampleId), path("${sampleId}_pileup.bed")

    script:
    def min_coverage = params.min_coverage ?: 5
    def prob_threshold = params.prob_threshold ?: 0.8
    """
    modkit pileup ${bam_file} ${sampleId}_pileup.bed \
        --ref ${reference} \
        --cpus ${task.cpus} \
        --log-filepath ${sampleId}_modkit.log \
        --filter-threshold ${prob_threshold} \
        --interval-size 1000 \
        --min-coverage ${min_coverage}
    """
}

process modkit_extract {
    container 'xiang2019/modkit:v1.0.0'
    cpus 4
    publishDir "${params.output_dir}/modkit_extract", mode: 'copy'

    input:
    tuple val(sampleId), path(bam_file)

    output:
    tuple val(sampleId), path("${sampleId}_mod_calls.tsv")

    script:
    """
    modkit extract ${bam_file} ${sampleId}_mod_calls.tsv \
        --cpus ${task.cpus} \
        --log-filepath ${sampleId}_extract.log
    """
}

process modkit_summary {
    container 'xiang2019/modkit:v1.0.0'
    cpus 2
    publishDir "${params.output_dir}/modkit_summary", mode: 'copy'

    input:
    tuple val(sampleId), path(bam_file)

    output:
    tuple val(sampleId), path("${sampleId}_summary.txt")

    script:
    """
    modkit summary ${bam_file} > ${sampleId}_summary.txt
    """
}

process filter_pseU {
    container 'xiang2019/modkit:v1.0.0'
    publishDir "${params.output_dir}/pseU_sites", mode: 'copy'

    input:
    tuple val(sampleId), path(pileup_bed)

    output:
    tuple val(sampleId), path("${sampleId}_pseU_sites.bed")

    script:
    """
    # Filter for pseudouridine (ψ) modifications
    # Modkit uses 'Y' for pseudouridine in the modification code
    grep -E "\\sY\\s|pseU" ${pileup_bed} > ${sampleId}_pseU_sites.bed || touch ${sampleId}_pseU_sites.bed
    """
}
