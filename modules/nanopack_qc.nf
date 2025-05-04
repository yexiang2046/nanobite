// BAM file QC using NanoPlot
process nanopack_plot {
    publishDir "${params.output_dir}/nanopack/plots", mode: 'copy'
    container "staphb/nanoplot:1.42.0"
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    path("${sample_id}_nanoplot"), emit: plots
    
    script:
    """
    NanoPlot \
        --bam ${bam} \
        --outdir ${sample_id}_nanoplot \
        --threads ${task.cpus} \
        --plots dot kde hex
    """
}

// Compare multiple BAM files using NanoComp
process nanopack_compare {
    publishDir "${params.output_dir}/nanopack/comparison", mode: 'copy'
    container "bikc/nanocomp:1.15.1"
    
    input:
    path bams
    
    output:
    path("nanocomp_results"), emit: comparison
    
    script:
    """
    NanoComp \
        --bam ${bams} \
        --outdir nanocomp_results \
        --threads ${task.cpus} \
        --plots dot kde hex \
        --minlength 0 \
        --maxlength 100000 \
        --minqual 0 \
        --maxqual 100
    """
}

// Generate BAM statistics using Cramino
process nanopack_stats {
    publishDir "${params.output_dir}/nanopack/stats", mode: 'copy'
    container "alexanrna/cramino:0.14.5"
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    path("${sample_id}_stats.txt"), emit: stats
    
    script:
    """
    cramino ${bam} > ${sample_id}_stats.txt
    """
}

// Generate read phasing visualization using Phasius
process nanopack_phasing {
    publishDir "${params.output_dir}/nanopack/phasing", mode: 'copy'
    container "xiang2019/phasius:v1.0.0"
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    path("${sample_id}_phasing.html"), emit: phasing
    
    script:
    """
    phasius \
        --output ${sample_id}_phasing.html \
        --region all \
        ${bam}
    """
}

// Generate BAM overview using Kyber
process nanopack_overview {
    publishDir "${params.output_dir}/nanopack/overview", mode: 'copy'
    container "alexanrna/kyber:0.1.0"
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    path("${sample_id}_overview"), emit: overview
    
    script:
    """
    kyber -o ${sample_id}_overview ${bam} 
    """
} 