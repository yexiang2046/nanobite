// Prepare input files for DESeq2 analysis
process prepare_deseq2 {
    publishDir "${params.output_dir}/deseq2/input", mode: 'copy'
    container "${params.container}"
    
    input:
    path count_files
    path sample_info
    
    output:
    path("count_matrix.txt"), emit: count_matrix
    path("design_matrix.txt"), emit: design_matrix
    
    script:
    """
    Rscript ${projectDir}/bin/prepare_deseq2.R \
        --count_files "${count_files.join(' ')}" \
        --sample_info "${sample_info}" \
        --count_matrix count_matrix.txt \
        --design_matrix design_matrix.txt
    """
}

// Run DESeq2 analysis
process deseq2_analysis {
    publishDir "${params.output_dir}/deseq2", mode: 'copy'
    container "${params.container}"
    
    input:
    path count_matrix
    path design_matrix
    
    output:
    path("results/"), emit: results
    path("plots/"), emit: plots
    
    script:
    """
    mkdir -p results plots
    
    Rscript ${projectDir}/bin/run_deseq2.R \
        --count_matrix ${count_matrix} \
        --design_matrix ${design_matrix} \
        --output_dir results \
        --plot_dir plots
    """
} 