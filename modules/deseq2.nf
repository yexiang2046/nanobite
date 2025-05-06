// Prepare input files for DESeq2 analysis
process prepare_deseq2 {
    container "${params.container}"
    publishDir "${params.output_dir}/input", mode: 'copy'
    
    input:
    path count_dir
    path sample_info
    path merge_script
    
    output:
    path "count_matrix.txt", emit: count_matrix
    path "design_matrix.txt", emit: design_matrix
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Make the R script executable
    chmod +x ${merge_script}
    
    # Run the merge script
    ${merge_script} "${sample_info}" "${count_dir}"
    """
}

// Run DESeq2 analysis
process deseq2_analysis {
    container "${params.container}"
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path count_matrix
    path design_matrix
    path run_script
    
    output:
    path "results/*.tsv", emit: results
    path "plots/*.pdf", emit: plots
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Make the R script executable
    chmod +x ${run_script}
    
    # Run DESeq2 analysis
    ${run_script} "${count_matrix}" "${design_matrix}"
    """
} 