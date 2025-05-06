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
    
    # Copy and make the R script executable
    cp ${merge_script} merge_counts.R
    chmod +x merge_counts.R
    
    # Run the merge script
    ./merge_counts.R "${sample_info}" "${count_dir}"
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
    
    # Copy and make the R script executable
    cp ${run_script} run_deseq2.R
    chmod +x run_deseq2.R
    
    # Run DESeq2 analysis
    ./run_deseq2.R "${count_matrix}" "${design_matrix}"
    """
} 