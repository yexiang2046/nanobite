nextflow.enable.dsl=2

// Import modules
include { prepare_deseq2, deseq2_analysis } from '../modules/deseq2.nf'

// Create a channel from sample info file
def create_sample_channel(sample_info) {
    Channel
        .fromPath(sample_info)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            // Validate required columns
            if (!row.containsKey('SampleId') || !row.containsKey('treatment')) {
                error "Sample info file must contain columns: SampleId and treatment"
            }
            
            def sample_id = row.SampleId?.trim()
            def treatment = row.treatment?.trim()
            
            // Validate values
            if (!sample_id) {
                error "SampleId cannot be empty in sample info file"
            }
            if (!treatment) {
                error "Treatment cannot be empty for sample ${sample_id}"
            }
            
            return tuple(sample_id, treatment)
        }
}

// Workflow definition
workflow DIFF_EXPRESSION {
    take:
        sample_info  // path to sample info file
        count_files  // path to count files directory

    main:
        // Create sample channel
        sample_ch = create_sample_channel(sample_info)

        // Prepare DESeq2 input files
        prepare_deseq2(
            Channel.fromPath("${count_files}/*_counts.txt").collect(),
            file(sample_info)
        )

        // Run DESeq2 analysis
        deseq2_analysis(
            prepare_deseq2.out.count_matrix,
            prepare_deseq2.out.design_matrix
        )

    emit:
        results = deseq2_analysis.out.results
        plots = deseq2_analysis.out.plots
}

// Entry point
workflow {
    // Parameter validation
    if (!params.sample_info) {
        error "Sample info file is required: --sample_info sample_info.txt"
    }
    if (!params.count_files) {
        error "Count files directory is required: --count_files path/to/counts"
    }

    // Set default parameters
    params.container = params.container ?: 'xiang2019/deseq2:v1.0.0'
    params.output_dir = params.output_dir ?: 'results/deseq2'

    DIFF_EXPRESSION(
        file(params.sample_info),
        params.count_files
    )
} 