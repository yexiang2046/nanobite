// Common utility functions for NanoBite pipeline

/**
 * Creates a sample channel from a sample info file
 *
 * @param sample_info_file Path to tab-separated file with columns: SampleId, treatment, pod5_path (or bam_path)
 * @param skip_basecalling If true, expects 'bam_path' column instead of 'pod5_path'
 * @param use_fastq If true, expects 'pod5_path' column for FASTQ basecalling
 * @return Channel of tuples [sample_id, treatment, file_path]
 */
def create_sample_channel(sample_info_file, skip_basecalling = false, use_fastq = false) {
    def sample_data = Channel
        .fromPath(sample_info_file)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            // Validate required columns
            if (!row.containsKey('SampleId')) {
                error "Sample info file must contain 'SampleId' column"
            }
            if (!row.containsKey('treatment')) {
                error "Sample info file must contain 'treatment' column"
            }

            def sample_id = row.SampleId?.trim()
            def treatment = row.treatment?.trim()

            // Validate common values
            if (!sample_id) {
                error "SampleId cannot be empty in sample info file"
            }
            if (!treatment) {
                error "Treatment cannot be empty for sample ${sample_id}"
            }

            if (skip_basecalling) {
                // When skipping basecalling, expect bam_path column
                if (!row.containsKey('bam_path')) {
                    error "Sample info file must contain 'bam_path' column when skip_basecalling is enabled"
                }

                def bam_path = row.bam_path?.trim()
                if (!bam_path) {
                    error "BAM path cannot be empty for sample ${sample_id}"
                }

                def bam_file = file(bam_path)
                if (!bam_file.exists()) {
                    error "BAM file not found for sample ${sample_id}: ${bam_file}"
                }

                return tuple(sample_id, treatment, bam_file)
            } else {
                // When running basecalling, expect pod5_path column
                if (!row.containsKey('pod5_path')) {
                    error "Sample info file must contain 'pod5_path' column when basecalling is enabled"
                }

                def pod5_path = row.pod5_path?.trim()
                if (!pod5_path) {
                    error "Pod5 path cannot be empty for sample ${sample_id}"
                }

                def pod5_dir = file(pod5_path)
                if (!pod5_dir.exists()) {
                    error "Pod5 directory not found for sample ${sample_id}: ${pod5_dir}"
                }

                return tuple(sample_id, treatment, pod5_dir)
            }
        }

    return sample_data
}

/**
 * Extracts sample ID from a BAM file name
 * Handles various naming patterns used in the pipeline
 *
 * @param bam_file File object or path
 * @return Sample ID string
 */
def extract_sample_id(bam_file) {
    def base_name = bam_file.getSimpleName()

    // Remove common suffixes
    base_name = base_name
        .replace('_aligned.srt', '')
        .replace('_aligned', '')
        .replace('.srt', '')
        .replace('.sorted', '')

    return base_name
}

/**
 * Creates a BAM channel with sample IDs from BAM files
 *
 * @param bam_files Channel of BAM files
 * @return Channel of tuples [sample_id, bam_file]
 */
def create_bam_channel(bam_files) {
    return bam_files.map { bam ->
        def sample_id = extract_sample_id(bam)
        tuple(sample_id, bam)
    }
}

/**
 * Creates a BAM + BAI index channel with sample IDs
 *
 * @param bam_bai_pairs Channel of [bam, bai] tuples
 * @return Channel of tuples [sample_id, bam, bai]
 */
def create_bam_bai_channel(bam_bai_pairs) {
    return bam_bai_pairs.map { bam, bai ->
        def sample_id = extract_sample_id(bam)
        tuple(sample_id, bam, bai)
    }
}
