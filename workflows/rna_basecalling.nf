nextflow.enable.dsl=2

// Import modules
include { basecalling_rna } from '../modules/basecalling.nf'

// Workflow definition
workflow RNA_BASECALLING {
    take:
        sample_id
        pod5_dir

    main:
        // Validate input
        if (!pod5_dir) {
            error "Input directory is required"
        }

        // Create channel from input directory
        pod5_ch = channel.fromPath(pod5_dir)

        // Run basecalling without modifications
        basecalling_rna(
            sample_id,
            pod5_ch,
        )

        align(
            sample_id,
            basecalling_rna.out.fastq,
        )

    emit:
        fastq = basecalling_rna.out.fastq
        bam = align.out.bam
}

// Entry point
workflow {
    RNA_BASECALLING(
        params.sample_id,
        params.input_dir,
        params.model
    )
} 