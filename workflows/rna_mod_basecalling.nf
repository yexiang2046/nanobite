nextflow.enable.dsl=2

// Import modules
include { mod_basecalling_rna } from '../modules/basecalling.nf'

// Workflow definition
workflow RNA_MOD_BASECALLING {
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

        // Run basecalling with modifications
        mod_basecalling_rna(
            sample_id,
            pod5_ch
        )

    emit:
        fastq = mod_basecalling_rna.out.fastq
}

// Entry point
workflow {
    RNA_MOD_BASECALLING(
        params.sample_id,
        params.input_dir
    )
} 