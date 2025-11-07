nextflow.enable.dsl=2

// Import modules
include { basecalling_rna } from '../modules/basecalling.nf'
include { align } from '../modules/align.nf'
include { process_sam } from '../modules/align.nf'
include { nanocount } from '../modules/nanocount.nf'
include {
    nanopack_plot;
    nanopack_compare;
    nanopack_stats;
    nanopack_phasing;
    nanopack_overview
} from '../modules/nanopack_qc.nf'
include { prepare_deseq2 } from '../modules/deseq2.nf'
include { deseq2_analysis } from '../modules/deseq2.nf'

// Import common utilities
include { create_sample_channel; extract_sample_id; create_bam_channel } from '../modules/common.nf'

// Workflow definition
workflow DRS_PIPELINE {
    take:
        sample_info  // path to sample info file
        reference    // path to reference genome

    main:
        // Create sample channel
        sample_ch = create_sample_channel(sample_info)

        // Run basecalling for each sample
        basecalling_rna(
            sample_ch.map { sample_id, treatment, pod5_dir -> 
                tuple(sample_id, pod5_dir)
            }
        )

        // Run alignment for each sample
        align(
            basecalling_rna.out
                .map { bam ->
                    tuple(params.mm2opts, bam.getSimpleName(), bam, file(params.reference))
                }
        )

        // Process BAM files (sort and index)
        process_sam(align.out)

        // Format BAM files for nanopack QC using utility function
        bam_bai_ch = create_bam_bai_channel(process_sam.out)

        // Run nanopack QC on aligned BAM files
        nanopack_plot(bam_bai_ch)
        nanopack_stats(bam_bai_ch)

        // Run nanocount for transcript quantification
        nanocount(
            process_sam.out[0].map { bam ->
                def sample_id = extract_sample_id(bam)
                tuple(bam, sample_id)
            }
        )

    emit:
        basecalled_bam = basecalling_rna.out
        bam = process_sam.out[0]      // The sorted BAM file
        bai = process_sam.out[1]      // The BAM index file
        counts = nanocount.out        // Transcript counts
        samples = sample_ch           // Sample information
        qc_plots = nanopack_plot.out.plots
        qc_stats = nanopack_stats.out.stats
}

// Entry point
workflow {
    // Parameter validation
    if (!params.sample_info) {
        error "Sample info file is required: --sample_info sample_info.txt"
    }
    if (!params.reference) {
        error "Reference genome is required: --reference reference.fasta"
    }

    DRS_PIPELINE(
        file(params.sample_info),
        file(params.reference)
    )
} 