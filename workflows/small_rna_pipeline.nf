#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { basecalling_small_rna } from '../modules/basecalling.nf'
include { basecalling_small_rna_fastq } from '../modules/basecalling.nf'
include { bam_to_fastq } from '../modules/align.nf'
include { nanofilt } from '../modules/align.nf'
include { bwa_index } from '../modules/align.nf'
include { bwa_align } from '../modules/align.nf'
include { minimap2_align_sr } from '../modules/align.nf'
include { sam_to_bam } from '../modules/align.nf'
include { process_sam } from '../modules/align.nf'
include { nanocount } from '../modules/nanocount.nf'
include {
    nanopack_plot;
    nanopack_stats
} from '../modules/nanopack_qc.nf'

// Import common utilities
include { create_sample_channel; extract_sample_id; create_bam_channel; create_bam_bai_channel } from '../modules/common.nf'

// Workflow definition
workflow SMALL_RNA_PIPELINE {
    take:
        sample_info  // path to sample info file
        reference    // path to reference genome

    main:
        // Create sample channel
        sample_ch = create_sample_channel(sample_info, params.skip_basecalling, params.use_fastq)

        if (params.use_fastq) {
            // Run basecalling with FASTQ output
            basecalling_small_rna_fastq(
                sample_ch.map { sample_id, treatment, pod5_dir ->
                    tuple(sample_id, pod5_dir)
                }
            )
            filtered_fastq_ch = basecalling_small_rna_fastq.out
        } else if (params.skip_basecalling) {
            // Skip basecalling - use provided BAM files
            basecalled_bam_ch = sample_ch.map { sample_id, treatment, bam_file ->
                bam_file
            }

            // Convert BAM to FASTQ
            bam_to_fastq(
                basecalled_bam_ch.map { bam ->
                    tuple(bam.getSimpleName(), bam)
                }
            )

            // Filter reads by Q score
            nanofilt(bam_to_fastq.out)
            filtered_fastq_ch = nanofilt.out
        } else {
            // Run standard basecalling (BAM output)
            basecalling_small_rna(
                sample_ch.map { sample_id, treatment, pod5_dir ->
                    tuple(sample_id, pod5_dir)
                }
            )

            // Convert BAM to FASTQ
            bam_to_fastq(
                basecalling_small_rna.out
                    .map { bam ->
                        tuple(bam.getSimpleName(), bam)
                    }
            )

            // Filter reads by Q score
            nanofilt(bam_to_fastq.out)
            filtered_fastq_ch = nanofilt.out
        }

        // Alignment: choose between BWA and minimap2
        if (params.use_minimap2) {
            // Use minimap2 for alignment (outputs SAM)
            minimap2_align_sr(
                filtered_fastq_ch
                    .map { sample_id, fastq ->
                        tuple(sample_id, fastq, file(params.reference))
                    }
            )

            // Convert SAM to BAM
            sam_to_bam(minimap2_align_sr.out)
            aligned_bam_ch = sam_to_bam.out
        } else {
            // Use BWA for alignment (default)
            // Build BWA index for reference genome
            bwa_index(file(params.reference))

            // Run BWA alignment for each sample
            bwa_align(
                filtered_fastq_ch
                    .combine(bwa_index.out)
                    .map { sample_id, fastq, reference, index_files ->
                        tuple(sample_id, fastq, reference, index_files)
                    }
            )

            // Convert SAM to BAM
            sam_to_bam(bwa_align.out)
            aligned_bam_ch = sam_to_bam.out
        }

        // Process BAM files (sort and index)
        process_sam(aligned_bam_ch)

        // Format BAM files for nanopack QC - process_sam.out is already a channel of tuples [bam, bai]
        bam_bai_ch = process_sam.out
            .map { bam, bai ->
                def sample_id = extract_sample_id(bam)
                tuple(sample_id, bam, bai)
            }

        // Run nanopack QC on aligned BAM files
        nanopack_plot(bam_bai_ch)
        nanopack_stats(bam_bai_ch)

        // Run nanocount for small RNA quantification
        nanocount(
            process_sam.out
                .map { bam, bai ->
                    def sample_id = extract_sample_id(bam)
                    tuple(bam, sample_id)
                }
        )

    emit:
        basecalled_bam = params.use_fastq ? Channel.empty() : (params.skip_basecalling ? sample_ch.map { s, t, bam -> bam } : basecalling_small_rna.out)
        basecalled_fastq = params.use_fastq ? basecalling_small_rna_fastq.out : Channel.empty()
        filtered_fastq = filtered_fastq_ch
        bam = process_sam.out[0]      // The sorted BAM file
        bai = process_sam.out[1]      // The BAM index file
        counts = nanocount.out        // Small RNA counts
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

    SMALL_RNA_PIPELINE(
        file(params.sample_info),
        file(params.reference)
    )
}
