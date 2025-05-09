nextflow.enable.dsl=2


// Import modules
include { mod_basecalling_rna } from '../modules/basecalling.nf'
include { 
    nanopack_plot;
    nanopack_stats;
} from '../modules/nanopack_qc.nf'

// Create a channel from sample info file
def create_sample_channel(sample_info) {
    Channel
        .fromPath(sample_info)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def sample_id = row.SampleId
            def treatment = row.treatment
            def pod5_dir = file(row.pod5_path)
            
            if (!pod5_dir.exists()) {
                error "Pod5 directory not found for sample ${sample_id}: ${pod5_dir}"
            }
            
            return tuple(sample_id, treatment, pod5_dir)
        }
}

// Workflow definition
workflow RNA_MOD_BASECALLING {
    take:
        sample_info  // path to sample info file
        reference    // path to reference genome

    main:
        // Create sample channel
        sample_ch = create_sample_channel(sample_info)

        // Run modification basecalling for each sample
        mod_basecalling_rna(
            sample_ch.map { sample_id, treatment, pod5_dir -> 
                tuple(sample_id, pod5_dir)
            }
        )


        // Run alignment for each sample
        align(
            mod_basecalling_rna.out.bam.map { sample_id, bam ->
                tuple(sample_id, bam, file(params.reference))
            }
        )

        // Run nanopack QC on aligned BAM files
        nanopack_plot(align.out.bam)
        nanopack_stats(align.out.bam)

    emit:
        bam = mod_basecalling_rna.out.bam
        align_bam = align.out.bam
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

    RNA_MOD_BASECALLING(
        file(params.sample_info),
        file(params.reference),
    )
} 