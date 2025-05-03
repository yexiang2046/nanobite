nextflow.enable.dsl=2

// Import modules
include { mod_basecalling_rna } from '../modules/basecalling.nf'
include { 
    modkit_rna_mods;
    modkit_inspect_rna_mods;
    modkit_extract_rna_mods;
    modkit_summary;
    modkit_prepare_regions 
} from '../modules/modkit.nf'
include { 
    nanopack_plot;
    nanopack_compare;
    nanopack_stats;
    nanopack_phasing;
    nanopack_overview 
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
        regions     // optional: path to regions bed file

    main:
        // Create sample channel
        sample_ch = create_sample_channel(sample_info)

        // Run modification basecalling for each sample
        mod_basecalling_rna(
            sample_ch.map { sample_id, treatment, pod5_dir -> 
                tuple(sample_id, pod5_dir)
            }
        )

        // Create default empty BED if none provided
        def regions_ch = regions ? Channel.fromPath(regions) : Channel.value(file('NO_FILE'))

        // Optional: Generate regions from motifs if requested
        if (params.generate_motif_regions) {
            modkit_prepare_regions(
                tuple(params.sample_id, reference),
                params.motif
            )
            regions_ch = modkit_prepare_regions.out.regions_bed
        }

        // Run alignment for each sample
        align(
            mod_basecalling_rna.out.bam.map { sample_id, bam ->
                tuple(sample_id, bam, file(params.reference))
            }
        )

        // Run nanopack QC on aligned BAM files
        nanopack_plot(align.out.bam)
        nanopack_stats(align.out.bam)
        nanopack_phasing(align.out.bam)
        nanopack_overview(align.out.bam)
        
        // Compare all BAM files
        nanopack_compare(align.out.bam.collect())

        // Run modkit analysis on basecalled files
        modkit_rna_mods(
            mod_basecalling_rna.out.bam.map { sample_id, bam -> 
                tuple(sample_id, bam)
            },
            reference,
            regions_ch
        )

        modkit_inspect_rna_mods(
            mod_basecalling_rna.out.bam,
            regions_ch
        )

        modkit_extract_rna_mods(
            mod_basecalling_rna.out.bam,
            regions_ch
        )

        modkit_summary(mod_basecalling_rna.out.bam)

    emit:
        bam = mod_basecalling_rna.out.bam
        align_bam = align.out.bam
        mod_beds = modkit_rna_mods.out.mod_beds
        mod_probs = modkit_inspect_rna_mods.out.inspection_results
        per_read_mods = modkit_extract_rna_mods.out.per_read_data
        summary = modkit_summary.out.summary
        samples = sample_ch
        qc_plots = nanopack_plot.out.plots
        qc_stats = nanopack_stats.out.stats
        qc_phasing = nanopack_phasing.out.phasing
        qc_overview = nanopack_overview.out.overview
        qc_comparison = nanopack_compare.out.comparison
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
        params.regions_bed
    )
} 