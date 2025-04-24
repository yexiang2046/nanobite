nextflow.enable.dsl=2

// Import modules
include { basecalling_rna } from '../modules/basecalling.nf'
include { align } from '../modules/align.nf'
include { nanocount } from '../modules/nanocount.nf'

// Create a channel from sample info file
def create_sample_channel(sample_info) {
    Channel
        .fromPath(sample_info)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            def sample_id = row.SampleId
            def treatment = row.treatment
            def pod5_dir = file(row.pod5_path.trim())
            
            if (!pod5_dir.exists()) {
                error "Pod5 directory not found for sample ${sample_id}: ${pod5_dir}"
            }
            
            return tuple(sample_id, treatment, pod5_dir)
        }
}

// Workflow definition
workflow RNA_BASECALLING {
    take:
        sample_info // path to sample info file

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
            basecalling_rna.out.map { bam ->
                tuple(params.mm2opts, sample_ch.map { id, treatment, _ -> id }.first(), bam, file(params.reference))
            }
        )

        // Run NanoCount for transcript abundance estimation
        nanocount(
            align.out[0].map { bam ->
                tuple(bam, sample_ch.map { id, _, _ -> id }.first())
            }
        )

    emit:
        bam = align.out[0]  // The sorted BAM file
        bai = align.out[1]  // The BAM index file
        counts = nanocount.out.counts  // Transcript counts
        log = nanocount.out.log  // NanoCount log files
        samples = sample_ch // emit sample channel for downstream processing
}

// Entry point
workflow {
    if (!params.sample_info) {
        error "Sample info file is required: --sample_info sample_info.txt"
    }

    RNA_BASECALLING(
        file(params.sample_info)
    )
} 