#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Mock processes for testing
process mock_basecalling {
    input:
    tuple val(sample_id), path(pod5_dir)

    output:
    path "${sample_id}.bam"

    script:
    """
    echo "Mock basecalling for ${sample_id}" > ${sample_id}.bam
    """
}

process mock_filter {
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_filtered.fastq.gz")

    script:
    """
    echo "Mock filtered reads" | gzip > ${sample_id}_filtered.fastq.gz
    """
}

// Test workflow
workflow test_pipeline_structure {
    main:
        // Create test sample channel
        Channel
            .of(
                tuple('sample1', 'control', file('tests/data/pod5/sample1')),
                tuple('sample2', 'treated', file('tests/data/pod5/sample2'))
            )
            .set { sample_ch }

        // Test sample channel creation
        sample_ch.view { sample_id, treatment, pod5_dir ->
            println "✓ Sample channel: [${sample_id}, ${treatment}, ${pod5_dir}]"
        }

        // Test basecalling mock
        mock_basecalling(
            sample_ch.map { sample_id, treatment, pod5_dir ->
                tuple(sample_id, pod5_dir)
            }
        )

        // Test filter mock
        mock_filter(
            mock_basecalling.out.map { bam ->
                tuple(bam.getSimpleName(), bam)
            }
        )

        mock_filter.out.view { sample_id, fastq ->
            println "✓ Filtered output: [${sample_id}, ${fastq}]"
        }
}

workflow {
    println "=== Testing Small RNA Pipeline Structure ==="
    println ""

    test_pipeline_structure()

    println ""
    println "=== Pipeline structure test completed ==="
}
