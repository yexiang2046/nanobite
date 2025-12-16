#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import functions to test
include { create_sample_channel; extract_sample_id; create_bam_channel } from '../modules/common.nf'

// Test: extract_sample_id function
workflow test_extract_sample_id {
    main:
        // Create test files
        def test_cases = [
            [name: 'sample1_aligned.srt.bam', expected: 'sample1'],
            [name: 'sample2_aligned.bam', expected: 'sample2'],
            [name: 'sample3.srt.bam', expected: 'sample3'],
            [name: 'sample4.sorted.bam', expected: 'sample4']
        ]

        test_cases.each { test ->
            def mock_file = [getSimpleName: { -> test.name.replaceAll(/\.bam$/, '') }]
            def result = extract_sample_id(mock_file)

            if (result == test.expected) {
                println "✓ extract_sample_id('${test.name}') = '${result}'"
            } else {
                println "✗ extract_sample_id('${test.name}') = '${result}' (expected: '${test.expected}')"
                System.exit(1)
            }
        }

        println "\nAll extract_sample_id tests passed!"
}

// Test: create_bam_channel function
workflow test_create_bam_channel {
    main:
        // Create mock BAM channel
        Channel
            .fromList([
                'sample1_aligned.srt.bam',
                'sample2_aligned.bam'
            ])
            .map { filename ->
                // Mock file object
                [getSimpleName: { -> filename.replaceAll(/\.bam$/, '') }]
            }
            .set { mock_bam_ch }

        // Test the function
        create_bam_channel(mock_bam_ch)
            .view { sample_id, bam ->
                println "✓ create_bam_channel: [${sample_id}, ${bam}]"
            }
}

// Main test workflow
workflow {
    println "=== Running Common Module Tests ==="
    println ""

    test_extract_sample_id()

    println ""
    test_create_bam_channel()

    println ""
    println "=== All tests completed ==="
}
