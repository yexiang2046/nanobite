# NanoBite: Direct RNA Sequencing Pipeline

A Nextflow pipeline for analyzing direct RNA sequencing data from Oxford Nanopore Technologies.

## Overview

This pipeline processes direct RNA sequencing data through the following steps:
1. Basecalling of raw POD5 files
2. Alignment to reference genome
3. BAM file processing (sorting and indexing)
4. Quality control using NanoPack tools
5. Transcript quantification using NanoCount
6. Differential expression analysis using DESeq2

## Requirements

- Nextflow (version 22.10.1 or later)
- Docker
- Sufficient disk space for POD5 files and analysis results

## Input

The pipeline requires two main inputs:

1. A sample information file (tab-separated) with the following columns:
   - `SampleId`: Unique identifier for each sample
   - `treatment`: Treatment group for differential expression analysis
   - `pod5_path`: Path to the directory containing POD5 files

2. Reference genome in FASTA format

Example sample information file:
```tsv
SampleId	treatment	pod5_path
sample1	control	/path/to/sample1/pod5
sample2	treatment	/path/to/sample2/pod5
```

## Usage

Run the pipeline using:

```bash
nextflow run workflows/DRS_transcriptomic_pipeline.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta
```

### Optional Parameters

- `--output_dir`: Directory for output files (default: "results")
- `--mm2opts`: Minimap2 alignment options (default: "-x map-ont -t 8")

## Output

The pipeline generates the following outputs in the specified output directory:

- `basecalling/`: Basecalled FASTQ files
- `alignment/`: Aligned BAM files and indices
- `nanopack/`: Quality control reports
  - `plots/`: NanoPlot visualizations
  - `stats/`: Cramino statistics
  - `phasing/`: Phasius phasing analysis
  - `overview/`: Kyber overview plots
  - `comparison/`: NanoComp comparison plots
- `nanocount/`: Transcript quantification results
- `deseq2/`: Differential expression analysis results
  - `results/`: DESeq2 analysis results
  - `plots/`: MA plots, PCA plots, and heatmaps

## Quality Control

The pipeline includes several quality control steps:
- NanoPlot for read quality visualization
- NanoComp for comparing multiple samples
- Cramino for BAM file statistics
- Phasius for read phasing analysis
- Kyber for BAM file overview

## Transcript Quantification

Transcript quantification is performed using NanoCount, which:
- Counts reads per transcript
- Generates count matrices for downstream analysis
- Provides transcript-level abundance estimates

## Differential Expression Analysis

The pipeline uses DESeq2 for differential expression analysis:
- Normalizes count data
- Performs statistical testing
- Generates visualizations (MA plots, PCA, heatmaps)
- Identifies differentially expressed transcripts

## Citation

If you use this pipeline in your research, please cite the following tools:
- Dorado (Oxford Nanopore Technologies)
- Minimap2 (in Dorado)
- NanoPack tools
- NanoCount
- DESeq2

## License

This pipeline is released under the MIT License.

