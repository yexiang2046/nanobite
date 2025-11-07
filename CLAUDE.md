# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NanoBite is a Nextflow-based pipeline for analyzing Oxford Nanopore direct RNA sequencing (DRS) data. The pipeline processes raw POD5 files through basecalling, alignment, quality control, and differential expression analysis.

## Key Technical Details

### Input Requirements
- **POD5 files**: Raw sequencing data input for basecalling (use reads with Q score > 10 for alignment)
- **Sample info file**: Tab-separated file with columns `SampleId`, `treatment`, and `pod5_path`
- **Reference genome**: FASTA format for alignment

### Workflow Architecture

The pipeline is organized into modular Nextflow workflows and processes:

1. **Main Workflows** (in `workflows/`):
   - `DRS_transcriptomic_pipeline.nf`: Full pipeline from basecalling to quantification
   - `differential_expression.nf`: Standalone DESeq2 analysis from count files
   - `featurecounts_analysis.nf`: Gene-level quantification using featureCounts
   - `rna_mod_basecalling.nf`: RNA modification basecalling and alignment
   - `small_rna_pipeline.nf`: Small RNA quantification using BWA and NanoCount

2. **Process Modules** (in `modules/`):
   - `basecalling.nf`: Four basecalling processes (mod_basecalling_rna, basecalling_rna, basecalling_dna, basecalling_small_rna)
   - `align.nf`: Alignment using Dorado aligner, BWA for small RNA, and BAM processing with samtools
   - `nanocount.nf`: Transcript-level quantification
   - `featurecounts.nf`: Gene-level quantification and aggregation
   - `deseq2.nf`: Differential expression analysis
   - `nanopack_qc.nf`: Quality control using NanoPack tools

### Data Flow Pattern

The pipeline uses Nextflow's DSL2 channel-based architecture:
- Sample info file → parsed into channels → split by process
- Basecalling output (BAM) → Alignment → Sorted BAM + BAI index
- Multiple samples collected → aggregated for comparative analysis

## Running the Pipeline

### Full Pipeline
```bash
nextflow run workflows/DRS_transcriptomic_pipeline.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta \
    --output_dir results
```

### Differential Expression Only
```bash
nextflow run workflows/differential_expression.nf \
    --sample_info sample_info.txt \
    --count_dir path/to/counts
```

### featureCounts Analysis
```bash
nextflow run workflows/featurecounts_analysis.nf \
    --bam_dir /path/to/bams \
    --gtf_file annotation.gtf \
    --sample_info sample_info.txt
```

### RNA Modification Basecalling
```bash
nextflow run workflows/rna_mod_basecalling.nf \
    --reference reference.fa \
    --sample_info sample_info.txt
```

### Small RNA Quantification
```bash
nextflow run workflows/small_rna_pipeline.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta \
    --output_dir results
```

### Configuration Profiles
- GPU profile: `nextflow run -profile gpu ...` (enables GPU for basecalling)
- CPU profile: `nextflow run -profile cpu ...` (CPU-only mode)

## Key Parameters (in nextflow.config)

- `--cpus`: Number of CPUs (default: 16)
- `--sample_info`: Path to sample information file
- `--reference`: Reference genome FASTA
- `--output_dir`: Output directory (default: "results")
- `--mm2opts`: Minimap2 alignment options for RNA-seq
- `--mm2opts_sr`: Minimap2 options for small RNA (default: "-ax sr")
- `--use_minimap2`: Use minimap2 instead of BWA for small RNA alignment (default: false)
- `--use_gpu`: Enable/disable GPU usage
- `--base_mods`: Modification types to analyze (default: "m,a,u")
- `--min_coverage`: Minimum coverage for modification calling (default: 5)
- `--prob_threshold`: Probability threshold for modification calls (default: 0.8)
- `--min_qscore`: Minimum Q score for read filtering (default: 10)

## Container Strategy

The pipeline uses Docker containers specified per process:
- `ontresearch/dorado:latest` or `staphb/dorado:0.9.0-cuda12.2.0`: Basecalling and alignment (GPU-enabled)
- `quay.io/biocontainers/nanofilt:2.8.0--py_0`: Quality filtering with NanoFilt
- `staphb/bwa:0.7.17`: BWA alignment for small RNA (default)
- `staphb/minimap2:2.28`: Minimap2 alignment (alternative for small RNA)
- `staphb/samtools:1.21`: BAM file processing
- `xiang2019/deseq2:v1.0.0`: Differential expression analysis
- `xiang2019/modkit:v1.0.0`: RNA modification analysis
- `biocontainers/subread:v1.6.3dfsg-1-deb_cv1`: featureCounts

## Output Structure

```
results/
├── basecalling_rna/          # Basecalled BAM files (RNA)
├── basecalling_small_rna/    # Basecalled BAM files (small RNA)
├── filtered/                 # Quality filtered FASTQ files
├── alignment/                # Aligned BAM files
├── bam/                      # Sorted BAM files and indices
├── nanopack/                 # QC reports and plots
├── nanocount/                # Transcript quantification
├── featurecounts/            # Gene-level counts
└── deseq2/                   # Differential expression results
    ├── results/              # TSV result files
    └── plots/                # PDF visualizations
```

## Important Implementation Notes

1. **Sample Channel Creation**: The `create_sample_channel()` function in workflows validates sample info files and creates tuples of (sample_id, treatment, pod5_path). Always validate that required columns exist and paths are valid.

2. **Process Dependencies**: The DRS pipeline has strict sequential dependencies:
   - Basecalling → Alignment → BAM sorting → QC/Quantification
   - DESeq2 requires aggregated count files from all samples

3. **GPU Usage**: Basecalling and alignment processes are GPU-enabled by default. Use `containerOptions = '--gpus all'` for GPU processes. The `use_gpu` parameter can be toggled via profiles.

4. **Error Handling**: Pipeline uses retry strategy (maxRetries = 3) defined in nextflow.config. Process errors should be handled with clear validation messages.

5. **R Scripts**: DESeq2 and featureCounts modules expect R scripts in `bin/` directory (merge_counts.R, run_deseq2.R).

6. **Quality Filtering**: The small RNA pipeline includes a NanoFilt filtering step:
   - Filters reads by Q score using the `--min_qscore` parameter (default: 10)
   - Converts BAM to FASTQ → filters with NanoFilt → outputs compressed FASTQ
   - Filtering occurs after basecalling and before alignment

7. **Small RNA Alignment Options**: The small RNA pipeline supports two aligners:
   - **BWA (default)**: The `bwa_align` process uses parameters optimized for Nanopore small RNA:
     - `-C`: Append FASTA/FASTQ comment to SAM output
     - `-W 13`: Band width for banded alignment
     - `-k 6`: Minimum seed length
     - `-x ont2d`: Nanopore 2D read alignment preset
     - Requires BWA indexing step before alignment
   - **Minimap2 (optional)**: Enable with `--use_minimap2 true`
     - Uses `-ax sr` preset for short read alignment (configurable via `--mm2opts_sr`)
     - Outputs BAM directly without requiring SAM-to-BAM conversion
     - Generally faster than BWA for Nanopore data

## Development Tips

- When adding new processes, follow the existing module pattern: define process in `modules/`, import in `workflows/`
- Use `publishDir` with `mode: 'copy'` for process outputs
- Channel operators commonly used: `map`, `collect`, `splitCsv`, `fromPath`
- Sample IDs are extracted using `.getSimpleName()` on file objects
- Always validate input parameters and provide clear error messages using `error "message"`
- avoid using multiple tools in one process, since most containers are tool specific