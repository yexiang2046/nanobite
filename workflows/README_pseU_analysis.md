# Pseudouridine (ψ) Modification Analysis Pipeline

## Overview

This workflow performs comprehensive pseudouridine modification analysis from Oxford Nanopore POD5 files. The pipeline includes:

1. **Basecalling with ψ detection** - Uses Dorado with pseudouridine modification models
2. **Quality filtering** - Filters reads with Q score > 10 using NanoFilt
3. **Reference alignment** - Aligns filtered reads to reference genome using minimap2
4. **Modification pileup** - Calculates per-site modification frequencies
5. **Site filtering** - Extracts high-confidence pseudouridine sites
6. **Read-level extraction** - Provides detailed modification calls per read
7. **Summary statistics** - Overall modification statistics

## Requirements

- Nextflow >= 20.07.1
- Docker enabled
- GPU recommended (CUDA-compatible for faster basecalling)
- Input: POD5 files from Oxford Nanopore direct RNA sequencing

## Quick Start

```bash
# Basic usage (full pipeline with basecalling)
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta

# Skip basecalling and start from existing BAM files
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta \
    --skip_basecalling \
    --bam_dir /path/to/basecalled_bams

# With custom parameters
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta \
    --output_dir pseU_results \
    --min_coverage 10 \
    --prob_threshold 0.9
```

## Input Files

### Sample Info File (sample_info.txt)

Tab-separated file with the following columns:

**For full pipeline (with basecalling):**

| Column | Description |
|--------|-------------|
| SampleId | Unique sample identifier |
| pod5_path | Full path to POD5 file or directory |

Example:
```
SampleId	pod5_path
control_rep1	/path/to/control_rep1.pod5
control_rep2	/path/to/control_rep2.pod5
treated_rep1	/path/to/treated_rep1.pod5
treated_rep2	/path/to/treated_rep2.pod5
```

**For skipping basecalling (starting from BAM files):**

| Column | Description |
|--------|-------------|
| SampleId | Unique sample identifier |
| bam_path | Full path to basecalled BAM file with modifications |

Example:
```
SampleId	bam_path
control_rep1	/path/to/bams/control_rep1.bam
control_rep2	/path/to/bams/control_rep2.bam
treated_rep1	/path/to/bams/treated_rep1.bam
treated_rep2	/path/to/bams/treated_rep2.bam
```

**Note:** When using `--skip_basecalling`, you can also use `--bam_dir` to specify the directory and omit the `bam_path` column if all BAM files are in the same directory.

### Reference Genome

FASTA format reference genome (e.g., transcriptome or genome assembly).

## Parameters

### Required Parameters

- `--reference`: Reference genome FASTA file
- `--sample_info`: Sample information file

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--output_dir` | "results" | Output directory |
| `--skip_basecalling` | false | Skip basecalling and start with existing BAM files |
| `--bam_dir` | null | Directory containing basecalled BAM files (required when `--skip_basecalling` is true) |
| `--min_coverage` | 5 | Minimum coverage for modification calling |
| `--prob_threshold` | 0.8 | Probability threshold for modification calls |
| `--min_qscore` | 10 | Minimum Q score for read filtering |
| `--mm2opts` | "-ax splice -uf -k14" | Minimap2 alignment options |

## Output Structure

```
results/
├── basecalling_pseU/         # Basecalled BAM files with modifications
│   ├── sample1.bam
│   └── sample2.bam
├── filtered/                 # Quality filtered FASTQ files (Q > 10)
│   ├── sample1_filtered.fastq.gz
│   └── sample2_filtered.fastq.gz
├── alignment/                # Minimap2 aligned BAM files
│   ├── sample1_aligned.bam
│   └── sample2_aligned.bam
├── bam/                      # Sorted BAM files with indices
│   ├── sample1.srt.bam
│   ├── sample1.srt.bam.bai
│   ├── sample2.srt.bam
│   └── sample2.srt.bam.bai
├── modkit_pileup/            # Per-site modification frequencies
│   ├── sample1_pileup.bed
│   └── sample2_pileup.bed
├── pseU_sites/               # Filtered pseudouridine sites
│   ├── sample1_pseU_sites.bed
│   └── sample2_pseU_sites.bed
├── modkit_extract/           # Read-level modification calls
│   ├── sample1_mod_calls.tsv
│   └── sample2_mod_calls.tsv
├── modkit_summary/           # Summary statistics
│   ├── sample1_summary.txt
│   └── sample2_summary.txt
├── pipeline_report.html      # Execution report
├── timeline.html             # Timeline of execution
└── trace.txt                 # Resource usage trace
```

## Output File Descriptions

### Pileup BED Files (`modkit_pileup/`)

Per-site modification frequencies in BED format:
- Chromosome/contig
- Start position
- End position
- Modification type
- Coverage
- Modification frequency
- Additional statistics

### Pseudouridine Sites BED (`pseU_sites/`)

Filtered sites containing only high-confidence pseudouridine modifications.

### Modification Calls TSV (`modkit_extract/`)

Read-level modification calls:
- Read ID
- Reference position
- Modification type
- Modification probability
- Quality scores

### Summary Statistics (`modkit_summary/`)

Overall statistics including:
- Total reads processed
- Number of modified sites
- Average modification frequency
- Coverage statistics

## Advanced Usage

### Skipping Basecalling

The `--skip_basecalling` option allows you to start the analysis from existing basecalled BAM files. This is useful when:

1. **Basecalling was already performed**: You've already run Dorado basecalling with pseudouridine modification detection
2. **Reanalyzing with different parameters**: You want to rerun alignment or modification analysis without re-basecalling
3. **Time/resource constraints**: Basecalling is computationally expensive; skipping it saves time and GPU resources
4. **Working with pre-processed data**: You received basecalled BAM files from a collaborator

**Important requirements:**
- BAM files must have been basecalled with pseudouridine modification detection enabled
- BAM files should contain modification tags (e.g., MM, ML tags from Dorado)
- All BAM files should be in a single directory specified by `--bam_dir`

Example workflow:

```bash
# Step 1: Initial basecalling (time-consuming)
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta \
    --output_dir results_run1

# Step 2: Reanalyze with different modification thresholds (skip basecalling)
nextflow run workflows/pseU_analysis.nf \
    --reference reference.fasta \
    --skip_basecalling \
    --bam_dir results_run1/basecalling_pseU \
    --output_dir results_run2 \
    --min_coverage 10 \
    --prob_threshold 0.95
```

### Using with GPU Profile

```bash
nextflow run workflows/pseU_analysis.nf \
    -profile gpu \
    --sample_info sample_info.txt \
    --reference reference.fasta
```

### High-Stringency Analysis

For higher confidence calls:

```bash
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta \
    --min_coverage 20 \
    --prob_threshold 0.95
```

### Custom Alignment Parameters

For different alignment strategies:

```bash
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference reference.fasta \
    --mm2opts "-ax map-ont -uf"
```

## Downstream Analysis

### Loading Pseudouridine Sites in Python

```python
import pandas as pd

# Load filtered pseudouridine sites
pseu_sites = pd.read_csv('results/pseU_sites/sample1_pseU_sites.bed',
                         sep='\t', header=None,
                         names=['chrom', 'start', 'end', 'mod_type',
                                'coverage', 'freq', 'stats'])

# Filter by frequency
high_freq = pseu_sites[pseu_sites['freq'] > 0.5]
```

### Comparing Between Samples

```bash
# Use bedtools to find common sites
bedtools intersect \
    -a results/pseU_sites/control_pseU_sites.bed \
    -b results/pseU_sites/treated_pseU_sites.bed \
    > common_sites.bed
```

## Troubleshooting

### GPU Not Available

If GPU is not available, the workflow will still run but basecalling will be slower. Ensure Docker has GPU access:

```bash
docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi
```

### Memory Issues

Increase memory allocation in `nextflow.config` for the `basecalling_pseU` process:

```groovy
withName: 'basecalling_pseU' {
    memory = '64 GB'
}
```

### Low Pseudouridine Site Counts

Try adjusting parameters:
- Lower `--min_coverage` (e.g., 3)
- Lower `--prob_threshold` (e.g., 0.7)
- Check input quality with NanoPlot

## Citations

If you use this pipeline, please cite:

- **Dorado**: Oxford Nanopore Technologies
- **modkit**: Oxford Nanopore Technologies modification analysis toolkit
- **Nextflow**: Di Tommaso et al. (2017) Nature Biotechnology

## Support

For issues specific to this workflow, please check:
- Input file format matches requirements
- Reference genome is properly indexed
- Docker containers are accessible
- Sufficient disk space and memory available
