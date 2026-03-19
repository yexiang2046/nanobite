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
   - `nerdseq_analysis.nf`: NERD-seq style analysis with minimap2 alignment preserving MM/ML tags and modkit processing
   - `pseU_analysis.nf`: Pseudouridine-specific modification analysis with comprehensive modkit tools
   - `genome_transcript_quantification.nf`: Transcript quantification from genome-aligned BAMs with corrected DRS parameters and configurable multi-mapper handling

2. **Process Modules** (in `modules/`):
   - `basecalling.nf`: Four basecalling processes (mod_basecalling_rna, basecalling_rna, basecalling_dna, basecalling_small_rna)
   - `align.nf`: Alignment using Dorado aligner, minimap2 with MM/ML tag preservation, BWA for small RNA, and BAM processing with samtools
   - `modkit.nf`: RNA modification analysis tools (pileup, extract, summary, filter_pseU)
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

### NERD-seq Style RNA Modification Analysis
```bash
nextflow run workflows/nerdseq_analysis.nf \
    --reference reference.fa \
    --bam_dir /path/to/basecalled_bams \
    --mm2opts_nerd "-ax sr" \
    --prob_threshold 0.8
```

### Pseudouridine Modification Analysis
```bash
nextflow run workflows/pseU_analysis.nf \
    --reference reference.fa \
    --sample_info sample_info.txt \
    --skip_basecalling true \
    --bam_dir /path/to/bams
```

### Genome Transcript Quantification
```bash
nextflow run workflows/genome_transcript_quantification.nf \
    --bam_dir /path/to/genome_aligned_bams \
    --gtf_file gencode.v44.annotation.gtf \
    --sample_info sample_info.txt \
    --output_dir results/transcript_counts
```

**With multi-mapper fractional counting (recommended for small RNAs):**
```bash
nextflow run workflows/genome_transcript_quantification.nf \
    --bam_dir /path/to/bams \
    --gtf_file gencode.v44.annotation.gtf \
    --sample_info sample_info.txt \
    --count_multimappers true \
    --use_fractional_counting true
```

### Configuration Profiles
- GPU profile: `nextflow run -profile gpu ...` (enables GPU for basecalling)
- CPU profile: `nextflow run -profile cpu ...` (CPU-only mode)

## Key Parameters (in nextflow.config)

### General Parameters
- `--cpus`: Number of CPUs (default: 16)
- `--sample_info`: Path to sample information file
- `--reference`: Reference genome FASTA
- `--output_dir`: Output directory (default: "results")
- `--mm2opts`: Minimap2 alignment options for RNA-seq
- `--mm2opts_sr`: Minimap2 options for small RNA (default: "-ax sr")
- `--mm2opts_nerd`: Minimap2 options for NERD-seq style analysis with MM/ML tag preservation (default: "-ax sr")
- `--use_minimap2`: Use minimap2 instead of BWA for small RNA alignment (default: false)
- `--use_gpu`: Enable/disable GPU usage
- `--base_mods`: Modification types to analyze (default: "m,a,u")
- `--min_coverage`: Minimum coverage for modification calling (default: 5)
- `--prob_threshold`: Probability threshold for modification calls (default: 0.8)
- `--min_qscore`: Minimum Q score for read filtering (default: 10)

### Transcript Quantification Parameters (genome_transcript_quantification.nf)
- `--stranded`: Strandedness (0=unstranded, 1=stranded for DRS, 2=reverse, default: 1)
- `--min_mapq`: Minimum mapping quality (default: 10)
- `--feature_type`: Feature type to count from GTF (default: "exon")
- `--gene_attribute`: GTF attribute for grouping (default: "gene_id")
- `--count_multimappers`: Enable counting of multi-mapping reads (default: false)
- `--use_fractional_counting`: Use fractional counting for multi-mappers (default: false, requires count_multimappers=true)
- `--count_overlapping`: Count reads overlapping multiple features (default: false)
- `--run_deseq2`: Run differential expression analysis (default: false)

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

## Genome Transcript Quantification Workflow Details

### Overview
The `genome_transcript_quantification.nf` workflow quantifies all transcripts from genome-aligned DRS BAM files using featureCounts with **corrected parameters** for single-end Nanopore long-read data.

### Key Improvements Over Original featureCounts Implementation
The original `featurecounts` process used **incorrect parameters** for DRS data:
- ❌ `-p`: Paired-end mode (DRS is single-end)
- ❌ `-B`: Require both ends mapped (not applicable to single-end)
- ❌ `-C`: Chimeric fragment filtering (for paired-end only)
- ❌ Missing `-L`: Long-read mode essential for Nanopore
- ❌ Missing `-s`: Strand specification (DRS is stranded)
- ❌ Missing `--primary`: To avoid double-counting secondary alignments

The new `featurecounts_genome_aligned` process uses:
- ✅ `-L`: Long-read mode for Nanopore error handling
- ✅ `-s 1`: Stranded mode (DRS preserves strand information)
- ✅ `--primary`: Only count primary alignments
- ✅ `-Q 10`: Mapping quality threshold
- ✅ Configurable multi-mapper handling

### Multi-mapper Handling Modes

**Mode 1: Exclude multi-mappers (Default, Conservative)**
```bash
--count_multimappers false
```
- Multi-mapping reads are NOT counted (reported as "Unassigned_MultiMapping")
- Best for: Quantifying uniquely-mapped genes
- Limitation: Underestimates gene families and small RNAs with multiple genomic copies

**Mode 2: Count all multi-mapper alignments (Permissive)**
```bash
--count_multimappers true --use_fractional_counting false
```
- Each alignment of a multi-mapper is counted separately
- A read mapping to 5 locations adds 5 to total counts
- Warning: Inflates counts for repetitive genes

**Mode 3: Fractional counting (Recommended for small RNAs)**
```bash
--count_multimappers true --use_fractional_counting true
```
- Multi-mapper counts are distributed proportionally
- A read mapping to 5 locations contributes 0.2 counts to each
- Best for: Small RNAs (U6 snRNA has ~50 genomic copies, tRNAs have 400-600 copies)
- Most accurate representation of read abundance across gene families

### Use Cases

**For mRNA/lncRNA quantification (mostly unique genes):**
```bash
nextflow run workflows/genome_transcript_quantification.nf \
    --bam_dir bams/ \
    --gtf_file gencode.v44.annotation.gtf \
    --sample_info samples.txt
```

**For small RNA quantification (miRNA, snoRNA, snRNA, tRNA):**
```bash
nextflow run workflows/genome_transcript_quantification.nf \
    --bam_dir bams/ \
    --gtf_file gencode.v44.annotation.gtf \
    --sample_info samples.txt \
    --count_multimappers true \
    --use_fractional_counting true
```

**For gene-level (not exon-level) counting:**
```bash
nextflow run workflows/genome_transcript_quantification.nf \
    --bam_dir bams/ \
    --gtf_file gencode.v44.annotation.gtf \
    --sample_info samples.txt \
    --feature_type gene \
    --gene_attribute gene_name
```

### Input Requirements
- **BAM files**: Genome-aligned DRS reads (sorted with .bai indices)
- **GTF file**: Annotation file (GENCODE, Ensembl, or custom)
- **Sample info**: Tab-separated file with `SampleId` and `treatment` columns

### Output Files
- `{sample}_featurecounts.txt`: Per-sample counts with gene annotations
- `{sample}_featurecounts.txt.summary`: Assignment statistics
- `count_matrix.txt`: Merged counts across all samples
- `count_matrix_clean.txt`: Counts only (ready for DESeq2)

### Quality Control Metrics
Check `.summary` files for:
- **Assigned**: Should be >50-70% for good quantification
- **Unassigned_MultiMapping**: High if multi-mappers excluded (expected for small RNAs)
- **Unassigned_NoFeatures**: Reads mapping outside annotated features
- **Unassigned_Ambiguity**: Reads overlapping multiple genes

## Development Tips

- When adding new processes, follow the existing module pattern: define process in `modules/`, import in `workflows/`
- Use `publishDir` with `mode: 'copy'` for process outputs
- Channel operators commonly used: `map`, `collect`, `splitCsv`, `fromPath`
- Sample IDs are extracted using `.getSimpleName()` on file objects
- Always validate input parameters and provide clear error messages using `error "message"`
- avoid using multiple tools in one process, since most containers are tool specific