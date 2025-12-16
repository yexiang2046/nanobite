# Small RNA Pipeline Flowchart

```mermaid
flowchart TD
    %% Input
    Start([POD5 files]) --> basecalling{Dorado}

    %% Basecalling paths
    basecalling --> BAM

    %% Standard mode
    BAM --> BAM2FASTQ2[bam_to_fastq]
    BAM2FASTQ2 --> NanoFilt2[nanofilt<br/>Q score filter]
    NanoFilt2 --> FilteredFastq

    %% Alignment decision
    FilteredFastq --> AlignMode{Aligner?}

    %% BWA alignment path
    AlignMode -->|BWA| BWA_ALIGN[sam file]
    BWA_ALIGN --> SAM2BAM1[bam file]

    %% Minimap2 alignment path
    AlignMode -->|minimap2| MM2_ALIGN[minimap2 short read mode]
    MM2_ALIGN --> ALIGNED_BAM[bam file]

    %% BAM processing
    AlignedBAM --> ProcessSAM[process_sam<br/>sort & index]
    ProcessSAM --> SortedBAM[(Sorted BAM + BAI)]

    %% Parallel QC and quantification
    SortedBAM --> NanoPack_Plot[nanopack_plot]
    SortedBAM --> NanoPack_Stats[nanopack_stats]
    SortedBAM --> NanoCount[nanocount]

    %% Outputs
    NanoPack_Plot --> QC_Plots([QC Plots])
    NanoPack_Stats --> QC_Stats([QC Statistics])
    NanoCount --> Counts([Small RNA Counts])

    %% Styling
    classDef inputStyle fill:#e1f5ff,stroke:#01579b,stroke-width:2px
    classDef processStyle fill:#fff9c4,stroke:#f57f17,stroke-width:2px
    classDef decisionStyle fill:#f3e5f5,stroke:#4a148c,stroke-width:2px
    classDef outputStyle fill:#c8e6c9,stroke:#1b5e20,stroke-width:2px
    classDef dataStyle fill:#ffe0b2,stroke:#e65100,stroke-width:2px

    class Start,BAM_Input inputStyle
    class BC_FASTQ,BC_BAM,BAM2FASTQ,BAM2FASTQ2,NanoFilt1,NanoFilt2,BWA_Index,BWA_Align,MM2_Align,SAM2BAM1,SAM2BAM2,ProcessSAM,NanoPack_Plot,NanoPack_Stats,NanoCount processStyle
    class CheckMode,AlignMode decisionStyle
    class QC_Plots,QC_Stats,Counts outputStyle
    class FilteredFastq,AlignedBAM,SortedBAM dataStyle
```

## Pipeline Overview

### Input Requirements
- **Sample Info File**: Tab-separated file with sample metadata
- **Reference Genome**: FASTA format for alignment

### Workflow Stages

#### 1. Basecalling (3 modes)
- **use_fastq mode**: Direct FASTQ output from basecalling
- **skip_basecalling mode**: Start from existing BAM files
- **Standard mode**: BAM output from basecalling

#### 2. Quality Filtering
- Converts BAM to FASTQ (if needed)
- Filters reads by Q score using NanoFilt (default: Q > 10)

#### 3. Alignment (2 options)
- **BWA** (default): Optimized for Nanopore small RNA
  - Uses preset: `-x ont2d -C -W 13 -k 6`
- **Minimap2** (optional): Faster alternative
  - Uses preset: `-ax sr`

#### 4. Post-processing
- **process_sam**: Sort BAM and create index
- **nanopack_plot**: Generate QC visualizations
- **nanopack_stats**: Compute QC statistics
- **nanocount**: Quantify small RNA expression

### Output Files
- Sorted BAM files with indices
- QC plots and statistics
- Small RNA count matrices

## Key Parameters

```bash
--use_fastq         # Output FASTQ directly from basecalling
--skip_basecalling  # Start from existing BAM files
--use_minimap2      # Use minimap2 instead of BWA
--min_qscore 10     # Minimum Q score for filtering
```

## Usage Example

```bash
# Standard run with BWA
nextflow run workflows/small_rna_pipeline.nf \
    --sample_info samples.txt \
    --reference genome.fa

# Use minimap2 for faster alignment
nextflow run workflows/small_rna_pipeline.nf \
    --sample_info samples.txt \
    --reference genome.fa \
    --use_minimap2 true

# Skip basecalling, start from BAM files
nextflow run workflows/small_rna_pipeline.nf \
    --sample_info samples.txt \
    --reference genome.fa \
    --skip_basecalling true
```
