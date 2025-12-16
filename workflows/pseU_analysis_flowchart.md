# Pseudouridine (ψ) Modification Analysis Pipeline Flowchart

```mermaid
flowchart TD
    %% Input
    Start([POD5 files]) --> basecalling{Dorado<br/>pseU basecalling}

    %% Basecalling with modification detection
    basecalling --> BAM[BAM with<br/>ψ modifications<br/>filtering Q > 10?]

    %% Alignment with Q score filtering
    BAM --> dorado_align[dorado aligner]
    dorado_align --> ALIGNED_BAM[Aligned BAM]

    %% Sort and index
    ALIGNED_BAM --> samtools_sort[samtools sort<br/>& index]
    samtools_sort --> SortedBAM[Sorted BAM + BAI]

    %% Parallel modkit analysis
    minimappedBAM --> modkit_pileup[modkit pileup<br/>genome-wide sites]
    minimappedBAM --> modkit_extract[modkit extract<br/>per-read calls]
    minimappedBAM --> modkit_summary[modkit summary<br/>statistics]

    %% Pileup processing
    modkit_pileup --> PILEUP[Pileup BED]
    PILEUP --> filter_pseU[filter pseU sites]

    %% Outputs
    filter_pseU --> PseU_Sites([ψ Sites BED])
    modkit_extract --> Mod_Calls([Modification TSV])
    modkit_summary --> Mod_Stats([Summary TXT])

    %% Styling
    classDef inputStyle fill:#e1f5ff,stroke:#01579b,stroke-width:2px
    classDef processStyle fill:#fff9c4,stroke:#f57f17,stroke-width:2px
    classDef dataStyle fill:#ffe0b2,stroke:#e65100,stroke-width:2px
    classDef outputStyle fill:#c8e6c9,stroke:#1b5e20,stroke-width:2px

    class Start inputStyle
    class basecalling,dorado_align,samtools_sort,modkit_pileup,modkit_extract,modkit_summary,filter_pseU processStyle
    class BAM,ALIGNED_BAM,SortedBAM,PILEUP dataStyle
    class PseU_Sites,Mod_Calls,Mod_Stats outputStyle
```

## Pipeline Overview

### Input Requirements
- **POD5 Files**: Raw sequencing data from Oxford Nanopore
- **Reference Genome**: FASTA format for alignment

### Workflow Stages

#### 1. Basecalling with Modification Detection
- **basecalling_pseU**: Dorado basecaller with `--modified-bases pseU`
- Detects pseudouridine (ψ) modifications during basecalling
- Output: BAM files with modification tags

#### 2. Alignment with Quality Filtering
- **dorado_align**: Aligns basecalled BAM to reference
- Filters reads by Q score using `--min-qscore` (default: 10)
- Preserves modification information throughout alignment

#### 3. Post-processing
- **samtools_sort**: Sort BAM by coordinates
- **samtools index**: Create BAI index

#### 4. Modification Analysis (3 parallel processes)

##### a. Genome-wide Modification Sites
- **modkit_pileup**: Generates per-site modification frequencies
  - Minimum coverage: 5 (default)
  - Probability threshold: 0.8 (default)
- **filter_pseU**: Extracts only pseudouridine sites from pileup

##### b. Per-read Modification Calls
- **modkit_extract**: Extracts read-level modification information
- Output: TSV with modification calls for each read

##### c. Summary Statistics
- **modkit_summary**: Generates overall modification statistics
- Output: Summary text file

### Output Files

```
results/
├── basecalling_pseU/     # Basecalled BAM files with ψ modifications
├── alignment/            # Dorado aligned BAM files (Q score filtered)
├── bam/                  # Sorted BAM files and indices
├── modkit_pileup/        # Per-site modification frequencies (BED)
├── pseU_sites/           # Filtered pseudouridine sites (BED)
├── modkit_extract/       # Read-level modification calls (TSV)
└── modkit_summary/       # Summary statistics (TXT)
```

## Key Parameters

```bash
--min_qscore 10         # Minimum Q score for read filtering during alignment
--min_coverage 5        # Minimum coverage for modkit pileup
--prob_threshold 0.8    # Probability threshold for modification calls
--mm2opts "-ax sr"      # Minimap2 options for dorado aligner
```

## Usage Example

```bash
# Standard run
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference genome.fa

# Custom parameters
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference genome.fa \
    --min_qscore 15 \
    --min_coverage 10 \
    --prob_threshold 0.9

# Skip basecalling (start from existing BAM files)
nextflow run workflows/pseU_analysis.nf \
    --sample_info sample_info.txt \
    --reference genome.fa \
    --skip_basecalling true \
    --bam_dir /path/to/bams
```

## Notes

- The pipeline uses dorado for both basecalling AND alignment
- Modification information is preserved in BAM tags throughout the workflow
- All three modkit processes run in parallel for efficiency
- Q score filtering happens during alignment (not as a separate step)
