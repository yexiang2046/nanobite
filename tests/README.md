# NanoBite Pipeline Tests

## Running Tests

### Unit Tests (Common Module)
```bash
nextflow run tests/test_common.nf
```

Tests utility functions:
- `extract_sample_id()` - sample ID extraction from BAM filenames
- `create_bam_channel()` - BAM channel creation with sample IDs

### Integration Tests (Small RNA Pipeline)
```bash
nextflow run tests/test_small_rna_pipeline.nf
```

Tests pipeline workflow structure:
- Sample channel creation
- Process input/output handling
- Data flow between processes

## Test Structure

```
tests/
├── README.md                      # This file
├── test_common.nf                 # Unit tests for common.nf
├── test_small_rna_pipeline.nf     # Integration test for pipeline
├── fixtures/                      # Test data fixtures
│   └── sample_info.txt           # Mock sample info file
└── data/                         # Mock data directory
    └── pod5/                     # Mock POD5 directories
```

## Adding New Tests

1. Create test workflow in `tests/test_*.nf`
2. Use mock processes to avoid heavy computation
3. Validate outputs with assertions or view commands
4. Document in this README
