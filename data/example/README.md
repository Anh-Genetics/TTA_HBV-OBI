# Example Data and Sample Sheet Templates

This directory contains templates and example files to help you structure your
input data for the TTA_HBV-OBI pipeline.

## Sample Sheet Format

The pipeline accepts a CSV file with the following columns:

| Column       | Required | Description                                     |
|--------------|----------|-------------------------------------------------|
| `sample_id`  | ✅ Yes   | Unique sample identifier (no spaces/special chars) |
| `ab1_forward`| ✅ Yes   | Path to the forward-direction `.ab1` file        |
| `ab1_reverse`| ❌ No    | Path to the reverse-direction `.ab1` file (leave blank for single-read) |

### Rules
- The header row is required.
- Column names are case-insensitive.
- Paths may be absolute or relative to the working directory.
- Duplicate `sample_id` values will cause a validation error.
- Files with extensions other than `.ab1` / `.abi` will generate a warning.

## Example Files

- `sample_sheet.csv` – Example sample sheet (with placeholder paths)
- `sample_sheet_template.csv` – Blank template to fill in

## Generating Synthetic Test Data

If you want to test the pipeline without real `.ab1` files, you can generate
synthetic FASTA files and skip the PARSE_ABI step. However, the MVP pipeline
expects `.ab1` files as input.

For a quick integration test using real data, obtain 2–3 `.ab1` files from a
Sanger run and fill in `sample_sheet.csv` with their paths.

## Expected Output Structure

```
results/
├── pipeline_info/
│   ├── timeline.html
│   ├── report.html
│   └── trace.txt
├── validation/
│   ├── validated_samplesheet.tsv
│   └── validation_report.txt
├── parsed_reads/
│   └── {sample_id}/
│       ├── {sample_id}_forward.fasta
│       ├── {sample_id}_forward_quality.tsv
│       ├── {sample_id}_reverse.fasta       [if paired]
│       └── {sample_id}_reverse_quality.tsv [if paired]
├── trimmed_reads/
│   └── {sample_id}/
│       ├── {sample_id}_forward_trimmed.fasta
│       └── {sample_id}_forward_trim_stats.tsv
├── consensus/
│   └── {sample_id}/
│       ├── {sample_id}_consensus.fasta
│       └── {sample_id}_consensus_stats.tsv
├── alignments/
│   └── {sample_id}/
│       └── {sample_id}_aligned.fasta
├── genotypes/
│   └── {sample_id}/
│       └── {sample_id}_genotype.tsv
├── annotations/
│   └── {sample_id}/
│       └── {sample_id}_variants.tsv
└── reports/
    └── {sample_id}/
        └── {sample_id}_report.html
```
