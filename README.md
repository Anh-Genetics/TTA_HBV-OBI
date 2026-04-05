# TTA_HBV-OBI — HBV pre-S1/pre-S2/S Sanger Analysis Pipeline

> **Status: MVP v0.1.0** — Correct, runnable skeleton. See [Limitations](#limitations) before use.

A Nextflow DSL2 pipeline for analysis of HBV Sanger sequencing data covering
the **pre-S1 / pre-S2 / S** region, targeting:

1. **HBV genotype assignment** (genotypes A–I) by BLAST-based alignment screen
2. **OBI-associated variant annotation** against a curated mutation catalogue
3. **Structured per-sample reports** (HTML / TSV / Markdown)

---

## Table of Contents

- [Background](#background)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Input Format](#input-format)
- [Running the Pipeline](#running-the-pipeline)
- [Expected Outputs](#expected-outputs)
- [Pipeline Architecture](#pipeline-architecture)
- [OBI Mutation Catalogue](#obi-mutation-catalogue)
- [Limitations](#limitations)
- [Reviewing Flagged Samples](#reviewing-flagged-samples)
- [Future Work](#future-work)
- [References](#references)

---

## Background

Occult hepatitis B infection (OBI) is defined by the presence of HBV DNA in serum
or liver tissue in individuals who test negative for HBsAg. Sanger sequencing of the
**pre-S1 / pre-S2 / S** region (~1200-1300 bp) is widely used to:

- Assign HBV genotype (A-I), which influences clinical behaviour and antiviral response
- Identify surface antigen variants that may cause **immune escape** or **diagnostic escape**
  (false-negative HBsAg tests), contributing to the OBI phenotype

This pipeline automates the bioinformatics steps from raw `.ab1` chromatogram files
through to a structured variant interpretation report.

---

## Prerequisites

### System
- Linux or WSL2 (Ubuntu 20.04 / 22.04 / 24.04)
- >= 4 GB RAM, >= 4 CPU cores recommended

### Software
| Tool | Version | Purpose |
|------|---------|---------|
| Java | >= 17 | Required by Nextflow |
| Nextflow | >= 23.04 | Workflow engine |
| micromamba / conda | any | Python environment management |
| MAFFT | >= 7.5 | Multiple sequence alignment |
| BLAST+ | >= 2.14 | Genotype pre-screen |
| Python | 3.10-3.12 | Helper scripts |
| Biopython | >= 1.81 | ABI parsing |

---

## Installation

### 1. Install Java (required for Nextflow)

```bash
sudo apt update && sudo apt install -y default-jre-headless
java -version   # should show >= 17
```

### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
nextflow -version
```

### 3. Install micromamba (recommended) or conda

```bash
# micromamba (fast, lightweight)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# Or use miniconda: https://docs.conda.io/en/latest/miniconda.html
```

### 4. Create the conda environment

```bash
# Clone this repository
git clone https://github.com/Anh-Genetics/TTA_HBV-OBI.git
cd TTA_HBV-OBI

# Create environment (takes 3-8 minutes)
micromamba create -f environment.yml -y
micromamba activate hbv-obi

# Verify tools
python --version                            # Python 3.11
python -c "import Bio; print(Bio.__version__)"  # >= 1.81
mafft --version
blastn -version
```

### 5. Prepare reference data

The pipeline requires a multi-FASTA reference file with representative HBV
genotype sequences. See `data/references/PLACEHOLDER.md` for full instructions.

**Quick setup with NCBI Entrez utilities:**

```bash
micromamba activate hbv-obi

# Download representative genomes (genotypes A through I)
efetch -db nuccore \
  -id NC_003977.2,AB073858.1,AB014381.1,AB090270.1,X75657.1,X69798.1,AF160501.1,AY090454.1,AB298362.1 \
  -format fasta > data/references/hbv_genotype_refs_raw.fasta

# Rename headers to include genotype tags (required for pipeline genotype parsing)
python - << 'PYEOF'
import sys
mapping = {
    "NC_003977": "HBV_GENOTYPE_A_NC_003977.2",
    "AB073858":  "HBV_GENOTYPE_B_AB073858",
    "AB014381":  "HBV_GENOTYPE_C_AB014381",
    "AB090270":  "HBV_GENOTYPE_D_AB090270",
    "X75657":    "HBV_GENOTYPE_E_X75657",
    "X69798":    "HBV_GENOTYPE_F_X69798",
    "AF160501":  "HBV_GENOTYPE_G_AF160501",
    "AY090454":  "HBV_GENOTYPE_H_AY090454",
    "AB298362":  "HBV_GENOTYPE_I_AB298362",
}
with open("data/references/hbv_genotype_refs_raw.fasta") as fin, \
     open("data/references/hbv_genotype_refs.fasta", "w") as fout:
    for line in fin:
        if line.startswith(">"):
            acc = line[1:].split()[0].split(".")[0]
            new_id = mapping.get(acc, acc)
            rest = " ".join(line[1:].split()[1:])
            fout.write(f">{new_id} {rest}\n")
        else:
            fout.write(line)
PYEOF

echo "Reference FASTA ready:"
grep "^>" data/references/hbv_genotype_refs.fasta
```

---

## Input Format

### Sample Sheet (CSV)

Provide a CSV file with the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | Yes | Unique sample identifier (no spaces) |
| `ab1_forward` | Yes | Absolute or relative path to forward `.ab1` file |
| `ab1_reverse` | No | Path to reverse `.ab1` file (leave blank for single-read) |

**Example:**

```csv
sample_id,ab1_forward,ab1_reverse
PATIENT001,/data/sanger/PATIENT001_F.ab1,/data/sanger/PATIENT001_R.ab1
PATIENT002,/data/sanger/PATIENT002_F.ab1,/data/sanger/PATIENT002_R.ab1
PATIENT003_single,/data/sanger/PATIENT003_F.ab1,
```

See `data/example/sample_sheet_template.csv` for a blank template.

### ABI Files (`.ab1`)

Standard Sanger sequencing output files from Applied Biosystems instruments
(3130, 3500, 3730 series). Files must contain embedded Phred quality scores
(standard for all modern ABI instruments).

---

## Running the Pipeline

### Basic run

```bash
micromamba activate hbv-obi

nextflow run main.nf \
    --sample_sheet data/example/sample_sheet.csv \
    --outdir results/ \
    -profile local
```

### All parameters

```
--sample_sheet       Path to CSV sample sheet (required)
--outdir             Output directory [results]
--reference_fasta    HBV genotype reference FASTA [data/references/hbv_genotype_refs.fasta]
--mutation_db        OBI mutation catalogue TSV [data/references/obi_mutation_catalogue.tsv]
--min_quality        Phred quality threshold for base trimming [20]
--min_length         Minimum sequence length after trimming (bp) [200]
--trim_ends_bases    Fixed bases to trim from each end [20]
--min_overlap        Minimum F/R overlap to attempt paired consensus (bp) [80]
--conflict_policy    How to resolve F/R conflicts: iupac|majority|n [iupac]
--blast_evalue       BLAST e-value cutoff [1e-10]
--min_genotype_pct   Minimum % identity to call genotype [90.0]
--report_format      Report format: html|tsv|markdown [html]
```

### Resume a failed run

```bash
nextflow run main.nf --sample_sheet sample_sheet.csv --outdir results/ -resume
```

---

## Expected Outputs

```
results/
+-- pipeline_info/
|   +-- timeline.html         Process execution timeline
|   +-- report.html           Nextflow resource report
|   +-- trace.txt             Per-task resource trace
+-- validation/
|   +-- validated_samplesheet.tsv
|   +-- validation_report.txt
+-- parsed_reads/{sample_id}/
|   +-- {sample_id}_forward.fasta
|   +-- {sample_id}_forward_quality.tsv
|   +-- (reverse equivalents if paired)
+-- trimmed_reads/{sample_id}/
|   +-- {sample_id}_forward_trimmed.fasta
|   +-- {sample_id}_forward_trim_stats.tsv
+-- consensus/{sample_id}/
|   +-- {sample_id}_consensus.fasta      << main output: per-sample consensus
|   +-- {sample_id}_consensus_stats.tsv
+-- alignments/{sample_id}/
|   +-- {sample_id}_aligned.fasta        MAFFT alignment vs references
+-- genotypes/{sample_id}/
|   +-- {sample_id}_genotype.tsv         Genotype assignment
+-- annotations/{sample_id}/
|   +-- {sample_id}_variants.tsv         OBI variant annotation
+-- reports/{sample_id}/
    +-- {sample_id}_report.html          Final report per sample
```

### Key output: genotype TSV

| Column | Description |
|--------|-------------|
| `genotype` | HBV genotype (A-I or UNKNOWN) |
| `subgenotype` | Subgenotype if determined |
| `pct_identity` | Best BLAST hit % identity |
| `blast_hit` | Best BLAST hit sequence ID |
| `confidence` | HIGH / MEDIUM / LOW |
| `note` | Warnings or flags |

### Key output: variant TSV

| Column | Description |
|--------|-------------|
| `aa_position` | S-protein amino acid position (1-based) |
| `ref_aa` | Reference amino acid |
| `alt_aa` | Variant amino acid |
| `notation` | Standard notation (e.g. sG145R) |
| `domain` | HBV domain (a_determinant, MHR, S, preS1, preS2) |
| `mechanism` | Proposed OBI mechanism |
| `evidence_level` | A (strong) / B (moderate) / C (weak/novel) |
| `notes` | Literature references |
| `catalogued` | yes / no |

---

## Pipeline Architecture

```
main.nf
|
+-- VALIDATE_SAMPLESHEET    (validate_samplesheet.py)
+-- PARSE_ABI               (parse_abi.py)              -- per read
+-- TRIM_READS              (trim_reads.py)             -- per read
+-- BUILD_CONSENSUS         (build_consensus.py)        -- per sample (grouped)
+-- ALIGN_CONSENSUS         (MAFFT)                     -- per sample
+-- ASSIGN_GENOTYPE         (BLAST + genotype_assign.py)-- per sample
+-- ANNOTATE_VARIANTS       (annotate_variants.py)      -- per sample
+-- GENERATE_REPORT         (generate_report.py)        -- per sample
```

**Key design decisions:**

- **No race conditions**: Each process writes only to its own Nextflow work directory.
  The per-sample `groupTuple` in `main.nf` ensures all reads for a sample are
  collected before consensus building begins.
- **Single-read graceful handling**: Samples with only a forward read produce a
  valid single-read consensus; no crash or silent skip.
- **IUPAC ambiguity**: Conflicting positions in F/R overlap are encoded as IUPAC
  codes by default (configurable via `--conflict_policy`).
- **Deterministic output**: All outputs are written to isolated directories per
  sample; no shared files are appended during parallel execution.

---

## OBI Mutation Catalogue

The bundled catalogue (`data/references/obi_mutation_catalogue.tsv`) contains
~50 well-documented surface antigen variants associated with OBI, grouped into:

| Evidence Level | Meaning |
|---------------|---------|
| **A** (Strong) | Reported in multiple OBI cohorts; functional or clinical validation available |
| **B** (Moderate) | Reported in OBI patients; limited functional data |
| **C** (Weak/Novel) | Single report, or variant not in catalogue (novel) |

**Key mutation groups covered:**

1. **MHR / a-determinant escape mutations** (S aa 99-169): G145R/A, D144A/E/N,
   T126I/S, Q129R/P, M133T, Y134C, P120L, C124R/Y, R122P, etc.
2. **Premature stop codons in S**: Q181* and similar
3. **Pre-S2 initiation codon mutations**: M1I, M1T, M1V (loss of pre-S2 protein)

> WARNING: The catalogue is a research reference only. All findings must be
> interpreted in clinical context (HBsAg serology, HBV DNA quantification,
> patient history). This tool does NOT diagnose OBI.

---

## Limitations

This is an MVP (v0.1.0). The following limitations apply:

1. **Genotype assignment**: BLAST-based pre-screen only. No maximum-likelihood
   phylogenetic tree is built automatically. For publication-quality genotype
   assignment, run IQ-TREE or MEGA on the `alignments/` output manually.

2. **Subgenotype reliability**: Subgenotype assignment from a partial sequence
   (~1245 bp) may be unreliable, especially for inter-genotype recombinants.
   Treat subgenotype calls with caution.

3. **OBI domain coordinates**: Pre-S1/pre-S2/S domain boundaries are approximate,
   based on genotype A (NC_003977.2). Alignment-based remapping is used for
   other genotypes.

4. **Quasispecies / minor variants**: Sanger sequencing cannot reliably detect
   minor variants below ~15-20% frequency. Dual-peak positions are encoded as
   IUPAC codes, but low-level mixed infections may be missed.

5. **OBI diagnosis**: This tool identifies variants *associated* with OBI in the
   literature. It does NOT diagnose OBI. OBI diagnosis requires:
   - HBsAg negativity (laboratory)
   - HBV DNA positivity (PCR)
   - Clinical correlation

6. **Biological completeness**: The variant catalogue covers published literature
   but is not exhaustive. Novel variants are flagged as Evidence Level C.

---

## Reviewing Flagged Samples

When the report flags a sample for review:

1. **Open the HTML report**: `results/reports/{sample_id}/{sample_id}_report.html`
2. **Check read QC**: Were both forward and reverse reads above quality threshold?
   Low-quality reads or FAIL status require re-sequencing.
3. **Verify genotype**: Is genotype assignment HIGH confidence? LOW confidence
   genotypes should be confirmed with IQ-TREE phylogenetic analysis.
4. **Review Evidence Level A variants**: These are most clinically significant
   (e.g. G145R, D144A, Q129R). Verify they appear in both F and R reads.
5. **Check for premature stop codons**: Strong indicators of HBsAg truncation.
6. **Examine the alignment**: Open `{sample_id}_aligned.fasta` in SeaView or
   Jalview to manually inspect variant positions.
7. **Cross-reference serology**: Correlate with HBsAg quantitative assay results.

---

## Development and Testing

```bash
# Run Python unit tests
micromamba activate hbv-obi
python -m pytest tests/ -v

# Run a specific test class
python -m pytest tests/test_pipeline_scripts.py::TestBuildConsensus -v

# Validate individual scripts
python workflow/bin/validate_samplesheet.py --help
python workflow/bin/parse_abi.py --help
python workflow/bin/build_consensus.py --help
```

---

## Future Work

- [ ] Full ML phylogenetic tree: integrate IQ-TREE automatically
- [ ] Expanded subgenotype reference panel
- [ ] Recombination detection (RDP5 / jpHMM integration)
- [ ] NGS mode: support Illumina amplicon paired-end FASTQ input
- [ ] Pre-S deletion detection (structural variant calling)
- [ ] N-glycosylation motif gain/loss screening
- [ ] Docker/Singularity container for fully reproducible runs
- [ ] Multi-sample aggregate summary report
- [ ] Polymerase drug resistance mutation annotation

---

## References

Key literature used to curate the OBI mutation catalogue:

1. Candotti D & Allain JP (2019). Occult Hepatitis B Virus Infection. PMC6784188
2. Gu C et al. (2021). OBI-associated HBsAg variants. PMC9325327
3. Mulrooney-Cousins PM & Michalak TI (2016). Pre-S mutations and OBI. PMC4868315
4. Yon M et al. (2020). HBV genotyping from partial sequences. PMC7416611
5. Wang C et al. (2021). Pre-S1/pre-S2/S Sanger sequencing for HBV OBI. PMC8537069
6. Mak LY et al. (2024). Hyperglycosylation insertion and HBsAg escape. PMC10144012

---

## Licence

This project is for research use only. No clinical diagnostic claims are made.
