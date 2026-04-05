# HBV Genotype Reference Sequences

## Purpose

This directory contains the reference FASTA used by the pipeline for:
1. **MAFFT alignment** of consensus sequences
2. **BLASTN genotype pre-screen** (makeblastdb + blastn)

## Required File: `hbv_genotype_refs.fasta`

The pipeline expects a multi-FASTA file containing representative sequences for
HBV genotypes A through I (and preferably key subgenotypes).

### Format requirements

Each FASTA header **must** contain the genotype identifier in a parseable form.
Any of the following naming schemes will be recognised automatically:

```
>HBV_GENOTYPE_A_NC_003977
>HBV_A_AB010291
>JN642165_A             ← genotype after last underscore
>NC_003977.2_GenotypeA
```

### Recommended sources (download yourself)

1. **NCBI Reference Sequences (recommended)**

   Obtain one representative complete genome per genotype/subgenotype from:
   https://www.ncbi.nlm.nih.gov/nuccore/

   Suggested accessions (complete genomes, widely cited):

   | Genotype | Accession   | Notes                    |
   |----------|-------------|--------------------------|
   | A        | NC_003977.2 | Genotype A2 (standard)   |
   | B        | AB073858    | Genotype B2              |
   | C        | AB014381    | Genotype C2              |
   | D        | AB090270    | Genotype D               |
   | E        | X75657      | Genotype E               |
   | F        | X69798      | Genotype F               |
   | G        | AF160501    | Genotype G               |
   | H        | AY090454    | Genotype H               |
   | I        | AB298362    | Genotype I               |

   **Note:** For pre-S/S focused analysis, you may optionally trim each genome to
   only the pre-S1–pre-S2–S region (~nt 2848–835 wrapping the origin).

2. **HBVdb / HBV-genotype reference panel**
   https://hbvdb.ibcp.fr/HBVdb/HBVdbDownload

3. **Prebuilt panel (if available)**
   A prebuilt `hbv_genotype_refs.fasta` with verified headers will be
   provided in a future release. Until then, build yours from NCBI.

### Quick download example (Linux/WSL2)

```bash
# Install NCBI Entrez utilities
conda install -c bioconda entrez-direct

# Fetch representative genomes (adjust accessions as needed)
efetch -db nuccore -id NC_003977.2,AB073858,AB014381,AB090270,X75657,X69798,AF160501,AY090454,AB298362 \
       -format fasta > hbv_genotype_refs.fasta

# Rename headers to include genotype labels (edit manually or use sed)
# Example:
sed -i 's/>NC_003977.2 />HBV_GENOTYPE_A_NC_003977.2 /' hbv_genotype_refs.fasta
sed -i 's/>AB073858\.1 />HBV_GENOTYPE_B_AB073858 /'  hbv_genotype_refs.fasta
# ... etc.
```

## Required File: `obi_mutation_catalogue.tsv`

A TSV file mapping known OBI-associated surface antigen variants.
See `obi_mutation_catalogue.tsv` in this directory for the bundled catalogue.
