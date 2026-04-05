// workflow/modules/genotype_assign.nf
// Assigns HBV genotype to a consensus sequence by:
//   1. BLASTN against a reference panel (fast pre-screen)
//   2. Parsing the alignment result to select the best genotype
//
// Inputs:  consensus_fasta [sample_id, fasta]
//          reference_fasta [path]
// Outputs: genotype_tsv    [sample_id, tsv]

process ASSIGN_GENOTYPE {
    label 'medium'
    tag   "${sample_id}"

    publishDir "${params.outdir}/genotypes/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_fasta)
    path  reference_fasta

    output:
    tuple val(sample_id), path("${sample_id}_genotype.tsv"), emit: genotype_tsv

    script:
    """
    # Check BLAST is available
    if ! command -v blastn &>/dev/null; then
        echo "[WARN] blastn not found – genotype pre-screen will be skipped." >&2
        echo -e "sample_id\tgenotype\tsubgenotype\tpct_identity\tblast_hit\tconfidence\tnote" \\
            > "${sample_id}_genotype.tsv"
        echo -e "${sample_id}\tUNKNOWN\tUNKNOWN\tNA\tNA\tLOW\tblastn_not_found" \\
            >> "${sample_id}_genotype.tsv"
        exit 0
    fi

    # Build a temporary BLAST database from the reference FASTA
    makeblastdb \\
        -in    "${reference_fasta}" \\
        -dbtype nucl \\
        -out   ref_db \\
        -title HBV_genotype_refs \\
        2>&1 | grep -v "^Building" || true

    # Run BLASTN
    blastn \\
        -query         "${consensus_fasta}" \\
        -db            ref_db \\
        -out           blast_raw.txt \\
        -outfmt        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
        -evalue        ${params.blast_evalue} \\
        -max_target_seqs 10 \\
        -num_threads   ${task.cpus}

    # Parse BLAST output to assign genotype
    genotype_assign.py \\
        --blast-result  blast_raw.txt \\
        --sample-id     "${sample_id}" \\
        --min-pct-id    ${params.min_genotype_pct} \\
        --out-tsv       "${sample_id}_genotype.tsv"
    """
}
