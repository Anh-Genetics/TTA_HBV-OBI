// workflow/modules/align_consensus.nf
// Aligns the per-sample consensus to the HBV genotype reference set using MAFFT.
// Inputs:  consensus_fasta [sample_id, fasta]
//          reference_fasta [path]
// Outputs: alignment_fasta [sample_id, fasta]

process ALIGN_CONSENSUS {
    label 'medium'
    tag   "${sample_id}"

    publishDir "${params.outdir}/alignments/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_fasta)
    path  reference_fasta

    output:
    tuple val(sample_id), path("${sample_id}_aligned.fasta"), emit: alignment_fasta

    script:
    """
    # Check MAFFT is available
    if ! command -v mafft &>/dev/null; then
        echo "[ERROR] mafft not found. Install with: conda install -c bioconda mafft" >&2
        exit 1
    fi

    # Combine consensus with reference sequences, then align
    cat "${reference_fasta}" "${consensus_fasta}" > combined_input.fasta

    mafft \\
        --auto \\
        --preservecase \\
        --thread ${task.cpus} \\
        --quiet \\
        combined_input.fasta \\
        > "${sample_id}_aligned.fasta"
    """
}
