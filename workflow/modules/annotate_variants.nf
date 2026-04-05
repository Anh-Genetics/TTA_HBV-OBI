// workflow/modules/annotate_variants.nf
// Compares the per-sample consensus sequence to the HBV reference (same genotype)
// and annotates variants against the OBI-associated mutation catalogue.
//
// Inputs:  consensus_fasta [sample_id, fasta]
//          alignment_fasta [sample_id, fasta]  (from ALIGN_CONSENSUS)
//          mutation_db     [path]
// Outputs: annotation_tsv  [sample_id, tsv]

process ANNOTATE_VARIANTS {
    label 'low'
    tag   "${sample_id}"

    publishDir "${params.outdir}/annotations/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_fasta)
    tuple val(sample_id), path(alignment_fasta)
    path  mutation_db

    output:
    tuple val(sample_id), path("${sample_id}_variants.tsv"), emit: annotation_tsv

    script:
    """
    annotate_variants.py \\
        --sample-id     "${sample_id}" \\
        --consensus     "${consensus_fasta}" \\
        --alignment     "${alignment_fasta}" \\
        --mutation-db   "${mutation_db}" \\
        --out-tsv       "${sample_id}_variants.tsv"
    """
}
