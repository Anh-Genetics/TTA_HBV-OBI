// workflow/modules/build_consensus.nf
// Builds a per-sample consensus FASTA from one (forward only) or two (F+R)
// trimmed reads. Uses a Python helper that:
//   - with a single read: uses that read as the consensus
//   - with two reads: attempts overlap assembly; falls back to longest read
//     with IUPAC ambiguity codes at conflicting positions
//
// Inputs:  [sample_id, [list_of_trimmed_fastas]]
//          (grouped by sample_id using groupTuple in main.nf)
// Outputs: consensus_fasta [sample_id, fasta]
//          consensus_stats [sample_id, tsv]

process BUILD_CONSENSUS {
    label 'low'
    tag   "${sample_id}"

    publishDir "${params.outdir}/consensus/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_fastas)

    output:
    tuple val(sample_id), path("${sample_id}_consensus.fasta"), emit: consensus_fasta
    tuple val(sample_id), path("${sample_id}_consensus_stats.tsv"), emit: consensus_stats

    script:
    // Pass all fasta files as space-separated list
    """
    build_consensus.py \\
        --sample-id      "${sample_id}" \\
        --reads          ${trimmed_fastas} \\
        --min-overlap    ${params.min_overlap} \\
        --conflict       ${params.conflict_policy} \\
        --out-fasta      "${sample_id}_consensus.fasta" \\
        --out-stats      "${sample_id}_consensus_stats.tsv"
    """
}
